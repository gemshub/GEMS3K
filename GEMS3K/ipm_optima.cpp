//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_optima.cpp
/// Chemical equilibrium via the Optima library's primal-dual interior-point NLP solver, as an
/// alternative to GEMS3K's own IPM/MBR loop. GEMS3K's G0[] and activity-coefficient models
/// provide the objective; Optima owns the Newton iteration and the box/linear-equality
/// handling.
///
/// Modes dispatched here from TNode::GEM_run(): AOP (cold start) and SOP (warm start), both
/// through CalculateEquilibriumStateOptima(), which reads pm.pNP; ROP (reference setup); and
/// HOP/SHP (native first, then Optima) through CalculateEquilibriumStateHOP(). See NODECODECH
/// in databr.h.
///
/// Optional control conditions (pH, Eh) are solved in the same Newton system as extra
/// unknowns with a fixed objective gradient (see ipm_optima.h). With none registered this is
/// a plain equilibrium solve.
//
// Copyright (c) 1992-2026
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//

#include "gems3k_timing.h"
#include <chrono>
#include "ms_multi.h"

#ifdef USE_OPTIMA_SOLVER

#include "node.h"
#include <Optima/Optima.hpp>
#include <Optima/Sensitivity.hpp>
#include <Optima/Trace.hpp>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <sstream>

// pa_OptimaCgSeed is a second attempt: the column-generation seed is used only while TNode
// arms it for the retry of a cold AOP call that failed from the ordinary LP-feasibility seed.
// Optima's solve(), timed when GEMS3K_TIMING_FILE is set.
template<class S, class... A>
static Optima::Result timed_solve( S& solver, A&&... args )
{
    gems3k_timing::Scope timing( gems3k_timing::OSolve );
    return solver.solve( std::forward<A>( args )... );
}

extern thread_local bool g_ipmNoProgFired;   // pa_IpmNoProgress (ipm_main.cpp)
extern thread_local bool g_optimaCgSeedArmed;   // defined in node.cpp (needed without USE_OPTIMA_SOLVER)
extern thread_local bool g_optimaTargetWarned;  // defined in node.cpp
extern thread_local long int g_optimaIterBudget;  // defined in node.cpp; -1 = no limit
extern thread_local long int g_optimaCapAllowance;   // defined in node.cpp
extern thread_local bool g_optimaDryFailed;   // defined in node.cpp
extern thread_local bool g_optimaDryOff;      // defined in node.cpp
extern thread_local bool g_optimaDryRedoCall; // defined in node.cpp
extern thread_local bool g_optimaPreSolveEmpty;   // defined in node.cpp; set when the pre-solve found nothing
// pa_OptimaPreSolveSeed: set when a seeded first pre-solve pass missed its trial; the call then skips the pre-solve.
static thread_local bool g_ps0GaveUp = false;
// pa_OptimaColdLevel: which leveled cold start this call used (0 = none).
static thread_local int g_coldLvlUsed = 0;
// pa_OptimaLevelStart: the column-generation LP's last round, and the certified start it gives (ipm_level.cpp).
extern thread_local bool g_cgConverged;
extern thread_local long int g_cgPseudo;
extern thread_local bool g_cgGap;
extern thread_local std::vector<double> g_cgDual;
static thread_local bool g_lvlShortcut = false;
// The leveled answer holds mixed phases: start from it, but solve it the ordinary way.
static thread_local bool g_lvlCgStart = false;
extern thread_local bool g_tpdStrongSearch;   // defined in ipm_main.cpp
extern thread_local bool g_aqueousAbsent;     // defined in ipm_chemical2.cpp; set when the aqueous phase has no water
extern thread_local long int g_lvlFullCalls;  // defined in node.cpp
extern thread_local bool g_optimaFixAccept;       // defined in node.cpp
// Set by the caller when PotentialSpaceFinish() starts from a call Optima reported OK. From
// such a start the finish may only improve: it must not raise G, worsen the mass balance or
// drop a present phase.
thread_local bool g_finishFromSuccess = false;


// pa_OptimaLineSearch: put Optima's merit line search on the unmasked error at the given
// trigger factor (needs the Optima fork's ErrorControl::execute). A negative value makes the
// line search a second attempt: the first Optima attempt runs without it, and
// TNode::GEM_run() re-runs a failed call once with g_optimaLineSearchRetry set, at |factor|.
extern thread_local bool g_optimaLineSearchRetry;   // defined in node.cpp (needed without USE_OPTIMA_SOLVER)
// Set by TNode::GEM_run() for the re-run of a failed Optima call with pa_OptimaLSStallEscape
// switched off (DECIDE escretry).
extern thread_local bool g_optimaLSEscapeOff;   // defined in node.cpp (needed without USE_OPTIMA_SOLVER)
extern thread_local bool g_optimaEscRetryNext;   // defined in node.cpp
// Objective memo for the line search. Optima's line-search trigger evaluates the objective at
// the new point, and the solver then evaluates it again at the same point. GEMS3K's objective
// is not a pure function: CalculateActivityCoefficients(LINK_UX_MODE) accumulates lnGmo and
// blends F0 through FitVar[3], so each extra call advances that history. With the line search
// on, a repeat call at a bit-identical x returns the stored result. Inactive when the line
// search is off.
static void memoize_objective_if_linesearch( Optima::Problem& problem, double factor )
{
    const bool lsOn = factor > 0. || ( factor < 0. && g_optimaLineSearchRetry );
    // g_optimaLSEscapeOff marks TNode::GEM_run()'s legacy retry: stall escape off and
    // memo off.
    if( !lsOn || g_optimaLSEscapeOff || !problem.f.initialized() )
        return;
    struct Memo { bool valid = false, hasFxx = false; Optima::Vector x; Optima::ObjectiveResult r; };
    auto memo = std::make_shared<Memo>();
    auto base = problem.f;
    problem.f = [base, memo]( Optima::ObjectiveResultRef res, Optima::VectorView x, Optima::VectorView p,
                              Optima::VectorView c, Optima::ObjectiveOptions opts )
    {
        if( memo->valid && memo->x.size() == x.size() && ( memo->x.array() == x.array() ).all()
            && ( !opts.eval.fxx || memo->hasFxx ) )
        {
            res.f = memo->r.f;
            res.fx = memo->r.fx;
            if( opts.eval.fxx ) res.fxx = memo->r.fxx;
            res.diagfxx = memo->r.diagfxx;
            res.fxx4basicvars = memo->r.fxx4basicvars;
            res.succeeded = memo->r.succeeded;
            return;
        }
        base( res, x, p, c, opts );
        memo->x = x;
        memo->r.f = res.f;
        memo->r.fx = res.fx;
        memo->hasFxx = opts.eval.fxx;
        if( opts.eval.fxx ) memo->r.fxx = res.fxx;
        memo->r.diagfxx = res.diagfxx;
        memo->r.fxx4basicvars = res.fxx4basicvars;
        memo->r.succeeded = res.succeeded;
        memo->valid = true;
    };
}

// The debug-trace setting GEMS3K_OPTIMA_STEPTRACE_<name>, or its older name REAKTORO_TRACE_<name>; null when unset.
static const char* optima_trace_env( const char* name )
{
    if( const char* e = std::getenv( ( std::string( "GEMS3K_OPTIMA_STEPTRACE_" ) + name ).c_str() ) )
        return e;
    return std::getenv( ( std::string( "REAKTORO_TRACE_" ) + name ).c_str() );
}

// Debug trace, off unless GEMS3K_OPTIMA_STEPTRACE_FILE is set. True when the Hessian matrices are to be
// written for the current Optima iteration (GEMS3K_OPTIMA_STEPTRACE_HESSIAN=<iteration>,
// GEMS3K_OPTIMA_STEPTRACE_HESSIAN_SOLVE=<solve>, default 1).
static bool optima_trace_hessian_dump_now()
{
    static const long it = []{ const char* e = optima_trace_env( "HESSIAN" ); return e ? std::atol(e) : -1L; }();
    static const long solve = []{ const char* e = optima_trace_env( "HESSIAN_SOLVE" ); return e ? std::atol(e) : 1L; }();
    const auto& ts = Optima::trace::state();
    static const int attempt = []{ const char* e = optima_trace_env( "HESSIAN_ATTEMPT" ); return e ? std::atoi(e) : 1; }();
    return Optima::trace::file() && it == ts.iter && solve == ts.solve && ts.attempt == attempt;
}

// Writes an n-by-n matrix as comma-separated rows next to the trace file.
template<class M>
static void optima_trace_write_matrix( const char* tag, const M& m, long n )
{
    const auto& ts = Optima::trace::state();
    const std::string path = std::string( optima_trace_env( "FILE" ) ) + ".hess." + std::to_string(ts.solve)
                           + "." + std::to_string(ts.iter) + "." + tag + ".csv";
    std::ofstream f( path );
    f.precision( 17 );
    for( long r = 0; r < n; r++ )
        for( long c = 0; c < n; c++ )
            f << m(r,c) << ( c + 1 < n ? "," : "\n" );
}

// pa_OptimaPartialPivLU: the smallest linear system Optima factorises with partial pivoting (0 = never).
static void apply_optima_lu( Optima::Options& o, long int minSize )
{
#ifdef OPTIMA_LU_PARTIAL_PIVOT   // needs an Optima build with this option
    o.newtonstep.linearsolver.lu_partial_min_size = minSize > 0 ? minSize : 0;
#else
    (void)o; (void)minSize;
#endif
}

// pa_OptimaSnapLeaving: how many limited steps in a row before a species heading to zero is put on its floor.
static void apply_optima_snap( Optima::Options& o, long int snapAfter )
{
#ifdef OPTIMA_BACKTRACK_SNAP_LEAVING   // needs an Optima build with this option
    o.backtracksearch.snap_leaving_after = snapAfter > 0 ? snapAfter : 0;
#else
    (void)o; (void)snapAfter;
#endif
}

static void apply_optima_linesearch( Optima::Options& o, double factor, long int stallEscape = 0,
                                     long int rejectWorse = 0 )
{
    if( factor < 0. )
    {
        if( !g_optimaLineSearchRetry ) return;   // first attempt: no line search
        factor = -factor;
    }
    if( !( factor > 0. ) ) return;
    o.linesearch.enabled = true;
    o.linesearch.use_unmasked_error = true;
#ifdef OPTIMA_LINESEARCH_STALL_ESCAPE   // pa_OptimaLSStallEscape needs an Optima build with this field
    o.linesearch.stall_escape_after = ( stallEscape > 0 && !g_optimaLSEscapeOff ) ? (std::size_t)stallEscape : 0;
#else
    (void)stallEscape;
#endif
#ifdef OPTIMA_LINESEARCH_REJECT_WORSE
    o.linesearch.reject_if_worse = rejectWorse > 0 && !g_optimaLSEscapeOff;
#else
    (void)rejectWorse;
#endif
    o.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = factor;
}

namespace {
// Symmetric eigenvalue floor for a small dense block (the end-members of one solution
// phase), by cyclic Jacobi rotations. Inside a miscibility gap the exact curvature of a
// solution phase is indefinite; Newton needs a positive-definite model, and flooring the
// eigenvalues at `ratio` times the block's own largest |eigenvalue| (but no higher than `absFloor`,
// if given) sets how far it steps
// along the unmixing direction, in the phase's own curvature units.
// Returns false and leaves `a` untouched if the rotations do not settle. If `indefinite` is
// given it is set to whether the block had a negative eigenvalue.
bool SymEigFloorInPlace( std::vector<double>& a, int n, double ratio, bool* indefinite = nullptr,
                         double absFloor = 0. )
{
    if( n < 1 ) return false;
    std::vector<double> v( (size_t)n*n, 0. );
    for( int i = 0; i < n; i++ ) v[(size_t)i*n+i] = 1.;
    const int maxSweeps = 60;
    bool converged = false;
    for( int sweep = 0; sweep < maxSweeps && !converged; sweep++ )
    {
        double off = 0.;
        for( int p = 0; p < n; p++ )
            for( int q = p+1; q < n; q++ )
                off += a[(size_t)p*n+q]*a[(size_t)p*n+q];
        double nrm = 0.;
        for( int i = 0; i < n*n; i++ ) nrm += a[(size_t)i]*a[(size_t)i];
        if( off <= 1e-24 * std::max( nrm, 1e-300 ) ) { converged = true; break; }
        for( int p = 0; p < n; p++ )
            for( int q = p+1; q < n; q++ )
            {
                const double apq = a[(size_t)p*n+q];
                if( std::fabs(apq) < 1e-300 ) continue;
                const double app = a[(size_t)p*n+p], aqq = a[(size_t)q*n+q];
                const double theta = ( aqq - app ) / ( 2.*apq );
                const double t = ( theta >= 0. ? 1. : -1. ) /
                                 ( std::fabs(theta) + std::sqrt( theta*theta + 1. ) );
                const double c = 1./std::sqrt( t*t + 1. ), sn = t*c;
                for( int k = 0; k < n; k++ )
                {
                    const double akp = a[(size_t)k*n+p], akq = a[(size_t)k*n+q];
                    a[(size_t)k*n+p] = c*akp - sn*akq;
                    a[(size_t)k*n+q] = sn*akp + c*akq;
                }
                for( int k = 0; k < n; k++ )
                {
                    const double apk = a[(size_t)p*n+k], aqk = a[(size_t)q*n+k];
                    a[(size_t)p*n+k] = c*apk - sn*aqk;
                    a[(size_t)q*n+k] = sn*apk + c*aqk;
                    const double vkp = v[(size_t)k*n+p], vkq = v[(size_t)k*n+q];
                    v[(size_t)k*n+p] = c*vkp - sn*vkq;
                    v[(size_t)k*n+q] = sn*vkp + c*vkq;
                }
            }
    }
    if( !converged ) return false;
    double lmax = 0.;
    for( int i = 0; i < n; i++ ) lmax = std::max( lmax, std::fabs( a[(size_t)i*n+i] ) );
    if( !( lmax > 0. ) ) return false;
    if( indefinite )
    {
        double lmin = a[0];
        for( int i = 1; i < n; i++ ) lmin = std::min( lmin, a[(size_t)i*n+i] );
        *indefinite = ( lmin < 0. );
    }
    const double headFloor = ratio * lmax;
    const double floorVal = ( absFloor > 0. ) ? std::min( headFloor, absFloor ) : headFloor;
    std::vector<double> lam( (size_t)n );
    for( int i = 0; i < n; i++ ) lam[(size_t)i] = std::max( a[(size_t)i*n+i], floorVal );
    for( int i = 0; i < n; i++ )
        for( int j = 0; j < n; j++ )
        {
            double sum = 0.;
            for( int k = 0; k < n; k++ ) sum += v[(size_t)i*n+k] * lam[(size_t)k] * v[(size_t)j*n+k];
            a[(size_t)i*n+j] = sum;
        }
    return true;
}

// Molality<->mole-fraction scale correction for the pH formula
// (ipm_chemical2.cpp's lnFmol = log(1000/MMC)), hardcoded for the common
// single/dominant-water-solvent case (MMC = water's molar mass)
const double kLnFmol = std::log(1000.0 / 18.01528);

// The system size Optima's default boxes are scaled by (species upper bounds with no DUL,
// the solvent-reseed ceiling, the control-condition titrant bound). With rescaling
// (pa_DG > 1e-5) that is pm.SMols = pa_DG, the internal total IC moles; without it pm.SMols
// is 0, so the actual total IC moles (SystemTotalMolesIC()'s sum) are used instead.
static double optima_default_box_moles( const MULTI& pm, double DG )
{
    if( DG > 1e-5 )
        return pm.SMols;
    double tot = 0.;
    for( long int i = 0; i < pm.N - pm.E; i++ )
        tot += pm.B[i];
    return tot;
}

} // namespace

// The two-phase simplex and the column-generation seed live in ipm_level.cpp (built without Optima too).
bool TwoPhaseSimplexMinSum( long int nRows, long int nCols,
                             const std::function<double(long int,long int)>& a,
                             const double* b, std::vector<double>& xOut,
                             const double* cost = nullptr, double* yOut = nullptr );

void TMultiBase::SetControlCondition_pH( double pH_target, double tolerance )
{
    if( node1 == nullptr )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_pH() requires an initialized TNode (call after GEM_init())" );

    long int xH  = node1->IC_name_to_xCH( "H" );
    long int xZz = node1->IC_name_to_xCH( "Zz" );
    long int xHplusDC = node1->DC_name_to_xCH( "H+" );
    if( xH < 0 || xZz < 0 || xHplusDC < 0 )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_pH(): could not resolve H/Zz IC or H+ DC index in this project" );

    optima_control_conditions.erase(
        std::remove_if( optima_control_conditions.begin(), optima_control_conditions.end(),
                         [](const EqControlCondition& c){ return c.name == "pH"; } ),
        optima_control_conditions.end() );

    EqControlCondition c;
    c.name = "pH";
    c.target = pH_target;
    c.tolerance = tolerance;
    c.stoich = { { xH, 1.0 }, { xZz, 1.0 } };
    // Default tolerance (used only when the caller passes tolerance<0):
    // pa_p->GAS is a chem.pot.-difference threshold (mol/mol); pH's own
    // formula below is linear in Muj with slope -ln_to_lg, so a Muj-space
    // tolerance of GAS maps to a pH-space tolerance of ln_to_lg*GAS.
    c.defaultToleranceFn = [this]() { return std::fabs( ln_to_lg ) * base_param()->GAS; };
    // pH = -ln_to_lg*(Muj - G0[H+] + lnFmol)  =>  Muj_target = G0[H+] - lnFmol - pH_target/ln_to_lg
    // (ipm_chemical2.cpp's own pH formula; Muj = sum_i U[i]*a(H+,i), which is
    // exactly sum_i U[i]*stoich[i] above, since stoich matches H+'s own
    // stoichiometry)
    c.fixedGradientFn = [this, xHplusDC]( double target )
    { return pm.G0[xHplusDC] - kLnFmol - target / ln_to_lg; };
    c.achievedValueFn = [this, xHplusDC]( double achievedMuj )
    { return -ln_to_lg * ( achievedMuj - pm.G0[xHplusDC] + kLnFmol ); };
    c.needsAqueousPhase = true;

    optima_control_conditions.push_back( c );
}

void TMultiBase::SetControlCondition_Eh( double Eh_target, double tolerance )
{
    if( node1 == nullptr )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_Eh() requires an initialized TNode (call after GEM_init())" );

    long int xZz = node1->IC_name_to_xCH( "Zz" );
    if( xZz < 0 )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_Eh(): could not resolve Zz IC index in this project" );

    optima_control_conditions.erase(
        std::remove_if( optima_control_conditions.begin(), optima_control_conditions.end(),
                         [](const EqControlCondition& c){ return c.name == "Eh"; } ),
        optima_control_conditions.end() );

    EqControlCondition c;
    c.name = "Eh";
    c.target = Eh_target;
    c.tolerance = tolerance;
    // Direct electron/charge (Zz-row) titrant: only a titrant acting on the charge row has a
    // closed-form relation to the Zz dual (and so to Eh), which makes the one-shot joint
    // solve possible.
    c.stoich = { { xZz, -1.0 } };
    // Default tolerance (used only when the caller passes tolerance<0):
    // same pa_p->GAS reuse as SetControlCondition_pH(), mapped through Eh's
    // own slope (0.000086*T) instead of pH's. Lazy (reads pm.T at solve
    // time, not here) for the same reason fixedGradientFn below is lazy.
    c.defaultToleranceFn = [this]() { return 0.000086 * pm.T * base_param()->GAS; };
    // Eh = 0.000086*U[Zz]*T  =>  U[Zz]_target = Eh_target/(0.000086*T);
    // achievedMuj = sum_i U[i]*stoich[i] = -U[Zz], so the fixed gradient
    // (pinned to achievedMuj, via the KKT stationarity of this slot) is
    // -U[Zz]_target.
    c.fixedGradientFn = [this]( double target )
    { return -( target / ( 0.000086 * pm.T ) ); };
    c.achievedValueFn = [this]( double achievedMuj )
    { return 0.000086 * ( -achievedMuj ) * pm.T; };
    c.needsAqueousPhase = true;

    optima_control_conditions.push_back( c );
}

void TMultiBase::ClearControlConditions()
{
    optima_control_conditions.clear();
}

double TMultiBase::OptimaMaxMassBalanceResidual()
{
    double maxres = 0.;
    for( long int i = 0; i < pm.N; i++ )
    {
        double sum = 0.;
        for( long int j = 0; j < pm.L; j++ )
            sum += pm.A[ i + j*pm.N ] * pm.Y[j];
        maxres = std::max( maxres, std::fabs( sum - pm.B[i] ) );
    }
    return maxres;
}

/// Does this system have an aqueous phase? Tested with the phase classifier (PH_AQUEL), not
/// with pm.LO: pm.LO is initialised to 0 and only reassigned when an aqueous phase exists,
/// so on an aqueous-free system it still points at species 0 and a range test on it cannot
/// tell "no solvent" from "the solvent is DC 0".
bool TMultiBase::HasAqueousPhase() const
{
    for( long int k = 0; k < pm.FIs; k++ )
        if( pm.PHC[k] == PH_AQUEL )
            return true;
    return false;
}

/// Lists the valences each element has across the allowed species (O -2, H +1; species of H and O only are
/// left out, since H2 and O2 are what carry the redox level). The bulk is redox-free if some choice of one
/// listed valence per element leaves it with no electrons to spare: then nothing but traces sets the redox level.
TMultiBase::RedoxFreeCheck TMultiBase::CheckRedoxFree() const
{
    RedoxFreeCheck r;
    const long int N = pm.N, L = pm.L;
    long int iO = -1, iH = -1, iZ = -1;
    for( long int i = 0; i < N; i++ )
    {
        if( pm.ICC[i] == IC_OXYGEN && iO < 0 ) iO = i;
        else if( pm.ICC[i] == IC_HYDROGEN && iH < 0 ) iH = i;
        else if( pm.ICC[i] == IC_CHARGE && iZ < 0 ) iZ = i;
    }
    if( iO < 0 || iH < 0 || iZ < 0 || !HasAqueousPhase() )
        return r;
    r.checked = true;
    auto a = [&]( long int j, long int i ) { return pm.A[ i + j*N ]; };
    auto isOther = [&]( long int i ) { return i != iO && i != iH && i != iZ; };
    // the valences of a species add up to its charge: sum_e v_e a_e - 2 nO + nH = z
    auto rest = [&]( long int j ) { return a( j, iZ ) + 2. * a( j, iO ) - a( j, iH ); };
    auto addVal = []( std::vector<double>& vs, double v ) {
        for( double w : vs ) if( std::fabs( w - v ) < 1e-6 ) return;
        vs.push_back( v );
    };

    std::vector<std::vector<double>> val( (size_t)N );
    std::vector<std::vector<long int>> multi;   // species with two or more other elements
    for( long int j = 0; j < L; j++ )
    {
        if( pm.DUL[j] <= 0. ) continue;   // switched off
        std::vector<long int> el;
        for( long int i = 0; i < N; i++ )
            if( isOther( i ) && a( j, i ) != 0. ) el.push_back( i );
        if( el.size() == 1 )
            addVal( val[(size_t)el[0]], rest( j ) / a( j, el[0] ) );
        else if( el.size() > 1 )
            multi.push_back( { j } );
    }
    // An element seen only together with others takes a valence from each choice of valences of the others there
    // (capped at 64 choices); a spare valence only widens the search below.
    for( bool changed = true; changed; )
    {
        changed = false;
        for( const auto& m : multi )
        {
            const long int j = m[0];
            long int unknown = -1, nUnknown = 0;
            std::vector<long int> known;
            size_t combos = 1;
            for( long int i = 0; i < N && nUnknown < 2; i++ )
            {
                if( !isOther( i ) || a( j, i ) == 0. ) continue;
                if( val[(size_t)i].empty() ) { unknown = i; nUnknown++; }
                else { known.push_back( i ); combos *= val[(size_t)i].size(); }
            }
            if( nUnknown != 1 || combos > 64 ) continue;
            const size_t before = val[(size_t)unknown].size();
            for( size_t c = 0; c < combos; c++ )
            {
                double s = rest( j );
                size_t q = c;
                for( long int i : known )
                {
                    const auto& vs = val[(size_t)i];
                    s -= vs[q % vs.size()] * a( j, i );
                    q /= vs.size();
                }
                addVal( val[(size_t)unknown], s / a( j, unknown ) );
            }
            if( val[(size_t)unknown].size() != before ) changed = true;
        }
    }

    // Search one valence per element with zero excess, largest amounts first, pruned by the reachable range.
    std::vector<long int> els;
    double scale = 2. * std::fabs( pm.B[iO] ) + std::fabs( pm.B[iH] ) + std::fabs( pm.B[iZ] );
    r.consistent = true;
    r.determined = true;
    for( long int i = 0; i < N; i++ )
    {
        if( !isOther( i ) || pm.B[i] == 0. ) continue;
        if( val[(size_t)i].empty() ) { r.determined = false; return r; }   // valence unknown
        if( val[(size_t)i].size() > 1 ) r.consistent = false;
        double vmax = 0.;
        for( double v : val[(size_t)i] ) vmax = std::max( vmax, std::fabs( v ) );
        scale += vmax * std::fabs( pm.B[i] );
        els.push_back( i );
    }
    std::sort( els.begin(), els.end(), [&]( long int x, long int y ) { return std::fabs( pm.B[x] ) > std::fabs( pm.B[y] ); } );
    const size_t ne = els.size();
    std::vector<double> loRest( ne + 1, 0. ), hiRest( ne + 1, 0. );
    for( size_t k = ne; k-- > 0; )
    {
        double lo = 1e300, hi = -1e300;
        for( double v : val[(size_t)els[k]] ) { lo = std::min( lo, v * pm.B[els[k]] ); hi = std::max( hi, v * pm.B[els[k]] ); }
        loRest[k] = loRest[k+1] + lo;
        hiRest[k] = hiRest[k+1] + hi;
    }
    const double tol = 1e-12 * std::max( scale, 1e-300 );   // rounding of the input amounts, not a real trace oxidant
    double best = 1e300, bestE = 0.;
    long int visits = 0;
    const long int kMaxVisits = 200000;
    std::function<bool( size_t, double )> search = [&]( size_t k, double e ) -> bool
    {
        if( ++visits > kMaxVisits ) return false;
        if( k == ne )
        {
            if( std::fabs( e ) < best ) { best = std::fabs( e ); bestE = e; }
            return std::fabs( e ) <= tol;
        }
        if( e + loRest[k] > tol || e + hiRest[k] < -tol )
            return false;
        for( double v : val[(size_t)els[k]] )
            if( search( k + 1, e + v * pm.B[els[k]] ) ) return true;
        return false;
    };
    r.free = search( 0, -2. * pm.B[iO] + pm.B[iH] - pm.B[iZ] );
    if( visits > kMaxVisits ) r.determined = false;
    r.excess = best < 1e300 ? bestE : 0.;
    r.rel = std::fabs( r.excess ) / std::max( scale, 1e-300 );
    return r;
}

/// Flags each phase that has a twin: another phase with the same end-members (same count and
/// same stoichiometry rows). Only such phases can show a miscibility gap; they are listed in the
/// trace and warned about if the floor never changed them.
std::vector<char> TMultiBase::PhaseHessianTwins() const
{
    std::vector<char> twin( (size_t)pm.FIs, 0 );
    const long int N = pm.N;
    std::vector<std::vector<std::vector<double>>> rows( (size_t)pm.FIs );
    long int p0 = 0;
    for( long int k = 0; k < pm.FIs; k++ )
    {
        const long int p1 = p0 + pm.L1[k];
        if( pm.PHC[k] != PH_AQUEL )
            for( long int j = p0; j < p1 && p1 <= pm.L; j++ )
            {
                std::vector<double> r( (size_t)N );
                for( long int i = 0; i < N; i++ ) r[(size_t)i] = pm.A[ i + j*N ];
                rows[(size_t)k].push_back( std::move(r) );
            }
        std::sort( rows[(size_t)k].begin(), rows[(size_t)k].end() );
        p0 = p1;
    }
    for( long int a = 0; a < pm.FIs; a++ )
        for( long int b = a+1; b < pm.FIs; b++ )
            if( rows[(size_t)a].size() > 1 && rows[(size_t)a] == rows[(size_t)b] )
                twin[(size_t)a] = twin[(size_t)b] = 1;
    return twin;
}

/// Starts the per-run record of the curvature floor: which phases are twins, none fired yet.
void TMultiBase::PhaseHessianBegin()
{
    phfTwin = ( base_param()->PhaseHessianFloor > 0. ) ? PhaseHessianTwins()
                                                      : std::vector<char>( (size_t)pm.FIs, 0 );
    phfFired.assign( (size_t)pm.FIs, 0 );
}

/// Ends the record: one trace line and, if twins never split, one warning.
void TMultiBase::PhaseHessianEnd()
{
    if( phfTwin.size() == (size_t)pm.FIs )
    {
        std::string names, fired;
        long int total = 0;
        std::vector<std::string> list;
        for( long int k = 0; k < pm.FIs; k++ )
            if( phfTwin[(size_t)k] )
            {
                const std::string nm = char_array_to_string( pm.SF[k], MAXPHNAME );
                names += ( names.empty() ? "" : "|" ) + nm;
                fired += ( fired.empty() ? "" : "|" ) + std::to_string( phfFired[(size_t)k] );
                total += phfFired[(size_t)k];
                list.push_back( nm );
            }
        if( !names.empty() )
        {
            native_trace_decide( "phfloor twins=%s fired=%s", names.c_str(), fired.c_str() );
            if( total == 0 )
            {
                std::string shown;
                for( size_t i = 0; i < list.size(); i++ )
                    shown += ( i ? ( i+1 == list.size() ? " and " : ", " ) : "" ) + list[i];
                ipm_logger->warn( "Solid solutions {} are copies and did not split. "
                                  "Try: if no miscibility gap is expected, remove one copy to speed up the calculation.",
                                  shown );
            }
        }
    }
    phfTwin.clear();
    phfFired.clear();
}

/// Warns when two copies of a solid solution end with the same composition: the split between them is
/// arbitrary, which can make the solve slow or unstable. One warning and one trace line per such pair.
void TMultiBase::WarnCollapsedTwins()
{
    const std::vector<char> twin = PhaseHessianTwins();
    if( twin.size() == (size_t)pm.FIs && pm.X != nullptr )
    {
        const long int N = pm.N;
        std::vector<long int> start( (size_t)pm.FIs + 1, 0 );
        for( long int k = 0; k < pm.FIs; k++ ) start[(size_t)k+1] = start[(size_t)k] + pm.L1[k];
        // A copy counts as present above a tiny fraction of all moles in the system.
        double sysTot = 0.;
        for( long int j = 0; j < pm.L; j++ ) sysTot += std::max( pm.X[j], 0. );
        const double presentMin = 1e-9 * sysTot;
        auto phName = [this]( long int k ) {
            std::string n = char_array_to_string( pm.SF[k] + MAXSYMB, MAXPHNAME );
            const auto e = n.find_last_not_of( ' ' );
            return e == std::string::npos ? std::string() : n.substr( 0, e + 1 );
        };
        auto sameColumn = [&]( long int i, long int j ) {
            for( long int r = 0; r < N; r++ )
                if( pm.A[ r + i*N ] != pm.A[ r + j*N ] ) return false;
            return true;
        };
        for( long int a = 0; a < pm.FIs; a++ )
            for( long int b = a+1; b < pm.FIs; b++ )
            {
                if( !twin[(size_t)a] || !twin[(size_t)b] || pm.L1[a] != pm.L1[b] || start[(size_t)b+1] > pm.L )
                    continue;
                double ta = 0., tb = 0.;
                for( long int j = start[(size_t)a]; j < start[(size_t)a+1]; j++ ) ta += std::max( pm.X[j], 0. );
                for( long int j = start[(size_t)b]; j < start[(size_t)b+1]; j++ ) tb += std::max( pm.X[j], 0. );
                if( !( ta > presentMin && tb > presentMin && ta > 1e-6 * ( ta + tb ) && tb > 1e-6 * ( ta + tb ) ) )
                    continue;
                double dmax = 0.;
                bool matched = true;
                for( long int i = start[(size_t)a]; i < start[(size_t)a+1] && matched; i++ )
                {
                    long int jm = -1;
                    for( long int j = start[(size_t)b]; j < start[(size_t)b+1]; j++ )
                        if( sameColumn( i, j ) ) { jm = j; break; }
                    if( jm < 0 ) { matched = false; break; }
                    dmax = std::max( dmax, std::fabs( std::max( pm.X[i], 0. ) / ta - std::max( pm.X[jm], 0. ) / tb ) );
                }
                if( !matched || dmax > 1e-3 )
                    continue;
                const std::string na = phName( a );
                const std::string nb = phName( b );
                native_trace_decide( "twincollapse pair=%s|%s amounts=%.6e|%.6e dfrac=%.3e",
                                     na.c_str(), nb.c_str(), ta, tb, dmax );
                // The log line is shown once per pair, so long runs with many calls are not flooded.
                static thread_local std::set<std::string> warnedPairs;
                if( !warnedPairs.insert( na + "|" + nb ).second )
                    continue;
                // All copy pairs of the system: removing only the collapsed one is not enough.
                std::string allPairs;
                long int nPairs = 0;
                for( long int c = 0; c < pm.FIs; c++ )
                    for( long int d = c+1; d < pm.FIs; d++ )
                        if( twin[(size_t)c] && twin[(size_t)d] && pm.L1[c] == pm.L1[d] )
                        {
                            bool same = start[(size_t)d+1] <= pm.L;
                            for( long int i = start[(size_t)c]; i < start[(size_t)c+1] && same; i++ )
                            {
                                bool found = false;
                                for( long int j = start[(size_t)d]; j < start[(size_t)d+1] && !found; j++ )
                                    found = sameColumn( i, j );
                                same = found;
                            }
                            if( same )
                            {
                                allPairs += ( allPairs.empty() ? "" : ", " ) + phName( c ) + "/" + phName( d );
                                nPairs++;
                            }
                        }
                ipm_logger->warn( "Solid solutions {} and {} are copies and ended with the same composition, so the split "
                                  "between them is arbitrary; this can make the calculation slow or unstable. "
                                  "Try: if no miscibility gap is expected here, remove one copy of each copied pair "
                                  "({} in this system: {}).", na, nb, nPairs, allPairs );
            }
    }
}

bool TMultiBase::DetectSolventCollapseAndReseed( const double* x, double upperBound, double& waterSeedOut )
{
    if( !HasAqueousPhase() || pm.LO < 0 || pm.LO >= pm.L )
        return false;
    long int j0 = 0;
    for( long int k = 0; k < pm.FIs; k++ )
    {
        const long int j1 = j0 + pm.L1[k];
        if( pm.LO >= j0 && pm.LO < j1 )
        {
            double otherTotal = 0.;
            for( long int j = j0; j < j1; j++ )
                if( j != pm.LO ) otherTotal += x[j];
            // An empty aqueous phase also needs water, not only one where the solutes outweigh it
            // (pa_OptimaEmptySolvent); otherwise the test hangs on round-off when both are zero.
            const bool waterAbsent = base_param()->OptimaEmptySolvent && x[pm.LO] <= std::max( pm.DSM, pm.XwMinM );
            if( x[pm.LO] >= otherTotal && !waterAbsent )
                return false; // solvent already dominates - no correction needed
            if( waterAbsent && x[pm.LO] >= otherTotal )
                native_trace_decide( "emptysolvent water=%.6g solutes=%.6g", x[pm.LO], otherTotal );
            // Reseed value: bounded by the most water the system's total H and O could form.
            double waterSeed = otherTotal;
            const long int xH = ( node1 != nullptr ) ? node1->IC_name_to_xCH( "H" ) : -1;
            const long int xO = ( node1 != nullptr ) ? node1->IC_name_to_xCH( "O" ) : -1;
            if( xH >= 0 && xO >= 0 )
            {
                const double bulkEstimate = std::min( pm.B[xH] / 2.0, pm.B[xO] );
                if( bulkEstimate > waterSeed )
                    waterSeed = bulkEstimate;
            }
            if( upperBound > 0. )
                waterSeed = std::min( waterSeed, upperBound );
            waterSeedOut = waterSeed;
            return true;
        }
        j0 = j1;
    }
    return false;
}

bool TMultiBase::DetectPhaseCollapseAndReseed( const double* x,
                                                std::vector<std::pair<long int,double>>& reseedsOut,
                                                long int excludePhaseIdx,
                                                bool stabilityGated )
{
    reseedsOut.clear();
    // Cold seed only: a phase is refilled only if the linearised-Gibbs LP says it should form,
    // and is filled in proportion to how favourable each end-member is.
    // Phases with a twin keep the plain equal split, so the twins stay symmetric.
    std::vector<double> redCost;
    std::vector<char> hasTwin;
    bool gated = false;
    if( stabilityGated && pm.G0 != nullptr )
    {
        std::vector<double> yd;
        if( LPGibbsDual( yd ) )
        {
            gated = true;
            hasTwin = PhaseHessianTwins();
            redCost.assign( (size_t)pm.L, 0. );
            for( long int j = 0; j < pm.L; j++ )
            {
                double z = 0.;
                for( long int i = 0; i < pm.N; i++ ) z += yd[(size_t)i] * pm.A[ i + j*pm.N ];
                redCost[(size_t)j] = pm.G0[j] - z;
            }
        }
    }
    long int j0 = 0;
    for( long int k = 0; k < pm.FIs; k++ ) // multicomponent phases only
    {
        const long int j1 = j0 + pm.L1[k];
        const long int nEnd = j1 - j0;
        const bool gatePhase = gated && !( hasTwin.size() == (size_t)pm.FIs && hasTwin[(size_t)k] );
        // A gas or fluid phase is not refilled: the bulk-limited amount would be far from mass-balanced.
        if( !gatePhase && ( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_FLUID ) ) { j0 = j1; continue; }
        if( k != excludePhaseIdx && nEnd > 1 ) // single-species phases have their own PhMinM guard
        {
            double phaseTotal = 0.;
            for( long int j = j0; j < j1; j++ )
                phaseTotal += x[j];
            // The aqueous phase has an extra, stricter threshold (Y[LO] <= pm.XwMinM) on top of
            // the generic pm.DSM one; both guard the same "phase reported absent" cliff in
            // PrimalChemicalPotentials().
            const double threshold = ( pm.PHC[k] == PH_AQUEL ) ? std::max( pm.DSM, pm.XwMinM ) : pm.DSM;
            // Safety margin: require the total to be comfortably clear of the cliff, not
            // merely non-zero.
            if( phaseTotal < threshold * 10. )
            {
                // Bound this phase's plausible total by the bulk composition available for the
                // ICs its end-members contain - the water formula min(bulk-H/2, bulk-O)
                // generalised to any stoichiometry - per end-member, by that end-member's own
                // limiting IC. (A single phase-wide minimum would let one trace IC in a minor
                // end-member cap the whole phase back onto the phase-absence cliff.)
                bool any = false;
                std::vector<double> emBound( (size_t)nEnd, -1. );
                for( long int j = j0; j < j1; j++ )
                {
                    double b = -1.;
                    for( long int i = 0; i < pm.N; i++ )
                    {
                        const double coef = pm.A[ i + j*pm.N ];
                        if( coef > 0. )
                        {
                            const double icBound = pm.B[i] / coef;
                            if( b < 0. || icBound < b ) b = icBound;
                        }
                    }
                    emBound[(size_t)(j-j0)] = b;
                    if( b > 0. ) any = true;
                }
                if( any )
                {
                    // Each end-member gets its own bound divided by the end-member count;
                    // Newton then corrects the actual mix.
                    // Stability gate and composition from the LP reduced costs.
                    std::vector<double> frac( (size_t)nEnd, 0. );
                    double gatedTotal = 0.;
                    if( gatePhase )
                    {
                        double smin = 1e300;
                        for( long int j = j0; j < j1; j++ )
                            if( emBound[(size_t)(j-j0)] > 0. ) smin = std::min( smin, redCost[(size_t)j] );
                        double sumW = 0.;
                        for( long int j = j0; j < j1; j++ )
                            if( emBound[(size_t)(j-j0)] > 0. )
                            {
                                frac[(size_t)(j-j0)] = std::exp( -( redCost[(size_t)j] - smin ) );
                                sumW += frac[(size_t)(j-j0)];
                            }
                        const double sigma = -smin + std::log( sumW );
                        for( long int j = j0; j < j1; j++ )
                        {
                            if( emBound[(size_t)(j-j0)] > 0. )
                                gatedTotal += std::max( emBound[(size_t)(j-j0)], threshold * 10. ) / (double)nEnd;
                            frac[(size_t)(j-j0)] /= sumW;
                        }
                        native_trace_decide( "reseedR phase=%s sigma=%.6g refilled=%d",
                                             char_array_to_string( pm.SF[k], MAXPHNAME ).c_str(), sigma, sigma > 0. ? 1 : 0 );
                        if( sigma <= 0. ) { j0 = j1; continue; }
                    }
                    for( long int j = j0; j < j1; j++ )
                    {
                        const double b = emBound[(size_t)(j-j0)];
                        if( b <= 0. ) continue;
                        double perSpecies = std::max( b, threshold * 10. ) / (double)nEnd;
                        if( gatePhase ) perSpecies = gatedTotal * frac[(size_t)(j-j0)];
                        if( perSpecies > x[j] )
                            reseedsOut.push_back( { j, perSpecies } );
                    }
                }
            }
        }
        j0 = j1;
    }
    return !reseedsOut.empty();
}


// pa_OptimaFinish: with the phase set fixed the problem is smooth, and an exchange Optima can
// cycle on (a pure compound against a solution end-member of the same composition, e.g.
// CaO(s) <-> CaO(l)) is one equation. Equality-constrained Newton on the species amounts of
// the present phases (pure phases linear, solution phases through a finite-difference Hessian
// of pm.F, scaled by the amount so trace and major end-members weigh alike) with the element
// multipliers as the dual, Levenberg-Marquardt damping, fraction-to-boundary and an Armijo
// line search on G. Set changes are rare and tabu'd: a pure phase whose amount reaches zero
// leaves; a pure phase (or an end-member of an already present solution phase) whose reduced
// gradient is negative at convergence enters. An absent solution phase is never entered here
// (that is the TPD search's job, pa_OptimaTpdAccept). On success pm.X, pm.Y and the
// multipliers of rows with a free carrier are overwritten, and the caller judges the state by
// the same KKT / mass-balance / stability checks as Optima's own. On failure nothing changes.
// In plain words: if Optima stops just short, finish the job with the phases it found.
bool TMultiBase::PotentialSpaceFinish( double dcFloor, const std::vector<double>& xlo, const std::vector<double>& xhi )
{
    gems3k_timing::Scope timing( gems3k_timing::Psf );
    const long int L = pm.L, N = pm.N;
    if( L <= 0 || N <= 0 || !pm.A || !pm.B || !pm.F || !pm.U )
        return false;
    const std::vector<double> Xsave( pm.X, pm.X + L ), Ysave( pm.Y, pm.Y + L ), Usave( pm.U, pm.U + N );
    auto restore = [&]() {
        for( long int j = 0; j < L; j++ ) { pm.X[j] = Xsave[(size_t)j]; pm.Y[j] = Ysave[(size_t)j]; }
        for( long int i = 0; i < N; i++ ) pm.U[i] = Usave[(size_t)i];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
    };
    std::vector<double> n( Xsave ), F( (size_t)L ), Fp( (size_t)L );
    double bScale = 0.;
    for( long int i = 0; i < N; i++ ) bScale = std::max( bScale, std::fabs( pm.B[i] ) );
    if( !( bScale > 0. ) ) return false;

    auto lowTol = [&]( long int j ) { return xlo[(size_t)j] + std::max( dcFloor, xlo[(size_t)j] * 1e-6 ); };
    auto degenerate = [&]( long int j ) { return xhi[(size_t)j] <= lowTol( j ); };
    // A start that does not hold the mass balance (a stalled Optima iterate can miss it by tens
    // of mol) leaves too few free species to repair it, so the finish then starts from a
    // leveled vertex: the column-generation vertex, else the plain feasibility vertex. The
    // start's G is then no reference for acceptance.
    bool seededStart = false;
    {
        double r0 = 0.;
        for( long int i = 0; i < N; i++ )
        {
            double sm = -pm.B[i];
            for( long int j = 0; j < L; j++ ) sm += pm.A[ i + j*N ] * n[(size_t)j];
            r0 = std::max( r0, std::fabs( sm ) );
        }
        if( r0 > 1e-6 * bScale )
        {
            std::vector<double> nSeed;
            if( ColumnGenerationSeed( nSeed, 1e-6 ) || LPFeasibilitySeed( nSeed ) )
            {
                for( long int j = 0; j < L; j++ )
                    n[(size_t)j] = std::min( std::max( nSeed[(size_t)j], xlo[(size_t)j] ), xhi[(size_t)j] );
                seededStart = true;
            }
        }
    }
    std::vector<char> freeSp( (size_t)L, 0 ), tabu( (size_t)L, 0 );
    for( long int j = 0; j < L; j++ )
        freeSp[(size_t)j] = ( n[(size_t)j] > std::max( lowTol( j ), 1e-11 * bScale ) && !degenerate( j ) && n[(size_t)j] < xhi[(size_t)j] * ( 1. - 1e-6 ) ) ? 1 : 0;

    // phase start index per species and phase kind
    std::vector<long int> phOf( (size_t)L, 0 ), phStart( (size_t)pm.FI + 1, 0 );
    {
        long int j0 = 0;
        for( long int k = 0; k < pm.FI; k++ )
        {
            phStart[(size_t)k] = j0;
            for( long int j = j0; j < j0 + pm.L1[k]; j++ ) phOf[(size_t)j] = k;
            j0 += pm.L1[k];
        }
        phStart[(size_t)pm.FI] = j0;
    }
    auto evalF = [&]( const std::vector<double>& v, std::vector<double>& f ) -> double {
        for( long int j = 0; j < L; j++ ) pm.X[j] = v[(size_t)j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );
        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
        double G = 0.;
        for( long int j = 0; j < L; j++ ) { f[(size_t)j] = pm.F[j]; G += v[(size_t)j] * pm.F[j]; }
        return G;
    };
    auto massResid = [&]( const std::vector<double>& v, std::vector<double>& r ) -> double {
        double m = 0.;
        for( long int i = 0; i < N; i++ )
        {
            double s = -pm.B[i];
            for( long int j = 0; j < L; j++ ) s += pm.A[ i + j*N ] * v[(size_t)j];
            r[(size_t)i] = s;
            m = std::max( m, std::fabs( s ) );
        }
        return m;
    };

    std::vector<double> rTmp( (size_t)N );
    const double rStart = massResid( Xsave, rTmp );          // the start's own mass-balance residual
    const std::vector<double> rStartRow( rTmp );              // per row: a trace IC's residual is tiny in absolute
                                                              // terms and large against its own b_i
    std::vector<char> presentAtStart( (size_t)pm.FI, 0 );     // phases Optima's answer had present
    for( long int k = 0; k < pm.FI; k++ )
    {
        double t = 0.;
        for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ ) t += Xsave[(size_t)j];
        presentAtStart[(size_t)k] = ( t > pm.DSM ) ? 1 : 0;
    }
    const bool keepPhases = g_finishFromSuccess && !seededStart;
    double G = evalF( n, F );
    const double G0 = seededStart ? 1e300 : G;
    std::vector<double> r( (size_t)N ), U( (size_t)N, 0. ), Un( Usave );
    long int changes = 0, it = 0, entered = 0, left = 0; bool vanished = false;
    bool converged = false;
    const long int kMaxIt = 60, kMaxChanges = 6, kMaxSize = 900;
    std::string why = "maxit";
    // pa_OptimaFinishStall: the lowest reduced gradient so far, the step it was reached at, and G at that step.
    const long int kStall = base_param()->OptimaFinishStall;
    double bestRed = 1e300, stallG = 0., stepMin = 0., stepMax = 0.; long int bestRedIt = 0, bestRedChanges = 0, stallSkips = 0;
    for( ; it < kMaxIt; it++ )
    {
        std::vector<long int> fr;
        for( long int j = 0; j < L; j++ ) if( freeSp[(size_t)j] ) fr.push_back( j );
        const long int nf = (long int)fr.size();
        if( nf == 0 ) { why = "nofree"; break; }
        // rows with at least one free carrier
        std::vector<long int> rows;
        std::vector<double> rho( (size_t)N, 1. );
        for( long int i = 0; i < N; i++ )
        {
            double s = 0.;
            for( long int a = 0; a < nf; a++ ) { const double v = pm.A[ i + fr[(size_t)a]*N ] * n[(size_t)fr[(size_t)a]]; s += v*v; }
            if( s > 0. ) { rows.push_back( i ); rho[(size_t)i] = 1. / std::sqrt( s ); }
        }
        const long int nr = (long int)rows.size();
        const long int M = nf + nr;
        if( M > kMaxSize ) { why = "size"; break; }
        massResid( n, r );

        // Hessian of the free species, finite differences of pm.F inside each phase that has >= 2 species
        std::vector<double> Hm( (size_t)nf * (size_t)nf, 0. );
        {
            std::vector<long int> pos( (size_t)L, -1 );
            for( long int a = 0; a < nf; a++ ) pos[(size_t)fr[(size_t)a]] = a;
            for( long int a = 0; a < nf; a++ )
            {
                const long int j = fr[(size_t)a];
                if( j >= pm.Ls ) continue;                 // pure phase: mu independent of its amount
                const long int k = phOf[(size_t)j];
                if( pm.L1[k] < 2 ) continue;
                std::vector<double> np( n );
                const double h = 1e-6 * n[(size_t)j];
                np[(size_t)j] += h;
                evalF( np, Fp );
                for( long int jb = phStart[(size_t)k]; jb < phStart[(size_t)k+1]; jb++ )
                {
                    const long int b = pos[(size_t)jb];
                    if( b < 0 ) continue;
                    Hm[(size_t)b*(size_t)nf + (size_t)a] = ( Fp[(size_t)jb] - F[(size_t)jb] ) / h;
                }
            }
            evalF( n, F );   // leave pm.* at the base point
            for( long int a = 0; a < nf; a++ )
                for( long int b = a + 1; b < nf; b++ )
                {
                    const double s = 0.5 * ( Hm[(size_t)a*(size_t)nf + (size_t)b] + Hm[(size_t)b*(size_t)nf + (size_t)a] );
                    Hm[(size_t)a*(size_t)nf + (size_t)b] = Hm[(size_t)b*(size_t)nf + (size_t)a] = s;
                }
        }
        // scaled KKT system: unknowns d (Dn = s.d) and w (U = rho.w)
        std::vector<double> s( (size_t)nf );
        for( long int a = 0; a < nf; a++ ) s[(size_t)a] = n[(size_t)fr[(size_t)a]];
        double mu = 1e-10;
        std::vector<double> gt( (size_t)nf ), dn( (size_t)nf ), Uw( (size_t)nr ), xsol, Uw0;
        for( long int a = 0; a < nf; a++ ) gt[(size_t)a] = s[(size_t)a] * F[(size_t)fr[(size_t)a]];
        double gDotDn = 0., dMaxRel = 0., lam = 0., dMerit = 0.; bool nullStep = false;
        const double rFloor = 1e-6 * bScale;   // rows below this are not penalised: the multipliers of an unrepairable row make lam*|r| noise
        auto rex = []( double v, double f ) { return std::max( std::fabs( v ) - f, 0. ); };
        double r1 = 0.; for( long int i = 0; i < N; i++ ) r1 += rex( r[(size_t)i], rFloor );
        bool solved = false;
        double dual = 1e-16;
        // LU of the KKT matrix with partial pivoting; one factorisation serves the solve and its refinement steps.
        std::vector<double> KF;
        std::vector<long int> piv;
        auto factorDense = [&]( const std::vector<double>& K ) -> bool
        {
            KF = K;
            piv.assign( (size_t)M, 0 );
            for( long int c = 0; c < M; c++ )
            {
                long int p = c; double pv = std::fabs( KF[(size_t)c*(size_t)M + (size_t)c] );
                for( long int q = c + 1; q < M; q++ )
                {
                    const double v = std::fabs( KF[(size_t)q*(size_t)M + (size_t)c] );
                    if( v > pv ) { pv = v; p = q; }
                }
                if( !( pv > 1e-300 ) ) return false;
                piv[(size_t)c] = p;
                if( p != c )
                    for( long int q = 0; q < M; q++ ) std::swap( KF[(size_t)c*(size_t)M + (size_t)q], KF[(size_t)p*(size_t)M + (size_t)q] );
                for( long int q = c + 1; q < M; q++ )
                {
                    const double f = KF[(size_t)q*(size_t)M + (size_t)c] / KF[(size_t)c*(size_t)M + (size_t)c];
                    KF[(size_t)q*(size_t)M + (size_t)c] = f;
                    if( f == 0. ) continue;
                    for( long int t = c + 1; t < M; t++ ) KF[(size_t)q*(size_t)M + (size_t)t] -= f * KF[(size_t)c*(size_t)M + (size_t)t];
                }
            }
            return true;
        };
        auto solveFactored = [&]( std::vector<double> x )
        {
            for( long int c = 0; c < M; c++ )
                if( piv[(size_t)c] != c ) std::swap( x[(size_t)c], x[(size_t)piv[(size_t)c]] );
            for( long int c = 0; c < M; c++ )
                for( long int q = c + 1; q < M; q++ )
                {
                    const double f = KF[(size_t)q*(size_t)M + (size_t)c];
                    if( f != 0. ) x[(size_t)q] -= f * x[(size_t)c];
                }
            for( long int c = M - 1; c >= 0; c-- )
            {
                double v = x[(size_t)c];
                for( long int t = c + 1; t < M; t++ ) v -= KF[(size_t)c*(size_t)M + (size_t)t] * x[(size_t)t];
                x[(size_t)c] = v / KF[(size_t)c*(size_t)M + (size_t)c];
            }
            xsol = x;
        };
        const double mu0 = mu;
        bool ignoreFloorRows = false;
        for( int pass = 0; pass < 2 && !solved; pass++ )
        {
        if( pass == 1 ) { ignoreFloorRows = true; mu = mu0; }   // chasing the rows first; a repair that only raises G is then left alone
        for( int tryMu = 0; tryMu < 10 && !solved; tryMu++ )
        {
            if( tryMu > 0 ) mu *= 100.;
            std::vector<double> K( (size_t)M * (size_t)M, 0. ), rhs( (size_t)M, 0. );
            double hmax = 0.;
            for( long int a = 0; a < nf; a++ )
                for( long int b = 0; b < nf; b++ )
                {
                    const double v = s[(size_t)a] * s[(size_t)b] * Hm[(size_t)a*(size_t)nf + (size_t)b];
                    K[(size_t)a*(size_t)M + (size_t)b] = v;
                    if( a == b ) hmax = std::max( hmax, std::fabs( v ) );
                }
            for( long int a = 0; a < nf; a++ ) K[(size_t)a*(size_t)M + (size_t)a] += mu * std::max( hmax, 1. );
            for( long int q = 0; q < nr; q++ )
            {
                const long int i = rows[(size_t)q];
                for( long int a = 0; a < nf; a++ )
                {
                    const double v = rho[(size_t)i] * pm.A[ i + fr[(size_t)a]*N ] * s[(size_t)a];
                    K[(size_t)( nf + q )*(size_t)M + (size_t)a] = v;
                    K[(size_t)a*(size_t)M + (size_t)( nf + q )] = -v;
                }
                K[(size_t)( nf + q )*(size_t)M + (size_t)( nf + q )] = -dual;
                rhs[(size_t)( nf + q )] = ( ignoreFloorRows && std::fabs( r[(size_t)i] ) <= rFloor ) ? 0. : -rho[(size_t)i] * r[(size_t)i];   // a row inside the floor is left alone: chasing it raises G for nothing
            }
            for( long int a = 0; a < nf; a++ ) rhs[(size_t)a] = -gt[(size_t)a];
            if( !factorDense( K ) ) { dual = std::min( dual * 100., 1e-10 ); continue; }
            solveFactored( rhs );
            std::vector<double> x = xsol;
            for( int ref = 0; ref < 2; ref++ )   // iterative refinement
            {
                std::vector<double> res( rhs );
                for( long int a = 0; a < M; a++ )
                    for( long int b = 0; b < M; b++ ) res[(size_t)a] -= K[(size_t)a*(size_t)M + (size_t)b] * x[(size_t)b];
                solveFactored( res );
                for( long int a = 0; a < M; a++ ) x[(size_t)a] += xsol[(size_t)a];
            }
            gDotDn = 0.; dMaxRel = 0.;
            for( long int a = 0; a < nf; a++ )
            {
                dn[(size_t)a] = s[(size_t)a] * x[(size_t)a];
                gDotDn += F[(size_t)fr[(size_t)a]] * dn[(size_t)a];
                dMaxRel = std::max( dMaxRel, std::fabs( x[(size_t)a] ) );
            }
            for( long int q = 0; q < nr; q++ ) Uw[(size_t)q] = x[(size_t)( nf + q )];
            if( tryMu == 0 ) Uw0 = Uw;
            {   // exact-penalty merit G + lam*|r|_1 (lam above the multipliers): a step that repairs the mass balance may raise G
                double um = 0.;
                for( long int q = 0; q < nr; q++ ) um = std::max( um, std::fabs( rho[(size_t)rows[(size_t)q]] * x[(size_t)( nf + q )] ) );
                lam = 2. * um + 1.;
                double rl1 = 0.;   // residual the LINEARISED step leaves (a row without a free carrier cannot be repaired)
                for( long int i = 0; i < N; i++ )
                {
                    double v = r[(size_t)i];
                    for( long int a = 0; a < nf; a++ ) v += pm.A[ i + fr[(size_t)a]*N ] * dn[(size_t)a];
                    rl1 += rex( v, rFloor );
                }
                dMerit = gDotDn + lam * ( rl1 - r1 );
            }
            if( dMerit <= 1e-14 * ( 1. + std::fabs( G ) ) || std::fabs( gDotDn ) <= 1e-10 * ( 1. + std::fabs( G ) ) ) solved = true;   // descent (r ~ 0), else damp harder
        }
        }
        if( !solved && r1 == 0. && Uw0.size() == (size_t)nr )
        {   // every row is inside the noise floor: the only step left is a repair that raises G, so the set is stationary here
            std::fill( dn.begin(), dn.end(), 0. ); gDotDn = 0.; dMaxRel = 0.; solved = true; nullStep = true;
        }
        if( !solved ) { why = "nodescent"; break; }
        if( nullStep ) Uw = Uw0;   // the damped solves' multipliers are biased; the undamped one belongs to this set
        std::fill( Un.begin(), Un.end(), 0. );
        for( long int q = 0; q < nr; q++ ) Un[(size_t)rows[(size_t)q]] = rho[(size_t)rows[(size_t)q]] * Uw[(size_t)q];

        // reduced gradient of the free species at the new multipliers, and convergence
        double redMax = 0.;
        for( long int a = 0; a < nf; a++ )
        {
            const long int j = fr[(size_t)a];
            double g = F[(size_t)j];
            for( long int i = 0; i < N; i++ ) g -= Un[(size_t)i] * pm.A[ i + j*N ];
            redMax = std::max( redMax, std::fabs( g ) );
        }
        const double rMax = massResid( n, r );
        const bool stationary = ( dMaxRel <= 1e-8 || std::fabs( gDotDn ) <= 1e-10 * ( 1. + std::fabs( G ) ) )
                                && redMax <= ( nullStep ? 1e-5 : 1e-7 ) && rMax <= rFloor;
        if( kStall > 0 && !stationary )
        {
            if( changes != bestRedChanges ) { bestRed = 1e300; bestRedIt = it; bestRedChanges = changes; }
            if( redMax < 0.5 * bestRed ) { bestRed = redMax; bestRedIt = it; stallG = G; stepMin = stepMax = dMaxRel; }
            else
            {
                stepMin = std::min( stepMin, dMaxRel ); stepMax = std::max( stepMax, dMaxRel );
                if( it - bestRedIt >= kStall )
                {
                    // a finish whose G still falls, or whose step still changes, is slow, not dead: the window starts again
                    if( G < stallG - 1e-12 * ( 1. + std::fabs( G ) ) || stepMax > 1.01 * stepMin )
                    { bestRedIt = it; stallG = G; stepMin = stepMax = dMaxRel; stallSkips++; }
                    else { why = "stall"; break; }
                }
            }
        }
        if( stationary )
        {
            // stationary on this set: does any species outside it want in?
            long int best = -1; double bestG = -1e-6;
            for( long int j = 0; j < L; j++ )
            {
                if( freeSp[(size_t)j] || tabu[(size_t)j] || degenerate( j ) ) continue;
                const long int k = phOf[(size_t)j];
                bool carrier = ( j >= pm.Ls );
                if( !carrier )
                    for( long int jb = phStart[(size_t)k]; jb < phStart[(size_t)k+1]; jb++ )
                        if( freeSp[(size_t)jb] ) { carrier = true; break; }
                double g = F[(size_t)j];
                for( long int i = 0; i < N; i++ ) g -= Un[(size_t)i] * pm.A[ i + j*N ];
                if( !carrier ) continue;
                if( g < bestG ) { bestG = g; best = j; }
            }
            long int tpdPhase = -1; std::vector<double> tpdY;
            if( best < 0 && changes < kMaxChanges )
            {
                // no species wants in: does an ABSENT non-ideal condensed phase lower the tangent plane at these potentials?
                for( long int i = 0; i < N; i++ ) pm.U[i] = Un[(size_t)i];
                double bestTpd = -1e-6;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int nk = pm.L1[k];
                    const char ph = pm.PHC[k];
                    if( nk < 2 || ph == PH_AQUEL || ph == PH_GASMIX || ph == PH_PLASMA || ph == PH_FLUID || ph == PH_SORPTION
                        || ph == PH_POLYEL || ph == PH_ADSORPT || ph == PH_IONEX ) continue;
                    bool absent = true, blocked = false;
                    for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ )
                    { if( freeSp[(size_t)j] ) absent = false; if( tabu[(size_t)j] || degenerate( j ) ) blocked = true; }
                    if( !absent || blocked ) continue;
                    std::vector<double> yb;
                    const double tpd = NativeTpdPhase( k, phStart[(size_t)k], yb );
                    if( tpd > -1e299 && tpd < bestTpd ) { bestTpd = tpd; tpdPhase = k; tpdY = yb; }
                }
                evalF( n, F );   // NativeTpdPhase restores pm.*, this leaves them at the base point regardless
            }
            if( tpdPhase >= 0 )
            {
                const double seed = 1e-6 * bScale;
                for( long int a = 0; a < pm.L1[tpdPhase]; a++ )
                {
                    const long int j = phStart[(size_t)tpdPhase] + a;
                    n[(size_t)j] = std::max( n[(size_t)j], seed * tpdY[(size_t)a] + dcFloor * 10. );
                    freeSp[(size_t)j] = 1; tabu[(size_t)j] = 1;
                }
                changes++; entered++;
                G = evalF( n, F );
                continue;
            }
            if( best < 0 || changes >= kMaxChanges ) { U = Un; converged = ( best < 0 ); why = converged ? "ok" : "changes"; break; }
            freeSp[(size_t)best] = 1; tabu[(size_t)best] = 1;
            double tot = 0.; for( long int j = 0; j < L; j++ ) tot += std::fabs( n[(size_t)j] );
            n[(size_t)best] = std::max( n[(size_t)best], 1e-6 * tot / (double)L + dcFloor * 10. );
            changes++; entered++;
            G = evalF( n, F );
            continue;
        }

        // step: fraction to the boundary, then Armijo on G
        double alpha = 1.;
        for( long int a = 0; a < nf; a++ )
        {
            const long int j = fr[(size_t)a];
            if( dn[(size_t)a] < 0. )
            {
                const double lim = ( j >= pm.Ls ) ? 1.0 : 0.9;
                alpha = std::min( alpha, lim * n[(size_t)j] / ( -dn[(size_t)a] ) );
            }
            else if( dn[(size_t)a] > 0. && xhi[(size_t)j] < 1e300 )
                alpha = std::min( alpha, 0.9 * ( xhi[(size_t)j] - n[(size_t)j] ) / dn[(size_t)a] );
        }
        std::vector<double> nn( n );
        double Gn = G; bool accepted = false;
        for( int ls = 0; ls < 30; ls++ )
        {
            for( long int a = 0; a < nf; a++ ) nn[(size_t)fr[(size_t)a]] = n[(size_t)fr[(size_t)a]] + alpha * dn[(size_t)a];
            Gn = evalF( nn, Fp );
            massResid( nn, r );
            double rn1 = 0.; for( long int i = 0; i < N; i++ ) rn1 += rex( r[(size_t)i], rFloor );
            if( Gn + lam * rn1 <= G + lam * r1 + 1e-4 * alpha * dMerit + 1e-12 * ( 1. + std::fabs( G ) ) ) { accepted = true; break; }
            alpha *= 0.5;
        }
        if( !accepted ) { why = "linesearch"; break; }
        n = nn; F = Fp; G = Gn;
        // a pure phase driven to zero leaves the set (one change, tabu on re-entry)
        for( long int a = 0; a < nf; a++ )
        {
            const long int j = fr[(size_t)a];
            if( n[(size_t)j] <= ( j >= pm.Ls ? 1e-9 : 1e-11 ) * bScale )
            {
                n[(size_t)j] = std::max( xlo[(size_t)j], dcFloor );
                freeSp[(size_t)j] = 0; tabu[(size_t)j] = 1; changes++; left++;
                G = evalF( n, F );
            }
        }
        // a solution phase whose free end-members have all fallen to a trace and are still falling is vanishing: its members
        // leave together (one at a time each would creep, because the phase curvature is large at a trace amount)
        for( long int k = 0; k < pm.FIs; k++ )
        {
            double tot = 0., totOld = 0.; long int cnt = 0;
            for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ )
                if( freeSp[(size_t)j] ) { tot += n[(size_t)j]; cnt++; }
            if( cnt < 1 || tot > 1e-8 * bScale ) continue;
            for( long int a = 0; a < nf; a++ ) if( phOf[(size_t)fr[(size_t)a]] == k ) totOld += n[(size_t)fr[(size_t)a]] - alpha * dn[(size_t)a];
            if( tot >= totOld ) continue;
            for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ )
                if( freeSp[(size_t)j] )
                {
                    n[(size_t)j] = std::max( xlo[(size_t)j], dcFloor );
                    freeSp[(size_t)j] = 0; tabu[(size_t)j] = 1; left++;
                }
            changes++; vanished = true; G = evalF( n, F );
        }
        if( changes > kMaxChanges ) { why = "changes"; break; }
    }

    if( converged )
    {
        // a trace species of a phase with no free species (its own mole fraction is set by the floors of its neighbours, so its
        // reduced gradient there says nothing) goes to the floor, or the post-solve KKT check reads it as interior; a trace
        // species of a PRESENT phase is interior at its own amount and is relaxed to stationarity instead (snapping it to the
        // floor makes its reduced gradient about -ln(n/floor), a violation the check rightly reports)
        std::vector<char> phHasFree( (size_t)pm.FI, 0 );
        for( long int j = 0; j < L; j++ ) if( freeSp[(size_t)j] ) phHasFree[(size_t)phOf[(size_t)j]] = 1;
        std::vector<long int> traces;
        for( long int j = 0; j < L; j++ )
            if( !freeSp[(size_t)j] && n[(size_t)j] <= 1e-11 * bScale && !degenerate( j ) )
            {
                if( phHasFree[(size_t)phOf[(size_t)j]] ) traces.push_back( j );
                else if( n[(size_t)j] > lowTol( j ) && !( keepPhases && presentAtStart[(size_t)phOf[(size_t)j]] ) )
                    n[(size_t)j] = std::max( xlo[(size_t)j], dcFloor );   // never removes a phase that was present
            }
        if( !traces.empty() )
        {
            for( int pass = 0; pass < 8; pass++ )
            {
                evalF( n, F );
                double maxg = 0.;
                for( long int j : traces )
                {
                    double g = F[(size_t)j];
                    for( long int i = 0; i < N; i++ ) g -= U[(size_t)i] * pm.A[ i + j*N ];
                    maxg = std::max( maxg, std::fabs( g ) );
                    double nn = n[(size_t)j] * std::exp( -std::max( -20., std::min( 20., g ) ) );
                    n[(size_t)j] = std::max( std::max( xlo[(size_t)j], dcFloor ), std::min( nn, 1e-10 * bScale ) );
                }
                if( maxg < 1e-6 ) break;
            }
        }
        G = evalF( n, F );
        // publish, then hold the result to the mass balance; a state that does not hold it is not returned
        for( long int j = 0; j < L; j++ ) { pm.Y[j] = n[(size_t)j]; pm.X[j] = n[(size_t)j]; }
        // From a reported-OK start G may not rise at all (rounding only).
        const double gTol = keepPhases ? 1e-13 : ( vanished ? 5e-8 : 1e-9 );
        if( !seededStart && G > G0 + gTol * ( 1. + std::fabs( G0 ) ) ) { converged = false; why = "Gup"; }
        else
        {
            for( long int i = 0; i < N; i++ ) if( U[(size_t)i] != 0. ) pm.U[i] = U[(size_t)i];
            // Repair the mass balance whenever it is worse than the start's (the absolute
            // cutoff of CheckMassBalanceResiduals() cannot see a ~1e-10 mol trace move).
            const double rNow = massResid( n, rTmp );
            if( rNow > std::max( rStart, 1e-15 * bScale ) || CheckMassBalanceResiduals( pm.Y ) >= 0 )
            {
                MassBalanceReproject( pm.Y );
                for( long int j = 0; j < L; j++ ) { pm.X[j] = pm.Y[j]; n[(size_t)j] = pm.Y[j]; }
                if( CheckMassBalanceResiduals( pm.Y ) >= 0 ) { converged = false; why = "massbalance"; }
            }
            if( converged && keepPhases )
            {
                massResid( n, rTmp );
                for( long int i = 0; converged && i < N; i++ )   // no ROW may end worse than it started
                    if( std::fabs( rTmp[(size_t)i] ) > std::max( std::fabs( rStartRow[(size_t)i] ) * ( 1. + 1e-6 ), 1e-16 * bScale ) )
                    { converged = false; why = "mbworse"; }
                for( long int k = 0; converged && k < pm.FI; k++ )
                {
                    if( !presentAtStart[(size_t)k] ) continue;
                    double t = 0.;
                    for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ ) t += n[(size_t)j];
                    if( !( t > pm.DSM ) ) { converged = false; why = "phaselost"; }
                }
                if( converged )
                {
                    // From a reported-OK start the finish is kept only if it strictly lowers G: a
                    // result at the same G changes nothing chemically but still overwrites the
                    // multipliers, which on a system with a free dual direction can fail the
                    // post-solve KKT check.
                    const double G2 = evalF( n, F );
                    if( !( G2 < G0 - 1e-12 * ( 1. + std::fabs( G0 ) ) ) ) { converged = false; why = "noimprove"; }
                    else G = G2;
                }
            }
        }
    }
    native_trace_decide( "finish ok=%d why=%s it=%ld changes=%ld entered=%ld left=%ld G0=%.12g G1=%.12g%s",
                         converged ? 1 : 0, why.c_str(), (long)it, (long)changes, (long)entered, (long)left, G0, G,
                         stallSkips ? ( " stallskips=" + std::to_string( stallSkips ) ).c_str() : "" );
    if( !converged ) { restore(); return false; }
    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
    CalculateActivityCoefficients( LINK_UX_MODE );
    return true;
}

bool TMultiBase::LPFeasibilitySeed( std::vector<double>& nOut )
{
    const long int N = pm.N;
    const long int L = pm.L;
    auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
    if( !TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nOut ) )
        return false;
    if( (long int)nOut.size() != L )
        return false;

    // Self-check before trusting the simplex output: A*n = b and n >= 0 within a tolerance
    // scaled by the bulk composition.
    double bScale = 1.;
    for( long int i = 0; i < N; i++ )
        bScale = std::max( bScale, std::fabs( pm.B[i] ) );
    const double tol = std::max( 1e-6, 1e-8 * bScale );
    double maxResid = 0.;
    for( long int i = 0; i < N; i++ )
    {
        double sum = 0.;
        for( long int j = 0; j < L; j++ )
            sum += pm.A[ i + j*N ] * nOut[j];
        maxResid = std::max( maxResid, std::fabs( sum - pm.B[i] ) );
    }
    double minVal = 0.;
    for( long int j = 0; j < L; j++ )
        minVal = std::min( minVal, nOut[j] );
    if( maxResid > tol || minVal < -tol )
    {
        ipm_logger->debug( "LPFeasibilitySeed: self-check failed (maxResid={}, minVal={}, tol={}) - discarding",
                           maxResid, minVal, tol );
        return false;
    }
    return true;
}

// Dual of the linearised-Gibbs LP: min sum_j G0[j]*n_j s.t. A n = b, n >= 0. Same simplex
// and rows as LPFeasibilitySeed(), different objective: pricing a species against this dual,
// s_j = G0[j] - sum_i y_i A[i,j], measures how far it is from being stable.
// Returns false without touching yOut if the LP or its dual is not trustworthy.
bool TMultiBase::LPGibbsDual( std::vector<double>& yOut, const double* cost )
{
    // `cost` defaults to pm.G0. Passing pm.G prices against the current chemical potentials
    // (G0 + fDQF + F0); used only by LpFillProbeReport()'s LPRELP record.
    const long int N = pm.N;
    const long int L = pm.L;
    if( N <= 0 || L <= 0 || pm.G0 == nullptr )
        return false;
    auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
    std::vector<double> nDummy;
    std::vector<double> y( (size_t)N, 0. );
    if( !TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nDummy, cost ? cost : pm.G0, y.data() ) )
        return false;
    for( long int i = 0; i < N; i++ )
        if( !std::isfinite( y[(size_t)i] ) )
            return false;
    yOut.swap( y );
    return true;
}

long int TMultiBase::WorstPhaseStabilityViolation( double presenceThreshold, double dcFloor,
                                                   const char* exemptSpecies,
                                                   double& violOut, bool& wasAbsentOut,
                                                   std::vector<PhStabViolation>* rankedOut,
                                                   PhStabCensus* censusOut, const char* probeStage )
{
    gems3k_timing::Scope timing( gems3k_timing::Wpsv );
    const BASE_PARAM *pa_p = base_param();
    const long int L = pm.L;
    std::vector<PhStabViolation> ranked;
    PhStabCensus census;
    census.phases = pm.FI;

    // Phase-assemblage stability check: a second correctness signal beside the per-species
    // KKT check, which cannot tell the current assemblage from a different, wrong local
    // optimum where a phase was never let in (or should have left). Uses
    // StabilityIndexes(), as native's PhaseSelect() does, computed from the converged dual.
    // Its inputs are pm.Gamma[], pm.Y_la[] (CalculateConcentrations() with the final pm.U[]),
    // pm.fDQF[]/pm.sMod[] and pm.YF[]/pm.YFA[], which are refreshed from pm.Y here.
    TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
    StabilityIndexes();

    double worstStabilityViol = 0.;
    long int worstStabilityPhase = -1;
    bool worstStabilityWasAbsent = false;
    long int j0stab = 0;
    for( long int k = 0; k < pm.FI; k++ )
    {
        const long int j1stab = j0stab + pm.L1[k];
        // A phase with any end-member under a kinetic restriction - DUL[j] < 1e6 (including
        // exactly 0) or DLL[j] > 0 - is exempt, as in native's KinConstrDC/KinConstrPh: such
        // a phase is expected to look "stable but absent" or "present but unstable".
        bool kinConstrPh = false;
        for( long int j = j0stab; j < j1stab; j++ )
            if( pm.DUL[j] < 1e6 || pm.DLL[j] > 0.0 )
            { kinConstrPh = true; break; }
        // Same exemption for a phase this call's phase-extinction retry deactivated: it is
        // interchangeable with a present twin and so always reports logSI ~ 0.
        if( !kinConstrPh && L > 0 && exemptSpecies != nullptr )
            for( long int j = j0stab; j < j1stab && j < L; j++ )
                if( exemptSpecies[(size_t)j] ) { kinConstrPh = true; break; }
        j0stab = j1stab;
        if( kinConstrPh ) { census.exempt++; continue; }
        census.scanned++;

        const double logSI = pm.Falp[k];
        // Is this phase's stability index a measurement or a guard value? StabilityIndexes()
        // clamps each species' dual activity to [-608, +609] before exponentiating. Detected
        // here from pm.NMU[j] = log(exp(ln_ax_dual)/gamma) with the effective gamma
        // pm.Gamma[j] (pm.K2 is 0 on this path; gamma outside [1e-33, 1e33] taken as 1).
        bool clampedPh = false;
        for( long int j = j1stab - pm.L1[k]; j < j1stab && j < L; j++ )
        {
            double g = pm.Gamma[j];
            if( g < 1e-33 || g > 1e33 ) g = 1.;
            if( g <= 0. || !std::isfinite( g ) ) continue;
            const double lnAxDual = pm.NMU[j] + std::log( g );
            if( lnAxDual >= 609. - 1e-6 || lnAxDual <= -608. + 1e-6 )
            { clampedPh = true; break; }
        }
        if( clampedPh ) census.clamped++;
        // "Present" is also decided structurally: a phase whose every species sits at its own
        // lower bound is absent, whatever its total (the magnitude rule alone grows with the
        // phase's species count while the threshold does not). Can only move a phase from
        // present to absent.
        bool anyOffFloor = false;
        for( long int j = j1stab - pm.L1[k]; j < j1stab && j < L; j++ )
        {
            const double lo = std::max( pm.DLL[j], dcFloor );
            if( pm.Y[j] > lo * ( 1. + 1e-6 ) ) { anyOffFloor = true; break; }
        }
        // ... and a third rule, OR-ed in: a phase whose every species holds a negligible
        // fraction (below kRelPresenceEps) of the most of that species the bulk composition
        // could supply is absent, whatever its absolute amount. The ceiling is
        // n_max(j) = min over ICs i with a(j,i) > 0 of B[i]/a(j,i), scale-invariant. The
        // charge IC is excluded (B[Zz] = 0 would give every cation a zero ceiling), and a
        // species with a genuine zero ceiling counts as absent. Like the structural rule, this
        // can only move a phase from present to absent.
        static const double kRelPresenceEps = 1e-12;
        bool anySignificant = false;
        for( long int j = j1stab - pm.L1[k]; j < j1stab && j < L; j++ )
        {
            double nMax = -1.;
            for( long int i = 0; i < pm.N; i++ )
            {
                if( pm.ICC != nullptr && pm.ICC[i] == IC_CHARGE ) continue;
                const double coef = pm.A[ i + j*pm.N ];
                if( coef <= 0. ) continue;
                const double icBound = pm.B[i] / coef;
                if( nMax < 0. || icBound < nMax ) nMax = icBound;
            }
            if( nMax <= 0. ) continue;      // trap (b): genuine zero ceiling -> absent
            if( pm.Y[j] >= kRelPresenceEps * nMax ) { anySignificant = true; break; }
        }
        // A water phase that has vanished (its water at or below XwMinM, as in PrimalChemicalPotentials() and the collapse
        // check, and its total below 1e-9 of all phases) counts as absent, so a water phase that should come back is judged by
        // its stability index (pa_OptimaCollapseCheck). A dried-out brine that keeps its solutes is not affected.
        double allPhases = 0.;
        for( long int kk = 0; kk < pm.FI; kk++ ) allPhases += std::max( pm.YF[kk], 0. );
        const bool aqDry = pa_p->OptimaCollapseCheck && pm.PHC[k] == PH_AQUEL && pm.LO >= 0 && pm.Y[pm.LO] <= pm.XwMinM
                           && pm.YF[k] <= 1e-9 * allPhases;
        const bool present = anyOffFloor && ( pm.YF[k] >= presenceThreshold ) && anySignificant && !aqDry;
        if( aqDry && logSI > pa_p->DF && anyOffFloor && pm.YF[k] >= presenceThreshold && anySignificant )
            native_trace_decide( "aqdryabsent water=%.3e logSI=%.3f stage=%s", pm.Y[pm.LO], logSI, probeStage ? probeStage : "check" );
        // GEMS3K_PHSTAB_PROBE=<path>: one line per phase per evaluation - amount, stability
        // index, the three presence clauses and the verdict. Zero cost when unset.
        {
            static FILE* dbg = []() -> FILE* {
                const char* fn = std::getenv( "GEMS3K_PHSTAB_PROBE" );
                return fn ? fopen( fn, "a" ) : nullptr; }();
            if( dbg )
                fprintf( dbg, "PH k=%ld %-16s YF=%.6e logSI=%+.6e present=%d offFloor=%d "
                              "signif=%d  DF=%.2e DFM=%.2e  verdict=%s\n",
                         (long)k, char_array_to_string( pm.SF[k], MAXPHNAME ).c_str(),
                         pm.YF[k], logSI, (int)present, (int)anyOffFloor, (int)anySignificant,
                         pa_p->DF, pa_p->DFM,
                         ( !present && logSI > pa_p->DF ) ? "ABSENT-BUT-STABLE"
                       : (  present && logSI < -pa_p->DFM ) ? "PRESENT-BUT-UNSTABLE" : "ok" );
        }
            // Same thresholds as native's PhaseSelect(): a stable phase (logSI > DF) not in
            // the assemblage, or an unstable one (logSI < -DFM) that is, is a wrong assemblage.
        if( !present && logSI > pa_p->DF )
        {
            const double viol = logSI - pa_p->DF;
            census.absentStable++;
            ranked.push_back( PhStabViolation{ k, viol, true, clampedPh } );
            if( viol > worstStabilityViol )
            { worstStabilityViol = viol; worstStabilityPhase = k; worstStabilityWasAbsent = true; }
        }
        else if( present && logSI < -pa_p->DFM )
        {
            const double viol = -pa_p->DFM - logSI;
            census.presentUnstable++;
            ranked.push_back( PhStabViolation{ k, viol, false, clampedPh } );
            if( viol > worstStabilityViol )
            { worstStabilityViol = viol; worstStabilityPhase = k; worstStabilityWasAbsent = false; }
        }
    }
    // Descending by size, ties broken by phase index (clamped entries share one viol).
    std::sort( ranked.begin(), ranked.end(),
               []( const PhStabViolation& a, const PhStabViolation& b )
               { return a.viol != b.viol ? a.viol > b.viol : a.k < b.k; } );
    {
        static FILE* dbg = []() -> FILE* {
            const char* fn = std::getenv( "GEMS3K_PHSTAB_PROBE" );
            return fn ? fopen( fn, "a" ) : nullptr; }();
        if( dbg )
            // stage: which step asked; try: whether the quick first try of a cold call is running
            fprintf( dbg, "CENSUS phases=%ld exempt=%ld scanned=%ld clamped=%ld "
                          "absent-but-stable=%ld present-but-unstable=%ld worst=%ld stage=%s try=%s\n",
                     (long)census.phases, (long)census.exempt, (long)census.scanned,
                     (long)census.clamped, (long)census.absentStable,
                     (long)census.presentUnstable, (long)worstStabilityPhase,
                     probeStage ? probeStage : "check", g_optimaFixAccept ? "quick" : "full" );
    }
    if( rankedOut != nullptr ) rankedOut->swap( ranked );
    if( censusOut != nullptr ) *censusOut = census;
    violOut = worstStabilityViol;
    wasAbsentOut = worstStabilityWasAbsent;
    return worstStabilityPhase;
}


// ---------------------------------------------------------------------------
// Species-level dimension-reduction pre-solve (pa_OptimaDimReduce)
// ---------------------------------------------------------------------------
// An active-set method, exact at its fixed point. Hold every omitted species j at its
// lower bound l_j; the mass-balance rows become
//     sum_{s in S} A[i, s] * x_s  =  b_i - sum_{j not in S} A[i, j] * l_j
// and the reduced problem minimises G over x_S subject to those rows and the remaining
// boxes. Its dual y is a full N-vector (every IC row survives), so the omitted columns are
// priced on it:
//     s_j = F[j] - sum_i U[i] A[i,j].
// s_j >= 0 is the full problem's KKT condition for a variable at its lower bound (Optima's
// is_lower_unstable test). When a pass readmits nothing, the reduced answer solves the full
// problem.
//  1. The initial set comes from two LPs over pm.A/pm.B, independent of any solver state:
//     the LP-feasibility seed's support, widened by every species priced below
//     pa_OptimaDimReduceTol against the LP-Gibbs dual.
//  2. Readmission is monotone (a species once active is never dropped), which bounds the
//     loop. Removal is left to the phase-selection loop and the phase-extinction tier.
// Never returns an answer on its own: it leaves its primal in pm.Y[] and its dual in pm.U[]
// for the full-dimension solve to warm-start from and verify. A pre-solve that fails or
// cannot omit anything is discarded.
// In plain words: solve with only the likely species first, add back any species the
// result says is missing, and repeat until nothing is missing.
bool TMultiBase::OptimaReducedPreSolve( long int maxPasses, double dcFloor, double dimTol,
                                        long int passBudget,
                                        long int& iterationsOut, long int& activeOut )
{
    gems3k_timing::Scope timing( gems3k_timing::PreSolve );
    iterationsOut = 0;
    activeOut = 0;
    const long int N = pm.N;
    const long int L = pm.L;
    if( L <= 0 || N <= 0 || maxPasses <= 0 )
        return false;

    const BASE_PARAM* pa_p = base_param();
    const bool kMoleFracHessian = ( pa_p->OptimaMoleFracHessian != 0 );
    const bool kFDHessian       = ( pa_p->OptimaFDHessian != 0 );
    const double kLogBarrierTau = pa_p->LogBarrierTau;
    const double kPhaseHessianFloor = pa_p->PhaseHessianFloor;
    if( phfTwin.size() != (size_t)pm.FIs ) PhaseHessianBegin();
    const std::vector<char> phTwin = phfTwin;
    const bool hasAq = HasAqueousPhase();

    // Box bounds, built as the full path builds them (pm.DUL[j] < 1e6 marks a real, possibly
    // zero, upper restriction).
    std::vector<double> xlo( (size_t)L ), xhi( (size_t)L );
    for( long int j = 0; j < L; j++ )
    {
        xlo[(size_t)j] = std::max( pm.DLL[j], dcFloor );
        xhi[(size_t)j] = ( pm.DUL[j] < 1e6 )
                          ? std::max( pm.DUL[j], dcFloor ) : std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 10.;
    }

    // Initial active set: the seed's own support. A species with a degenerate box (a hard
    // kinetic exclusion) is never active.
    std::vector<char> act( (size_t)L, 0 );
    for( long int j = 0; j < L; j++ )
    {
        if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;
        if( pm.Y[j] > xlo[(size_t)j] * ( 1. + 1e-6 ) ) act[(size_t)j] = 1;
    }
    // The solvent is never omitted.
    if( hasAq && pm.LO >= 0 && pm.LO < L && xhi[(size_t)pm.LO] > xlo[(size_t)pm.LO] )
        act[(size_t)pm.LO] = 1;

    // ... and widen that set by pricing every species against the dual of the linearised-
    // Gibbs LP (pa_OptimaDimReduceTol). The seed's support alone is a vertex (N species),
    // too sparse for large systems. Pass 0 starts from a zero dual; later passes inherit
    // their predecessor's.
    {
        std::vector<double> lpDual;
        if( dimTol != 0. && pm.G0 != nullptr && LPGibbsDual( lpDual ) )
        {
            // Price every candidate once; the sign of the rule then selects how prices
            // become a set: > 0 an absolute threshold in RT, < 0 the rank rule -|tol| x N
            // (the cheapest-priced species until the set reaches that multiple of the IC count).
            std::vector<std::pair<double,long int>> priced;
            for( long int j = 0; j < L; j++ )
            {
                if( act[(size_t)j] ) continue;
                if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;
                double z = 0.;
                for( long int i = 0; i < N; i++ ) z += lpDual[(size_t)i] * pm.A[ i + j*N ];
                priced.push_back( std::make_pair( pm.G0[j] - z, j ) );
            }
            long int added = 0;
            if( dimTol > 0. )
            {
                for( size_t k = 0; k < priced.size(); k++ )
                    if( priced[k].first < dimTol ) { act[(size_t)priced[k].second] = 1; added++; }
                ipm_logger->debug( "OptimaReducedPreSolve: LP-Gibbs pricing admitted {} extra species "
                                   "(tol={} RT)", added, dimTol );
            }
            else
            {
                long int already = 0;
                for( long int j = 0; j < L; j++ ) if( act[(size_t)j] ) already++;
                const long int target = (long int)( -dimTol * (double)N + 0.5 );
                std::sort( priced.begin(), priced.end() );
                for( size_t k = 0; k < priced.size() && already + added < target; k++ )
                    { act[(size_t)priced[k].second] = 1; added++; }
                ipm_logger->debug( "OptimaReducedPreSolve: LP-Gibbs rank rule admitted {} extra species "
                                   "(target {} = {} x N, seed support {})",
                                   added, target, -dimTol, already );
            }
        }
    }

    std::vector<long int> nxToJ, jToNx( (size_t)L, -1 );
    std::vector<double>   Fbase( (size_t)L, 0. );

    // Only a fixed point is handed over. An intermediate pass's answer solves a different
    // (smaller) problem, so its dual is wrong about every species not yet active; on any
    // exit other than "priced nothing back in", restore what the caller had.
    std::vector<double> Y0( pm.Y, pm.Y + L );
    std::vector<double> U0( pm.U, pm.U + N );
    auto discard = [&]() -> bool {
        for( long int j = 0; j < L; j++ ) pm.Y[j] = Y0[(size_t)j];
        for( long int i = 0; i < N; i++ ) pm.U[i] = U0[(size_t)i];
        activeOut = 0;
        return false;
    };

    Optima::Options options;
    // max(2000, pa_IIM), or less on the first attempt (pa_OptimaPreSolveFirstIters).
    options.maxiters = (unsigned)passBudget;
    options.convergence.tolerance = pa_p->OptimaTol;
    apply_optima_linesearch( options, pa_p->OptimaLineSearch, pa_p->OptimaLSStallEscape,
                             pa_p->OptimaLSRejectWorse );
    apply_optima_lu( options, pa_p->OptimaPartialPivLU );

        // ---- Stall / wall-clock guard for the pre-solve ----
        // A pass that has stopped improving is discarded early instead of running to its full
        // budget (it would be discarded anyway). The wall-clock budget is shared across all
        // passes. As on the full path, `check` returning true means converged to Optima, so the
        // flag is folded back to succeeded = false, which lands on the discard path.
    struct PreStallWatch {
        long int window = 0;
        double bestErr = 0., bestComp = 0.;
        // Best-so-far one whole window ago, plus the live error's range within the current
        // window (both are needed - see the test).
        double refErr = 0., refComp = 0., winLo = 0., winHi = 0.;
        long int run = 0;              // iterations elapsed in the current window
        bool stalled = false, timedOut = false;
        double maxSeconds = 0.;
        std::chrono::steady_clock::time_point started;
        void reset() {
            bestErr = bestComp = std::numeric_limits<double>::infinity();
            refErr  = refComp  = std::numeric_limits<double>::infinity();
            winLo = std::numeric_limits<double>::infinity(); winHi = 0.;
            run = 0; stalled = false;
        }
        bool overBudget() const {
            if( g_optimaIterBudget == 0 ) return true;   // the call's iteration budget is used up
            if( maxSeconds <= 0. ) return false;
            return std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - started ).count() > maxSeconds;
        }
    };
    auto preWatch = std::make_shared<PreStallWatch>();
    preWatch->window     = pa_p->OptimaStallWindow;
    preWatch->maxSeconds = pa_p->OptimaMaxSeconds;
    preWatch->started    = std::chrono::steady_clock::now();
    preWatch->reset();
    if( preWatch->window > 0 || preWatch->maxSeconds > 0. || g_optimaIterBudget >= 0 )
        options.convergence.check =
            [preWatch]( Optima::ConvergenceCheckArgs const& args ) -> bool
            {
                if( g_optimaIterBudget > 0 ) --g_optimaIterBudget;
                if( preWatch->overBudget() ) { preWatch->timedOut = true; return true; }
                if( preWatch->window <= 0 ) return false;
                const double e = args.E.errorx();
                double comp = 0.;
                const auto ex = args.E.ex();
                const auto xv = args.u.x;
                const Optima::Index nx = ex.size() < xv.size() ? ex.size() : xv.size();
                for( Optima::Index j = 0; j < nx; ++j )
                {
                    const double c = std::fabs( ex[j] ) * std::fabs( xv[j] );
                    if( c > comp ) comp = c;
                }
                // Asks whether the pass is deadlocked, not whether it converges fast enough.
                // Over one whole window it declares a stall only when both
                //   (a) best-so-far has not fallen by a relative kStallRel against its value at
                //       the window's start, and
                //   (b) the live error has not moved: its range within the window is below the
                //       same relative amount.
                // (a) alone would kill a pass that converges through a long excursion with
                // best-so-far frozen; (b) is measured as a range over the window, so a slow
                // but smooth descent never looks static. A pass that is converging but too
                // slowly for its budget is a budget question, not a stall. `best` is updated on
                // every fall, however small.
                static const double kPreStallRelProgress = 1e-8;
                if( e < preWatch->bestErr )  preWatch->bestErr  = e;
                if( comp < preWatch->bestComp ) preWatch->bestComp = comp;
                if( e < preWatch->winLo ) preWatch->winLo = e;
                if( e > preWatch->winHi ) preWatch->winHi = e;
                if( ++preWatch->run < preWatch->window ) return false;

                const bool firstWindow = !std::isfinite( preWatch->refErr )
                                      && !std::isfinite( preWatch->refComp );
                const bool dropped =
                    ( std::isfinite( preWatch->refErr )
                      && preWatch->bestErr  < preWatch->refErr  * ( 1. - kPreStallRelProgress ) )
                 || ( std::isfinite( preWatch->refComp )
                      && preWatch->bestComp < preWatch->refComp * ( 1. - kPreStallRelProgress ) );
                const double lo = preWatch->winLo;
                const bool moved = ( preWatch->winHi - lo )
                                   > kPreStallRelProgress * std::max( lo, 1e-300 );
                if( !firstWindow && !dropped && !moved )
                {
                    preWatch->stalled = true;
                    return true;
                }
                // Window survived: roll the reference forward and start a new one.
                preWatch->refErr  = preWatch->bestErr;
                preWatch->refComp = preWatch->bestComp;
                preWatch->winLo   = std::numeric_limits<double>::infinity();
                preWatch->winHi   = 0.;
                preWatch->run     = 0;
                return false;
            };

    bool haveAnswer = false;
    for( long int pass = 0; pass < maxPasses; pass++ )
    {
        nxToJ.clear();
        std::fill( jToNx.begin(), jToNx.end(), -1 );
        for( long int j = 0; j < L; j++ )
            if( act[(size_t)j] ) { jToNx[(size_t)j] = (long int)nxToJ.size(); nxToJ.push_back( j ); }
        const long int nS = (long int)nxToJ.size();
        if( nS <= 0 )
            return false;
        if( nS >= L )
        {
            // Nothing left to omit: the reduced problem is the full one. Whatever the previous
            // pass produced (if any) stands as the warm start.
            ipm_logger->debug( "OptimaReducedPreSolve: active set reached the full species "
                               "count ({}) at pass {} - no reduction available", L, pass );
            // The previous pass converged and readmission has brought every species back: that
            // answer (readmitted species at their lower bound) is the warm start for the full
            // solve.
            if( haveAnswer )
            {
                native_trace_decide( "dimreducefull pass=%ld active=%ld of=%ld", (long)pass,
                                     (long)activeOut, (long)L );
                return true;
            }
            return discard();
        }

        Optima::Dims dims;
        dims.x  = nS;
        dims.be = N;
        Optima::Problem problem( dims );

        for( long int i = 0; i < N; i++ )
        {
            for( long int s = 0; s < nS; s++ )
                problem.Aex(i, s) = pm.A[ i + nxToJ[(size_t)s]*N ];
            // Omitted species are held at their lower bound; their constant contribution to
            // each IC row moves to the right-hand side.
            double rhs = pm.B[i];
            for( long int j = 0; j < L; j++ )
                if( !act[(size_t)j] )
                    rhs -= pm.A[ i + j*N ] * xlo[(size_t)j];
            problem.be[i] = rhs;
        }
        for( long int s = 0; s < nS; s++ )
        {
            const long int j = nxToJ[(size_t)s];
            problem.xlower[s] = xlo[(size_t)j];
            problem.xupper[s] = xhi[(size_t)j];
        }

        // Objective: as on the full path. The gradient of the full Gibbs energy with respect
        // to an active species is still its chemical potential. The phase-decay counters are
        // not maintained here (the full solve keeps its own).
        problem.f = [this, L, nS, dcFloor, kLogBarrierTau, kPhaseHessianFloor, phTwin,
                     kFDHessian, kMoleFracHessian, hasAq, &nxToJ, &jToNx, &xlo, &act, &Fbase]
                    ( Optima::ObjectiveResultRef res, Optima::VectorView x,
                      Optima::VectorView /*p*/, Optima::VectorView /*c*/,
                      Optima::ObjectiveOptions opts )
        {
            gems3k_timing::Scope timingObj( gems3k_timing::Obj );
            for( long int j = 0; j < L; j++ )
                pm.X[j] = act[(size_t)j] ? x[ jToNx[(size_t)j] ] : xlo[(size_t)j];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );

            // f is the full Gibbs energy (omitted species included), so the value stays
            // comparable with the full path's; the omitted terms are a constant shift.
            double fval = 0.;
            for( long int j = 0; j < L; j++ )
                fval += pm.X[j] * pm.F[j];
            for( long int j = pm.Ls; j < L; j++ )
                fval -= kLogBarrierTau * std::log( std::max( pm.X[j], dcFloor ) );
            for( long int s = 0; s < nS; s++ )
            {
                const long int j = nxToJ[(size_t)s];
                res.fx[s] = pm.F[j];
                if( j >= pm.Ls )
                    res.fx[s] -= kLogBarrierTau / std::max( pm.X[j], dcFloor );
            }
            res.f = fval;

            if( opts.eval.fxx )
            {
                res.fxx.setZero();
                res.diagfxx = false;
                // Per-multicomponent-phase ideal-mixing curvature, gathered into the reduced
                // slots; same expressions as the full path's.
                long int j0 = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int j1 = j0 + pm.L1[k];
                    const bool isAqueousLike = ( hasAq && pm.LO >= j0 && pm.LO < j1 );
                    double Xf = 0.;
                    for( long int j = j0; j < j1; j++ )
                        Xf += std::max( pm.X[j], dcFloor );
                    if( isAqueousLike )
                    {
                        const long int w = pm.LO;
                        const long int sw = jToNx[(size_t)w];
                        const double Xw = std::max( pm.X[w], dcFloor );
                        for( long int j = j0; j < j1; j++ )
                        {
                            if( j == w ) continue;
                            const long int sj = jToNx[(size_t)j];
                            if( sj < 0 ) continue;
                            const double Xj = std::max( pm.X[j], dcFloor );
                            res.fxx(sj,sj) = 1.0 / Xj;
                            if( sw >= 0 )
                            {
                                res.fxx(sj,sw) = -1.0 / Xw;
                                res.fxx(sw,sj) = -1.0 / Xw;
                            }
                        }
                        if( sw >= 0 )
                            res.fxx(sw,sw) = ( Xf - Xw ) / ( Xw * Xw );
                    }
                    else
                    {
                        // Same present-only restriction as the main objective.
                        double phTot = 0.;
                        for( long int j = j0; j < j1; j++ ) phTot += std::max( pm.X[j], 0. );
                        const double presThr = std::max( dcFloor * 1e3, phTot * 1e-6 );
                        for( long int j = j0; j < j1; j++ )
                        {
                            const long int sj = jToNx[(size_t)j];
                            if( sj < 0 ) continue;
                            const double Xj = std::max( pm.X[j], dcFloor );
                            if( kMoleFracHessian && pm.X[j] > presThr )
                            {
                                for( long int i = j0; i < j1; i++ )
                                {
                                    const long int si = jToNx[(size_t)i];
                                    if( si >= 0 && pm.X[i] > presThr ) res.fxx(sj,si) = -1.0 / Xf;
                                }
                                res.fxx(sj,sj) = 1.0 / Xj - 1.0 / Xf;
                            }
                            else
                                res.fxx(sj,sj) = 1.0 / Xj;
                        }
                    }
                    j0 = j1;
                }
                for( long int j = pm.Ls; j < L; j++ )
                {
                    const long int sj = jToNx[(size_t)j];
                    if( sj < 0 ) continue;
                    const double Xj = std::max( pm.X[j], dcFloor );
                    res.fxx(sj,sj) += kLogBarrierTau / ( Xj * Xj );
                }

                for( long int j = 0; j < L; j++ ) Fbase[(size_t)j] = pm.F[j];

                // pa_OptimaFDHessian: exact columns for Optima's own basic
                // variables. opts.ibasicvars indexes the REDUCED space.
                for( Optima::Index bk = 0; kFDHessian && bk < opts.ibasicvars.size(); bk++ )
                {
                    const long int si = (long int)opts.ibasicvars[bk];
                    if( si < 0 || si >= nS ) continue;
                    const long int i = nxToJ[(size_t)si];
                    const double Xi = pm.X[i];
                    const double h = std::max( std::fabs(Xi) * 1e-7, dcFloor * 10. );
                    pm.X[i] = Xi + h;
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                    for( long int sr = 0; sr < nS; sr++ )
                    {
                        const long int r = nxToJ[(size_t)sr];
                        res.fxx(sr,si) = ( pm.F[r] - Fbase[(size_t)r] ) / h;
                    }
                    pm.X[i] = Xi;
                }

                // pa_PhaseHessianFloor: exact, eigenvalue-floored curvature for the non-aqueous
                // multicomponent phases; only end-members both present and active take part.
                const double regRatio = kPhaseHessianFloor;
                if( regRatio > 0. )
                {
                    long int p0 = 0;
                    for( long int k = 0; k < pm.FIs; k++ )
                    {
                        const long int p1 = p0 + pm.L1[k];
                        const long int nEnd = p1 - p0;
                        const bool isAq = ( pm.LO >= p0 && pm.LO < p1 );
                        if( !isAq && nEnd > 1 && p1 <= L )
                        {
                            double phTot = 0.;
                            for( long int j = p0; j < p1; j++ ) phTot += std::max( pm.X[j], 0. );
                            std::vector<long int> pres;   // reduced indices
                            for( long int j = p0; j < p1; j++ )
                                if( jToNx[(size_t)j] >= 0
                                    && pm.X[j] > std::max( dcFloor * 1e3, phTot * 1e-6 ) )
                                    pres.push_back( jToNx[(size_t)j] );
                            const int nP = (int)pres.size();
                            if( nP > 1 )
                            {
                                for( int c = 0; c < nP; c++ )
                                {
                                    const long int i = nxToJ[(size_t)pres[(size_t)c]];
                                    const double Xi = pm.X[i];
                                    const double h = std::fabs(Xi) * 1e-7;
                                    if( !( h > 0. ) ) continue;
                                    pm.X[i] = Xi + h;
                                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                                    CalculateActivityCoefficients( LINK_UX_MODE );
                                    PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                                    for( int r = 0; r < nP; r++ )
                                    {
                                        const long int rj = nxToJ[(size_t)pres[(size_t)r]];
                                        res.fxx(pres[(size_t)r],pres[(size_t)c])
                                            = ( pm.F[rj] - Fbase[(size_t)rj] ) / h;
                                    }
                                    pm.X[i] = Xi;
                                }
                                std::vector<double> blk( (size_t)nP*nP );
                                for( int a = 0; a < nP; a++ )
                                    for( int b = 0; b < nP; b++ )
                                        blk[(size_t)a*nP+b] = 0.5 * ( res.fxx(pres[(size_t)a],pres[(size_t)b])
                                                                    + res.fxx(pres[(size_t)b],pres[(size_t)a]) );
                                bool indef = false;
                                double presTot = 0.;
                                for( int a = 0; a < nP; a++ ) presTot += pm.X[nxToJ[(size_t)pres[(size_t)a]]];
                                const double absFloor = ( presTot > 0. ) ? regRatio * nP / presTot : 0.;
                                const bool okEig = SymEigFloorInPlace( blk, nP, regRatio, &indef, absFloor );
                                if( okEig )
                                {
                                    if( indef && phTwin[(size_t)k] ) phfFired[(size_t)k]++;
                                    for( int a = 0; a < nP; a++ )
                                        for( int b = 0; b < nP; b++ )
                                            res.fxx(pres[(size_t)a],pres[(size_t)b]) = blk[(size_t)a*nP+b];
                                }
                            }
                        }
                        p0 = p1;
                    }
                }

                // Restore the base state the gradient above was computed from.
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                for( long int j = 0; j < L; j++ )
                    pm.F[j] = Fbase[(size_t)j];
            }
            res.succeeded = true;
        };
        memoize_objective_if_linesearch( problem, pa_p->OptimaLineSearch );

        Optima::State state( dims );
        for( long int s = 0; s < nS; s++ )
        {
            const long int j = nxToJ[(size_t)s];
            state.x[s] = std::min( std::max( pm.Y[j], problem.xlower[s] ), problem.xupper[s] );
        }
        // Seed the dual from the previous pass's answer, all or nothing (every IC row
        // survives the reduction, so this dual is complete).
        if( haveAnswer )
        {
            bool duals_usable = true;
            for( long int i = 0; i < N; i++ )
                if( !std::isfinite( pm.U[i] ) ) { duals_usable = false; break; }
            if( duals_usable )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -pm.U[i];
        }
        // pa_OptimaPreSolveSeed: the first pass starts from the linear-Gibbs solution and its element
        // potentials, each amount capped at what its elements allow, and gets 300 steps to settle.
        int ps0Used = 0;
        long int ps0Trial = 0;
        if( !haveAnswer )
            if( pa_p->OptimaPreSolveSeed > 0 )
        {
            // 1 = dual always, 2 = dual where water is at least a fifth of the system, 3 = dual and amounts
            const int mode = 3;
            const bool capOn = true;
            int gate = 1;
            if( mode == 2 )
            {
                double sumB = 0., bH = 0.;
                for( long int i = 0; i < N; i++ )
                {
                    sumB += pm.B[i];
                    std::string icName = char_array_to_string( pm.SB[i], MAXICNAME );
                    icName.erase( icName.find_last_not_of( " \t" ) + 1 );
                    if( icName == "H" ) bH = pm.B[i];
                }
                gate = ( sumB > 0. && 0.5 * bH / sumB >= 0.2 ) ? 1 : 0;
            }
            std::vector<double> yd;
            long int nVert = 0, nCapped = 0;
            if( mode > 0 && gate && LPGibbsDual( yd ) )
            {
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -yd[(size_t)i];
                ps0Used = 1;
                std::vector<double> nv;
                auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
                if( mode == 3 && TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nv, pm.G0 ) && (long int)nv.size() == L )
                {
                    ps0Used = 2;
                    long int j0 = 0;
                    for( long int k = 0; k < pm.FI; k++ )
                    {
                        const long int j1 = j0 + pm.L1[k];
                        double nPh = 0.;
                        for( long int j = j0; j < j1; j++ ) nPh += std::max( nv[(size_t)j], 0. );
                        for( long int j = j0; j < j1; j++ )
                        {
                            const long int sx = jToNx[(size_t)j];
                            if( sx < 0 ) continue;
                            double v = std::max( nv[(size_t)j], 0. );
                            if( v > 0. ) nVert++;
                            else if( k < pm.FIs && nPh > 0. )
                            {
                                double z = 0.;
                                for( long int i = 0; i < N; i++ ) z += yd[(size_t)i] * pm.A[ i + j*N ];
                                v = nPh * std::exp( -std::min( std::max( pm.G0[j] - z, 0. ), 700. ) );
                                if( capOn )
                                {
                                    const double cap = PhaseInsertionCeiling( j );
                                    if( cap >= 0. && v > cap ) { v = cap; nCapped++; }
                                }
                            }
                            state.x[sx] = std::min( std::max( v, problem.xlower[sx] ), problem.xupper[sx] );
                        }
                        j0 = j1;
                    }
                }
            }
            if( ps0Used )
                ps0Trial = 300;
            native_trace_decide( "ps0dual pass=%ld mode=%d gate=%d used=%d vertex=%ld capped=%ld trial=%ld",
                                 (long)pass, mode, gate, ps0Used, nVert, nCapped, ps0Trial );
        }


        // GEMS3K_OPTIMA_TRACE_FILE: Optima's per-iteration table for the reduced problem, as
        // for the full solve. The names are the active set's, so they are rebuilt each pass;
        // the file is appended to, one block per pass.
        if( const char* preTraceFn = std::getenv("GEMS3K_OPTIMA_TRACE_FILE") )
        {
            options.output.active = true;
            options.output.filename = preTraceFn;
            options.output.xnames.clear();
            options.output.ynames.clear();
            for( long int s = 0; s < nS; s++ )
                options.output.xnames.push_back(
                    char_array_to_string( pm.SM[ nxToJ[(size_t)s] ], MAXDCNAME ) );
            for( long int i = 0; i < N; i++ )
                options.output.ynames.push_back( char_array_to_string(pm.SB[i],3) );
        }

        Optima::Solver solver;
        solver.setOptions( options );
        preWatch->reset();          // per-pass stall state; the wall-clock budget is not reset
        if( ps0Trial > 0 && ps0Trial < passBudget )
        {
            Optima::Options trialOpts = options;
            trialOpts.maxiters = (unsigned)ps0Trial;
            solver.setOptions( trialOpts );
        }
        Optima::Result result = timed_solve( solver, problem, state );
        if( ps0Trial > 0 && ps0Trial < passBudget )
        {
            const bool trialOk = result.succeeded && !preWatch->stalled && !preWatch->timedOut;
            native_trace_decide( "ps0trial pass=%ld ok=%d iters=%ld trial=%ld",
                                 (long)pass, trialOk ? 1 : 0, (long)result.iterations, ps0Trial );
            const bool giveUp = !trialOk;
            if( giveUp )
            {
                ipm_logger->debug( "Optima pre-solve: the seeded first pass did not settle in {} steps; pre-solve dropped", ps0Trial );
                // The seeded start did not settle in time: drop the pre-solve for this call.
                g_ps0GaveUp = true;
                result.succeeded = false;
                native_trace_decide( "ps0giveup pass=%ld iters=%ld trial=%ld", (long)pass, (long)result.iterations, ps0Trial );
            }
            else if( !trialOk )
            {
                // The seeded start did not settle in time: start the pass again the ordinary way.
                iterationsOut += (long int)result.iterations;
                state = Optima::State( dims );
                for( long int s = 0; s < nS; s++ )
                {
                    const long int j = nxToJ[(size_t)s];
                    state.x[s] = std::min( std::max( pm.Y[j], problem.xlower[s] ), problem.xupper[s] );
                }
                Optima::Solver solver2;
                solver2.setOptions( options );
                preWatch->reset();
                result = timed_solve( solver2, problem, state );
            }
            else
                solver.setOptions( options );
        }
        iterationsOut += (long int)result.iterations;
        if( preWatch->stalled || preWatch->timedOut ) result.succeeded = false;

        if( !result.succeeded )
        {
            // GEMS3K_PRESOLVE_RESID_PROBE=<file>: on a failing pass, append the ranking of
            // the residual (same definition as the post-solve KKT check) of the state Optima
            // stopped at, to show which species the pass is stuck on. Zero cost when unset.
            if( const char* rp = std::getenv( "GEMS3K_PRESOLVE_RESID_PROBE" ) )
            {
                std::vector<double> Fsave( pm.F, pm.F + L ), Xsave( pm.X, pm.X + L );
                for( long int j = 0; j < L; j++ )
                    pm.X[j] = act[(size_t)j] ? state.x[ jToNx[(size_t)j] ] : xlo[(size_t)j];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                std::vector<std::pair<double,long int>> rk;
                double tot = 0.;
                for( long int j = 0; j < L; j++ )
                {
                    if( !act[(size_t)j] ) continue;
                    double g = pm.F[j];
                    for( long int i = 0; i < N; i++ ) g -= ( -state.ye[i] ) * pm.A[ i + j*N ];
                    const long int s2 = jToNx[(size_t)j];
                    double r;
                    if( problem.xupper[s2] <= problem.xlower[s2] * ( 1. + 1e-12 ) ) r = 0.;
                    else if( state.x[s2] <= problem.xlower[s2] * ( 1. + 1e-12 ) )   r = std::max( -g, 0. );
                    else if( state.x[s2] >= problem.xupper[s2] * ( 1. - 1e-12 ) )   r = std::max(  g, 0. );
                    else                                                            r = std::fabs( g );
                    rk.push_back( { r, j } );
                    tot += r;
                }
                std::sort( rk.begin(), rk.end(), []( const std::pair<double,long int>& a,
                                                     const std::pair<double,long int>& b )
                                                  { return a.first > b.first; } );
                FILE* fp = fopen( rp, "a" );
                if( fp )
                {
                    fprintf( fp, "# PASS %ld  nS=%ld L=%ld iters=%ld sum=%.8e top1frac=%.6f\n",
                             (long)pass, (long)nS, (long)L, (long)result.iterations, tot,
                             tot > 0. ? rk[0].first/tot : 0. );
                    for( size_t k = 0; k < rk.size() && k < 25; k++ )
                    {
                        const long int j = rk[k].second, s2 = jToNx[(size_t)j];
                        fprintf( fp, "%3zu %-22s resid=%.6e x=%.6e xlo=%.6e xhi=%.6e\n",
                                 k, char_array_to_string(pm.SM[j],MAXDCNAME).c_str(),
                                 rk[k].first, state.x[s2], problem.xlower[s2], problem.xupper[s2] );
                    }
                    fclose( fp );
                }
                for( long int j = 0; j < L; j++ ) { pm.F[j] = Fsave[(size_t)j]; pm.X[j] = Xsave[(size_t)j]; }
            }
            ipm_logger->debug( "OptimaReducedPreSolve: pass {} did not converge on {} of {} "
                               "species{} - discarding the reduced pre-solve", pass, nS, L,
                               preWatch->timedOut ? " (time or iteration budget)"
                                                  : ( preWatch->stalled ? " (stalled)" : "" ) );
            // DECIDE record per failed pass (the caller's `dimreduce` record carries only the
            // total). `tol` identifies which attempt the pass belongs to.
            native_trace_decide( "dimreducepass pass=%ld tol=%.6g ns=%ld of=%ld iters=%ld ok=0 "
                                 "readmit=-1 stop=%s budget=%ld",
                                 (long)pass, dimTol, (long)nS, (long)L,
                                 (long)result.iterations,
                                 preWatch->timedOut ? "timeout"
                                                    : ( preWatch->stalled ? "stall" : "nonconv" ),
                                 (long)passBudget );
            return discard();
        }

        // Scatter back to the full state and commit primal + dual.
        for( long int j = 0; j < L; j++ )
        {
            const double v = act[(size_t)j] ? state.x[ jToNx[(size_t)j] ] : xlo[(size_t)j];
            pm.X[j] = v;
            pm.Y[j] = v;
        }
        for( long int i = 0; i < N; i++ )
            pm.U[i] = -state.ye[i];
        haveAnswer = true;
        activeOut = nS;

        // Price the omitted columns on this dual.
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );
        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );

        // Readmit every omitted column whose reduced gradient is negative - the full
        // problem's KKT test for a variable at its lower bound (Optima's is_lower_unstable).
        // Readmitted species start at their lower bound, unless pa_OptimaReadmitSeed is set.
        long int readmitted = 0, seeded = 0;
        for( long int j = 0; j < L; j++ )
        {
            if( act[(size_t)j] ) continue;
            if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;  // fixed by its box
            double dual = 0.;
            for( long int i = 0; i < N; i++ )
                dual += pm.U[i] * pm.A[ i + j*N ];
            const double sj = pm.F[j] - dual;
            if( sj >= 0. ) continue;
            act[(size_t)j] = 1; readmitted++;

        }

        ipm_logger->debug( "OptimaReducedPreSolve: pass {} - {} of {} species active, "
                           "{} Optima iterations, {} readmitted ({} seeded above the floor)",
                           pass, nS, L, result.iterations, readmitted, seeded );
        // DECIDE record per converged pass; readmit=0 marks the fixed point that ends the loop.
        native_trace_decide( "dimreducepass pass=%ld tol=%.6g ns=%ld of=%ld iters=%ld ok=1 "
                             "readmit=%ld seeded=%ld budget=%ld",
                             (long)pass, dimTol, (long)nS, (long)L, (long)result.iterations,
                             (long)readmitted, (long)seeded, (long)passBudget );

        if( readmitted == 0 )
            return true;   // fixed point: the reduced answer satisfies the full KKT conditions
    }

    ipm_logger->debug( "OptimaReducedPreSolve: readmission did not settle within {} passes - "
                       "discarding (an unsettled active set is not a fixed point, and its dual is "
                       "wrong about everything still omitted)", maxPasses );
    return discard();
}

double TMultiBase::CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool referenceMode,
                                                    bool runKinetics )
{
    gems3k_timing::Scope timing( gems3k_timing::OptimaCall );
    phase_floor_parked.clear();   // set again below only for an accepted answer
    zeroed_floor_amount.clear();
    // Disable the IPM-2 chemical-potential smoothing for this call (see
    // optima_disable_smoothing). RAII, so the flag is restored on every exit.
    struct SmoothingGuard {
        TMultiBase* m;
        bool prev;
        explicit SmoothingGuard( TMultiBase* mm )
            : m(mm), prev(mm->optima_disable_smoothing)
{ m->optima_disable_smoothing = true; }
        ~SmoothingGuard() { m->optima_disable_smoothing = prev; }
    } smoothingGuard( this );

    PhaseHessianBegin();
    struct FloorRecordGuard {
        TMultiBase* m;
        ~FloorRecordGuard() { m->PhaseHessianEnd(); }
    } floorRecordGuard{ this };

#ifdef OPTIMA_LU_PARTIAL_PIVOT
    // pa_OptimaPartialPivLU: how many linear solves of this call used partial pivoting, and how many of them fell back.
    struct LuRecord {
        long p0 = Optima::trace::state().lupartial, f0 = Optima::trace::state().lufallbacks;
        ~LuRecord() { const long p = Optima::trace::state().lupartial - p0, f = Optima::trace::state().lufallbacks - f0;
                      if( p > 0 ) native_trace_decide( "lupartial n=%ld fallback=%ld", p, f ); }
    } luRecord;
#endif

#ifdef OPTIMA_BACKTRACK_SNAP_LEAVING
    // pa_OptimaSnapLeaving: how many species this call put on their floor, for the trace.
    struct SnapRecord {
        long n0 = Optima::trace::state().snaps;
        ~SnapRecord() { const long n = Optima::trace::state().snaps - n0; if( n > 0 ) native_trace_decide( "snapleaving n=%ld", n ); }
    } snapRecord;
#endif

    double ScFact = 1.;
    const BASE_PARAM* pa_p = base_param();
    // dcFloor, the species lower bound: pa_OptimaDcFloor if set, else pa_DHB.
    const double dcFloor = pa_p->OptimaDcFloor > 0. ? pa_p->OptimaDcFloor
                                                    : std::max( pa_p->DHB, 1e-300 );
    // pa_OptimaDryRetry: the start for the dry-state retry. Same LP as LPFeasibilitySeed(),
    // but water counts against the total, so H and O go into liquid water first.
    auto waterFirstSeed = [this]( std::vector<double>& nOut ) -> bool
    {
        if( !HasAqueousPhase() || pm.LO < 0 || pm.LO >= pm.L ) return false;
        const long int Nn = pm.N, Ll = pm.L;
        std::vector<double> cost( (size_t)Ll, 1. );
        cost[(size_t)pm.LO] = -1.;
        auto aFn = [this,Nn]( long int i, long int j ) { return pm.A[ i + j*Nn ]; };
        if( !TwoPhaseSimplexMinSum( Nn, Ll, aFn, pm.B, nOut, cost.data() ) || (long int)nOut.size() != Ll )
            return false;
        double bScale = 1., maxResid = 0., minVal = 0.;
        for( long int i = 0; i < Nn; i++ ) bScale = std::max( bScale, std::fabs( pm.B[i] ) );
        for( long int i = 0; i < Nn; i++ )
        {
            double sum = 0.;
            for( long int j = 0; j < Ll; j++ ) sum += pm.A[ i + j*Nn ] * nOut[(size_t)j];
            maxResid = std::max( maxResid, std::fabs( sum - pm.B[i] ) );
        }
        for( long int j = 0; j < Ll; j++ ) minVal = std::min( minVal, nOut[(size_t)j] );
        const double tol = std::max( 1e-6, 1e-8 * bScale );
        return maxResid <= tol && minVal >= -tol;
    };

    InitalizeGEM_IPM_Data();

    // Redundant species are held at zero for this call only; the guard restores their
    // metastability settings on every exit, exceptions included.
    const std::vector<RedundantDCHold> redundantHeld = ExcludeRedundantDCs();
    struct RedundantRestore {
        TMultiBase* self; const std::vector<RedundantDCHold>& held;
        ~RedundantRestore() { self->RestoreRedundantDCs( held ); }
    } redundantRestore{ this, redundantHeld };

    pm.t_start = clock();
    pm.t_end = pm.t_start;
    pm.t_elap_sec = 0.0;
    pm.ITF = pm.ITG = 0;
    pm.Ec = pm.MK = pm.PZ = 0;
    setErrorMessage( 0, "", "" );

    // One kinetics/metastability time step, at the same position native runs it: after
    // InitalizeGEM_IPM_Data() and ExcludeRedundantDCs(), before the internal rescaling, so
    // TKinMet sees the caller's real units. runKinetics is false only for the HOP Optima
    // leg, whose native leg has already advanced the step.
    if( runKinetics )
        RunKineticsStep();

    // pa_OptimaRedoxFree: if nothing in the bulk fixes the redox level, say so, and at 2 aim this call at pe 4.
    // The target added here is removed again when the call returns.
    bool redoxTargetAdded = false;
    if( pa_p->OptimaRedoxFree > 0 )
    {
        const RedoxFreeCheck rf = CheckRedoxFree();
        const bool apply = pa_p->OptimaRedoxFree >= 2 && rf.free && !referenceMode
                           && optima_control_conditions.empty();
        if( apply )
        {
            const double Eh = 4. * std::log( 10. ) * 0.000086 * pm.T;   // same Eh-to-dual factor as the Eh target
            SetControlCondition_Eh( Eh );
            redoxTargetAdded = true;
        }
        if( rf.checked )
            native_trace_decide( "redoxfree free=%d consistent=%d determined=%d excess=%.3e rel=%.3e action=%s",
                                 (int)rf.free, (int)rf.consistent, (int)rf.determined, rf.excess, rf.rel,
                                 apply ? "pe4" : ( rf.free ? ( optima_control_conditions.empty() ? "log" : "target_set" ) : "none" ) );
        if( rf.free )
            ipm_logger->info( "Nothing in the bulk fixes the redox level (electron excess {:.3e}, relative {:.3e}){}",
                              rf.excess, rf.rel, apply ? "; solving at pe 4" : "" );
    }
    struct RedoxTargetRemove {
        TMultiBase* m; bool added;
        ~RedoxTargetRemove()
        {
            if( !added ) return;
            auto& cc = m->optima_control_conditions;
            for( const auto& c : cc )
                if( c.name == "Eh" )
                    native_trace_decide( "redoxfree_done titrant=%.3e achieved=%.4g target=%.4g met=%d",
                                         c.titrantAmount, c.achievedValue, c.target, (int)c.targetMet );
            cc.erase( std::remove_if( cc.begin(), cc.end(), []( const EqControlCondition& c ){ return c.name == "Eh"; } ),
                      cc.end() );
        }
    } redoxTargetRemove{ this, redoxTargetAdded };

    if( pa_p->DG > 1e-5 )
    {
        ScFact = SystemTotalMolesIC();
        ScaleSystemToInternal( ScFact );
    }
    // Warn about elements with less material than their species' floor amounts can hold.
    // After rescaling, so pm.B, pm.DLL and dcFloor are in the same internal units.
    SubFloorElementCheck( dcFloor );

    // Scale of this node's previous solve, so a warm call can tell which species it left at the floor.
    const double prevScFactIn = prev_scfact;
    struct PrevScFactGuard { double& r; const double& v; ~PrevScFactGuard() { r = v; } } prevScFactGuard{ prev_scfact, ScFact };

    // Trace bookkeeping: one traced solve per call; each Optima solve inside it is an attempt.
    static thread_local std::vector<double> traceXPrev;
    static thread_local std::vector<char> tracePhPrev;
    static thread_local long traceCumIt = 0;
    if( Optima::trace::file() )
    {
        Optima::trace::beginSolve();
        Optima::trace::state().auto_attempt = true;
    }

    try
    {
        // The GUI loads the mixing model parameters here when they are not loaded yet, as GEM_IPM_Init()
        // does for the native path. Standalone GEMS3K reads them from the input files, and the call does nothing.
        {
            bool multiComponent = false;
            for( long int k = 0; k < pm.FI; k++ )
                if( pm.L1[k] > 1 )
                    multiComponent = true;
            if( pm.FIs && multiComponent )
                GEM_IPM_Init_gui1();
        }
        // Allocates/parametrises each multicomponent phase's TSolMod at the current T,P (the
        // native path does this in GEM_IPM_Init()).
        CalculateActivityCoefficients( LINK_TP_MODE );

        // pm.pNP is set by TNode::GEM_run(): 0 = cold start (AOP, like AIA), 1 = warm start
        // from the existing pm.Y[] (SOP, like SIA). No native refinement seeds or corrects
        // this solve. The cold start uses the LP-feasibility seed (both AOP and ROP).
        //
        // A warm start needs a previous solution on this node. pm.U[] (the dual) is non-zero
        // after any real solve by any solver; all zero (or non-finite) means nothing has been
        // solved here, and the project file's stored speciation is not a warm start. Then the
        // call warns and falls back to the cold seed. The dual seed below uses the same flag,
        // so the primal and dual always come from the same place.
        bool warmStateUsable = false;
        if( pm.pNP != 0 )
        {
            for( long int i = 0; i < pm.N; i++ )
            {
                if( !std::isfinite( pm.U[i] ) ) { warmStateUsable = false; break; }
                if( pm.U[i] != 0. ) warmStateUsable = true;
            }
            if( !warmStateUsable )
                ipm_logger->warn( "Warm start (SOP) requested, but this node has no earlier solution "
                                  "(the speciation stored in the .dbr file is not one). Solving cold. "
                                  "Try: use AOP for a first solve, and SOP only after a solve." );
        }

        if( pm.pNP == 0 || !warmStateUsable )
        {
            // Cold seed: the column-generation seed (pa_OptimaCgSeed, when armed) or the
            // LP-feasibility seed, min sum(n) s.t. A*n = b, n >= 0 (LPFeasibilitySeed()). An LP
            // vertex puts the mass in as few species as possible, so the result is passed once,
            // up front, through the solvent- and phase-collapse redistribution. If the LP fails,
            // a uniform tiny seed is used and the retries below remain the safety net.
            std::vector<double> lpSeed;
            const double cgSeedTol = pa_p->OptimaCgSeed;      // pa_OptimaCgSeed; 0 = off
            bool cgUsed = false;
            if( cgSeedTol > 0. && g_optimaCgSeedArmed && ColumnGenerationSeed( lpSeed, cgSeedTol ) )
                cgUsed = true;
            g_coldLvlUsed = 0;
            g_lvlShortcut = false;
            g_lvlCgStart = false;
            // The leveled starts serve the quick try and the first full attempt only; a retry starts the ordinary way.
            const bool lvlAllowed = g_optimaFixAccept || g_lvlFullCalls == 0;
            if( !g_optimaFixAccept ) g_lvlFullCalls++;
            if( !lvlAllowed && ( pa_p->OptimaLevelStart > 0 || pa_p->OptimaColdLevel > 0 ) )
                native_trace_decide( "lvlretry ordinary=1 call=%ld", g_lvlFullCalls );
            // pa_OptimaLevelStart: when the linear solve with mixed phases priced by their tangent plane
            // ends with pure phases only and no mixed phase wanting to form, that is the answer; start Optima
            // from it and its element potentials, so it only has to confirm it.
            bool lvlGap = false;
            if( lvlAllowed && !cgUsed && pa_p->OptimaLevelStart > 0 )
            {
                bool checkable = true;
                for( long int k = 0; k < pm.FIs && checkable; k++ )
                {
                    const char ph = pm.PHC[k];
                    if( pm.L1[k] > 1 && ( ph == PH_AQUEL || ph == PH_GASMIX || ph == PH_PLASMA || ph == PH_FLUID || ph == PH_SORPTION
                                          || ph == PH_POLYEL || ph == PH_ADSORPT || ph == PH_IONEX ) )
                        checkable = false;
                }
                std::vector<double> nCg;
                int why = checkable ? 0 : 1;
                g_tpdStrongSearch = true;
                const bool cgOk = checkable && ColumnGenerationSeed( nCg, 1e-6 );
                g_tpdStrongSearch = false;
                if( cgOk )
                {
                    long int inSoln = 0;
                    for( long int j = 0; j < pm.Ls && j < pm.L; j++ ) if( nCg[(size_t)j] > 0. ) inSoln++;
                    if( !g_cgConverged ) why = 3;
                    else if( g_cgPseudo > 0 || inSoln > 0 )
                    {
                        why = 4;
                        // Start each mixed phase at the composition the linear solve found for it.
                        // Two compositions of one phase far apart need two copies; start the ordinary way then.
                        if( g_cgGap )
                        {
                            lvlGap = true;
                            native_trace_decide( "lvlcgstart-gap pseudo=%ld", (long)g_cgPseudo );
                        }
                        else
                        {
                            lpSeed = nCg;
                            g_lvlCgStart = true;
                            native_trace_decide( "lvlcgstart pseudo=%ld insoln=%ld", (long)g_cgPseudo, (long)inSoln );
                        }
                    }
                    else
                    {
                        lpSeed = nCg;
                        g_lvlShortcut = true;
                    }
                }
                else if( checkable ) why = 2;
                native_trace_decide( "lvlshortcut certified=%d why=%d", g_lvlShortcut ? 1 : 0, why );
                ipm_logger->debug( "Optima leveled start: certified={} why={} composition-start={}", g_lvlShortcut, why, g_lvlCgStart );
            }
            if( g_lvlShortcut || g_lvlCgStart || cgUsed || LPFeasibilitySeed( lpSeed ) )
            {
                for( long int j = 0; j < pm.L; j++ )
                    pm.Y[j] = lpSeed[j];
                // pa_OptimaColdLevel: where no pre-solve runs, start from the linear-Gibbs vertex instead,
                // other members of a mixed phase at the phase amount times exp(-price), each capped at what its
                // elements allow, with the linear-Gibbs dual; copies of one solid solution start apart.
                // Not when the leveling found a gap: that start drops one side of it.
                if( lvlAllowed && !cgUsed && !g_lvlShortcut && !g_lvlCgStart && !lvlGap && pa_p->OptimaColdLevel > 0
                    && optima_dimreduce_passes( pa_p->OptimaDimReduce, pm.L ) <= 0 )
                {
                    const long int N = pm.N, L = pm.L;
                    std::vector<double> yd, nv;
                    auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
                    long int nVert = 0, nCapped = 0;
                    if( LPGibbsDual( yd ) && TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nv, pm.G0 ) && (long int)nv.size() == L )
                    {
                        long int j0 = 0;
                        for( long int k = 0; k < pm.FI; k++ )
                        {
                            const long int j1 = j0 + pm.L1[k];
                            double nPh = 0.;
                            for( long int j = j0; j < j1; j++ ) nPh += std::max( nv[(size_t)j], 0. );
                            for( long int j = j0; j < j1; j++ )
                            {
                                double v = std::max( nv[(size_t)j], 0. );
                                if( v > 0. ) nVert++;
                                else if( k < pm.FIs && nPh > 0. )
                                {
                                    double z = 0.;
                                    for( long int i = 0; i < N; i++ ) z += yd[(size_t)i] * pm.A[ i + j*N ];
                                    v = nPh * std::exp( -std::min( std::max( pm.G0[j] - z, 0. ), 700. ) );
                                    const double cap = PhaseInsertionCeiling( j );
                                    if( cap >= 0. && v > cap ) { v = cap; nCapped++; }
                                }
                                pm.Y[j] = v;
                            }
                            j0 = j1;
                        }
                        g_coldLvlUsed = 3;
                        // Start each pair of copies of one solid solution at two different
                        // compositions, so that a miscibility gap can open from the start.
                        if( g_coldLvlUsed >= 3 )
                        {
                            const std::vector<char> twin = PhaseHessianTwins();
                            std::vector<long int> st( (size_t)pm.FIs + 1, 0 );
                            for( long int k = 0; k < pm.FIs; k++ ) st[(size_t)k+1] = st[(size_t)k] + pm.L1[k];
                            auto sameCol = [&]( long int a, long int b ) {
                                for( long int r = 0; r < N; r++ ) if( pm.A[ r + a*N ] != pm.A[ r + b*N ] ) return false;
                                return true; };
                            long int nPairs = 0;
                            std::vector<char> used( (size_t)pm.FIs, 0 );
                            for( long int a = 0; a < pm.FIs && twin.size() == (size_t)pm.FIs; a++ )
                                for( long int b = a+1; b < pm.FIs; b++ )
                                {
                                    if( !twin[(size_t)a] || !twin[(size_t)b] || used[(size_t)a] || used[(size_t)b] || pm.L1[a] != pm.L1[b] || pm.L1[a] < 2 )
                                        continue;
                                    const long int n = pm.L1[a];
                                    std::vector<long int> mb( (size_t)n, -1 );
                                    bool ok = true;
                                    for( long int i = 0; i < n && ok; i++ )
                                    {
                                        for( long int j = 0; j < n; j++ )
                                            if( sameCol( st[(size_t)a]+i, st[(size_t)b]+j ) ) { mb[(size_t)i] = st[(size_t)b]+j; break; }
                                        ok = mb[(size_t)i] >= 0;
                                    }
                                    // Strict: copies of one phase also share the phase type and each member's standard energy.
                                    if( ok )
                                    {
                                        if( pm.PHC[a] != pm.PHC[b] ) ok = false;
                                        for( long int i = 0; i < n && ok; i++ )
                                        {
                                            const double g1 = pm.G0[st[(size_t)a]+i], g2 = pm.G0[mb[(size_t)i]];
                                            if( std::fabs( g1 - g2 ) > 1e-9 * ( 1. + std::fabs( g1 ) ) ) ok = false;
                                        }
                                        if( !ok ) native_trace_decide( "coldlvltwin-skip a=%ld b=%ld", (long)a, (long)b );
                                    }
                                    if( !ok ) continue;
                                    // Amount to share: both copies' seeded amounts, else half the matching pure phases' amounts.
                                    double T = 0.;
                                    std::vector<double> share( (size_t)n, 0. );
                                    for( long int i = 0; i < n; i++ )
                                    {
                                        share[(size_t)i] = pm.Y[st[(size_t)a]+i] + pm.Y[mb[(size_t)i]];
                                        T += share[(size_t)i];
                                    }
                                    if( !( T > 1e-12 ) )
                                    {
                                        T = 0.;
                                        for( long int i = 0; i < n; i++ )
                                        {
                                            share[(size_t)i] = 0.;
                                            for( long int j = st[(size_t)pm.FIs]; j < L; j++ )
                                                if( sameCol( j, st[(size_t)a]+i ) ) share[(size_t)i] += 0.5 * pm.Y[j];
                                            T += share[(size_t)i];
                                        }
                                    }
                                    if( !( T > 1e-12 ) ) continue;
                                    // Copy a leans to the end-member with the largest share, copy b to the next one.
                                    long int m1 = 0, m2 = -1;
                                    for( long int i = 1; i < n; i++ ) if( share[(size_t)i] > share[(size_t)m1] ) m1 = i;
                                    for( long int i = 0; i < n; i++ ) if( i != m1 && ( m2 < 0 || share[(size_t)i] > share[(size_t)m2] ) ) m2 = i;
                                    const double half = 0.5 * T, rest = ( n > 1 ) ? 0.1 / (double)( n - 1 ) : 0.;
                                    for( long int i = 0; i < n; i++ )
                                    {
                                        double va = half * ( i == m1 ? 0.9 : rest ), vb = half * ( i == m2 ? 0.9 : rest );
                                        const double ca = PhaseInsertionCeiling( st[(size_t)a]+i );
                                        if( ca >= 0. ) { va = std::min( va, ca ); vb = std::min( vb, ca ); }
                                        pm.Y[st[(size_t)a]+i] = va;
                                        pm.Y[mb[(size_t)i]] = vb;
                                    }
                                    used[(size_t)a] = used[(size_t)b] = 1;
                                    nPairs++;
                                    native_trace_decide( "coldlvltwin a=%ld b=%ld total=%.6g lean=%ld/%ld", (long)a, (long)b, T, (long)m1, (long)m2 );
                                }
                        }
                    }
                    native_trace_decide( "coldlvl used=%d vertex=%ld capped=%ld", g_coldLvlUsed, nVert, nCapped );
                    ipm_logger->debug( "Optima cold start from the linear-Gibbs solution: {} vertex species, {} amounts capped", nVert, nCapped );
                }
                if( HasAqueousPhase() && pm.LO >= 0 && pm.LO < pm.L )
                {
                    // pm.DUL[j] < 1e6 marks a real (possibly zero) upper restriction, as
                    // in ipm_chemical.cpp.
                    const double loUpper = ( pm.DUL[pm.LO] < 1e6 )
                                            ? pm.DUL[pm.LO] : std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 10.;
                    double waterSeed = 0.;
                    if( DetectSolventCollapseAndReseed( pm.Y, loUpper, waterSeed ) )
                        pm.Y[pm.LO] = waterSeed;
                }
                long int aqueousPhaseIdx = -1;
                {
                    long int j0seed = 0;
                    for( long int k = 0; k < pm.FIs; k++ )
                    {
                        if( pm.LO >= j0seed && pm.LO < j0seed + pm.L1[k] ) { aqueousPhaseIdx = k; break; }
                        j0seed += pm.L1[k];
                    }
                }
                std::vector<std::pair<long int,double>> seedReseeds;
                // A cold-level start already places each mixed phase; refilling it again only costs steps.
                const bool seedReseedFound = !g_lvlShortcut && !g_lvlCgStart && g_coldLvlUsed == 0 && DetectPhaseCollapseAndReseed( pm.Y, seedReseeds, aqueousPhaseIdx, true );
                if( seedReseedFound )
                    for( const auto& rs : seedReseeds )
                        pm.Y[ rs.first ] = rs.second;
            }
            else
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: the LP-feasibility seed is unavailable "
                                   "(the LP was infeasible or failed numerically, or its self-check "
                                   "rejected the result) - falling back to the uniform seed" );
                const double uniformSeed = std::max( dcFloor, 1e-16 * ScFact );
                for( long int j = 0; j < pm.L; j++ )
                    pm.Y[j] = uniformSeed;
            }
        }
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        for( long int j = 0; j < pm.L; j++ )
            pm.X[j] = pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );

        const long int N = pm.N;
        const long int L = pm.L;
        std::vector<EqControlCondition>& conditions = optima_control_conditions;
        const long int R = (long int)conditions.size();

        // ---- Species-level dimension-reduction pre-solve (pa_OptimaDimReduce) ----
        // Solves the equilibrium over a reduced set of species first and leaves its primal in
        // pm.Y[] and its dual in pm.U[]; the full-dimension solve below then warm-starts from
        // both, which makes it a cheap verification with every check at full dimension.
        // Skipped with control conditions and in ROP. Cold starts only (pm.pNP == 0): a warm
        // start already carries a consistent (primal, dual) pair, which a reduction would
        // overwrite. Exception: the HOP leg under an explicit positive setting (optima_hop_leg).
        // The pass count, including AUTO's size gate, comes from optima_dimreduce_passes().

        // Fallback initial-set rule, used once if the configured rule's pre-solve is discarded:
        // the threshold rule and the rank rule tend to fail on different projects. A discarded
        // pre-solve has restored pm.Y/pm.U, so the second attempt costs the wasted pass, which
        // can be a full pre-solve budget, max(2000, pa_IIM).
        static const double kDimReduceDefaultTol  =  10.;   // the shipped threshold
        static const double kDimReduceFallbackRank = -3.;   // 3 x N: the rank rule fallback

        const long int dimReducePasses = optima_dimreduce_passes( pa_p->OptimaDimReduce, L );

        // On the HOP leg the reduction runs only for an explicit pa_OptimaDimReduce > 0; AUTO
        // does not reach it (there the warm verification is usually already cheap).
        const bool hopReduce = optima_hop_leg && pa_p->OptimaDimReduce > 0;

        long int dimReduceIters = 0;
        bool dimReduceDone = false;
        // pa_OptimaEmptyPreSolveSkip: a retry of this call does not repeat a pre-solve that found nothing.
        const bool skipEmptyPreSolve = pa_p->OptimaEmptyPreSolveSkip && g_optimaPreSolveEmpty;
        if( skipEmptyPreSolve && !referenceMode && R == 0 && ( pm.pNP == 0 || hopReduce ) && dimReducePasses > 0 )
            native_trace_decide( "dimreduceskip reason=empty_earlier passes=%ld of=%ld", (long)dimReducePasses, (long)L );
        if( !skipEmptyPreSolve && !referenceMode && R == 0 && ( pm.pNP == 0 || hopReduce ) && dimReducePasses > 0 )
        {
            // Re-establish a consistent (Y, X, XF/XFA, activity coefficients) state after
            // each attempt: the pre-solve's objective changes these, and on the discard path
            // pm.Y[] is still the seed. Also run between the two attempts.
            auto reestablish = [&]() {
                TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                for( long int j = 0; j < L; j++ )
                    pm.X[j] = pm.Y[j];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
            };

            const double configuredTol = pa_p->OptimaDimReduceTol;
            long int nActive = 0;
            // First attempt at the lowered per-pass budget (pa_OptimaPreSolveFirstIters); the
            // fallback keeps max(2000, pa_IIM).
            g_ps0GaveUp = false;
            dimReduceDone = OptimaReducedPreSolve( dimReducePasses, dcFloor, configuredTol,
                                                   optima_presolve_pass_budget(
                                                       pa_p->OptimaPreSolveFirstIters,
                                                       (long int)pa_p->IIM, true ),
                                                   dimReduceIters, nActive );
            reestablish();

            // Iterations of a discarded attempt, tracked separately from the total.
            long int dimReduceWasted = 0;
            long int dimReduceAttempts = 1;

            if( !dimReduceDone && !g_ps0GaveUp )
            {
                // The configured rule was discarded: try the other one once. A negative
                // configured value (rank rule) falls back to the threshold rule; anything else
                // falls back to the rank rule.
                const double fallbackTol = ( configuredTol < 0. )
                                            ? kDimReduceDefaultTol : kDimReduceFallbackRank;
                long int fallbackIters = 0;
                ipm_logger->debug( "CalculateEquilibriumStateOptima: dimension-reduction pre-solve "
                                   "discarded at tol={} after {} iterations - retrying once at "
                                   "tol={}", configuredTol, dimReduceIters, fallbackTol );
                // DECIDE record of the discard: a discarded first attempt can be most of a
                // run's cost while giving the same answer.
                native_trace_decide( "dimreducediscard tol=%.6g iters=%ld active=%ld of=%ld "
                                     "fallbacktol=%.6g",
                                     configuredTol, (long)dimReduceIters, (long)nActive,
                                     (long)L, fallbackTol );
                dimReduceWasted = dimReduceIters;
                dimReduceAttempts = 2;
                dimReduceDone = OptimaReducedPreSolve( dimReducePasses, dcFloor, fallbackTol,
                                                       optima_presolve_pass_budget(
                                                           pa_p->OptimaPreSolveFirstIters,
                                                           (long int)pa_p->IIM, false ),
                                                       fallbackIters, nActive );
                dimReduceIters += fallbackIters;
                reestablish();
                native_trace_decide( "dimreduceretry tol=%.6g done=%d iters=%ld active=%ld of=%ld",
                                     fallbackTol, dimReduceDone ? 1 : 0, (long)fallbackIters,
                                     (long)nActive, (long)L );
            }

            if( dimReduceDone )
                ipm_logger->debug( "CalculateEquilibriumStateOptima: dimension-reduction pre-solve "
                                   "produced a warm start over {} of {} species in {} iterations",
                                   nActive, L, dimReduceIters );
            // DECIDE dimreduce: iters = wasted + productive, where wasted is a discarded
            // attempt's iterations (all of them if no attempt succeeded).
            if( !dimReduceDone )
            {
                dimReduceWasted = dimReduceIters;
                g_optimaPreSolveEmpty = true;
            }
            native_trace_decide( "dimreduce done=%d active=%ld of=%ld iters=%ld passes=%ld "
                                 "attempts=%ld wasted=%ld",
                                 dimReduceDone ? 1 : 0, (long)nActive, (long)L,
                                 (long)dimReduceIters, (long)dimReducePasses,
                                 (long)dimReduceAttempts, (long)dimReduceWasted );
        }

        // GEMS3K_PRESOLVE_HANDOVER_PROBE=<path>: dump, element by element, every piece of state
        // the main solve is about to start from (Y, X, lnGam, Gamma, F, F0, DUL, DLL, U, B,
        // XF, XFA, YF), so runs with and without the pre-solve can be diffed. Placed outside
        // the pre-solve block so both arms write at the same point.
        if( const char* hp = std::getenv( "GEMS3K_PRESOLVE_HANDOVER_PROBE" ) )
        {
            if( FILE* fh = fopen( hp, "a" ) )
            {
                // The header identifies the call (pNP, HOP leg, reference mode), since every
                // Optima call that reaches this point appends a block.
                fprintf( fh, "# HANDOVER L=%ld N=%ld FI=%ld pNP=%ld hop=%d reaktoro=%d "
                             "IT=%ld ITG=%ld ITF=%ld K2=%ld FitVar3=%.17g FitVar4=%.17g\n",
                         (long)L, (long)N, (long)pm.FI, (long)pm.pNP,
                         optima_hop_leg ? 1 : 0, referenceMode ? 1 : 0, (long)pm.IT,
                         (long)pm.ITG, (long)pm.ITF, (long)pm.K2,
                         pm.FitVar[3], pm.FitVar[4] );
                for( long int j = 0; j < L; j++ )
                    fprintf( fh, "DC %5ld %-22s Y=%.17g X=%.17g lnGam=%.17g Gamma=%.17g "
                                 "F=%.17g F0=%.17g DUL=%.17g DLL=%.17g\n",
                             (long)j, char_array_to_string(pm.SM[j],MAXDCNAME).c_str(),
                             pm.Y[j], pm.X[j], pm.lnGam[j], pm.Gamma[j],
                             pm.F[j], pm.F0[j], pm.DUL[j], pm.DLL[j] );
                for( long int i = 0; i < N; i++ )
                    fprintf( fh, "IC %5ld U=%.17g B=%.17g\n", (long)i, pm.U[i], pm.B[i] );
                for( long int k = 0; k < pm.FI; k++ )
                    fprintf( fh, "PH %5ld XF=%.17g XFA=%.17g YF=%.17g\n",
                             (long)k, pm.XF[k], pm.XFA[k], pm.YF[k] );
                fclose( fh );
            }
        }

        // Resolve each control condition's fixed objective gradient once (needs this call's
        // G0[] and T) and assign its unknown slot.
        std::vector<double> fixedGrad( R, 0. );
        for( long int k = 0; k < R; k++ )
        {
            conditions[k].slot = L + k;
            fixedGrad[k] = conditions[k].fixedGradientFn( conditions[k].target );
        }

        Optima::Dims dims;
        dims.x  = L + R;
        dims.be = N;
        // Sensitivity parameters c := the bulk composition b, so Sensitivity::xc is dn/db
        // (optima_want_sensitivity). With dims.c = 0 Optima skips the sensitivity solve.
        const bool wantSens = optima_want_sensitivity;
        if( wantSens ) dims.c = N;

        Optima::Problem problem( dims );

        for( long int i = 0; i < N; i++ )
            for( long int j = 0; j < L; j++ )
                problem.Aex(i, j) = pm.A[ i + j*N ];

        if( wantSens )
        {
            // be[i] is b[i], so d(be)/dc is the identity; c carries the current bulk composition.
            for( long int i = 0; i < N; i++ )
            {
                problem.c[i] = pm.B[i];
                for( long int k = 0; k < N; k++ )
                    problem.bec(i, k) = ( i == k ) ? 1. : 0.;
            }
        }
        for( long int k = 0; k < R; k++ )
        {
            for( long int i = 0; i < N; i++ )
                problem.Aex(i, L+k) = 0.;
            for( const auto& rc : conditions[k].stoich )
                problem.Aex( rc.first, L+k ) = rc.second;
        }
        for( long int i = 0; i < N; i++ )
            problem.be[i] = pm.B[i];

        // Box constraints: the species' kinetic/metastability bounds, floored by dcFloor.
        // pm.DUL[j] < 1e6 marks a real (possibly zero) upper restriction, as in
        // ipm_chemical.cpp; a DUL = DLL = 0 species gets the degenerate box [dcFloor, dcFloor].
        for( long int j = 0; j < L; j++ )
        {
            problem.xlower[j] = std::max( pm.DLL[j], dcFloor );
            problem.xupper[j] = ( pm.DUL[j] < 1e6 )
                                 ? std::max( pm.DUL[j], dcFloor ) : std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 10.;
        }
        // pa_OptimaLevelStart: a mixed phase the certified start leaves empty was shown not to form (its tangent
        // plane distance is positive), so its members stay at the floor; one by one their potentials are
        // undefined at zero amount and the solver would chase them.
        if( g_lvlShortcut && pm.pNP == 0 )
        {
            long int nPin = 0, j0 = 0;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                const long int j1 = j0 + pm.L1[k];
                bool empty = pm.L1[k] > 1;
                for( long int j = j0; j < j1 && empty; j++ ) if( pm.Y[j] > dcFloor ) empty = false;
                if( empty )
                    for( long int j = j0; j < j1; j++ ) { problem.xlower[j] = dcFloor; problem.xupper[j] = dcFloor; nPin++; }
                j0 = j1;
            }
            native_trace_decide( "lvlpin species=%ld", nPin );
        }
        // pa_OptimaLevelFinish: a leveled start with mixed phases is first finished by Newton on its own phases;
        // Optima then starts from that point and its potentials. A finish that fails leaves the start as it was.
        if( pa_p->OptimaLevelFinish > 0 && pm.pNP == 0 && g_lvlCgStart && g_cgDual.size() == (size_t)pm.N && R == 0 )
        {
            const std::vector<double> sx( pm.X, pm.X + pm.L ), sy( pm.Y, pm.Y + pm.L ), su( pm.U, pm.U + pm.N );
            for( long int j = 0; j < pm.L; j++ ) pm.X[j] = std::max( pm.Y[j], problem.xlower[j] );
            for( long int i = 0; i < pm.N; i++ ) pm.U[i] = g_cgDual[(size_t)i];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            std::vector<double> xlo( (size_t)pm.L ), xhi( (size_t)pm.L );
            for( long int j = 0; j < pm.L; j++ ) { xlo[(size_t)j] = problem.xlower[j]; xhi[(size_t)j] = problem.xupper[j]; }
            const bool ok = PotentialSpaceFinish( dcFloor, xlo, xhi );
            if( ok )
            {
                for( long int j = 0; j < pm.L; j++ ) pm.Y[j] = pm.X[j];
                for( long int i = 0; i < pm.N; i++ ) g_cgDual[(size_t)i] = pm.U[i];
            }
            else
            {
                for( long int j = 0; j < pm.L; j++ ) { pm.X[j] = sx[(size_t)j]; pm.Y[j] = sy[(size_t)j]; }
                for( long int i = 0; i < pm.N; i++ ) pm.U[i] = su[(size_t)i];
            }
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            native_trace_decide( "lvlfinish ok=%d", ok ? 1 : 0 );
            ipm_logger->debug( "Optima leveled start finished by Newton first: {}", ok ? "converged" : "did not converge" );
        }
        // Titrant unknowns are free in sign, bounded only generously as a numerical safety net.
        const double titrantBound = std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 2.0;
        for( long int k = 0; k < R; k++ )
        {
            problem.xlower[L+k] = -titrantBound;
            problem.xupper[L+k] =  titrantBound;
        }

        // Objective: Gibbs energy in RT units. fx = chemical potentials, recomputed with the
        // activity-coefficient models at every Optima iteration. Control-condition slots get
        // a fixed gradient equal to their target's chemical potential (see ipm_optima.h).
        // Pure single-species phases (pm.Ls..pm.L) get a logarithmic-barrier term,
        // -tau*ln(X_j), tau = pa_LogBarrierTau.
        const double kLogBarrierTau = pa_p->LogBarrierTau;
        // The objective, gradient and Hessian are the same for AOP/SOP and ROP; the modes
        // differ only in Optima::Options and in their retry chains.
        const double kPhaseHessianFloor = pa_p->PhaseHessianFloor;
        if( phfTwin.size() != (size_t)pm.FIs ) PhaseHessianBegin();
        const std::vector<char> phTwin = phfTwin;
        // Captured by value like the constants above - pa_p is not in scope inside the lambda.
        const bool kFDHessian = ( pa_p->OptimaFDHessian != 0 );
        const bool hasAq = HasAqueousPhase();
        // Last objective value seen by Optima, for the trace's ANSWER record.
        auto traceLastF = std::make_shared<double>( 0. );

        // Per-phase trend counters (last total, consecutive decreases, peak), kept here
        // because only the objective sees every iterate.
        auto phLast = std::make_shared<std::vector<double>>( pm.FIs, -1. );
        auto phDec  = std::make_shared<std::vector<long int>>( pm.FIs, 0 );
        auto phMax  = std::make_shared<std::vector<double>>( pm.FIs, 0. );
        // Phases already reported as getting the near-zero mixing block (pa_OptimaFloorCurvature 3).
        auto nzSeen = std::make_shared<std::vector<char>>( pm.FIs, 0 );
        // Exact twin copies for pa_OptimaFloorCurvature 3: phases with the same end-member formulas and the
        // same standard energies, i.e. one solid solution entered twice so it can unmix.
        std::vector<long int> phStart( (size_t)pm.FIs + 1, 0 );
        for( long int k = 0; k < pm.FIs; k++ ) phStart[(size_t)k+1] = phStart[(size_t)k] + pm.L1[k];
        std::vector<std::vector<long int>> twinCopies( (size_t)pm.FIs );
        if( pa_p->OptimaFloorCurvature == 3 )
        {
            const std::vector<char> sameFormulas = PhaseHessianTwins();
            for( long int a = 0; a < pm.FIs; a++ )
                for( long int b = a + 1; b < pm.FIs; b++ )
                {
                    if( !sameFormulas[(size_t)a] || !sameFormulas[(size_t)b] || pm.L1[a] != pm.L1[b] ) continue;
                    bool same = true;
                    for( long int m = 0; m < pm.L1[a] && same; m++ )
                    {
                        for( long int i = 0; i < pm.N; i++ )
                            if( pm.A[ i + (phStart[(size_t)a]+m)*pm.N ] != pm.A[ i + (phStart[(size_t)b]+m)*pm.N ] ) { same = false; break; }
                        if( std::fabs( pm.G0[phStart[(size_t)a]+m] - pm.G0[phStart[(size_t)b]+m] ) > 1e-6 ) same = false;
                    }
                    if( same ) { twinCopies[(size_t)a].push_back( b ); twinCopies[(size_t)b].push_back( a ); }
                }
        }
        // Tiered Hessian (pa_OptimaFDHessianDelay): set while the short cheap-Hessian attempt
        // runs, so the objective behaves as with pa_OptimaFDHessian = 0. An attempt-and-
        // restart rather than an in-place switch: switching the FD columns on mid-trajectory
        // does not recover a point the cheap iterations moved off course. shared_ptr because
        // the objective lambda outlives this scope inside Optima.
        auto fdSuppress = std::make_shared<bool>(false);
        // AOP/SOP only: ROP leaves Optima::Options at the library defaults (200 iterations per
        // attempt), so a delay would keep it from ever reaching the exact columns.
        // Cold start only: a warm call already carries a consistent (primal, dual) pair.
        const long int kFDDelay = ( !referenceMode && pm.pNP == 0 && pa_p->OptimaFDHessianDelay > 0 )
                                  ? pa_p->OptimaFDHessianDelay : 0;
        const bool kMoleFracHessian = ( pa_p->OptimaMoleFracHessian != 0 );
        // pa_OptimaFDPureSkip; not with electrostatic phases, which couple potentials across phases.
        bool fdPureSkip = ( pa_p->OptimaFDPureSkip > 0 );
        for( long int k = 0; k < pm.FI; k++ )
            if( pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL ) fdPureSkip = false;
        if( fdPureSkip ) native_trace_decide( "fdpureskip on=1" );
        const long int kFloorCurvature = pa_p->OptimaFloorCurvature;   // 0 off, 1 floor diagonal, 2 phase-sized step and full block, 3 = 1 + near-zero block
        problem.f = [this, L, R, dcFloor, &fixedGrad, fdPureSkip, kLogBarrierTau, kPhaseHessianFloor, phTwin, kFDHessian, kMoleFracHessian, kFloorCurvature, fdSuppress, hasAq, phLast, phDec, phMax, nzSeen, phStart, twinCopies, traceLastF]
                    ( Optima::ObjectiveResultRef res, Optima::VectorView x,
                      Optima::VectorView /*p*/, Optima::VectorView /*c*/,
                      Optima::ObjectiveOptions opts )
        {
            gems3k_timing::Scope timingObj( gems3k_timing::Obj );
            for( long int j = 0; j < L; j++ )
                pm.X[j] = x[j];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );

            double fval = 0.;
            for( long int j = 0; j < L; j++ )
            {
                res.fx[j] = pm.F[j];
                fval += pm.X[j] * pm.F[j];
            }
            for( long int j = pm.Ls; j < L; j++ )
            {
                const double Xj = std::max( pm.X[j], dcFloor );
                res.fx[j] -= kLogBarrierTau / Xj;
                fval -= kLogBarrierTau * std::log( Xj );
            }
            for( long int k = 0; k < R; k++ )
            {
                res.fx[L+k] = fixedGrad[k];
                fval += x[L+k] * fixedGrad[k];
            }
            res.f = fval;
            *traceLastF = fval;
            if( opts.eval.fxx )
            {
                gems3k_timing::Scope timingFxx( gems3k_timing::ObjFxx );
                res.fxx.setZero();
                // Ideal-mixing curvature per multicomponent phase (j < pm.Ls). Single-species
                // (pure) phases get no 1/X[j] term: their chemical potential does not depend on
                // their own amount, so their only curvature is the log-barrier term (below).
                //
                // Aqueous phase (the one containing pm.LO, found structurally): solute rows are
                // the molality Jacobian, d ln m_i/d n = delta_ij/n_j - delta_jw/n_w, which
                // couples each solute to the solvent amount; the solvent row uses the ideal
                // water-activity convention
                //     ln a_w = -(1 - x_w)/x_w = -(nSum - n_w)/n_w
                //  => d ln a_w/d n_i = -1/n_w (i != w), (nSum-n_w)/n_w^2 (i = w).
                // This is deliberately not d ln x_w/d n_i (that form converges worse).
                //
                // Non-aqueous solution phases: pa_p->OptimaMoleFracHessian selects the form.
                //   0 (default) - diag(1/X[j]) only.
                //   1           - the full ideal mole-fraction Jacobian
                //                 d ln x_j/d n_i = delta_ij/X[j] - 1/Xf, the exact derivative
                //                 of F for DC_SYMMETRIC species. The extra rank-1 term is zero
                //                 along the unmixing direction and affects only how freely a
                //                 phase's total amount moves.
                res.diagfxx = false;
                long int j0 = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int j1 = j0 + pm.L1[k];
                    const bool isAqueousLike =
                        ( hasAq && pm.LO >= j0 && pm.LO < j1 );
                    double Xf = 0.;
                    for( long int j = j0; j < j1; j++ )
                        Xf += std::max( pm.X[j], dcFloor );
                    // Consecutive-decrease counter per phase (the early-stability trend trigger).
                    {
                        double& prev = (*phLast)[k];
                        if( Xf > (*phMax)[k] ) (*phMax)[k] = Xf;
                        if( prev >= 0. )
                        {
                            if( Xf < prev * 0.999999 ) (*phDec)[k]++;
                            else if( Xf > prev * 1.000001 ) (*phDec)[k] = 0;
                        }
                        prev = Xf;
                    }
                    if( isAqueousLike )
                    {
                        const long int w = pm.LO;
                        const double Xw = std::max( pm.X[w], dcFloor );
                        for( long int j = j0; j < j1; j++ )
                        {
                            if( j == w ) continue;
                            const double Xj = std::max( pm.X[j], dcFloor );
                            res.fxx(j,j) = 1.0 / Xj;
                            res.fxx(j,w) = -1.0 / Xw;
                            res.fxx(w,j) = -1.0 / Xw;
                        }
                        res.fxx(w,w) = ( Xf - Xw ) / ( Xw * Xw );
                    }
                    else
                    {
                        // The rank-1 -(1/Xf) coupling applies only over the present
                        // end-members (same criterion as the pa_PhaseHessianFloor block). An
                        // absent end-member keeps diag(1/X_j); otherwise a phase driven to the
                        // floor would be left with a singular, unregularised block and could
                        // not recover.
                        double phTot = 0.;
                        for( long int j = j0; j < j1; j++ ) phTot += std::max( pm.X[j], 0. );
                        const double presThr = std::max( dcFloor * 1e3, phTot * 1e-6 );
                        // pa_OptimaFloorCurvature 3: the full block also on a nearly vanished gas or fluid phase, or
                        // on a vanishing copy of a twin pair whose other copy is present, so it settles in a few steps
                        // instead of shrinking by a small fraction per step. So small a phase cannot be moved far by
                        // one step. Other solids keep the diagonal: a solid that should form again can be held near
                        // zero by the full block.
                        bool twinHeld = false;
                        if( kFloorCurvature == 3 && (*phMax)[k] > 1e-6 )
                            for( long int b : twinCopies[(size_t)k] )
                            {
                                double tb = 0.;
                                for( long int j = phStart[(size_t)b]; j < phStart[(size_t)b+1]; j++ ) tb += std::max( pm.X[j], 0. );
                                if( tb > 1e-6 ) { twinHeld = true; break; }
                            }
                        const bool nearZero = kFloorCurvature == 3 && Xf < 1e7 * kLogBarrierTau
                                              && ( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_FLUID || twinHeld );
                        if( nearZero && !(*nzSeen)[k] )
                        {
                            (*nzSeen)[k] = 1;
                            native_trace_decide( "nearzeroblock phase=%s amount=%.3e peak=%.3e",
                                                 char_array_to_string( pm.SF[k], MAXPHNAME ).c_str(), Xf, (*phMax)[k] );
                            ipm_logger->info( "Optima: full mixing curvature for the nearly vanished phase {}",
                                               char_array_to_string( pm.SF[k], MAXPHNAME ) );
                        }
                        for( long int j = j0; j < j1; j++ )
                        {
                            const double Xj = std::max( pm.X[j], dcFloor );
                            // pa_OptimaFloorCurvature 2: the ideal-mixing block over all members.
                            if( ( kFloorCurvature == 2 || nearZero ) && j1 - j0 > 1 )
                            {
                                for( long int i = j0; i < j1; i++ ) res.fxx(j,i) = -1.0 / Xf;
                                res.fxx(j,j) = 1.0 / Xj - 1.0 / Xf;
                            }
                            else if( kMoleFracHessian && pm.X[j] > presThr )
                            {
                                for( long int i = j0; i < j1; i++ )
                                    if( pm.X[i] > presThr ) res.fxx(j,i) = -1.0 / Xf;
                                res.fxx(j,j) = 1.0 / Xj - 1.0 / Xf;
                            }
                            else
                                res.fxx(j,j) = 1.0 / Xj;
                        }
                    }
                    j0 = j1;
                }
                for( long int j = pm.Ls; j < L; j++ )
                {
                    const double Xj = std::max( pm.X[j], dcFloor );
                    res.fxx(j,j) += kLogBarrierTau / ( Xj * Xj );
                }
                // Control-condition slots get no direct curvature; they couple to the system
                // through the mass-balance rows.

                // PartiallyExact Hessian: starting from the per-phase approximate block above
                // (ideal mixing plus the log-barrier term for pure phases), the full column of
                // each Optima basic variable j (opts.ibasicvars) is replaced by the derivative
                // d(chem.potential_i)/dX_j, here by forward finite differences. Applied in all
                // Optima modes: the ideal-mixing block is convex in composition, so without
                // these columns the solver's model has no miscibility gap and an unmixed state
                // slides back to a homogeneous one.
                {
                    // pa_OptimaFDHessian = 0 skips this loop; pa_OptimaFDHessianDelay
                    // suppresses it during the cheap-Hessian attempt (fdSuppress).
                    const bool skipFD = !kFDHessian || *fdSuppress;
                    std::vector<double> Fbase( L );
                    for( long int j = 0; j < L; j++ ) Fbase[j] = pm.F[j];
                    const bool traceDump = optima_trace_hessian_dump_now();
                    std::vector<double> traceIdeal;
                    if( traceDump )
                        for( long int r = 0; r < L; r++ )
                            for( long int c = 0; c < L; c++ ) traceIdeal.push_back( res.fxx(r,c) );
                    for( Optima::Index bk = 0; !skipFD && bk < opts.ibasicvars.size(); bk++ )
                    {
                        const long int i = (long int)opts.ibasicvars[bk];
                        if( i >= L ) continue; // a control-condition virtual slot, not a species
                        const double Xi = pm.X[i];
                        gems3k_timing::Scope timingCol( i < pm.Ls ? gems3k_timing::FdMix : gems3k_timing::FdPure );
                        // pa_OptimaFloorCurvature 2: for a species of a solution phase the step is also at least
                        // 1e-7 of the phase amount, so the other potentials change well above round-off.
                        double hPhase = 0.;
                        if( kFloorCurvature == 2 && i < pm.Ls )
                        {
                            long int jj = 0;
                            for( long int kk = 0; kk < pm.FIs; jj += pm.L1[kk], kk++ )
                                if( i >= jj && i < jj + pm.L1[kk] ) { hPhase = pm.XF[kk] * 1e-7; break; }
                        }
                        const double h = std::max( { std::fabs(Xi) * 1e-7, dcFloor * 10., hPhase } );
                        if( fdPureSkip && i >= pm.Ls )
                        {   // a pure phase's amount changes no potential: its column is zero
                            for( long int j = 0; j < L; j++ ) res.fxx(j,i) = 0.;
                            if( kFloorCurvature == 2 )
                            {
                                const double Xb = std::max( Xi, dcFloor );
                                res.fxx(i,i) += kLogBarrierTau / ( Xb * Xb );
                            }
                            continue;
                        }
                        pm.X[i] = Xi + h;
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                        for( long int j = 0; j < L; j++ )
                            res.fxx(j,i) = ( pm.F[j] - Fbase[j] ) / h;
                        // The difference over h gives ln(1 + h/X)/h for the ln X term, too small when h is
                        // not small next to X (near the floor); put the exact 1/X back (pa_OptimaFloorCurvature).
                        if( kFloorCurvature && i < pm.Ls && Xi > 0. )
                            res.fxx(i,i) += 1.0 / Xi - std::log1p( h / Xi ) / h;
                        // Value 2 also keeps the log-barrier curvature of a pure phase.
                        if( kFloorCurvature == 2 && i >= pm.Ls )
                        {
                            const double Xb = std::max( Xi, dcFloor );
                            res.fxx(i,i) += kLogBarrierTau / ( Xb * Xb );
                        }
                        pm.X[i] = Xi;
                    }
                    // Exact, regularised curvature for non-aqueous multicomponent phases
                    // (pa_PhaseHessianFloor). Optima's basic set has only pm.N variables, so
                    // some phases get no exact columns from the loop above and keep the
                    // ideal-mixing 1/X_j, which is wrong near a critical point. The block of
                    // present end-members is finite-differenced, symmetrised, and its
                    // eigenvalues are floored at a fraction of its own largest (an exact
                    // block is indefinite inside a miscibility gap, and an unmodified Newton
                    // step on it would be unbounded).
                    //
                    // The aqueous phase is excluded: its analytic ideal-molality block has the
                    // right shape, and FD columns of its near-floor trace species are noise.
                    const double regRatio = kPhaseHessianFloor;
                    if( regRatio > 0. )
                    {
                        long int p0 = 0;
                        gems3k_timing::Scope timingBlock( gems3k_timing::FdBlock );
                        for( long int k = 0; k < pm.FIs; k++ )
                        {
                            const long int p1 = p0 + pm.L1[k];
                            const long int nEnd = p1 - p0;
                            const bool isAq = ( pm.LO >= p0 && pm.LO < p1 );
                            if( !isAq && nEnd > 1 && p1 <= L )
                            {
                                // Only end-members that are actually present: a floor
                                // end-member's 1/X_j would set the block's largest eigenvalue
                                // (and so the floor), and its FD column is a secant over a
                                // step larger than the species itself.
                                double phTot = 0.;
                                for( long int j = p0; j < p1; j++ ) phTot += std::max( pm.X[j], 0. );
                                std::vector<long int> pres;
                                for( long int j = p0; j < p1; j++ )
                                    if( pm.X[j] > std::max( dcFloor * 1e3, phTot * 1e-6 ) )
                                        pres.push_back( j );
                                const int nP = (int)pres.size();
                                if( nP > 1 )
                                {
                                    for( int c = 0; c < nP; c++ )
                                    {
                                        const long int i = pres[(size_t)c];
                                        const double Xi = pm.X[i];
                                        const double h = std::fabs(Xi) * 1e-7;
                                        if( !( h > 0. ) ) continue;
                                        pm.X[i] = Xi + h;
                                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                                        CalculateActivityCoefficients( LINK_UX_MODE );
                                        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                                        for( int r = 0; r < nP; r++ )
                                            res.fxx(pres[(size_t)r],i) = ( pm.F[pres[(size_t)r]] - Fbase[pres[(size_t)r]] ) / h;
                                        pm.X[i] = Xi;
                                    }
                                    std::vector<double> blk( (size_t)nP*nP );
                                    for( int a = 0; a < nP; a++ )
                                        for( int b = 0; b < nP; b++ )
                                            blk[(size_t)a*nP+b] = 0.5 * ( res.fxx(pres[(size_t)a],pres[(size_t)b])
                                                                        + res.fxx(pres[(size_t)b],pres[(size_t)a]) );
                                    bool indef = false;
                                    double presTot = 0.;
                                    for( int a = 0; a < nP; a++ ) presTot += pm.X[pres[(size_t)a]];
                                    const double absFloor = ( presTot > 0. ) ? regRatio * nP / presTot : 0.;
                                    const bool okEig = SymEigFloorInPlace( blk, nP, regRatio, &indef, absFloor );
                                    if( okEig )
                                    {
                                        if( indef && phTwin[(size_t)k] ) phfFired[(size_t)k]++;
                                        for( int a = 0; a < nP; a++ )
                                            for( int b = 0; b < nP; b++ )
                                                res.fxx(pres[(size_t)a],pres[(size_t)b]) = blk[(size_t)a*nP+b];
                                    }
                                }
                            }
                            p0 = p1;
                        }
                    }

                    // Restore the base state: XF/XFA, activity coefficients and pm.F were left
                    // at the last perturbation's values.
                    {
                        gems3k_timing::Scope timingRestore( gems3k_timing::ObjRestore );
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                    }
                    for( long int j = 0; j < L; j++ )
                        pm.F[j] = Fbase[j];

                    if( auto* tfp = Optima::trace::file() )
                    {
                        // Comparison trace: which variables got exact columns, and the matrix used.
                        auto& ts = Optima::trace::state();
                        std::vector<long> cur;
                        long nq = 0;
                        std::string names;
                        for( Optima::Index bk = 0; bk < opts.ibasicvars.size(); bk++ )
                        {
                            const long int i = (long int)opts.ibasicvars[bk];
                            cur.push_back( i );
                            if( i >= L ) { nq++; continue; }
                            names += ( names.empty() ? "" : "," ) + char_array_to_string( pm.SM[i], MAXDCNAME );
                        }
                        auto a = cur, b = ts.basis;
                        std::sort( a.begin(), a.end() );
                        std::sort( b.begin(), b.end() );
                        if( ts.basis.empty() || a != b )
                        {
                            fprintf( tfp, "HESSBASIS solve=%ld attempt=%d it=%ld nbasic=%ld nq_skipped=%ld names=%s\n",
                                     ts.solve, ts.attempt, ts.iter, (long)cur.size(), nq, names.c_str() );
                            ts.basis = cur;
                        }
                        const long nexact = skipFD ? 0 : (long)cur.size() - nq;
                        ts.fxevals += nexact;
                        double dmin = 1e300, dmax = -1e300, taumax = 0.;
                        for( long int j = 0; j < L; j++ ) { dmin = std::min( dmin, (double)res.fxx(j,j) ); dmax = std::max( dmax, (double)res.fxx(j,j) ); }
                        long nnz = 0;
                        if( !skipFD )
                            for( long i : cur )
                                if( i < L )
                                    for( long int r = 0; r < L; r++ )
                                        if( r != i && res.fxx(r,i) != 0. ) nnz++;
                        for( long int j = pm.Ls; j < L; j++ )
                        {
                            const double Xj = std::max( pm.X[j], dcFloor );
                            taumax = std::max( taumax, kLogBarrierTau / ( Xj * Xj ) );
                        }
                        fprintf( tfp, "HESS solve=%ld attempt=%d it=%ld hessian=%s nexact_cols=%ld nx=%ld nn_approx_cols=%ld diag_min=%.6e diag_max=%.6e"
                                      " offdiag_nnz_exact_cols=%ld tau_term_max=%.6e fx_evals=%ld src=pa_OptimaFDHessian\n",
                                 ts.solve, ts.attempt, ts.iter, skipFD ? "Approx" : "PartiallyExact", nexact, (long)res.fxx.rows(),
                                 (long)L - nexact, dmin, dmax, nnz, taumax, ts.fxevals );
                        fflush( tfp );
                        if( traceDump )
                        {
                            optima_trace_write_matrix( "used", res.fxx, L );
                            std::vector<double> ex( (size_t)L * L );
                            for( long int i = 0; i < L; i++ )
                            {
                                const double Xi = pm.X[i];
                                const double h = std::max( std::fabs(Xi) * 1e-7, dcFloor * 10. );
                                pm.X[i] = Xi + h;
                                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                                CalculateActivityCoefficients( LINK_UX_MODE );
                                PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                                for( long int r = 0; r < L; r++ ) ex[(size_t)r*L+i] = ( pm.F[r] - Fbase[r] ) / h;
                                pm.X[i] = Xi;
                            }
                            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                            CalculateActivityCoefficients( LINK_UX_MODE );
                            for( long int j = 0; j < L; j++ ) pm.F[j] = Fbase[j];
                            optima_trace_write_matrix( "exact", [&]( long r, long c ) { return ex[(size_t)r*L+c]; }, L );
                            optima_trace_write_matrix( "ideal", [&]( long r, long c ) { return traceIdeal[(size_t)r*L+c]; }, L );
                            std::ofstream f( std::string( optima_trace_env( "FILE" ) ) + ".hess." + std::to_string(ts.solve)
                                             + "." + std::to_string(ts.iter) + ".state.csv" );
                            f.precision( 17 );
                            f << "index,name,n,basic\n";
                            for( long int i = 0; i < L; i++ )
                                f << i << "," << char_array_to_string( pm.SM[i], MAXDCNAME ) << "," << pm.X[i] << ","
                                  << ( std::find( cur.begin(), cur.end(), i ) != cur.end() ? 1 : 0 ) << "\n";
                        }
                    }
                }
            }
            res.succeeded = true;
        };
        memoize_objective_if_linesearch( problem, pa_p->OptimaLineSearch );

        Optima::State state( dims );
        for( long int j = 0; j < L; j++ )
            state.x[j] = std::max( pm.Y[j], problem.xlower[j] );
        for( long int k = 0; k < R; k++ )
            state.x[L+k] = 0.0;
        // On a warm call, a species the previous answer left at the floor starts exactly at the
        // floor; the change of scale between calls would otherwise lift it slightly.
        if( pm.pNP != 0 && warmStateUsable && prevScFactIn > 0. )
            for( long int j = 0; j < L; j++ )
                if( pm.Y[j] * prevScFactIn / ScFact <= 10. * dcFloor && problem.xlower[j] <= dcFloor )
                    state.x[j] = problem.xlower[j];
        // pa_OptimaDryOut: a warm start whose previous answer had no water starts already dried out
        const bool dryWarmStart = pm.pNP != 0 && warmStateUsable && pm.LO >= 0 && pm.LO < L
                                  && state.x[pm.LO] <= problem.xlower[pm.LO];

        // A warm start must seed the dual too, not just the primal: Optima is primal-dual, and
        // with y = 0 a correct x still shows the full gradient as its stationarity error, so
        // the solver rebuilds the dual from zero while the primal waits. pm.U[] holds the
        // previous dual (written by native's IPM and by this function); Optima's ye has the
        // opposite sign. Used only on a warm call with a usable state, or after the reduced
        // pre-solve (which writes both pm.Y[] and pm.U[]). On a cold call pm.U[] is left over
        // from an unrelated calculation and is not used.
        if( ( pm.pNP != 0 && warmStateUsable ) || dimReduceDone )
        {
            bool duals_usable = true;
            for( long int i = 0; i < N; i++ )
                if( !std::isfinite( pm.U[i] ) ) { duals_usable = false; break; }
            // All or nothing: a dual right for some ICs and zero for others is inconsistent
            // and costs more than a zero dual.
            if( duals_usable )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -pm.U[i];
        }
        else if( pm.pNP == 0 )
        {
            // Cold start: take the multipliers from the linear-programming dual, unless water is a
            // minor share of the system, where that dual is not a good guide.
            double sumB = 0., bH = 0.;
            for( long int i = 0; i < N; i++ )
            {
                sumB += pm.B[i];
                std::string icName = char_array_to_string( pm.SB[i], MAXICNAME );
                icName.erase( icName.find_last_not_of( " \t" ) + 1 );
                if( icName == "H" ) bH = pm.B[i];
            }
            std::vector<double> ydual;
            if( ( g_lvlShortcut || g_lvlCgStart ) && g_cgDual.size() == (size_t)N )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -g_cgDual[(size_t)i];
            else
            if( sumB > 0. && ( 0.5 * bH / sumB >= 0.2 || g_coldLvlUsed >= 2 ) && LPGibbsDual( ydual ) )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -ydual[(size_t)i];
        }

        Optima::Options options;
        if( !referenceMode )
        {
            // IIM is the native IPM loop's iteration cap and OptimaTol the convergence
            // tolerance (pa_DK is a different quantity). maxiters is floored at 2000: raising
            // a cap can only turn a budget-exhausted failure into a converged result.
            options.maxiters = (unsigned)std::max( 2000L, (long int)pa_p->IIM );
            // pa_OptimaEarlyStabilityAt's cap is enforced in the convergence hook, not as a
            // maxiters clamp, so maxiters keeps the full budget.
            options.convergence.tolerance = pa_p->OptimaTol;
            apply_optima_linesearch( options, pa_p->OptimaLineSearch, pa_p->OptimaLSStallEscape,
                             pa_p->OptimaLSRejectWorse );
            apply_optima_lu( options, pa_p->OptimaPartialPivLU );
            // pa_OptimaLSWarmTrials: a warm call starts near the answer, where the full line search costs many
            // evaluations per step for little gain; cap its trials there.
            if( pm.pNP != 0 && pa_p->OptimaLSWarmTrials > 0 && options.linesearch.enabled )
            {
                options.linesearch.maxiterations = (unsigned)pa_p->OptimaLSWarmTrials;
                native_trace_decide( "lswarmcap trials=%ld", (long)pa_p->OptimaLSWarmTrials );
            }
            // Plain warm calls only: on a cold start, or on a leg that starts from native's answer, a species
            // far from its answer can limit the step many times on its way to a value that is not zero.
            if( !referenceMode && pm.pNP != 0 && !optima_hop_leg )
                apply_optima_snap( options, pa_p->OptimaSnapLeaving );
            // The quick first try of a cold call (pa_OptimaFastFirst): full step, each amount clipped to its bounds.
            // Not from a certified leveled start (pa_OptimaLevelStart): there the ordinary step only has to confirm it.
            if( g_optimaFixAccept && !g_lvlShortcut )
                options.backtracksearch.apply_min_max_fix_and_accept = true;
            // A certified leveled start makes the quick try pointless: one step, then the full solve confirms it.
            if( g_optimaFixAccept && g_lvlShortcut )
                options.maxiters = 1;
        }
        // else (referenceMode): Optima::Options stay at the library defaults; GEMS3K's own
        // tolerance and iteration overrides are not applied in the reference mode.

        // GEMS3K_OPTIMA_TRACE_FILE=<path>: Optima's own per-iteration output table, off by
        // default (zero cost when unset).
        if( const char* fn = std::getenv("GEMS3K_OPTIMA_TRACE_FILE") )
        {
            options.output.active = true;
            options.output.filename = fn;
            for( long int j = 0; j < L; j++ )
                options.output.xnames.push_back( char_array_to_string(pm.SM[j],MAXDCNAME) );
            for( long int k = 0; k < R; k++ )
                options.output.xnames.push_back( conditions[k].name );
            for( long int i = 0; i < N; i++ )
                options.output.ynames.push_back( char_array_to_string(pm.SB[i],3) );
        }

        // ---- Stall/freeze limit (pa_OptimaStallWindow) and wall-clock budget ----
        // Turns a frozen iterate (Optima error bit-identical from iteration 0) into an early
        // failure, so the retry tiers can start. Implemented through
        // Optima::ConvergenceOptions::check, a caller-supplied hook. Two points:
        //  1. `check` returning true means converged, and every retry tier is gated on
        //     !result.succeeded, so a stall is folded back to succeeded = false.
        //  2. The test is on the best-so-far error (with complementarity), not on the
        //     objective or the iterate displacement, which do not discriminate.
        // Reset before every solve() call: the retries reuse this `options` object.
        struct StallWatch {
            long int window = 0;
            double bestErr = 0.;    // best-so-far ||ex||inf (Optima's own criterion)
            double bestComp = 0.;   // best-so-far max_j |ex_j| * x_j (complementarity)
            // Best-so-far one whole window ago (the test is cumulative over a window).
            double refErr = 0., refComp = 0.;
            long int run = 0;       // iterations elapsed in the current window
            bool stalled = false;
            // Wall-clock budget (pa_OptimaMaxSeconds), for the whole call: not cleared by reset().
            double maxSeconds = 0.;
            std::chrono::steady_clock::time_point started;
            bool timedOut = false;
            // pa_OptimaDryRetry: water index, the most water the bulk could form, dry iterations so far.
            long int dryLO = -1;
            double dryCap = 0.;
            long int dryRun = 0;
            bool dryStop = false;
            // pa_OptimaFrozenCut: steps in a row with energy and error unchanged, and their last values.
            long int fzRun = 0;
            double fzErr = -1., fzF = 0.;
            // pa_OptimaTailStop: energy at the start of the current flat run, and its length.
            double tailRef = 0.; long int tailRun = 0; bool tailStop = false;
            std::vector<long int> tailPhases;   // pa_OptimaTailDryOnly 2: the vanished phases that held the error
            double tailErrPrev = -1., tailRatioMin = 1e300, tailRatioMax = 0.;   // step-to-step error ratio over the flat run
            bool vanishedOff = false;   // the scoped tail stop has been used in this call (not cleared by reset())
            // pa_OptimaDryOut: the solve stopped because the water ran out, and where.
            bool dryOut = false; long int dryOutAt = 0; double dryOutW = 0.;
            long int wLO = -1;   // pa_OptimaDryOut: the water index, -1 = off (not cleared by reset())
            bool dryDone = false;   // pa_OptimaDryOut: the water phase was taken out; do not look again
            bool wSeenWet = false;  // pa_OptimaDryOut: water was above the threshold at some step of this call
            long int stallAt = 0;   // pa_OptimaStallNetFactor: iterations of the first solve that the watch stopped (not cleared by reset())
            void reset() {
                bestErr = bestComp = std::numeric_limits<double>::infinity();
                refErr  = refComp  = std::numeric_limits<double>::infinity();
                run = 0; stalled = false; dryRun = 0; dryStop = false;
                fzRun = 0; fzErr = -1.; fzF = 0.;
                tailRef = 0.; tailRun = 0; tailStop = false; tailPhases.clear(); tailErrPrev = -1.; tailRatioMin = 1e300; tailRatioMax = 0.;
                dryOut = false; dryOutAt = 0; dryOutW = 0.;
            }
            bool overBudget() const {
                if( g_optimaIterBudget == 0 ) return true;   // the call's iteration budget is used up
                if( maxSeconds <= 0. ) return false;
                return std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - started ).count() > maxSeconds;
            }
        };
        auto stallWatch = std::make_shared<StallWatch>();
        stallWatch->window = pa_p->OptimaStallWindow;
        stallWatch->maxSeconds = pa_p->OptimaMaxSeconds;
        stallWatch->started = std::chrono::steady_clock::now();
        if( pa_p->OptimaDryOut > 0 && HasAqueousPhase() )   // pa_OptimaDryOut
        {
            // Water cannot vanish if an element in the bulk has no carrier outside the water phase.
            long int ja0 = 0, ja1 = 0;
            for( long int k = 0, jb = 0; k < pm.FIs; jb += pm.L1[k], k++ )
                if( pm.LO >= jb && pm.LO < jb + pm.L1[k] ) { ja0 = jb; ja1 = jb + pm.L1[k]; break; }
            std::string onlyAq;
            for( long int i = 0; i < pm.N; i++ )
            {
                if( !( pm.B[i] > 0. ) ) continue;
                if( pm.ICC && ( pm.ICC[i] == IC_CHARGE || pm.ICC[i] == IC_VOLUME ) ) continue;
                bool other = false;
                for( long int j = 0; j < pm.L && !other; j++ )
                    if( ( j < ja0 || j >= ja1 ) && pm.A[i + j*pm.N] != 0.
                        && !( pm.DUL && pm.DUL[j] < 1e6 && !( pm.DUL[j] > 0. ) ) ) other = true;
                if( !other )
                {
                    std::string s = char_array_to_string( pm.SB[i], MAXICNAME );
                    s.erase( s.find_last_not_of( " \t" ) + 1 );
                    onlyAq += ( onlyAq.empty() ? "" : "," ) + s;
                }
            }
            if( onlyAq.empty() ) { stallWatch->wLO = pm.LO; stallWatch->wSeenWet = dryWarmStart; }
            else native_trace_decide( "dryoff onlyaq=%s", onlyAq.c_str() );
            if( !onlyAq.empty() )
                ipm_logger->info( "Optima: water-only elements ({}), the dry-out is off (pa_OptimaDryOut)", onlyAq );
        }
        stallWatch->reset();
        if( pa_p->OptimaDryRetry && !g_optimaDryOff && pm.LO >= 0 && pm.LO < pm.L && node1 != nullptr )
        {
            const long int xHd = node1->IC_name_to_xCH( "H" );
            const long int xOd = node1->IC_name_to_xCH( "O" );
            if( xHd >= 0 && xOd >= 0 )
            {
                stallWatch->dryLO = pm.LO;
                stallWatch->dryCap = std::min( pm.B[xHd] / 2., pm.B[xOd] );
            }
        }

        // pa_OptimaEarlyStabilityAt, resolved by optima_earlystability_at() (0 = AUTO), so the
        // solver and the trace's EFF line agree.
        //   < 0  trend trigger: -N ends the first attempt once some multicomponent phase has
        //        fallen monotonically for N consecutive objective evaluations (counters kept in
        //        the objective callback). Spent only on a run
        //        that shows the signature.
        //   > 0  iteration cap (below).
        // A false trigger costs the abandoned attempt only: the safety net re-solves at the
        // full budget when the early look found nothing to repair.
        const long int earlyStabilityAt =
            optima_earlystability_at( pa_p->OptimaEarlyStabilityAt,
                                      optima_multisite_phase_count( pm.sMod, pm.FIs ),
                                      pm.pNP != 0 );
        const long int earlyTrendN = earlyStabilityAt < 0 ? -earlyStabilityAt : 0;
        // Only a guard that the phase really is below its own running peak, not a magnitude
        // test; the consecutive-fall count is what discriminates.
        const double kEarlyTrendDropRatio = 0.999;
        // Dual-settled gate: the repair loop's verdict is a stability index computed from the
        // dual, so the trend trigger fires only when max|dw|/max|w| <= pa_OptimaTol.
        const double kEarlyTrendDualSettled = pa_p->OptimaTol;
        // Rate clause: the average fractional fall per evaluation, (1 - frac)/falls, must be
        // at least kEarlyTrendMinRate. Separates a dissolving phase from one that is only
        // settling slowly (where, on a warm leg, the dual gate cannot help). Scale-free in the
        // length of the fall.
        const double kEarlyTrendMinRate = 1e-4;
        auto earlyTrend      = std::make_shared<bool>( false );
        auto earlyTrendArmed = std::make_shared<bool>( earlyTrendN > 0 );
        // The cap form (> 0): stop the first attempt at iteration N, provided the dual has
        // settled. Same gate and threshold as the trend form - the early look's verdict is a
        // stability index computed from the dual, so it is worth what the dual is worth.
        const long int earlyCapN = earlyStabilityAt > 0 ? earlyStabilityAt : 0;
        auto earlyCap      = std::make_shared<bool>( false );
        // Armed immediately before the primary solve, not here: the compaction probe and the
        // cheap-Hessian attempt copy `options` (and this hook), and a cap firing inside either
        // would mark the primary solve as failed.
        auto earlyCapArmed = std::make_shared<bool>( false );
        // The dual movement at the moment either form fired, for the DECIDE record.
        auto earlyStopDwRel = std::make_shared<double>( -1. );
        // GEMS3K_EARLYTREND_PROBE=<path>: one line the first time each phase satisfies the
        // trend trigger's count-and-drop clause, whether or not the field is armed and whether
        // or not the dual gate would let it fire. Prints the phase's fraction of its own peak,
        // the average fractional fall per evaluation and the dual's movement. Zero cost when
        // unset. GEMS3K_EARLYTREND_PROBE_N sets the fall count surveyed when the field is off
        // (default 200); when it is armed, the field's own -N is used.
        static FILE* earlyTrendDbg = []() -> FILE* {
            const char* fn = std::getenv( "GEMS3K_EARLYTREND_PROBE" );
            return fn ? fopen( fn, "a" ) : nullptr; }();
        static const long int kEarlyTrendProbeN = []() -> long int {
            const char* v = std::getenv( "GEMS3K_EARLYTREND_PROBE_N" );
            return v ? std::atol( v ) : 200; }();
        const long int earlyTrendProbeN = earlyTrendN > 0 ? earlyTrendN : kEarlyTrendProbeN;
        auto earlyTrendSeen = std::make_shared<std::vector<char>>( pm.FIs, (char)0 );
        // The solvent phase is excluded: the phase-selection loop never acts on it, so a
        // trigger on it can only waste a probe.
        long int aqPhaseIdxTrend = -1;
        {
            long int j0aq = 0;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                if( hasAq && pm.LO >= j0aq && pm.LO < j0aq + pm.L1[k] ) { aqPhaseIdxTrend = k; break; }
                j0aq += pm.L1[k];
            }
        }
        const long int ffFrozenN = pa_p->OptimaFastFrozen;   // pa_OptimaFastFrozen
        const double ffDivR = pa_p->OptimaFastDiverge;        // pa_OptimaFastDiverge
        const long int ffRisesN = pa_p->OptimaFastRises;      // pa_OptimaFastRises
        const double ffCollapseR = pa_p->OptimaFastCollapse;  // pa_OptimaFastCollapse
        // pa_OptimaFrozenCut; off when Optima starts from a native answer (HOP/SHP, native finish),
        // where a long still stretch can still end well.
        const long int frozenCutN = runKinetics ? pa_p->OptimaFrozenCut : 0;
        const double kTailTol = pa_p->OptimaTol;   // pa_OptimaTailStop: the error that still counts
        // pa_OptimaTailStop / pa_OptimaTailStopWarm: flat steps needed (a warm call starts next to an answer)
        const long int kTailN = pm.pNP == 1 ? pa_p->OptimaTailStopWarm : pa_p->OptimaTailStop;
        // pa_OptimaTailDryOnly: 0 = every solve, 1 = only the dry solve, 2 = also a warm call whose left-over
        // error sits only on members of nearly vanished non-aqueous solution phases
        const long int kTailMode = pa_p->OptimaTailDryOnly;
        const bool kTailWarm = pm.pNP != 0;
        const double kVanished = 1e7 * pa_p->LogBarrierTau;   // the nearly-vanished amount of pa_OptimaFloorCurvature 3
        std::vector<long int> solidPhaseOf( (size_t)L, -1 );
        for( long int k = 0; k < pm.FIs; k++ )
            if( pm.L1[k] > 1 && !( hasAq && pm.LO >= phStart[(size_t)k] && pm.LO < phStart[(size_t)k+1] ) )
                for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ ) solidPhaseOf[(size_t)j] = k;
        if( stallWatch->window > 0 || stallWatch->maxSeconds > 0. || earlyTrendN > 0
            || earlyCapN > 0 || g_optimaIterBudget >= 0 )
            options.convergence.check =
                [stallWatch, earlyTrendN, kEarlyTrendDropRatio, kEarlyTrendDualSettled,
                 kEarlyTrendMinRate,
                 earlyTrend, earlyTrendArmed, earlyTrendProbeN, earlyTrendSeen,
                 earlyCapN, earlyCap, earlyCapArmed, earlyStopDwRel,
                 aqPhaseIdxTrend, phDec, phMax, phLast, ffFrozenN, ffDivR, ffRisesN, ffCollapseR, frozenCutN, kTailTol, kTailN, kTailMode, kTailWarm, kVanished, solidPhaseOf]
                ( Optima::ConvergenceCheckArgs const& args ) -> bool
                {
                    // pa_OptimaFastCollapse: the quick first try stops at once when its energy has grown to more
                    // than this many times its start; such a try does not come back.
                    if( ffCollapseR > 0. && g_optimaFixAccept )
                    {
                        OptimaQuickTryWatch& cw = g_optimaQuickTry;
                        const double fNow = args.F.result().f.f;
                        if( cw.cstep == 0 ) cw.cf0 = fNow;
                        else if( cw.cf0 != 0. && std::fabs( fNow ) > ffCollapseR * std::fabs( cw.cf0 ) )
                        {
                            native_trace_decide( "fastcollapse step=%ld ratio=%.3e", cw.cstep, std::fabs( fNow / cw.cf0 ) );
                            ipm_logger->info( "Optima: the quick first try stopped at step {}: energy grew to {:.3g} times its start",
                                               cw.cstep, std::fabs( fNow / cw.cf0 ) );
                            cw.cstep++;
                            g_optimaIterBudget = 0; stallWatch->timedOut = true; return true;
                        }
                        cw.cstep++;
                    }
                    // pa_OptimaFastDiverge: once the energy of the quick first try has jumped above 10 times its
                    // start, the try stops 10 steps later if its lowest error since the jump is still above
                    // this many times the lowest before it.
                    if( ffDivR > 0. && g_optimaFixAccept )
                    {
                        OptimaQuickTryWatch& dv = g_optimaQuickTry;
                        const double fNow = args.F.result().f.f, eNow = args.E.error();
                        if( dv.step == 0 ) { dv.f0 = fNow; dv.preMin = eNow; }
                        else if( dv.jumpAt < 0 )
                        {
                            if( std::fabs( fNow ) > 10. * std::fabs( dv.f0 ) ) { dv.jumpAt = dv.step; dv.postMin = eNow; }
                            else dv.preMin = std::min( dv.preMin, eNow );
                        }
                        else
                        {
                            dv.postMin = std::min( dv.postMin, eNow );
                            if( dv.step == dv.jumpAt + 10 && dv.postMin > ffDivR * dv.preMin )
                            {
                                native_trace_decide( "fastdiverge jump_at=%ld pre=%.3e post=%.3e", dv.jumpAt, dv.preMin, dv.postMin );
                                dv.step++;
                                g_optimaIterBudget = 0; stallWatch->timedOut = true; return true;
                            }
                        }
                        dv.step++;
                    }
                    // pa_OptimaFastFrozen: the quick first try stops after this many steps in a row that moved
                    // the state and the error only at rounding level; such a step repeats itself.
                    if( ffFrozenN > 0 && g_optimaFixAccept )
                    {
                        OptimaQuickTryWatch& fz = g_optimaQuickTry;
                        double dx = 0., xm = 0., dw = 0., wm = 0.;
                        const bool sizes = args.u.x.size() == args.uo.x.size() && args.u.w.size() == args.uo.w.size();
                        if( sizes )
                        {
                            for( long int j = 0; j < (long int)args.u.x.size(); j++ ) { dx = std::max( dx, std::fabs( args.u.x[j] - args.uo.x[j] ) ); xm = std::max( xm, std::fabs( args.u.x[j] ) ); }
                            for( long int i = 0; i < (long int)args.u.w.size(); i++ ) { dw = std::max( dw, std::fabs( args.u.w[i] - args.uo.w[i] ) ); wm = std::max( wm, std::fabs( args.u.w[i] ) ); }
                        }
                        const double err = args.E.error();
                        const double dxRel = xm > 0. ? dx / xm : 0., dwRel = wm > 0. ? dw / wm : 0.;
                        const bool still = sizes && dxRel <= 1e-14 && dwRel <= 1e-13 && fz.errPrev >= 0. && std::fabs( err - fz.errPrev ) <= 1e-12 * err;
                        fz.errPrev = err;
                        fz.frozenRun = still ? fz.frozenRun + 1 : 0;
                        if( fz.frozenRun >= ffFrozenN )
                        {
                            native_trace_decide( "fastfrozen steps=%d", fz.frozenRun );
                            g_optimaIterBudget = 0; stallWatch->timedOut = true; return true;
                        }
                    }
                    // pa_OptimaFastRises: the quick first try stops once its energy has gone up by more than 1 %
                    // on this many steps; a try that keeps climbing back rarely finishes.
                    if( ffRisesN > 0 && g_optimaFixAccept )
                    {
                        OptimaQuickTryWatch& rs = g_optimaQuickTry;
                        const double fNow = args.F.result().f.f;
                        if( rs.riseSteps > 0 && fNow - rs.fPrev > 0.01 * std::fabs( rs.fPrev ) ) rs.rises++;
                        rs.fPrev = fNow; rs.riseSteps++;
                        if( rs.rises >= ffRisesN )
                        {
                            native_trace_decide( "fastrises rises=%ld step=%ld", rs.rises, rs.riseSteps );
                            ipm_logger->info( "Optima: the quick first try stopped after its energy went up on {} steps (step {})",
                                               rs.rises, rs.riseSteps );
                            g_optimaIterBudget = 0; stallWatch->timedOut = true; return true;
                        }
                    }
                    if( g_optimaIterBudget > 0 ) --g_optimaIterBudget;
                    if( g_optimaIterBudget == 0 ) { stallWatch->timedOut = true; return true; }   // budget used up
                    // The dual's relative movement over the last step, max|dw|/max|w| on
                    // u = (x, p, w); the hook receives both iterates. One definition for the
                    // probe, the trend trigger and the cap.
                    auto dualRelMove = [&args]() -> double
                    {
                        double dwAbs = 0., wAbs = 0.;
                        for( long int iw = 0; iw < (long int)args.u.w.size(); iw++ )
                        {
                            const double d = std::fabs( args.u.w[iw] - args.uo.w[iw] );
                            if( d > dwAbs ) dwAbs = d;
                            const double a = std::fabs( args.u.w[iw] );
                            if( a > wAbs ) wAbs = a;
                        }
                        return ( wAbs > 0. ? dwAbs / wAbs : 0. );
                    };
                    // Diagnostic: record what the trigger would see, once per phase,
                    // independently of arming and of the gate.
                    if( earlyTrendDbg != nullptr )
                        for( size_t k = 0; k < phDec->size(); ++k )
                            if( !(*earlyTrendSeen)[k]
                                && (*phDec)[k] >= earlyTrendProbeN
                                && (*phMax)[k] > 0.
                                && (*phLast)[k] >= 0.
                                && (*phLast)[k] < (*phMax)[k] * kEarlyTrendDropRatio )
                            {
                                (*earlyTrendSeen)[k] = 1;
                                const double dwRelProbe = dualRelMove();
                                const double frac = (*phLast)[k] / (*phMax)[k];
                                fprintf( earlyTrendDbg,
                                         "TREND k=%ld falls=%ld last=%.6e peak=%.6e frac=%.6e "
                                         "rate=%.6e dwRel=%.6e aq=%d armed=%d\n",
                                         (long)k, (*phDec)[k], (*phLast)[k], (*phMax)[k], frac,
                                         ( 1. - frac ) / (double)(*phDec)[k],
                                         dwRelProbe,
                                         (int)( (long int)k == aqPhaseIdxTrend ),
                                         (int)( earlyTrendN > 0 ) );
                                fflush( earlyTrendDbg );
                            }
                    // Trend trigger first: cheap, independent of the stall signals, armed only
                    // for the first attempt.
                    if( *earlyTrendArmed && !*earlyTrend )
                    {
                        // Fire only once the dual has settled enough for the repair loop's
                        // stability index to be trusted; otherwise re-test next evaluation.
                        const double dwRel = dualRelMove();
                        if( dwRel <= kEarlyTrendDualSettled )
                        for( size_t k = 0; k < phDec->size(); ++k )
                            if( (long int)k != aqPhaseIdxTrend
                                && (*phDec)[k] >= earlyTrendN
                                && (*phMax)[k] > 0.
                                && (*phLast)[k] >= 0.
                                && (*phLast)[k] < (*phMax)[k] * kEarlyTrendDropRatio
                                && ( 1. - (*phLast)[k] / (*phMax)[k] )
                                       / (double)(*phDec)[k] >= kEarlyTrendMinRate )
                            {
                                *earlyTrend = true;
                                *earlyStopDwRel = dwRel;
                                ipm_logger->debug( "CalculateEquilibriumStateOptima: phase {} has "
                                                  "fallen for {} consecutive evaluations to {:.3e} "
                                                  "from a peak of {:.3e} - ending the first attempt "
                                                  "so the phase-selection loop can look "
                                                  "(pa_OptimaEarlyStabilityAt = -{}; dual settled, "
                                                  "max|dw|/max|w| = {:.3e} <= {:.3e}; fall rate "
                                                  "{:.3e}/eval >= {:.3e})",
                                                  k, (*phDec)[k], (*phLast)[k], (*phMax)[k],
                                                  earlyTrendN, dwRel, kEarlyTrendDualSettled,
                                                  ( 1. - (*phLast)[k] / (*phMax)[k] )
                                                      / (double)(*phDec)[k], kEarlyTrendMinRate );
                                break;
                            }
                        if( *earlyTrend ) return true;   // folded back to a failure below
                    }
                    // The cap form (pa_OptimaEarlyStabilityAt > 0): at iteration N, if the dual
                    // has settled, end the first attempt so the repair loop can look. Asked
                    // once and then disarmed, whatever the answer; waiting at N for the dual to
                    // settle would turn a bounded cost into nearly a whole solve on ordinary
                    // runs. args.result.iterations is current here (Optima checks convergence
                    // before its budget test). A `>=` latch, so a hook not called at exactly
                    // iteration N cannot turn the field off.
                    if( *earlyCapArmed
                        && (long int)args.result.iterations >= earlyCapN )
                    {
                        *earlyCapArmed = false;   // asked and answered, either way
                        const double dwRelCap = dualRelMove();
                        // Recorded either way, so the DECIDE record can report a blocked cap.
                        *earlyStopDwRel = dwRelCap;
                        if( dwRelCap <= kEarlyTrendDualSettled )
                        {
                            *earlyCap = true;
                            ipm_logger->debug( "CalculateEquilibriumStateOptima: reached iteration {} "
                                              "with the dual settled - ending the first attempt so "
                                              "the phase-selection loop can look "
                                              "(pa_OptimaEarlyStabilityAt = {}; max|dw|/max|w| = "
                                              "{:.3e} <= {:.3e})",
                                              (long)args.result.iterations, earlyCapN,
                                              dwRelCap, kEarlyTrendDualSettled );
                            return true;   // folded back to a failure below
                        }
                        ipm_logger->debug( "CalculateEquilibriumStateOptima: reached iteration {} but "
                                          "the dual is still moving (max|dw|/max|w| = {:.3e} > "
                                          "{:.3e}) - NOT ending the first attempt; a stability "
                                          "verdict taken here would be read off noise "
                                          "(pa_OptimaEarlyStabilityAt = {})",
                                          (long)args.result.iterations, dwRelCap,
                                          kEarlyTrendDualSettled, earlyCapN );
                    }
                    // Two signals (optimality error and complementarity), and both must
                    // stagnate: each alone can stagnate on a run that does converge. The
                    // wall-clock guard comes first, since it must fire even while the iterate
                    // is still improving.
                    if( stallWatch->overBudget() )
                    {
                        stallWatch->timedOut = true;
                        return true;   // stop now; folded back to a failure below
                    }
                    // pa_OptimaDryRetry: water below 1e-3 of what the bulk could form for 100
                    // iterations ends the attempt, so the dry retry runs. Not in the safety net
                    // (window 0), which must run to the end.
                    if( stallWatch->window > 0 && stallWatch->dryLO >= 0
                        && stallWatch->dryLO < (long int)args.u.x.size() )
                    {
                        if( args.u.x[stallWatch->dryLO] < 1e-3 * stallWatch->dryCap )
                        {
                            if( ++stallWatch->dryRun >= 100 )
                            {
                                stallWatch->stalled = true;
                                stallWatch->dryStop = true;
                                return true;   // folded back to a failure below
                            }
                        }
                        else stallWatch->dryRun = 0;
                    }
                    // pa_OptimaFrozenCut: energy and error unchanged for this many steps means the solve
                    // repeats itself; stop and let the finish and retries take over. Not in the quick try
                    // or the safety net.
                    if( frozenCutN > 0 && stallWatch->window > 0 && !g_optimaFixAccept )
                    {
                        const double fNow = args.F.result().f.f, eNow = args.E.error();
                        const bool same = stallWatch->fzErr >= 0.
                            && std::fabs( fNow - stallWatch->fzF ) <= 1e-15 * std::fabs( fNow )
                            && std::fabs( eNow - stallWatch->fzErr ) <= 1e-12 * eNow;
                        stallWatch->fzF = fNow; stallWatch->fzErr = eNow;
                        stallWatch->fzRun = same ? stallWatch->fzRun + 1 : 0;
                        if( stallWatch->fzRun >= frozenCutN )
                        {
                            native_trace_decide( "frozencut steps=%ld it=%ld err=%.6e", stallWatch->fzRun,
                                                 (long)args.result.iterations, eNow );
                            ipm_logger->info( "Optima: solve stopped at step {}: energy and error unchanged for {} steps",
                                               (long)args.result.iterations, stallWatch->fzRun );
                            stallWatch->stalled = true;
                            return true;   // folded back to a failure below
                        }
                    }
                    // pa_OptimaDryOut: the water has run out during the solve; stop so the water
                    // phase can be taken out (its species fixed at the floor) and the solids solved alone.
                    {
                        if( !stallWatch->dryDone && stallWatch->window > 0 && !g_optimaFixAccept && stallWatch->wLO >= 0
                            && stallWatch->wLO < (long int)args.u.x.size() )
                        {
                            const auto xv = args.u.x;
                            double xSum = 0.;
                            for( long int j = 0; j < (long int)xv.size(); j++ ) xSum += std::max( (double)xv[j], 0. );
                            // only water that ran out during this call counts, not a start that is already nearly dry
                            if( xv[stallWatch->wLO] >= 1e-6 * xSum ) stallWatch->wSeenWet = true;
                            else if( stallWatch->wSeenWet )
                            {
                                stallWatch->dryOut = true;
                                stallWatch->dryOutAt = (long int)args.result.iterations;
                                stallWatch->dryOutW = xv[stallWatch->wLO];
                                stallWatch->stalled = true;
                                return true;   // folded back to a failure below
                            }
                        }
                    }
                    // pa_OptimaTailStop: the energy has not moved for M steps and all the error left
                    // sits on species holding a negligible amount; stop and let the finish and the answer
                    // check decide.
                    {
                        const long int kTail = kTailN;
                        // pa_OptimaTailDryOnly: where the tail stop acts (see above)
                        const bool vanishedOnly = kTailMode == 2 && !stallWatch->dryDone && !stallWatch->vanishedOff;
                        if( kTail > 0 && stallWatch->window > 0 && !g_optimaFixAccept
                            && ( kTailMode == 0 || stallWatch->dryDone || ( vanishedOnly && kTailWarm ) ) )
                        {
                            const double fNow = args.F.result().f.f;
                            const double eNow = args.E.errorx();
                            if( stallWatch->tailRun > 0 && std::fabs( fNow - stallWatch->tailRef ) <= 1e-12 * std::fabs( fNow ) )
                            {
                                stallWatch->tailRun++;
                                if( stallWatch->tailErrPrev > 0. )
                                {
                                    const double q = eNow / stallWatch->tailErrPrev;
                                    stallWatch->tailRatioMin = std::min( stallWatch->tailRatioMin, q );
                                    stallWatch->tailRatioMax = std::max( stallWatch->tailRatioMax, q );
                                }
                            }
                            else { stallWatch->tailRef = fNow; stallWatch->tailRun = 1; stallWatch->tailRatioMin = 1e300; stallWatch->tailRatioMax = 0.; }
                            stallWatch->tailErrPrev = eNow;
                            if( stallWatch->tailRun >= kTail )
                            {
                                const auto exv = args.E.ex();
                                const auto xv = args.u.x;
                                const long int n = std::min( (long int)exv.size(), (long int)xv.size() );
                                double xSum = 0.;
                                for( long int j = 0; j < n; j++ ) xSum += std::max( (double)xv[j], 0. );
                                const double tol = kTailTol;
                                bool onlyTrace = true; long int nCarry = 0; double xCarryMax = 0.;
                                std::vector<long int> carrierPhases;
                                for( long int j = 0; j < n; j++ )
                                    if( std::fabs( (double)exv[j] ) > tol )
                                    {
                                        nCarry++; xCarryMax = std::max( xCarryMax, (double)xv[j] );
                                        if( xv[j] > 1e-9 * xSum ) { onlyTrace = false; break; }
                                        if( vanishedOnly )
                                        {
                                            const long int kp = j < (long int)solidPhaseOf.size() ? solidPhaseOf[(size_t)j] : -1;
                                            if( kp < 0 || !( (*phLast)[(size_t)kp] >= 0. && (*phLast)[(size_t)kp] < kVanished ) )
                                            { onlyTrace = false; break; }
                                            if( std::find( carrierPhases.begin(), carrierPhases.end(), kp ) == carrierPhases.end() )
                                                carrierPhases.push_back( kp );
                                        }
                                    }
                                if( onlyTrace && nCarry > 0 )
                                {
                                    char rr[96] = "";
                                    if( vanishedOnly ) std::snprintf( rr, sizeof rr, " scope=vanished ratio=%.3f..%.3f", stallWatch->tailRatioMin, stallWatch->tailRatioMax );
                                    native_trace_decide( "tailstop it=%ld flat=%ld carriers=%ld xmax=%.3e xsum=%.3e err=%.3e%s",
                                                         (long)args.result.iterations, stallWatch->tailRun, nCarry,
                                                         xCarryMax, xSum, args.E.errorx(), rr );
                                    stallWatch->tailStop = true;
                                    stallWatch->tailPhases = carrierPhases;
                                    stallWatch->stalled = true;
                                    return true;   // folded back to a failure below
                                }
                            }
                        }
                    }
                    if( stallWatch->window <= 0 ) return false;
                    const double e = args.E.errorx();
                    double comp = 0.;
                    const auto ex = args.E.ex();
                    const auto xv = args.u.x;
                    const Optima::Index nx = ex.size() < xv.size() ? ex.size() : xv.size();
                    for( Optima::Index j = 0; j < nx; ++j )
                    {
                        const double c = std::fabs( ex[j] ) * std::fabs( xv[j] );
                        if( c > comp ) comp = c;
                    }
                    // The test is cumulative over a whole window: fire only when best-so-far
                    // has failed to fall by a relative kStallRelProgress against its value one
                    // window ago. A slow plateau that still makes real progress survives; a
                    // frozen best-so-far fires. Unlike the pre-solve watch in
                    // OptimaReducedPreSolve(), there is no "has the live error moved" clause:
                    // a live error running a limit cycle must not count as progress here.
                    // `best` is updated on every fall, however small.
                    static const double kStallRelProgress = 1e-8;
                    if( e    < stallWatch->bestErr  ) stallWatch->bestErr  = e;
                    if( comp < stallWatch->bestComp ) stallWatch->bestComp = comp;
                    if( ++stallWatch->run < stallWatch->window ) return false;

                    const bool firstWindow = !std::isfinite( stallWatch->refErr )
                                          && !std::isfinite( stallWatch->refComp );
                    const bool dropped =
                        ( std::isfinite( stallWatch->refErr )
                          && stallWatch->bestErr  < stallWatch->refErr  * ( 1. - kStallRelProgress ) )
                     || ( std::isfinite( stallWatch->refComp )
                          && stallWatch->bestComp < stallWatch->refComp * ( 1. - kStallRelProgress ) );
                    if( !firstWindow && !dropped )
                    {
                        stallWatch->stalled = true;
                        return true;   // stop now; folded back to a failure below
                    }
                    // Window survived: roll the reference forward, start a new one.
                    stallWatch->refErr  = stallWatch->bestErr;
                    stallWatch->refComp = stallWatch->bestComp;
                    stallWatch->run     = 0;
                    return false;
                };
        // Applied after each solve(): a stall or timeout is turned into a failure, logged and
        // recorded in the trace.
        auto applyStall = [stallWatch]( Optima::Result& r ) {
            // Does not fold the early-trend stop: that is done once, at the primary solve,
            // since the flag stays set during the retries.
            if( stallWatch->stalled || stallWatch->timedOut )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: full solve abandoned after {}"
                                  " - pa_OptimaStallWindow={}",
                                  stallWatch->timedOut ? "the time or iteration budget"
                                                       : "a window with no meaningful progress",
                                  stallWatch->window );
                // DECIDE record: this site abandons a solve and hands it to the retry tiers.
                // `iters` is Optima's count for the abandoned solve, `run` how far into the
                // current window it got, `kind` separates a stall from a wall-clock timeout.
                native_trace_decide( "optimastall kind=%s window=%ld run=%ld iters=%ld"
                                     " besterr=%.6e bestcomp=%.6e",
                                     stallWatch->timedOut ? "timeout" : "stall",
                                     (long)stallWatch->window, (long)stallWatch->run,
                                     (long)r.iterations, stallWatch->bestErr, stallWatch->bestComp );
                if( stallWatch->stallAt == 0 && !stallWatch->timedOut )
                    stallWatch->stallAt = (long int)r.iterations;
                r.succeeded = false;
            }
        };

        Optima::Solver solver;
        solver.setOptions( options );
        // The original state, for retries that re-solve from it (ROP's fallback, the safety
        // nets). Saved before the first solve modifies `state`.
        const Optima::State initialState = state;

        // Debug trace (GEMS3K_OPTIMA_STEPTRACE_FILE, off unless set): the per-solve setup records.
        long int traceFirstIt = 0, traceFupFirst = 0, traceFxFirst = 0;
        bool traceFirstOk = false;
        std::vector<double> traceX0;
        if( auto* tfp = Optima::trace::file() )
        {
            const auto& ts = Optima::trace::state();
            const bool coldStart = !( ( pm.pNP != 0 && warmStateUsable ) || dimReduceDone );
            long nipps = 0;
            for( long int k = pm.FIs; k < pm.FI; k++ ) nipps++;
            fprintf( tfp, "RUN   solve=%ld T=%.6f P=%.6e Nx=%ld Nn=%ld Nq=%ld Np=%ld Nc=%ld Nw=%ld nspecies=%ld nphase=%ld nipps=%ld\n",
                     ts.solve, pm.Tc, pm.P * 1e5, (long)dims.x, (long)L, (long)R, (long)dims.p, (long)N, 0L,
                     (long)L, (long)pm.FI, nipps );
            const char* lsm = options.newtonstep.linearsolver.method == Optima::LinearSolverMethod::Nullspace ? "nullspace"
                            : options.newtonstep.linearsolver.method == Optima::LinearSolverMethod::Rangespace ? "rangespace" : "fullspace";
            fprintf( tfp, "SET   solve=%ld epsilon=%.3e logarithm_barrier_factor=%.3e tau=%.3e warmstart=%d hessian=%s use_ideal_activity_models=0"
                          " tolerance=%.3e maxiters=%u check_used=%d linsolver=%s fix_and_accept=%d linesearch=%s errorstatus_inc=%.3g errorstatus_dec=%.3g"
                          " epsilon_src=%s tau_src=pa_LogBarrierTau fd_hessian_src=pa_OptimaFDHessian=%ld\n",
                     ts.solve, dcFloor, dcFloor > 0. ? kLogBarrierTau / dcFloor : 0., kLogBarrierTau, coldStart ? 0 : 1,
                     kFDHessian ? "PartiallyExact" : "Approx", options.convergence.tolerance, options.maxiters,
                     options.convergence.check ? 1 : 0, lsm, options.backtracksearch.apply_min_max_fix_and_accept ? 1 : 0,
                     options.linesearch.enabled ? "on" : "off", options.errorstatus.significantly_increased,
                     options.errorstatus.significantly_decreased,
                     pa_p->OptimaDcFloor > 0. ? "pa_OptimaDcFloor" : "pa_DHB", (long)pa_p->OptimaFDHessian );
            long nloweps = 0, nlowabove = 0, nupfinite = 0;
            double lmin = 1e300, lmax = -1e300;
            for( long int j = 0; j < L; j++ )
            {
                const double lo = problem.xlower[j];
                lmin = std::min( lmin, lo ); lmax = std::max( lmax, lo );
                if( lo <= dcFloor ) nloweps++; else nlowabove++;
                if( std::isfinite( problem.xupper[j] ) ) nupfinite++;
            }
            fprintf( tfp, "EFF   solve=%ld ScFact=%.9g xlower_min=%.3e xlower_max=%.3e n_lower_at_eps=%ld n_lower_above_eps=%ld n_upper_finite=%ld"
                          " nplower_finite=0 npupper_finite=0\n", ts.solve, ScFact, lmin, lmax, nloweps, nlowabove, nupfinite );
            traceX0.assign( state.x.data(), state.x.data() + L );
            auto inf = []( const auto& v ) { return v.size() ? (double)v.cwiseAbs().maxCoeff() : 0.0; };
            auto nnz = []( const auto& v ) { long k = 0; for( Eigen::Index i = 0; i < v.size(); ++i ) if( v[i] != 0.0 ) k++; return k; };
            long nzero = 0, nateps = 0, nclamp = 0;
            double xmin = 1e300, xmax = -1e300, dprev = -1.0;
            for( long int j = 0; j < L; j++ )
            {
                const double x = traceX0[j];
                xmin = std::min( xmin, x ); xmax = std::max( xmax, x );
                if( x == 0.0 ) nzero++;
                if( x <= dcFloor ) nateps++;
                if( x < problem.xlower[j] ) nclamp++;
            }
            if( (long)traceXPrev.size() == L )
            {
                dprev = 0.0;
                for( long int j = 0; j < L; j++ )
                    dprev = std::max( dprev, std::abs( std::log( std::max( traceX0[j], dcFloor ) ) - std::log( std::max( traceXPrev[j], dcFloor ) ) ) );
            }
            fprintf( tfp, "INIT  solve=%ld from=%s x_min=%.6e x_max=%.6e n_x_zero=%ld n_x_at_or_below_eps=%ld n_x_clamped_up=%ld"
                          " ye_inf=%.6e yg_inf=%.6e ze_inf=%.6e zg_inf=%.6e s_inf=%.6e n_y_nonzero=%ld n_s_nonzero=%ld stored_js=%ld stored_ju=%ld stored_jb=%ld"
                          " p_inf=%.6e dist_prev_solve=%.6e\n",
                     ts.solve, dimReduceDone ? "dimreduce" : coldStart ? "cold" : "warm", xmin, xmax, nzero, nateps, nclamp,
                     inf( state.ye ), inf( state.yg ), inf( state.ze ), inf( state.zg ), inf( state.s ),
                     nnz( state.ye ) + nnz( state.yg ), nnz( state.s ), (long)state.js.size(), (long)state.ju.size(), (long)state.jb.size(),
                     inf( state.p ), dprev );
            fflush( tfp );
        }

        // Presence threshold: not bare pm.DSM. Optima never reaches exactly 0 (every species
        // has a dcFloor lower bound), so a phase held exactly at dcFloor is Optima's form of
        // "absent"; with a project's pa_DS far below the floor it would read as present.
        const double presenceThreshold = std::max( pm.DSM, dcFloor * 10. );

        // Phase-selection state, shared by the compaction probe and the phase-selection loop
        // (a phase pinned by the probe must be readmissible by that loop):
        //   phSelState[k]: 0 = untouched, 1 = deactivated, 2 = readmitted (terminal - never
        //                  dropped again, which bounds the loop).
        std::vector<char> phSelState( (size_t)std::max( pm.FI, 1L ), 0 );
        std::vector<double> savedLo( (size_t)std::max( L, 1L ), 0. );
        std::vector<double> savedHi( (size_t)std::max( L, 1L ), 0. );


        // ---- Tiered Hessian: a short attempt with the cheap Hessian first ----
        // pa_OptimaFDHessianDelay = N runs up to N iterations with the FD PartiallyExact loop
        // suppressed. If that converges, its state is adopted and the primary solve below
        // re-runs from it, still cheap. If not, the state is left untouched and the primary
        // solve starts from the original seed with FD on (the cheap trajectory is discarded,
        // not corrected).
        long int fdCheapIterations = 0;
        bool cheapAttemptWon = false;
        // Structural skip: multisite / reciprocal solid solutions (Berman, CEF, Modified
        // Bragg-Williams) and fluid EoS phases need the FD columns from the first iteration
        // (in a multisite model an end-member at a low amount still carries curvature through
        // its site fractions, which the exact per-phase block does not cover). Decided once
        // from the phase models, before spending anything.
        bool fdRequiredByModel = false;
        for( long int k = 0; k < pm.FIs && !fdRequiredByModel; k++ )
        {
            const char mc = pm.sMod[k][SPHAS_TYP];
            if( mc == SM_BERMAN || mc == SM_CEF || mc == SM_MBW ||
                mc == SM_CGFLUID || mc == SM_PRFLUID || mc == SM_PCFLUID ||
                mc == SM_STFLUID || mc == SM_PR78FL || mc == SM_CORKFL ||
                mc == SM_REFLUID || mc == SM_SRFLUID )
                fdRequiredByModel = true;
        }
        if( kFDDelay > 0 && kFDHessian && fdRequiredByModel )
            ipm_logger->warn( "pa_OptimaFDHessianDelay ignored: this system has a multisite or "
                              "fluid-EoS phase, which needs the full Hessian from the start. "
                              "Set it to 0 to silence this (in the project's -ipm file; not in GEMS)." );
        if( kFDDelay > 0 && kFDHessian && !fdRequiredByModel )
        {
            Optima::Options cheapOpts = options;
            cheapOpts.maxiters = kFDDelay;
            Optima::Solver cheapSolver;
            cheapSolver.setOptions( cheapOpts );
            Optima::State cheapState = state;
            *fdSuppress = true;
            stallWatch->reset();
            // The early-trend trigger must not fire inside this attempt: it would mark the
            // primary solve below as failed. The cheap attempt is itself a discarded probe.
            const bool trendWasArmed = *earlyTrendArmed;
            *earlyTrendArmed = false;
            Optima::Result cheapResult = timed_solve( cheapSolver, problem, cheapState );
            *earlyTrendArmed = trendWasArmed;
            *earlyTrend = false;
            *fdSuppress = false;
            fdCheapIterations = cheapResult.iterations;
            cheapAttemptWon = ( cheapResult.succeeded && !stallWatch->stalled && !stallWatch->timedOut );
            if( cheapAttemptWon )
                state = cheapState;
            else
                ipm_logger->debug( "CalculateEquilibriumStateOptima: cheap-Hessian attempt did not "
                                  "converge in {} iterations - discarding it and restarting with the "
                                  "finite-difference Hessian (pa_OptimaFDHessianDelay = {})",
                                  cheapResult.iterations, kFDDelay );
            // If the cheap attempt won, the primary solve re-runs from its state with the FD
            // loop still suppressed, only to confirm an already converged point. Confirming it
            // with the FD columns instead fails: those columns are near-singular at a
            // converged point.
            *fdSuppress = cheapAttemptWon;
        }

        // ---- Diagnostic: phase-stability scan of the starting state ------------
        // GEMS3K_PRESOLVE_PHSTAB_PROBE=<path>: run the phase-assemblage stability check once
        // on the state this call starts from, before the primary solve, and write what it
        // finds (plus a DECIDE presolvephstab record). Does not feed the solve; not reached
        // when the variable is unset.
        {
            static FILE* preDbg = []() -> FILE* {
                const char* fn = std::getenv( "GEMS3K_PRESOLVE_PHSTAB_PROBE" );
                return fn ? fopen( fn, "a" ) : nullptr; }();
            if( preDbg && L > 0 )
            {
                // The scan reads derived arrays, so refresh them from the starting state;
                // they are recomputed after the solve on every path.
                for( long int j = 0; j < L; j++ )
                { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
                for( long int i = 0; i < N; i++ )
                    pm.U[i] = -state.ye[i];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                double preViol = 0.;
                bool preAbsent = false;
                std::vector<PhStabViolation> preRanked;
                PhStabCensus preCensus;
                const long int preWorst = WorstPhaseStabilityViolation(
                            presenceThreshold, dcFloor, nullptr,
                            preViol, preAbsent, &preRanked, &preCensus, "presolve" );
                fprintf( preDbg,
                         "PRESOLVE pNP=%ld warm=%d capN=%ld trendN=%ld "
                         "worst=%ld viol=%.6e absent=%d "
                         "absent-but-stable=%ld present-but-unstable=%ld "
                         "scanned=%ld clamped=%ld ranked=%ld\n",
                         (long)pm.pNP, (int)warmStateUsable, (long)earlyCapN,
                         (long)earlyTrendN, (long)preWorst, preViol, (int)preAbsent,
                         (long)preCensus.absentStable, (long)preCensus.presentUnstable,
                         (long)preCensus.scanned, (long)preCensus.clamped,
                         (long)preRanked.size() );
                for( size_t r = 0; r < preRanked.size() && r < 8; r++ )
                    fprintf( preDbg, "PRESOLVE-RANK %zu k=%ld %-16s viol=%.6e absent=%d clamped=%d\n",
                             r, (long)preRanked[r].k,
                             char_array_to_string( pm.SF[preRanked[r].k], MAXPHNAME ).c_str(),
                             preRanked[r].viol, (int)preRanked[r].wasAbsent,
                             (int)preRanked[r].clamped );
                fflush( preDbg );
                native_trace_decide( "presolvephstab warm=%d worst=%ld viol=%.6e "
                                     "absentstable=%ld presentunstable=%ld",
                                     (int)warmStateUsable, (long)preWorst, preViol,
                                     (long)preCensus.absentStable,
                                     (long)preCensus.presentUnstable );
            }
        }

        std::vector<char> extinctFixed( (size_t)std::max(L,1L), 0 );
        bool tailAccepted = false;   // pa_OptimaTailStop: the stop point was taken as the answer
        double dryZzSaved = 0.; bool dryZzSet = false;   // pa_OptimaDryOut: the charge-row target before the dry solve
        std::unique_ptr<Optima::State> dryStopSaved;      // pa_OptimaDryOut: where the dry-out stopped the solve
        long int dryJ0 = 0, dryJ1 = 0;                     // pa_OptimaDryOut: the water phase's species
        std::vector<double> extLo( (size_t)std::max(L,1L), 0. ), extHi( (size_t)std::max(L,1L), 0. );
        stallWatch->reset();
        Optima::Sensitivity sensitivity( dims );
        // Arm the cap here and nowhere earlier: late arming scopes it to the first attempt.
        *earlyCapArmed = ( earlyCapN > 0 );
        Optima::Result result = wantSens ? timed_solve( solver, problem, state, sensitivity )
                                         : timed_solve( solver, problem, state );
        // Returning true from convergence.check means "stop", and Optima reports that as
        // success, so an early stop is folded back to a failure here, before anything reads
        // result.succeeded.
        if( *earlyTrend || *earlyCap ) result.succeeded = false;
        if( Optima::trace::file() )
        {
            // First pass = reduced pre-solve and cheap-Hessian attempt (if any) plus the first full solve.
            traceFirstIt = (long int)result.iterations + fdCheapIterations + dimReduceIters;
            traceFirstOk = result.succeeded;
            traceFupFirst = Optima::trace::state().fupdates;
            traceFxFirst = Optima::trace::state().fxevals;
        }
        // Did pa_OptimaEarlyStabilityAt, rather than the problem, stop the first attempt?
        // Covers the cap and the trend trigger, each of which sets its own flag in the
        // convergence hook. Consumed by the early-probe safety net.
        const bool earlyCapHit = !result.succeeded && ( *earlyCap || *earlyTrend );
        // ---- Early-probe safety net: restart or resume ------------------------
        // The net re-solves after a probe that found nothing. Restarting from `initialState`
        // costs N + full; resuming from the probe's end state can cost ~full. The probe state
        // is captured here, from the primary solve on the original boxes, not from `state`
        // at the net (retry tiers in between reassign it and may pin species).
        // Mode (optima_net_resume_mode(); GEMS3K_OPTIMA_NET_RESUME overrides it):
        //   0  restart;
        //   1  resume (can lose answers);
        //   2  resume, and if that does not converge, restart once;
        //   3  resume, and restart only at the end of the retry ladder if the call has
        //      still failed (default).
        const int netResumeMode = optima_net_resume_mode();
        const long int earlyProbeIters = earlyCapHit ? (long int)result.iterations : 0;
        std::unique_ptr<Optima::State> earlyProbeState;
        if( earlyCapHit && netResumeMode >= 1 && netResumeMode <= 3 )
            earlyProbeState.reset( new Optima::State( state ) );
        // Did the problem fail, or did the early stop end the attempt? Every failure-only
        // retry tier below is gated on genuineFailure(): a state truncated on purpose must
        // not be handed to repair machinery that assumes the solve failed (the extinction
        // tier would read a merely shrinking phase as vanishing and return a wrong answer).
        // After an early stop only the phase-selection loop acts; if it finds nothing, the
        // net re-solves at the full budget. A lambda, not a bool: `result` is reassigned by
        // each tier, so each test asks about the state reached by the tiers above it.
        // earlyCapHit is a snapshot of how the primary attempt ended.
        auto genuineFailure = [&]() -> bool { return !result.succeeded && !earlyCapHit; };
        // Disarm both for every retry: the early look is given once.
        *earlyTrendArmed = false;
        *earlyCapArmed   = false;
        // DECIDE record of which early-stop form fired, when, and at what dual movement.
        if( earlyCapHit )
            native_trace_decide( "earlystop form=%s fired=1 at=%ld dwrel=%.6e",
                                 ( *earlyCap ? "cap" : "trend" ),
                                 (long int)result.iterations, *earlyStopDwRel );
        else if( earlyCapN > 0 && *earlyStopDwRel >= 0. )
            // The cap reached N and the dual gate refused. Traced, so "did not fire" can be
            // told apart from "never reached".
            native_trace_decide( "earlystop form=cap fired=0 at=%ld dwrel=%.6e",
                                 earlyCapN, *earlyStopDwRel );
        if( wantSens )
        {
            // Captured from the primary solve only: a retry solves a different problem
            // (reseeded, or with phases pinned).
            optima_dndb_rows = L + R;
            optima_dndb_cols = N;
            optima_dndb.assign( (size_t)(( L + R ) * N), 0. );
            for( long int j = 0; j < L + R; j++ )
                for( long int i = 0; i < N; i++ )
                    optima_dndb[ (size_t)( j * N + i ) ] = sensitivity.xc( j, i );
        }
        // Every retry below runs with the exact (FD) columns: reaching a retry is itself
        // the escalation trigger.
        *fdSuppress = false;
        applyStall( result );
        // Optima's iteration count, summed over every solve() this function makes, goes
        // into pm.ITG below; pm.ITF stays 0.
        long int optimaIterTotal = result.iterations + fdCheapIterations + dimReduceIters;
        // pa_OptimaDryOut: the water ran out - fix every water-phase species at the floor, set the
        // charge row to what those fixed amounts carry (native drops that row once the water is gone),
        // and solve the rest from where the solve stopped.
        if( stallWatch->dryOut && !result.succeeded && pm.LO >= 0 )
        {
            long int j0 = 0, j1 = 0;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                if( pm.LO >= j0 && pm.LO < j0 + pm.L1[k] ) { j1 = j0 + pm.L1[k]; break; }
                j0 += pm.L1[k];
            }
            const long int at = stallWatch->dryOutAt; const double wAt = stallWatch->dryOutW;
            dryStopSaved.reset( new Optima::State( state ) );
            for( long int j = j0; j < j1; j++ )
            {
                if( !extinctFixed[(size_t)j] ) { extLo[(size_t)j] = problem.xlower[j]; extHi[(size_t)j] = problem.xupper[j]; }
                problem.xlower[j] = dcFloor; problem.xupper[j] = dcFloor; state.x[j] = dcFloor;
                extinctFixed[(size_t)j] = 1;
            }
            double zz = 0.;
            if( pm.E > 0 )
            {
                for( long int j = j0; j < j1; j++ ) zz += pm.A[ (N-1) + j*N ] * dcFloor;
                if( !dryZzSet ) { dryZzSaved = problem.be[N-1]; dryZzSet = true; }
                problem.be[N-1] = zz;
            }
            stallWatch->dryDone = true;
            Optima::Solver solverDryOut;
            solverDryOut.setOptions( options );
            stallWatch->reset();
            result = timed_solve( solverDryOut, problem, state );
            applyStall( result );
            optimaIterTotal += result.iterations;
            native_trace_decide( "dryout at=%ld water=%.3e fixed=%ld zz=%.3e ok=%d it=%ld", at, wAt, (long)( j1 - j0 ), zz,
                                 result.succeeded ? 1 : 0, (long)result.iterations );
            ipm_logger->info( "Optima: the water ran out at step {}; the solids were solved alone, ok={} (pa_OptimaDryOut)",
                               at, result.succeeded ? 1 : 0 );
            dryJ0 = j0; dryJ1 = j1;
        }
        // pa_OptimaTailDryOnly 2: a warm call stopped on the slow tail of nearly vanished solid solutions is solved on
        // from the stop point with those phases held at the floor, up to three times (each round may find another such
        // phase); the answer is then checked like any other. If a round fails, or a held phase wants to form, the call goes on from the first stop
        // point with no phase held, as if there had been no stop. DECIDE tailpin.
        if( stallWatch->tailStop && !result.succeeded && R == 0 && !stallWatch->tailPhases.empty() )
        {
            const Optima::State stopState = state;
            std::vector<long int> pinned;
            std::vector<double> pinLo, pinHi;
            long int rounds = 0, nPh = 0;
            bool pinOk = false;
            std::string names;
            while( stallWatch->tailStop && !result.succeeded && !stallWatch->tailPhases.empty() && rounds < 3 )
            {
                for( long int k : stallWatch->tailPhases )
                {
                    nPh++;
                    names += ( names.empty() ? "" : "," ) + char_array_to_string( pm.SF[k] + MAXSYMB, MAXPHNAME );
                    for( long int j = phStart[(size_t)k]; j < phStart[(size_t)k+1]; j++ )
                        if( !extinctFixed[(size_t)j] && std::find( pinned.begin(), pinned.end(), j ) == pinned.end() )
                        {
                            pinned.push_back( j ); pinLo.push_back( problem.xlower[j] ); pinHi.push_back( problem.xupper[j] );
                            problem.xlower[j] = dcFloor; problem.xupper[j] = dcFloor; state.x[j] = dcFloor;
                        }
                }
                rounds++;
                stallWatch->vanishedOff = rounds >= 3;
                Optima::Solver solverPin;
                solverPin.setOptions( options );
                stallWatch->reset();
                result = timed_solve( solverPin, problem, state );
                applyStall( result );
                optimaIterTotal += result.iterations;
                pinOk = result.succeeded;
            }
            stallWatch->vanishedOff = true;
            // A held phase stays held for the rest of the call only if, at the converged potentials, neither it nor any
            // other absent solid mixture wants to form; otherwise the call goes on from the first stop point.
            double pinTpd = 1e300;
            if( pinOk )
            {
                for( long int j = 0; j < L; j++ ) { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
                for( long int i = 0; i < N; i++ ) pm.U[i] = -state.ye[i];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                std::vector<long int> heldPh;
                for( long int j : pinned )
                {
                    long int k = 0;
                    while( k + 1 < (long int)phStart.size() && phStart[(size_t)k+1] <= j ) k++;
                    if( std::find( heldPh.begin(), heldPh.end(), k ) == heldPh.end() ) heldPh.push_back( k );
                }
                // the full composition search: holding a phase for the rest of the call needs more than the local one
                const bool strongWas = g_tpdStrongSearch;
                g_tpdStrongSearch = true;
                for( long int k : heldPh )
                {
                    std::vector<double> yb;
                    pinTpd = std::min( pinTpd, NativeTpdPhase( k, phStart[(size_t)k], yb ) );
                }
                // a solid mixture that was already absent is checked too: the held solve never looked at it
                for( long int k = 0; k < pm.FIs; k++ )
                    if( phStart[(size_t)k] < L && solidPhaseOf[(size_t)phStart[(size_t)k]] == k && pm.XF[k] < kVanished
                        && std::find( heldPh.begin(), heldPh.end(), k ) == heldPh.end() )
                    {
                        std::vector<double> yb;
                        pinTpd = std::min( pinTpd, NativeTpdPhase( k, phStart[(size_t)k], yb ) );
                    }
                g_tpdStrongSearch = strongWas;
            }
            const bool keepHeld = pinOk && !( pinTpd < -std::max( pa_p->OptimaTpdAccept, 1e-6 ) );
            if( !keepHeld )
                for( size_t q = 0; q < pinned.size(); q++ ) { problem.xlower[pinned[q]] = pinLo[q]; problem.xupper[pinned[q]] = pinHi[q]; }
            long int itBack = 0;
            if( !keepHeld )
            {   // a failed round, or a held phase that wants to form: a solve from the floor cannot grow it back, so go on
                // from the first stop point
                state = stopState;
                Optima::Solver solverOn;
                solverOn.setOptions( options );
                stallWatch->reset();
                result = timed_solve( solverOn, problem, state );
                applyStall( result );
                optimaIterTotal += result.iterations;
                itBack = (long int)result.iterations;
            }
            native_trace_decide( "tailpin rounds=%ld phases=%ld species=%ld ok=%d held=%d tpd=%.3e back=%d it_back=%ld list=%s", rounds, nPh,
                                 (long)pinned.size(), pinOk ? 1 : 0, keepHeld ? 1 : 0, pinTpd, keepHeld ? 0 : 1, itBack, names.c_str() );
        }
        // pa_OptimaTailStop: take the stop point as it stands when no phase wants to form or go.
        // The finish is not run first: it can drop a genuine trace phase (a carrier of a trace element).
        if( stallWatch->tailStop && !result.succeeded && R == 0 )
        {
            for( long int j = 0; j < L; j++ ) { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
            for( long int i = 0; i < N; i++ ) pm.U[i] = -state.ye[i];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            double tViol = 0.; bool tAbsent = false;
            const long int tWorst = WorstPhaseStabilityViolation( presenceThreshold, dcFloor,
                        extinctFixed.empty() ? nullptr : extinctFixed.data(), tViol, tAbsent, nullptr, nullptr, "accept" );
            // Every absent mixed solid must also pass the composition search: a near-degenerate swap
            // (a solid solution replaced by a pure phase of almost the same energy) passes the check above.
            double tTpd = 1e300; long int tTpdPh = -1, tNtpd = 0;
            if( tWorst < 0 )
            {
                long int jb = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int n = pm.L1[k];
                    const char ph = pm.PHC[k];
                    bool allLo = n > 1;
                    for( long int j = jb; j < jb + n && allLo; j++ )
                        if( pm.X[j] > problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) ) allLo = false;
                    if( allLo && ph != PH_AQUEL && ph != PH_GASMIX && ph != PH_PLASMA && ph != PH_FLUID
                        && ph != PH_SORPTION && ph != PH_POLYEL && ph != PH_ADSORPT && ph != PH_IONEX )
                    {
                        std::vector<double> yb;
                        const double tpd = NativeTpdPhase( k, jb, yb );
                        tNtpd++;
                        if( tpd < tTpd ) { tTpd = tpd; tTpdPh = k; }
                    }
                    jb += n;
                }
            }
            // pa_OptimaTailDryOnly 2: the vanished phases that still held the error are searched too, whatever their
            // members' amounts; one that wants to form makes the stop point no answer.
            if( tWorst < 0 )
                for( long int k : stallWatch->tailPhases )
                {
                    long int jb = 0;
                    for( long int q = 0; q < k; q++ ) jb += pm.L1[q];
                    std::vector<double> yb;
                    const double tpd = NativeTpdPhase( k, jb, yb );
                    tNtpd++;
                    if( tpd < tTpd ) { tTpd = tpd; tTpdPh = k; }
                }
            const bool tTpdOk = !( tTpd < -std::max( pa_p->OptimaTpdAccept, 1e-6 ) );
            native_trace_decide( "tailaccept stable=%d worst=%ld viol=%.3e absent=%d ntpd=%ld min_tpd=%.3e tpd_phase=%s",
                                 tWorst < 0 ? 1 : 0, tWorst, tViol, tAbsent ? 1 : 0, tNtpd, tTpd,
                                 tTpdPh >= 0 ? char_array_to_string( pm.SF[tTpdPh], MAXPHNAME ).c_str() : "-" );
            if( tWorst < 0 && tTpdOk ) { result.succeeded = true; tailAccepted = true; }
            ipm_logger->info( "Optima: stopped on a flat tail; the stop point was {} (pa_OptimaTailStop)",
                               tailAccepted ? "checked and taken as the answer" : "rejected by the stability check" );
        }
        // pa_OptimaDryOut: a dry solve that failed (the tail-stop check included) must not cost the
        // answer: put the water phase back and go on with the ordinary solve from where the dry-out stopped it.
        if( stallWatch->dryDone && !result.succeeded && dryStopSaved )
        {
            for( long int j = dryJ0; j < dryJ1; j++ )
            {
                problem.xlower[j] = extLo[(size_t)j]; problem.xupper[j] = extHi[(size_t)j];
                extinctFixed[(size_t)j] = 0;
            }
            if( dryZzSet ) problem.be[N-1] = dryZzSaved;
            state = *dryStopSaved;
            Optima::Solver solverBack;
            solverBack.setOptions( options );
            stallWatch->reset();
            result = timed_solve( solverBack, problem, state );
            applyStall( result );
            optimaIterTotal += result.iterations;
            native_trace_decide( "dryback ok=%d it=%ld", result.succeeded ? 1 : 0, (long)result.iterations );
            ipm_logger->info( "Optima: the dry solve failed; solved again with the water back, ok={} (pa_OptimaDryOut)",
                               result.succeeded ? 1 : 0 );
        }
        // pa_OptimaDryOut: after a dry solve, check the dry answer against a solve that keeps the water.
        if( stallWatch->dryDone && result.succeeded && pm.LO >= 0 )
        {
            long int j0 = 0, j1 = 0, kAq = -1;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                if( pm.LO >= j0 && pm.LO < j0 + pm.L1[k] ) { j1 = j0 + pm.L1[k]; kAq = k; break; }
                j0 += pm.L1[k];
            }
            for( long int j = 0; j < L; j++ ) { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
            for( long int i = 0; i < N; i++ ) pm.U[i] = -state.ye[i];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            std::vector<char> exempt( extinctFixed );
            for( long int j = j0; j < j1; j++ ) exempt[(size_t)j] = 0;
            double dViol = 0.; bool dAbsent = false;
            std::vector<PhStabViolation> dRanked;
            WorstPhaseStabilityViolation( presenceThreshold, dcFloor, exempt.data(), dViol, dAbsent, &dRanked, nullptr, "dryout" );
            double aqViol = 0.; bool aqPhaseWants = false;
            for( const auto& v : dRanked )
                if( v.k == kAq && v.wasAbsent && !v.clamped ) { aqPhaseWants = true; aqViol = v.viol; }
            // Recorded only: whether pure water would form at the dry answer's element potentials. Neither this nor the
            // phase check above decides: a brine can form where pure water would not, and with no water the activity
            // model is not meaningful.
            double gWater = pm.G0[pm.LO] + ( pm.fDQF ? pm.fDQF[pm.LO] : 0. );
            for( long int i = 0; i < N; i++ ) gWater -= pm.U[i] * pm.A[ i + pm.LO*N ];
            native_trace_decide( "drycheck gwater=%.3e phase_wanted=%d viol=%.3e", gWater, aqPhaseWants ? 1 : 0, aqViol );
            // The dry answer is compared with a solve that keeps the water, started where the water ran out;
            // the lower energy wins, and a failed solve with water leaves the dry answer.
            auto energyOf = [&]( const Optima::State& st ) -> double {
                for( long int j = 0; j < L; j++ ) pm.X[j] = st.x[j];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                double g = 0.;
                for( long int j = 0; j < L; j++ ) g += st.x[j] * pm.F[j];
                return g;
            };
            if( dryStopSaved )
            {
                const Optima::State dryState = state;
                const Optima::Result dryResult = result;
                const bool dryTail = tailAccepted;
                const double dryBe = problem.be[N-1];
                const double gDry = energyOf( dryState );
                std::vector<double> dryLo( problem.xlower.data() + j0, problem.xlower.data() + j1 );
                std::vector<double> dryHi( problem.xupper.data() + j0, problem.xupper.data() + j1 );
                std::vector<char> dryExt( extinctFixed.begin() + j0, extinctFixed.begin() + j1 );
                for( long int j = j0; j < j1; j++ )
                {
                    problem.xlower[j] = extLo[(size_t)j]; problem.xupper[j] = extHi[(size_t)j];
                    extinctFixed[(size_t)j] = 0;
                }
                if( dryZzSet ) problem.be[N-1] = dryZzSaved;
                state = *dryStopSaved;
                Optima::Solver solverWet;
                solverWet.setOptions( options );
                stallWatch->reset();
                result = timed_solve( solverWet, problem, state );
                applyStall( result );
                optimaIterTotal += result.iterations;
                tailAccepted = false;
                const double gWet = result.succeeded ? energyOf( state ) : 1e300;
                const bool wetWins = result.succeeded && gWet < gDry - 1e-9 * ( 1. + std::fabs( gDry ) );
                native_trace_decide( "drywet ok=%d it=%ld g_dry=%.12g g_wet=%.12g chosen=%s", result.succeeded ? 1 : 0,
                                     (long)result.iterations, gDry, result.succeeded ? gWet : 0., wetWins ? "wet" : "dry" );
                ipm_logger->info( "Optima: dry answer compared with a solve keeping the water: the {} answer is kept (pa_OptimaDryOut)",
                                   wetWins ? "wet" : "dry" );
                if( !wetWins )
                {
                    for( long int j = j0; j < j1; j++ )
                    {
                        problem.xlower[j] = dryLo[(size_t)(j-j0)]; problem.xupper[j] = dryHi[(size_t)(j-j0)];
                        extinctFixed[(size_t)j] = dryExt[(size_t)(j-j0)];
                    }
                    problem.be[N-1] = dryBe;
                    state = dryState; result = dryResult; tailAccepted = dryTail;
                }
            }
        }

        // The phase-extinction tier, as a re-runnable unit: assembled below for AOP/SOP and
        // left empty for ROP. It must be callable a second time, by the early-probe safety
        // net, on the full-budget re-solve's state - the net is the last step, so without
        // this nothing would be left to rescue a re-solve that stalls the same way the
        // unarmed primary solve does.
        std::function<void(bool)> runExtinctionTier;

        if( referenceMode )
        {
            // ROP retry order: toggle retry first, then the solvent retry. The problem is
            // non-convex, so a different retry order means different starting points and can
            // reach a different (possibly wrong) stationary point; do not change the order
            // without re-validating the full test set. Each retry uses a fresh Optima::Solver.
            auto tryToggleRetry = [&]()
            {
                // Reference fallback: toggle backtracksearch.apply_min_max_fix_and_accept and
                // re-solve from the original state.
                Optima::Options toggled = options;
                toggled.backtracksearch.apply_min_max_fix_and_accept =
                    !toggled.backtracksearch.apply_min_max_fix_and_accept;
                Optima::Solver solverToggle;
                solverToggle.setOptions( toggled );
                state = initialState;
                stallWatch->reset();
                result = timed_solve( solverToggle, problem, state );
                applyStall( result );
                optimaIterTotal += result.iterations;
            };
            // Detects the "aqueous solvent collapsed to its floor" signature on the current
            // state. Inline rather than DetectSolventCollapseAndReseed(), which computes no
            // candidate once the solvent dominates.
            auto solventTrapped = [&]() -> bool
            {
                if( pm.LO < 0 || pm.LO >= L ) return false;
                long int j0 = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int j1 = j0 + pm.L1[k];
                    if( pm.LO >= j0 && pm.LO < j1 )
                    {
                        double otherTotal = 0.;
                        for( long int j = j0; j < j1; j++ )
                            if( j != pm.LO ) otherTotal += state.x[j];
                        return state.x[pm.LO] < otherTotal;
                    }
                    j0 = j1;
                }
                return false;
            };
            // Solvent retry for ROP. Returns whether it ran a solve (false when the current
            // state does not warrant one, e.g. water already exceeds the computed bound).
            auto tryWaterRetry = [&]() -> bool
            {
                if( pm.LO < 0 || pm.LO >= L ) return false;
                long int j0 = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    const long int j1 = j0 + pm.L1[k];
                    if( pm.LO >= j0 && pm.LO < j1 )
                    {
                        double otherTotal = 0.;
                        for( long int j = j0; j < j1; j++ )
                            if( j != pm.LO ) otherTotal += state.x[j];
                        double waterSeed = otherTotal;
                        const long int xH = ( node1 != nullptr ) ? node1->IC_name_to_xCH( "H" ) : -1;
                        const long int xO = ( node1 != nullptr ) ? node1->IC_name_to_xCH( "O" ) : -1;
                        if( xH >= 0 && xO >= 0 )
                        {
                            const double bulkEstimate = std::min( pm.B[xH] / 2.0, pm.B[xO] );
                            if( bulkEstimate > waterSeed )
                                waterSeed = bulkEstimate;
                        }
                        if( waterSeed <= initialState.x[pm.LO] )
                            return false;
                        Optima::State waterState = initialState;
                        waterState.x[pm.LO] = waterSeed;
                        Optima::Solver solverWater;
                        solverWater.setOptions( options );
                        stallWatch->reset();
                        Optima::Result waterResult = timed_solve( solverWater, problem, waterState );
                        applyStall( waterResult );
                        optimaIterTotal += waterResult.iterations;
                        state = waterState;
                        result = waterResult;
                        return true;
                    }
                    j0 = j1;
                }
                return false;
            };

            if( genuineFailure() )
                tryToggleRetry();
            // Trigger on either `!result.succeeded` or the trap signature. The trap-signature
            // arm also runs after an early stop: a collapsed solvent is a direct observation of
            // the state, not an inference from a failure.
            if( genuineFailure() || solventTrapped() )
                tryWaterRetry();

            // General phase-collapse check for every multicomponent phase, in addition to the
            // aqueous-specific retry above.
            std::vector<std::pair<long int,double>> reseeds;
            long int aqueousPhaseIdx = -1;
            {
                long int j0b = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    if( pm.LO >= j0b && pm.LO < j0b + pm.L1[k] ) { aqueousPhaseIdx = k; break; }
                    j0b += pm.L1[k];
                }
            }
            if( DetectPhaseCollapseAndReseed( state.x.data(), reseeds, aqueousPhaseIdx ) )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: general phase-collapse retry (ROP) - "
                                   "{} species reseeded", reseeds.size() );
                Optima::State retryState = state; // build on whatever the prior attempt(s) left, not the original seed
                for( const auto& rs : reseeds )
                    retryState.x[rs.first] = rs.second;
                Optima::Solver solverPhase; // fresh instance - see the note above tryToggleRetry's own solver
                solverPhase.setOptions( options );
                stallWatch->reset();
                Optima::Result retryResult = timed_solve( solverPhase, problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                state = retryState;
                result = retryResult;
            }
        }
        else
        {
        // Solvent-seed retry: when the aqueous solvent (pm.LO) does not dominate its own
        // phase by mass in the final state. An LP-vertex cold start can put the solvent at its
        // floor and the H/O mass into trace redox species; the solvent's ideal-mixing Hessian
        // diagonal (Xf-Xw)/Xw^2 then diverges and Newton cannot grow it back. Such a trapped
        // solve can converge cleanly at the wrong answer, so the retry is gated on the
        // solvent's dominance, not on result.succeeded; as a retry rather than an up-front
        // reseed, projects that recover on their own pay nothing.
        // In plain words: if the first guess left almost no water, put the water back and
        // solve again.
        {
            double waterSeed = 0.;
            // Guard against a no-op retry: if the reseed value does not exceed the original
            // seed's solvent amount, a re-run would reproduce the same trajectory.
            // pa_OptimaWaterKeep = 2: a dry answer just certified by the tail-stop check (phase
            // stability and composition search) is not re-solved with water put back.
            const bool dryCertified = pa_p->OptimaWaterKeep >= 2 && tailAccepted && result.succeeded;
            if( dryCertified ) native_trace_decide( "waterskip certified=1 water=%.3e", state.x[pm.LO] );
            if( !dryCertified && DetectSolventCollapseAndReseed( state.x.data(), problem.xupper[pm.LO], waterSeed )
                && waterSeed > initialState.x[pm.LO] )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: solvent-collapse retry - "
                                   "Xw={}, reseeding to {}", state.x[pm.LO], waterSeed );
                Optima::State retryState( dims );
                for( long int j = 0; j < L; j++ )
                    retryState.x[j] = std::max( pm.Y[j], problem.xlower[j] );
                for( long int rk = 0; rk < R; rk++ )
                    retryState.x[L+rk] = 0.0;
                retryState.x[pm.LO] = waterSeed;
                Optima::Solver solverWater; // fresh instance - see the note above tryToggleRetry's own solver
                solverWater.setOptions( options );
                stallWatch->reset();
                Optima::Result retryResult = timed_solve( solverWater, problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                // pa_OptimaWaterKeep: a converged answer is kept when the water retry fails;
                // a dry answer can be the right one.
                if( pa_p->OptimaWaterKeep > 0 && result.succeeded && !retryResult.succeeded )
                {
                    native_trace_decide( "waterkeep water=%.3e retry_it=%ld", state.x[pm.LO], (long)retryResult.iterations );
                    ipm_logger->info( "Optima: the retry with water put back failed; the converged answer is kept (pa_OptimaWaterKeep)" );
                }
                else
                {
                    state = retryState;
                    result = retryResult;
                }
            }

            // General phase-collapse retry: the "phase reported absent" cliff in
            // PrimalChemicalPotentials() (YF[k] <= pm.DSM) applies to every multicomponent
            // phase. Checked on the state the solvent retry left, excluding the aqueous phase.
            std::vector<std::pair<long int,double>> reseeds;
            long int aqueousPhaseIdx = -1;
            {
                long int j0 = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    if( pm.LO >= j0 && pm.LO < j0 + pm.L1[k] ) { aqueousPhaseIdx = k; break; }
                    j0 += pm.L1[k];
                }
            }
            if( DetectPhaseCollapseAndReseed( state.x.data(), reseeds, aqueousPhaseIdx ) )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: general phase-collapse retry - "
                                   "{} species reseeded", reseeds.size() );
                // Base this retry on the current state (including the solvent retry's fix), not
                // the original seed.
                Optima::State retryState = state;
                for( const auto& rs : reseeds )
                    retryState.x[rs.first] = rs.second;
                Optima::Solver solverPhase; // fresh instance - see the note above tryToggleRetry's own solver
                solverPhase.setOptions( options );
                stallWatch->reset();
                Optima::Result retryResult = timed_solve( solverPhase, problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                state = retryState;
                result = retryResult;
            }

            // pa_OptimaDryRetry: a failed solve that left almost no liquid water, while the bulk
            // could form much more, is solved once more from the water-first start; kept only if
            // it converges.
            if( pa_p->OptimaDryRetry && !g_optimaDryOff && !result.succeeded && pm.LO >= 0 && pm.LO < L && node1 != nullptr )
            {
                const long int xH = node1->IC_name_to_xCH( "H" );
                const long int xO = node1->IC_name_to_xCH( "O" );
                const double wCap = ( xH >= 0 && xO >= 0 ) ? std::min( pm.B[xH] / 2., pm.B[xO] ) : 0.;
                const double xw = state.x[pm.LO];
                std::vector<double> seed;
                if( wCap > 0. && xw < 1e-3 * wCap && waterFirstSeed( seed ) )
                {
                    long int aqIdx = -1;
                    {
                        long int j0d = 0;
                        for( long int k = 0; k < pm.FIs; k++ )
                        {
                            if( pm.LO >= j0d && pm.LO < j0d + pm.L1[k] ) { aqIdx = k; break; }
                            j0d += pm.L1[k];
                        }
                    }
                    std::vector<std::pair<long int,double>> dryReseeds;
                    if( DetectPhaseCollapseAndReseed( seed.data(), dryReseeds, aqIdx, true ) )
                        for( const auto& rs : dryReseeds ) seed[(size_t)rs.first] = rs.second;
                    Optima::State retryState( dims );
                    for( long int j = 0; j < L; j++ )
                        retryState.x[j] = std::min( std::max( seed[(size_t)j], problem.xlower[j] ), problem.xupper[j] );
                    for( long int rk = 0; rk < R; rk++ )
                        retryState.x[L+rk] = 0.0;
                    const bool dryStopped = stallWatch->dryStop;
                    Optima::Solver solverDry;   // fresh instance, as for the retries above
                    solverDry.setOptions( options );
                    stallWatch->reset();
                    Optima::Result dryResult = timed_solve( solverDry, problem, retryState );
                    applyStall( dryResult );
                    optimaIterTotal += dryResult.iterations;
                    native_trace_decide( "dryretry water=%.6g cap=%.6g drystop=%d seedwater=%.6g ok=%d it=%ld",
                                         xw, wCap, dryStopped ? 1 : 0, seed[(size_t)pm.LO],
                                         dryResult.succeeded ? 1 : 0, (long)dryResult.iterations );
                    if( dryResult.succeeded )
                    {
                        state = retryState;
                        result = dryResult;
                    }
                    else if( dryStopped )
                    {   // end this call now; GEM_run() runs it again with the dry stop off.
                        // Only inside that rerun wrapper, which puts the budget back afterwards.
                        g_optimaDryFailed = true;
                        if( g_optimaDryRedoCall )
                        {
                            g_optimaIterBudget = 0;
                            stallWatch->timedOut = true;
                        }
                    }
                }
            }

            // Fourth tier: phase-extinction deactivation. Fires only when everything above has
            // not converged. Target: a system with two solution phases of identical
            // end-members and G0, one of which is redundant. The losing one decays only
            // geometrically toward its bound, and its floor end-members carry large reduced
            // gradients that change sign between iterations; Optima masks a variable's
            // optimality error only when x == xlower and s > 0, so the error toggles and never
            // falls below tolerance while the composition is frozen. Fixing the box
            // (xlower == xupper) masks the variable for either sign of s, which ends the
            // toggle. Only problem.xlower/xupper are fixed, not pm.DUL/pm.DLL, so the final
            // phase-stability check still runs on the deactivated phase and rejects the result
            // if the phase should in fact be present.
            // In plain words: when two identical mixed phases compete, remove the one that is
            // disappearing anyway.
            // This tier also runs after an early stop, with the extra "still falling"
            // clause below.
            runExtinctionTier = [&]( bool onEarlyStopPath )
            {
                // The aqueous phase index is recomputed here, because this lambda is also
                // called after the enclosing block's variable has gone out of scope.
                long int aqueousPhaseIdx = -1;
                {
                    long int j0aq = 0;
                    for( long int k = 0; k < pm.FIs; k++ )
                    {
                        if( pm.LO >= j0aq && pm.LO < j0aq + pm.L1[k] ) { aqueousPhaseIdx = k; break; }
                        j0aq += pm.L1[k];
                    }
                }
                // Detection is structural, not a magnitude threshold: two solution phases with
                // identical stoichiometry and identical standard potentials are
                // interchangeable, so a state with one vanishing is the same state as one with
                // it removed. A redundant twin cannot be found by its total (it is heading for
                // the floor, not there yet) or by its logSI (twins share a tangent plane, so
                // logSI ~ 0). To keep a genuine two-limb split of such a pair, one twin must
                // also be smaller than the other by kTwinRatio.
                const double kTwinRatio = 1.0e3;
                std::vector<long int> deactivated;
                std::vector<long int> ph0, ph1;
                std::vector<double> phTot;
                {
                    long int j0 = 0;
                    for( long int k = 0; k < pm.FIs; k++ )
                    {
                        const long int j1 = j0 + pm.L1[k];
                        double tot = 0.;
                        if( j1 <= L )
                            for( long int j = j0; j < j1; j++ ) tot += std::max( state.x[j], 0. );
                        ph0.push_back( j0 ); ph1.push_back( j1 ); phTot.push_back( tot );
                        j0 = j1;
                    }
                }
                auto interchangeable = [&]( long int ka, long int kb ) -> bool
                {
                    if( pm.L1[ka] != pm.L1[kb] ) return false;
                    if( ph1[(size_t)ka] > L || ph1[(size_t)kb] > L ) return false;
                    for( long int d = 0; d < pm.L1[ka]; d++ )
                    {
                        const long int ja = ph0[(size_t)ka] + d, jb = ph0[(size_t)kb] + d;
                        if( std::fabs( pm.G0[ja] - pm.G0[jb] ) >
                            1e-9 * std::max( 1., std::fabs( pm.G0[ja] ) ) ) return false;
                        for( long int i = 0; i < N; i++ )
                            if( pm.A[ i + ja*N ] != pm.A[ i + jb*N ] ) return false;
                    }
                    return true;
                };
                for( long int ka = 0; ka < pm.FIs; ka++ )
                {
                    if( ka == aqueousPhaseIdx || pm.L1[ka] <= 1 || ph1[(size_t)ka] > L ) continue;
                    bool alreadyFixed = true;
                    for( long int j = ph0[(size_t)ka]; j < ph1[(size_t)ka]; j++ )
                        if( problem.xupper[j] > problem.xlower[j] ) { alreadyFixed = false; break; }
                    if( alreadyFixed ) continue;
                    for( long int kb = 0; kb < pm.FIs; kb++ )
                    {
                        if( kb == ka || kb == aqueousPhaseIdx || pm.L1[kb] <= 1 ) continue;
                        if( phTot[(size_t)kb] <= phTot[(size_t)ka] * kTwinRatio ) continue;
                        if( !interchangeable( ka, kb ) ) continue;
                        // "Still falling" clause, on the early-stop path only: the twin must
                        // have decreased at the last evaluation. On a state truncated by
                        // pa_OptimaEarlyStabilityAt the ratio alone cannot tell a vanishing
                        // twin from a genuine two-limb split that is still rebounding.
                        if( onEarlyStopPath && (*phDec)[(size_t)ka] <= 0 )
                        {
                            ipm_logger->debug( "CalculateEquilibriumStateOptima: phase {} looks like a "
                                              "vanishing twin of {} (total {:.3e} against {:.3e}) but "
                                              "is NOT falling - declining to deactivate it on a state "
                                              "the early-stability probe truncated (peak {:.3e})",
                                              ka, kb, phTot[(size_t)ka], phTot[(size_t)kb],
                                              (*phMax)[(size_t)ka] );
                            continue;
                        }
                        for( long int j = ph0[(size_t)ka]; j < ph1[(size_t)ka]; j++ )
                            deactivated.push_back( j );
                        break;
                    }
                }
                if( !deactivated.empty() )
                {
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-extinction retry - "
                                       "fixing {} species of a vanishing interchangeable phase at the floor",
                                       deactivated.size() );
                    // The kTwinRatio clause firing.
                    native_trace_decide( "extinction fixedsp=%ld",
                                         (long)deactivated.size() );
                    Optima::State retryState = state;
                    for( long int j : deactivated )
                    {
                        if( !extinctFixed[(size_t)j] ) { extLo[(size_t)j] = problem.xlower[j]; extHi[(size_t)j] = problem.xupper[j]; }
                        problem.xlower[j] = dcFloor;
                        problem.xupper[j] = dcFloor;
                        retryState.x[j]   = dcFloor;
                        extinctFixed[(size_t)j] = 1;
                    }
                    Optima::Solver solverExt; // fresh instance - see the note above tryToggleRetry's own solver
                    solverExt.setOptions( options );
                    stallWatch->reset();
                    // pa_OptimaCallCap: when the cap is used up, this repair still gets a small allowance,
                    // since it often finishes in a few iterations.
                    if( g_optimaIterBudget == 0 && g_optimaCapAllowance > 0 && !g_optimaDryFailed )
                    {
                        g_optimaIterBudget = g_optimaCapAllowance;
                        g_optimaCapAllowance = -1;
                        stallWatch->timedOut = false;
                        native_trace_decide( "extallow budget=%ld", g_optimaIterBudget );
                    }
                    Optima::Result retryResult = timed_solve( solverExt, problem, retryState );
                    applyStall( retryResult );
                    optimaIterTotal += retryResult.iterations;
                    state = retryState;
                    result = retryResult;
                }
            };
            if( genuineFailure() || earlyCapHit )
                runExtinctionTier( earlyCapHit );
            // Fifth tier: phase-selection loop - an action attached to the phase-stability
            // detector. Native drops absent phases on every call (PSSC, pa_PC = 2); this
            // solver otherwise carries every absent phase as a live unknown, and a floor
            // species with a negative reduced gradient keeps the optimality error up.
            // Iterative and post-solve (a stability index needs a dual), bounded at
            // kMaxPhaseSelectLoops = 5 like native's PSSC. Each phase moves at most
            // untouched -> deactivated -> readmitted, so the loop cannot oscillate; a wrongly
            // dropped phase is seen next pass as "absent but stable" and put back, and a
            // readmitted phase that still disagrees is left to the final check.
            // Deactivation fixes only problem.xlower/xupper at dcFloor, never pm.DUL/pm.DLL:
            // Optima then masks the variable, while the final stability check still sees it.
            // In plain words: like the original solver, remove phases that should not be
            // there and put back ones that were removed by mistake.
            {
                const int kMaxPhaseSelectLoops = 5;
                // Same threshold as the final check, so this loop never chases a violation
                // that check would not report.
                for( int psLoop = 0; psLoop < kMaxPhaseSelectLoops; psLoop++ )
                {
                    // Evaluate the assemblage on the current state through the same
                    // pipeline the final check uses (pm.Gamma, pm.Y_la, pm.YF/YFA).
                    for( long int j = 0; j < L; j++ )
                    {
                        pm.Y[j] = state.x[j];
                        pm.X[j] = state.x[j];
                    }
                    for( long int i = 0; i < N; i++ )
                        pm.U[i] = -state.ye[i];
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    CalculateConcentrations( pm.X, pm.XF, pm.XFA );

                    double psViol = 0.;
                    bool psAbsent = false;
                    std::vector<PhStabViolation> psRanked;
                    PhStabCensus psCensus;
                    WorstPhaseStabilityViolation(
                                presenceThreshold, dcFloor,
                                extinctFixed.empty() ? nullptr : extinctFixed.data(),
                                psViol, psAbsent, &psRanked, &psCensus, "phasesel" );
                    if( psRanked.empty() )
                        break;                       // assemblage is self-consistent

                    // Walk the ranking rather than stopping at its head: an unactionable
                    // violation (the solvent phase, or a phSelState this tier may not touch)
                    // says nothing about those ranked behind it. Acting on any entry
                    // invalidates the state, so the loop re-solves and re-evaluates from the
                    // top. A warning is issued when most scanned phases have a saturated
                    // stability index, since the ranking is then not informative.
                    if( psCensus.clamped * 2 > psCensus.scanned && psCensus.scanned > 0 )
                        ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                           "{} of {} scanned phases have a SATURATED stability index "
                                           "(overflow guard, not a measured driving force); the "
                                           "ranking of {} violation(s) is not ordered by anything",
                                           psLoop, psCensus.clamped, psCensus.scanned,
                                           (long)psRanked.size() );

                    bool acted = false;
                    long int psSkipped = 0;
                    for( size_t psIdx = 0; psIdx < psRanked.size() && !acted; psIdx++ )
                    {
                        const long int kBad = psRanked[psIdx].k;
                        psAbsent = psRanked[psIdx].wasAbsent;
                        psViol = psRanked[psIdx].viol;
                        if( kBad == aqueousPhaseIdx )
                        { psSkipped++; continue; }       // never remove or reseed the solvent phase here

                        long int jb = 0;
                        for( long int k = 0; k < kBad; k++ )
                            jb += pm.L1[k];
                        const long int je = jb + pm.L1[kBad];
                        if( je > L )
                        { psSkipped++; continue; }

                        if( !psAbsent && phSelState[(size_t)kBad] == 0 )
                        {
                            // Present but unstable: drop it from the assemblage.
                            for( long int j = jb; j < je; j++ )
                            {
                                savedLo[(size_t)j] = problem.xlower[j];
                                savedHi[(size_t)j] = problem.xupper[j];
                                problem.xlower[j] = dcFloor;
                                problem.xupper[j] = dcFloor;
                                state.x[j] = dcFloor;
                            }
                            phSelState[(size_t)kBad] = 1;
                            acted = true;
                            ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                               "deactivating phase {} (present but unstable, logSI gap {})",
                                               psLoop, kBad, psViol );
                            native_trace_decide( "phasesel loop=%ld deactivate=%ld logsigap=%.6e",
                                                 (long)psLoop, (long)kBad, psViol );
                        }
                        else if( psAbsent && phSelState[(size_t)kBad] == 1 )
                        {
                            // Absent but stable, and deactivated by the compaction probe or by
                            // this loop: that removal was wrong. Restore the original box and
                            // seed the phase off the floor, each end-member bounded by its own
                            // limiting IC (as in DetectPhaseCollapseAndReseed()). Every phase
                            // that now looks wrongly removed is readmitted at once;
                            // deactivation stays one at a time.
                            long int nReadmit = 0, jbR = 0;
                            for( long int k = 0; k < pm.FI; k++ )
                            {
                                const long int jeR = jbR + pm.L1[k];
                                if( phSelState[(size_t)k] == 1 && jeR <= L
                                    && pm.YF[k] < presenceThreshold && pm.Falp[k] > pa_p->DF )
                                {
                                    const double nEnd = (double)( jeR - jbR );
                                    for( long int jj = jbR; jj < jeR; jj++ )
                                    {
                                        problem.xlower[jj] = savedLo[(size_t)jj];
                                        problem.xupper[jj] = savedHi[(size_t)jj];
                                        double b = -1.;
                                        for( long int i = 0; i < N; i++ )
                                        {
                                            const double coef = pm.A[ i + jj*N ];
                                            if( coef > 0. )
                                            {
                                                const double icBound = pm.B[i] / coef;
                                                if( b < 0. || icBound < b ) b = icBound;
                                            }
                                        }
                                        const double seed = ( b > 0. ? b : presenceThreshold * 10. ) / nEnd;
                                        state.x[jj] = std::min( std::max( seed, problem.xlower[jj] ), problem.xupper[jj] );
                                    }
                                    phSelState[(size_t)k] = 2;   // terminal - never dropped again
                                    nReadmit++;
                                }
                                jbR = jeR;
                            }
                            if( nReadmit > 0 )
                            {
                                acted = true;
                                ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                                   "READMITTING {} wrongly removed phase(s) (worst: phase {}, "
                                                   "absent but stable, logSI gap {})",
                                                   psLoop, nReadmit, kBad, psViol );
                            }
                        }
                        if( !acted )
                            psSkipped++;
                    }   // for psIdx - fall through to the next-worst violation
                    if( acted && psSkipped > 0 )
                        ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                           "skipped {} unactionable violation(s) ranked ahead of the one "
                                           "acted on (worst was phase {}, gap {})",
                                           psLoop, psSkipped, psRanked[0].k, psRanked[0].viol );
                    if( !acted )
                    {
                        // Nothing in the whole ranking is actionable.
                        ipm_logger->debug( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                           "none of {} violation(s) is actionable ({} saturated at the "
                                           "overflow guard); worst is phase {}, gap {} - leaving it to "
                                           "the final check",
                                           psLoop, (long)psRanked.size(), psCensus.clamped,
                                           psRanked[0].k, psRanked[0].viol );
                        break;   // nothing this loop is allowed to do - let the final check report it
                    }

                    Optima::Solver solverPS; // fresh instance, per the retry-ordering note above
                    solverPS.setOptions( options );
                    stallWatch->reset();
                    Optima::State psState = state;
                    Optima::Result psResult = timed_solve( solverPS, problem, psState );
                    applyStall( psResult );
                    optimaIterTotal += psResult.iterations;
                    state = psState;
                    result = psResult;
                }
            }

        }
        } // else (!referenceMode) - AOP/SOP's own solvent-collapse retry

        // pa_OptimaTwinReadmit: a phase the extinction tier switched off is put back when the composition
        // search says it should form: a twin can be a genuine second limb of a miscibility gap.
        if( pa_p->OptimaTwinReadmit && !referenceMode && result.succeeded
            && std::find( extinctFixed.begin(), extinctFixed.end(), 1 ) != extinctFixed.end() )
        {
            for( long int j = 0; j < L; j++ ) { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
            for( long int i = 0; i < N; i++ ) pm.U[i] = -state.ye[i];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            Optima::State trial = state;
            std::vector<long int> readm;
            long int jb = 0;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                const long int je = jb + pm.L1[k];
                bool ext = false;
                for( long int j = jb; j < je && j < L; j++ ) if( extinctFixed[(size_t)j] ) ext = true;
                if( ext )
                {
                    std::vector<double> yb;
                    const double tpd = NativeTpdPhase( k, jb, yb );
                    double amt = 0.;
                    long int kb0 = 0;
                    for( long int kb = 0; kb < pm.FIs; kb++ )
                    {
                        const long int kbe = kb0 + pm.L1[kb];
                        if( kb != k && pm.L1[kb] == pm.L1[k] )
                        {
                            double t = 0.;
                            for( long int j = kb0; j < kbe; j++ ) t += state.x[j];
                            amt = std::max( amt, 0.5 * t );
                        }
                        kb0 = kbe;
                    }
                    const bool form = tpd < -1e-6 && amt > 0. && tpd < 1e299;
                    native_trace_decide( "twinreadmit phase=%s tpd=%.6e amount=%.6e form=%d",
                                         char_array_to_string( pm.SF[k], MAXPHNAME ).c_str(), tpd, amt, form ? 1 : 0 );
                    if( form )
                        for( long int j = jb; j < je; j++ )
                        {
                            problem.xlower[j] = extLo[(size_t)j];
                            problem.xupper[j] = extHi[(size_t)j];
                            trial.x[j] = std::min( std::max( yb[(size_t)(j-jb)] * amt, problem.xlower[j] ), problem.xupper[j] );
                            readm.push_back( j );
                        }
                }
                jb = je;
            }
            if( !readm.empty() )
            {
                Optima::Solver solverTR;
                solverTR.setOptions( options );
                stallWatch->reset();
                Optima::Result trResult = timed_solve( solverTR, problem, trial );
                applyStall( trResult );
                optimaIterTotal += trResult.iterations;
                native_trace_decide( "twinreadmit solved=%d it=%ld", trResult.succeeded ? 1 : 0, (long)trResult.iterations );
                if( trResult.succeeded )
                {
                    for( long int j : readm ) extinctFixed[(size_t)j] = 0;
                    state = trial;
                    result = trResult;
                }
                else
                    for( long int j : readm ) { problem.xlower[j] = dcFloor; problem.xupper[j] = dcFloor; }
            }
        }

        // ---- A net must restore the problem, not only the state ---------------
        // The safety nets below re-solve from `initialState`, which restores the state but
        // not the boxes. Phases pinned at the floor by the compaction probe or the
        // phase-selection loop (phSelState == 1, originals in savedLo/savedHi) are un-pinned
        // for the re-solve; the extinction tier's pins are kept (they are the rescue).
        // If the net does not win, the pins are restored, because the KKT check reads the
        // boxes afterwards.
        std::vector<long int> netUnpinned;
        auto unpinForNet = [&]() -> long int
        {
            netUnpinned.clear();
            long int j0u = 0;
            for( long int k = 0; k < pm.FI; k++ )
            {
                const long int j1u = j0u + pm.L1[k];
                if( phSelState[(size_t)k] == 1 && j1u <= L )
                    for( long int j = j0u; j < j1u; j++ )
                    {
                        problem.xlower[j] = savedLo[(size_t)j];
                        problem.xupper[j] = savedHi[(size_t)j];
                        netUnpinned.push_back( j );
                    }
                j0u = j1u;
            }
            if( !netUnpinned.empty() )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: safety net - restoring the "
                                   "original box of {} species pinned out by a probe-length look",
                                   netUnpinned.size() );
                native_trace_decide( "netunpin sp=%ld", (long)netUnpinned.size() );
            }
            return (long int)netUnpinned.size();
        };
        auto repinAfterNet = [&]()
        {
            for( long int j : netUnpinned )
            {
                problem.xlower[j] = dcFloor;
                problem.xupper[j] = dcFloor;
            }
            netUnpinned.clear();
        };


        // pa_OptimaFinishFirst: a stalled failed solve first tries the Newton finish on its own end point.
        // The finish keeps its phase set, so its answer is taken only if no phase then wants to form;
        // otherwise the re-solves below run as before.
        bool finishFirstFailed = false;
        if( pa_p->OptimaFinishFirst > 0 && pa_p->OptimaFinish > 0 && !result.succeeded && stallWatch->stalled && R == 0 )
        {
            for( long int j = 0; j < L; j++ ) { pm.Y[j] = state.x[j]; pm.X[j] = state.x[j]; }
            for( long int i = 0; i < N; i++ ) pm.U[i] = -state.ye[i];
            std::vector<double> xlo( (size_t)L ), xhi( (size_t)L );
            for( long int j = 0; j < L; j++ ) { xlo[(size_t)j] = problem.xlower[j]; xhi[(size_t)j] = problem.xupper[j]; }
            g_finishFromSuccess = false;
            const bool ok = PotentialSpaceFinish( dcFloor, xlo, xhi );
            long int ffWorst = -1; double ffViol = 0.; bool ffAbsent = false;
            if( ok )
            {
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                ffWorst = WorstPhaseStabilityViolation( presenceThreshold, dcFloor,
                            extinctFixed.empty() ? nullptr : extinctFixed.data(), ffViol, ffAbsent, nullptr, nullptr, "finishfirst" );
            }
            native_trace_decide( "finishfirst ok=%d stable=%d worst=%ld viol=%.3e absent=%d",
                                 ok ? 1 : 0, ffWorst < 0 ? 1 : 0, ffWorst, ffViol, ffAbsent ? 1 : 0 );
            if( ok && ffWorst < 0 )
            {
                for( long int j = 0; j < L; j++ ) state.x[j] = pm.X[j];
                for( long int i = 0; i < N; i++ ) state.ye[i] = -pm.U[i];
                result.succeeded = true;
            }
            else if( !ok )
                finishFirstFailed = true;
        }

        // ---- Last resort: the stall watch's own safety net ----------------------
        // pa_OptimaStallWindow hands control to the retry tiers sooner; it does not decide
        // the verdict. If it fired and no tier recovered, re-solve once with it disarmed,
        // from the primary solve's starting state. Bounded by one ordinary solve, paid only
        // when the watch fired and every retry failed.
        // Gated on `stalled`, not on `timedOut`: pa_OptimaMaxSeconds is a budget the caller
        // set on purpose. Not mode-gated (ROP gets the same net).
        // pa_OptimaColdRetry = 2: a warm call (pm.pNP == 1) skips the full-budget re-solves
        // below and fails fast, so TNode::GEM_run_optima_cold_retry() takes over at once.
        const bool warmFailFast = ( pa_p->OptimaColdRetry == 2 );   // pa_OptimaColdRetry = 2
        const bool skipFullBudgetResolves = warmFailFast && pm.pNP == 1;
        if( skipFullBudgetResolves && !result.succeeded )
            native_trace_decide( "warmfailfast stalled=%d", stallWatch->stalled ? 1 : 0 );
        // When the call will be re-run with the stall escape off (escretry), that attempt carries the net instead.
        const bool netLeftToRetry = g_optimaEscRetryNext && !g_optimaLSEscapeOff;
        if( netLeftToRetry && !result.succeeded && stallWatch->stalled && stallWatch->window > 0 )
            native_trace_decide( "stallnet skipped=escretry" );
        if( !skipFullBudgetResolves && !netLeftToRetry && !result.succeeded && stallWatch->stalled && stallWatch->window > 0 )
        {
            ipm_logger->info( "Optima: all retries failed after a stall. Solving once more with the "
                               "stall guard off (pa_OptimaStallWindow, in the project's -ipm file)." );
            const long int savedWindow = stallWatch->window;
            stallWatch->window = 0;      // the check lambda reads this on every call
            stallWatch->reset();
            unpinForNet();               // the PROBLEM too, not only the state - see above
            Optima::Solver solverNS;     // fresh instance, per the retry-ordering note above
            // pa_OptimaStallNetFactor: the net gets a multiple of the iterations after which the first solve stalled
            Optima::Options netOptions = options;
            if( pa_p->OptimaStallNetFactor > 0 && stallWatch->stallAt > 0 )
                netOptions.maxiters = (unsigned)std::min( (long int)options.maxiters,
                                          std::max( 2000L, pa_p->OptimaStallNetFactor * stallWatch->stallAt ) );
            solverNS.setOptions( netOptions );
            Optima::State  nsState  = initialState;
            Optima::Result nsResult = timed_solve( solverNS, problem, nsState );
            optimaIterTotal += nsResult.iterations;
            native_trace_decide( "stallnet outcome=%s it=%ld cap=%ld stallat=%ld", nsResult.succeeded ? "converged" : "failed",
                                 (long int)nsResult.iterations, (long int)netOptions.maxiters, stallWatch->stallAt );
            stallWatch->window = savedWindow;
            if( !nsResult.succeeded )
                repinAfterNet();         // the net did not win - leave the boxes as they were
            if( nsResult.succeeded )
            {
                netUnpinned.clear();     // this state was solved on the restored boxes; keep them
                state  = nsState;
                result = nsResult;
                ipm_logger->debug( "CalculateEquilibriumStateOptima: the disarmed re-solve converged"
                                   " in {} iterations - the stall was a false positive",
                                   nsResult.iterations );
            }
        }

        // ---- Last resort: pa_OptimaEarlyStabilityAt's safety net ----------------
        // The capped first attempt is a probe. If the cap is what stopped it and nothing
        // downstream recovered, re-solve once at the full budget, so a probe that found
        // nothing costs only itself.
        if( !skipFullBudgetResolves && !result.succeeded && earlyCapHit )
        {
            ipm_logger->debug( "CalculateEquilibriumStateOptima: the pa_OptimaEarlyStabilityAt probe"
                               " found nothing to repair - re-solving once at the full budget" );
            Optima::Solver solverEP;     // fresh instance, per the retry-ordering note above
            solverEP.setOptions( options );   // the full budget; the cap is disarmed by now
            // The stall watch stays armed here. Optima reports a stop requested by the
            // convergence hook as success, so a stalled re-solve is folded to a failure below.
            stallWatch->reset();
            unpinForNet();               // the PROBLEM too, not only the state - see above
            // Resume from the probe's end state (*earlyProbeState) or restart from
            // `initialState`, per optima_net_resume_mode(); recorded in the trace.
            const bool epResumed = (bool)earlyProbeState;
            Optima::State  epState  = epResumed ? *earlyProbeState : initialState;
            // dx = max|x_probe - x_0| over the primal; 0 would mean the start point never
            // changed.
            double epStartDx = 0.;
            for( long int j = 0; j < L; j++ )
                epStartDx = std::max( epStartDx, std::fabs( epState.x[j] - initialState.x[j] ) );
            native_trace_decide( "netresume from=%s probeat=%ld dx=%.6e",
                                 epResumed ? "probe" : "initial", earlyProbeIters, epStartDx );
            Optima::Result epResult = timed_solve( solverEP, problem, epState );
            optimaIterTotal += epResult.iterations;
            const bool epStalled = ( stallWatch->stalled || stallWatch->timedOut );
            if( epStalled ) epResult.succeeded = false;
            // What stopped the re-solve (converged / stalled / failed) is traced: a
            // stalled re-solve spent a budget, so its start point does not change its cost.
            native_trace_decide( "netresolve outcome=%s it=%ld",
                                 epResult.succeeded ? "converged"
                                                    : ( epStalled ? "stalled" : "failed" ),
                                 (long int)epResult.iterations );
            // Mode 2 fallback: a resume that did not converge is followed by one ordinary
            // restart from `initialState`, with a fresh solver and a reset stall watch.
            if( netResumeMode == 2 && epResumed && !epResult.succeeded )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: the resumed re-solve did not"
                                   " converge in {} iterations - falling back to the restart",
                                   epResult.iterations );
                Optima::Solver solverFB;
                solverFB.setOptions( options );
                stallWatch->reset();
                Optima::State  fbState  = initialState;
                Optima::Result fbResult = timed_solve( solverFB, problem, fbState );
                optimaIterTotal += fbResult.iterations;
                const bool fbStalled = ( stallWatch->stalled || stallWatch->timedOut );
                if( fbStalled ) fbResult.succeeded = false;
                native_trace_decide( "netfallback outcome=%s it=%ld",
                                     fbResult.succeeded ? "converged"
                                                        : ( fbStalled ? "stalled" : "failed" ),
                                     (long int)fbResult.iterations );
                // Taken unconditionally: everything below is written against the restart's state.
                epState  = fbState;
                epResult = fbResult;
            }
            if( epResult.succeeded )
            {
                netUnpinned.clear();     // this state was solved on the restored boxes; keep them
                state  = epState;
                result = epResult;
                ipm_logger->debug( "CalculateEquilibriumStateOptima: the full-budget re-solve"
                                   " converged in {} iterations ({} the probe's {} iterations)",
                                   epResult.iterations,
                                   epResumed ? "resumed from" : "restarted, discarding",
                                   earlyProbeIters );
            }
            else if( runExtinctionTier )
            {
                // The re-solve reproduced the unarmed primary solve, including its failure,
                // so it gets the phase-extinction rescue an unarmed call gets, judged on the
                // ordinary path (onEarlyStopPath = false). The pre-net state and result are
                // restored unless the tier converges.
                const Optima::State  savedState  = state;
                const Optima::Result savedResult = result;
                state  = epState;
                result = epResult;
                runExtinctionTier( false );
                if( result.succeeded )
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: the full-budget re-solve"
                                       " stalled at {} iterations and the phase-extinction tier"
                                       " converged it - the probe cost {} iterations and nothing else",
                                       epResult.iterations, earlyCapN > 0 ? earlyCapN : 0 );
                else
                {
                    state  = savedState;
                    result = savedResult;
                    repinAfterNet();     // nothing won - the boxes go back with the state
                }
                if( result.succeeded )
                    netUnpinned.clear();
            }
            else
                repinAfterNet();         // no rescue available and the re-solve failed

            // ---- Mode 3: the deferred fallback --------------------------------
            // Asked only here, after the resumed re-solve and every rescue below it have had
            // their turn, so `!result.succeeded` means the call has lost its answer.
            if( netResumeMode == 3 && epResumed && !result.succeeded )
            {
                ipm_logger->debug( "CalculateEquilibriumStateOptima: the resumed re-solve and its"
                                   " rescue both failed - falling back to the restart the unarmed"
                                   " call would have used" );
                unpinForNet();           // repinAfterNet() may have put the boxes back
                stallWatch->reset();
                Optima::Solver solverFB;
                solverFB.setOptions( options );
                Optima::State  fbState  = initialState;
                Optima::Result fbResult = timed_solve( solverFB, problem, fbState );
                optimaIterTotal += fbResult.iterations;
                const bool fbStalled = ( stallWatch->stalled || stallWatch->timedOut );
                if( fbStalled ) fbResult.succeeded = false;
                native_trace_decide( "netfallback outcome=%s it=%ld",
                                     fbResult.succeeded ? "converged"
                                                        : ( fbStalled ? "stalled" : "failed" ),
                                     (long int)fbResult.iterations );
                if( fbResult.succeeded )
                {
                    netUnpinned.clear();
                    state  = fbState;
                    result = fbResult;
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: the restart converged in {}"
                                       " iterations where the resume did not", fbResult.iterations );
                }
                else if( runExtinctionTier )
                {
                    const Optima::State  savedState2  = state;
                    const Optima::Result savedResult2 = result;
                    state  = fbState;
                    result = fbResult;
                    runExtinctionTier( false );
                    if( result.succeeded )
                        netUnpinned.clear();
                    else
                    {
                        state  = savedState2;
                        result = savedResult2;
                        repinAfterNet();
                    }
                }
                else
                    repinAfterNet();
            }
        }

        ipm_logger->debug( "CalculateEquilibriumStateOptima: pNP={} referenceMode={} succeeded={} iterations={} nConditions={}",
                           pm.pNP, referenceMode, result.succeeded, result.iterations, R );

        for( long int j = 0; j < L; j++ )
        {
            pm.Y[j] = state.x[j];
            pm.X[j] = state.x[j];
        }
        // Optima's multipliers ye have the opposite sign convention from U[].
        for( long int i = 0; i < N; i++ )
            pm.U[i] = -state.ye[i];
        // pa_OptimaTraceRespread: an element all of whose species sit at their floors, and whose floors add up to more than
        // its amount, gets its true amount spread over its solution species, x_j * exp(-s_j + a_ij * d), s_j the reduced gradient, d set so it balances.
        if( pa_p->OptimaTraceRespread > 0 && result.succeeded )
        {
            bool changed = false, potentials = false;
            for( long int i = 0; i < N - pm.E; i++ )
            {
                if( !( pm.B[i] > 0. ) ) continue;
                double tot = 0., fixedPart = 0.;
                bool present = false;   // a species of this element above its floor fixes its potential: leave it alone
                std::vector<long int> fl;
                for( long int j = 0; j < L; j++ )
                {
                    const double a = pm.A[ i + j*N ];
                    if( !( a > 0. ) ) continue;
                    tot += a * pm.X[j];
                    if( pm.X[j] <= problem.xlower[j] * 10. && pm.DLL[j] <= 0. ) fl.push_back( j );
                    else { fixedPart += a * pm.X[j]; present = true; }
                }
                if( present ) continue;
                for( long int k = 0; k < R; k++ ) { const double t = problem.Aex( i, L+k ) * state.x[L+k]; tot += t; fixedPart += t; }
                if( !( tot - pm.B[i] > 1e-6 * pm.B[i] ) || fl.empty() ) continue;   // the balance holds
                const double target = pm.B[i] - fixedPart;
                if( !( target > 0. ) ) continue;
                if( !potentials )   // only when some element needs it
                {
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                    potentials = true;
                }
                std::vector<double> lnx0( fl.size() ), aa( fl.size() );
                std::vector<char> sol( fl.size() );
                for( size_t q = 0; q < fl.size(); q++ )
                {
                    const long int j = fl[q];
                    double mu = 0.;
                    for( long int k = 0; k < N; k++ ) mu += pm.A[ k + j*N ] * pm.U[k];
                    sol[q] = j < pm.Ls;
                    lnx0[q] = std::log( pm.X[j] ) - ( pm.F[j] - mu );
                    aa[q] = pm.A[ i + j*N ];
                }
                auto mass = [&]( double d ) {
                    double m = 0.;
                    for( size_t q = 0; q < fl.size(); q++ )
                        if( sol[q] ) m += aa[q] * std::exp( std::min( lnx0[q] + aa[q] * d, 700. ) );
                    return m; };
                double lo = -2000., hi = 2000.;
                if( !( mass( lo ) < target && mass( hi ) > target ) ) continue;
                for( int it = 0; it < 200; it++ ) { const double mid = 0.5 * ( lo + hi ); ( mass( mid ) < target ? lo : hi ) = mid; }
                const double d = 0.5 * ( lo + hi );
                for( size_t q = 0; q < fl.size(); q++ )
                {
                    const long int j = fl[q];
                    pm.X[j] = sol[q] ? std::exp( lnx0[q] + aa[q] * d ) : 0.;
                    state.x[j] = pm.X[j];
                }
                pm.U[i] += d;
                state.ye[i] = -pm.U[i];
                changed = true;
                native_trace_decide( "tracerespread ic=%ld species=%ld delta=%.6g before=%.3e target=%.3e", (long)i,
                                     (long)fl.size(), d, tot - fixedPart, target );
            }
            if( changed )
                for( long int j = 0; j < L; j++ ) pm.Y[j] = pm.X[j];
        }

        // pa_OptimaFinish: Newton finish on the fixed phase set from Optima's last point; the
        // checks below judge its result like Optima's own.
        bool finishTrace = false;   // a reported success with a free species at trace level
        if( pa_p->OptimaFinish > 0 && result.succeeded && R == 0 )
        {
            double bSc = 0.;
            for( long int i = 0; i < N; i++ ) bSc = std::max( bSc, std::fabs( pm.B[i] ) );
            for( long int j = 0; j < L && !finishTrace; j++ )
                if( pm.X[j] > problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) && pm.X[j] <= 1e-11 * bSc &&
                    problem.xupper[j] > problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) )
                    finishTrace = true;
        }
        // A finish already tried on this same failed end point is not repeated.
        if( pa_p->OptimaFinish > 0 && ( !result.succeeded || finishTrace ) && R == 0
            && !( finishFirstFailed && !result.succeeded ) )
        {
            std::vector<double> xlo( (size_t)L ), xhi( (size_t)L );
            for( long int j = 0; j < L; j++ ) { xlo[(size_t)j] = problem.xlower[j]; xhi[(size_t)j] = problem.xupper[j]; }
            g_finishFromSuccess = result.succeeded;   // the trace case starts from a reported-OK call
            const bool finishOk = PotentialSpaceFinish( dcFloor, xlo, xhi );
            g_finishFromSuccess = false;
            if( finishOk )
            {
                for( long int j = 0; j < L; j++ ) state.x[j] = pm.X[j];
                for( long int i = 0; i < N; i++ ) state.ye[i] = -pm.U[i];
                result.succeeded = true;
            }
        }

        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );

        // Post-solve concentration/output pipeline, as native's CalculateEquilibriumState()
        // calls it: packDataBr() reads pmm->VXc/FVOL[]/FWGT[]/IC, which only
        // CalculateConcentrations() fills (it also sets pm.pH/pe/Eh). Uses the final pm.X and pm.U.
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );

        // General KKT (complementary slackness) check for every ordinary species. The mass
        // balance alone cannot tell a Gibbs-energy minimum from any feasible point. An
        // interior species needs g[j] = F[j] - sum_i U[i]*A[i,j] ~ 0; a species at its lower
        // bound needs g[j] >= 0, at its upper bound g[j] <= 0. A mineral wrongly left at its
        // lower bound shows up as a strongly negative g[j].
        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
        // A phase below its minimum amount gets no chemical potentials above (they stay 0), so it would
        // always look correctly absent. Judge it by its pure species instead (aqueous phase: the water).
        // A single-species phase below that amount gets the potential it has when present (pm.G).
        if( pa_p->OptimaCollapseCheck )
        {
            long int jc = 0, nset = 0, nsingle = 0;
            for( long int k = 0; k < pm.FI; k++ )
            {
                const long int jn = jc + pm.L1[k];
                const bool vanished = pm.L1[k] > 1 && ( pm.XF[k] <= pm.DSM
                    || ( pm.PHC[k] == PH_AQUEL && pm.LO >= 0 && pm.X[pm.LO] <= pm.XwMinM ) );
                if( vanished )
                    for( long int j = jc; j < jn; j++ )
                    {
                        if( pm.PHC[k] == PH_AQUEL && j != pm.LO ) continue;
                        pm.F[j] = pm.G0[j] + ( pm.fDQF ? pm.fDQF[j] : 0. );
                        nset++;
                    }
                // the same skip as in PrimalChemicalPotentials()
                const bool singleSkipped = pm.L1[k] == 1 && ( pm.XF[k] < pm.PhMinM || pm.XF[k] <= pm.DSM );
                if( singleSkipped && pm.DCCW[jc] == DC_SINGLE )
                {
                    pm.F[jc] = pm.G[jc];
                    nsingle++;
                }
                jc = jn;
            }
            if( nset || nsingle ) native_trace_decide( "collapsecheck species=%ld single=%ld", nset, nsingle );
        }
        double maxKKTResidual = 0.;
        long int worstKKTSpecies = -1;
        std::vector<double> gradSaved( (size_t)L, 0. );   // for pa_OptimaZeroAbsent below
        const double kktTol = std::max( pa_p->GAS, 1e-300 );
        long int nBarrierPassed = 0;   // species the barrier term alone would have failed
        for( long int j = 0; j < L; j++ )
        {
            // The gradient Optima minimised includes, for pure-phase species, the
            // log-barrier term -tau/X[j]; stationarity is checked for that objective.
            double gradJ = pm.F[j];
            // A pure-phase species whose barrier term exceeds the tolerance is held there by the barrier:
            // it is judged as sitting on its floor, on its gradient without that term (a tie then passes).
            const double barrier = ( j >= pm.Ls ) ? kLogBarrierTau / std::max( pm.X[j], dcFloor ) : 0.;
            const bool barrierHeld = barrier > kktTol;
            if( !barrierHeld )
                gradJ -= barrier;
            for( long int i = 0; i < N; i++ )
                gradJ -= pm.U[i] * pm.A[ i + j*N ];
            gradSaved[(size_t)j] = gradJ;
            double resid;
            // A kinetically fixed species (degenerate box, xlower == xupper) is an equality
            // constraint: its multiplier is unrestricted in sign, so no one-sided test applies.
            if( problem.xupper[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) )
                resid = 0.;                        // fixed (degenerate box): no sign test applies
            else if( barrierHeld )
            {
                resid = std::max( -gradJ, 0. );   // held at its floor by the barrier: gradJ should be >= 0
                if( resid <= kktTol && std::max( barrier - gradJ, 0. ) > kktTol ) nBarrierPassed++;
            }
            else if( pm.X[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) )
                resid = std::max( -gradJ, 0. );   // at lower bound: gradJ should be >= 0
            else if( pm.X[j] >= problem.xupper[j] * ( 1. - 1e-6 ) )
                resid = std::max(  gradJ, 0. );   // at upper bound: gradJ should be <= 0
            else
                resid = std::fabs( gradJ );        // interior: gradJ should be ~0
            if( resid > maxKKTResidual )
            {
                maxKKTResidual = resid;
                worstKKTSpecies = j;
            }
        }
        if( nBarrierPassed > 0 )
            native_trace_decide( "kktbarrier passed=%ld tau=%.1e tol=%.1e", (long)nBarrierPassed, kLogBarrierTau, kktTol );
        // pa_OptimaTailStop: after a tail-stop answer, a species holding under 1e-9 of the total amount is not
        // held to the stationarity test, as the stop itself did not; recorded with the largest skipped residual.
        if( tailAccepted && maxKKTResidual > kktTol )
        {
            double xSum = 0.;
            for( long int j = 0; j < L; j++ ) xSum += std::max( pm.X[j], 0. );
            double kept = 0., skipped = 0.; long int nSkip = 0, worstKept = -1;
            for( long int j = 0; j < L; j++ )
            {
                double gradJ = gradSaved[(size_t)j];
                double resid;
                if( problem.xupper[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) ) resid = 0.;
                else if( pm.X[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) ) resid = std::max( -gradJ, 0. );
                else if( pm.X[j] >= problem.xupper[j] * ( 1. - 1e-6 ) ) resid = std::max( gradJ, 0. );
                else resid = std::fabs( gradJ );
                if( resid <= kktTol ) continue;
                if( pm.X[j] <= 1e-9 * xSum ) { nSkip++; skipped = std::max( skipped, resid ); }
                else if( resid > kept ) { kept = resid; worstKept = j; }
            }
            native_trace_decide( "tailkkt skipped=%ld worst_skipped=%.3e worst_kept=%.3e", nSkip, skipped, kept );
            maxKKTResidual = kept;
            worstKKTSpecies = worstKept;
        }
        bool kktOk = maxKKTResidual <= kktTol;

        // Dual determinacy. The sign test above reads each reduced gradient at the dual Optima
        // returned. When the interior species (whose g = 0 fixes the dual) span fewer than N
        // directions, the dual is free along the rest, and every bound-active species'
        // gradient depends on where along them the solver stopped. The question is then
        // whether some dual in that free set satisfies every bound condition. Runs only when
        // the plain test failed; one free direction (typically redox) is handled exactly,
        // more are recorded and left failing.
        // In plain words: when some element potentials are not pinned down by the answer,
        // check whether any admissible choice of them makes the answer valid.
        if( !kktOk )
        {
            auto boxFixed = [&]( long int j ) {
                return problem.xupper[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ); };
            auto atLower = [&]( long int j ) {
                return pm.X[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ); };
            auto atUpper = [&]( long int j ) { return pm.X[j] >= problem.xupper[j] * ( 1. - 1e-6 ); };

            // Orthonormal basis of span{a_j : j interior}, by modified Gram-Schmidt with the same
            // rank tolerance MassBalanceReproject() uses.
            std::vector<double> Q, r( (size_t)N );
            long int rank = 0;
            auto orthogonalize = [&]( const std::vector<double>& basis, long int nb ) {
                for( long int c = 0; c < nb; c++ )
                {
                    double d = 0.;
                    for( long int i = 0; i < N; i++ ) d += basis[(size_t)(c*N + i)] * r[(size_t)i];
                    for( long int i = 0; i < N; i++ ) r[(size_t)i] -= d * basis[(size_t)(c*N + i)];
                }
                double nrm = 0.;
                for( long int i = 0; i < N; i++ ) nrm += r[(size_t)i] * r[(size_t)i];
                return std::sqrt( nrm );
            };
            for( long int j = 0; j < L && rank < N; j++ )
            {
                if( boxFixed( j ) || atLower( j ) || atUpper( j ) ) continue;
                double nrm0 = 0.;
                for( long int i = 0; i < N; i++ ) { r[(size_t)i] = pm.A[ i + j*N ]; nrm0 += r[(size_t)i]*r[(size_t)i]; }
                nrm0 = std::sqrt( nrm0 );
                if( !( nrm0 > 0. ) ) continue;
                const double nrm = orthogonalize( Q, rank );
                if( nrm < 1e-8 * nrm0 ) continue;
                Q.resize( (size_t)((rank+1)*N) );
                for( long int i = 0; i < N; i++ ) Q[(size_t)(rank*N + i)] = r[(size_t)i] / nrm;
                rank++;
            }
            // Free directions: unit vectors orthogonalized against the span and each other.
            std::vector<double> D; long int nDirs = 0;
            for( long int e = 0; e < N && rank + nDirs < N; e++ )
            {
                std::fill( r.begin(), r.end(), 0. ); r[(size_t)e] = 1.;
                double nrm = orthogonalize( Q, rank );
                if( nrm < 1e-8 ) continue;
                for( long int i = 0; i < N; i++ ) r[(size_t)i] /= nrm;
                nrm = orthogonalize( D, nDirs );
                if( nrm < 1e-8 ) continue;
                D.resize( (size_t)((nDirs+1)*N) );
                for( long int i = 0; i < N; i++ ) D[(size_t)(nDirs*N + i)] = r[(size_t)i] / nrm;
                nDirs++;
            }

            double tLo = -std::numeric_limits<double>::infinity(), tHi = std::numeric_limits<double>::infinity();
            // The same interval without the tolerance slack. The point is taken from it when
            // it is not empty: an end of the slack interval puts the bounding species at a
            // residual of exactly kktTol, where rounding would decide the verdict.
            double tLo0 = -std::numeric_limits<double>::infinity(), tHi0 = std::numeric_limits<double>::infinity();
            long int jLo = -1, jHi = -1;
            double tStar = 0., residAfter = maxKKTResidual;
            int resolved = 0;
            const double before = maxKKTResidual;
            if( nDirs == 1 )
            {
                std::vector<double> cdir( (size_t)L, 0. );
                bool infeasible = false;
                for( long int j = 0; j < L; j++ )
                {
                    double c = 0.;
                    for( long int i = 0; i < N; i++ ) c += pm.A[ i + j*N ] * D[(size_t)i];
                    cdir[(size_t)j] = c;
                    if( boxFixed( j ) ) continue;
                    const double g = gradSaved[(size_t)j];
                    const bool lower = atLower( j ), upper = !lower && atUpper( j );
                    if( !lower && !upper ) continue;             // interior: c = 0 by construction
                    // lower: g - t c >= -kktTol ; upper: g - t c <= kktTol
                    const double rhs = lower ? ( g + kktTol ) : ( g - kktTol );
                    if( std::fabs( c ) < 1e-12 )
                    {
                        if( lower ? ( g < -kktTol ) : ( g > kktTol ) ) infeasible = true;
                        continue;
                    }
                    const double bound = rhs / c, bound0 = g / c;
                    const bool isUpperOnT = ( lower == ( c > 0. ) );
                    if( isUpperOnT ) { if( bound < tHi ) { tHi = bound; jHi = j; } tHi0 = std::min( tHi0, bound0 ); }
                    else             { if( bound > tLo ) { tLo = bound; jLo = j; } tLo0 = std::max( tLo0, bound0 ); }
                }
                if( !infeasible && tLo <= tHi )
                {
                    tStar = ( tLo0 <= tHi0 ) ? std::min( std::max( 0., tLo0 ), tHi0 )
                                             : std::min( std::max( 0., tLo ), tHi );
                    residAfter = 0.; long int worstAfter = -1;
                    for( long int j = 0; j < L; j++ )
                    {
                        if( boxFixed( j ) ) continue;
                        const double g = gradSaved[(size_t)j] - tStar * cdir[(size_t)j];
                        const double res = atLower( j ) ? std::max( -g, 0. )
                                         : atUpper( j ) ? std::max(  g, 0. ) : std::fabs( g );
                        if( res > residAfter ) { residAfter = res; worstAfter = j; }
                    }
                    if( residAfter <= kktTol )
                    {
                        for( long int j = 0; j < L; j++ ) gradSaved[(size_t)j] -= tStar * cdir[(size_t)j];
                        maxKKTResidual = residAfter;
                        worstKKTSpecies = worstAfter;
                        kktOk = true;
                        resolved = 1;
                    }
                }
            }
            auto dcName = [&]( long int j ) {
                if( j < 0 ) return std::string( "-" );
                std::string s = char_array_to_string( pm.SM[j], MAXDCNAME );
                while( !s.empty() && ( s.back() == ' ' || s.back() == '\0' ) ) s.pop_back();
                return s;
            };
            // Report the verdict on the CERT line. dual_resolved stays -1 ("not searched") on
            // every row whose plain sign test passed.
            native_cert_dualfree_set( resolved );
            native_trace_decide( "dualfree rank=%ld of=%ld dirs=%ld tlo=%.4g thi=%.4g tstar=%.4g resolved=%d "
                                 "before=%.3e after=%.3e lo=%s hi=%s",
                                 (long)rank, (long)N, (long)nDirs, tLo, tHi, tStar, resolved,
                                 before, residAfter, dcName( jLo ).c_str(), dcName( jHi ).c_str() );
        }

        double worstStabilityViol = 0.;
        bool worstStabilityWasAbsent = false;
        std::vector<PhStabViolation> stabRanked;
        long int worstStabilityPhase = WorstPhaseStabilityViolation(
                    presenceThreshold, dcFloor,
                    extinctFixed.empty() ? nullptr : extinctFixed.data(),
                    worstStabilityViol, worstStabilityWasAbsent, &stabRanked, nullptr, "final" );
        bool stabilityOk = ( worstStabilityPhase < 0 );

        // Commit each active condition's titrant into pm.B[], so that the mass-balance check
        // and every consumer of pm.B[] (including CNode->bIC[]) see the titrated bulk
        // composition. Optima's equality is A*Y + sum_k stoich_k*xi_k = B, so the effective
        // bulk is B - sum_k stoich_k*xi_k. The system is still in pa_DG's internal scale and
        // RescaleSystemFromInternal() does not know the titrant: pm.B is committed in
        // internal units, titrantAmount is stored in real moles (xi / ScFact).
        double bulkTotal = 0.;   // all element amounts before titration, internal units
        for( long int i = 0; i < N; i++ )
            if( !pm.ICC || ( pm.ICC[i] != IC_CHARGE && pm.ICC[i] != IC_VOLUME ) )
                bulkTotal += fabs( pm.B[i] );
        for( long int k = 0; k < R; k++ )
        {
            for( const auto& rc : conditions[k].stoich )
                pm.B[ rc.first ] -= rc.second * state.x[L+k];
            conditions[k].titrantAmount = state.x[L+k] / ScFact;
        }
        // An unusual target, or one that needs more titrant than the whole system, is only reported:
        // the answer may be right for the model but mean little physically.
        for( long int k = 0; k < R; k++ )
        {
            const auto& c = conditions[k];
            const bool farTarget = ( c.name == "pH" && ( c.target < -2. || c.target > 16. ) )
                                || ( c.name == "Eh" && fabs( c.target ) > 2. );
            const bool hugeTitrant = bulkTotal > 0. && fabs( state.x[L+k] ) > bulkTotal;
            if( ( farTarget || hugeTitrant ) && !g_optimaTargetWarned )
            {
                g_optimaTargetWarned = true;
                ipm_logger->warn( "{} target {} {}. Check that the target is meant; the result may have little physical meaning.",
                                  c.name, c.target,
                                  hugeTitrant ? fmt::format( "needs {:.3g} mol of titrant, more than the whole system holds",
                                                             fabs( c.titrantAmount ) )
                                              : std::string( "is outside the usual range" ) );
            }
        }

        // Verify each active condition actually reached its target. At a bound-active
        // titrant slot the box constraint's dual absorbs the KKT residual, so Optima can
        // report success while the condition is unsatisfied.
        bool allTargetsMet = true;
        std::string targetMissBuf;
        for( long int k = 0; k < R; k++ )
        {
            double achievedMuj = 0.;
            for( const auto& rc : conditions[k].stoich )
                achievedMuj += pm.U[ rc.first ] * rc.second;
            conditions[k].achievedValue = conditions[k].achievedValueFn( achievedMuj );
            double tol = conditions[k].tolerance;
            if( tol < 0. && conditions[k].defaultToleranceFn )
                tol = conditions[k].defaultToleranceFn();
            conditions[k].targetMet = std::fabs( conditions[k].achievedValue - conditions[k].target ) <= tol;
            // A pH or Eh with no water left is not an answer: the value above then describes
            // nothing, and the reported pH/Eh would be left over from before. Same water test
            // as the molality calculation.
            const bool aqueousGone = conditions[k].needsAqueousPhase
                && !( pm.FIs > 0 && pm.PHC[0] == PH_AQUEL && pm.XF[0] > pm.DSM && pm.XFA[0] > pm.XwMinM );
            if( aqueousGone )
                conditions[k].targetMet = false;
            if( !conditions[k].targetMet )
            {
                allTargetsMet = false;
                targetMissBuf += conditions[k].name + ": target=" + std::to_string(conditions[k].target)
                               + " achieved=" + std::to_string(conditions[k].achievedValue)
                               + ( aqueousGone ? " (no water left)" : "" ) + "; ";
            }
        }

        // Mass-balance check with GEMS3K's own per-IC tolerance (pm.DHBM), as native's
        // testMulti() path does, against the post-titration bulk composition.
        long int massBalanceBadIC = CheckMassBalanceResiduals( pm.Y );

        // Report Optima's total iteration count in pm.ITG (the NumIterIPM channel);
        // pm.ITF (MBR-equivalent) stays 0 - this solver has no such phase.
        pm.ITG = optimaIterTotal;

        if( auto* tfp = Optima::trace::file() )
        {
            // Answer and per-step records (f in J, err is Optima's final error).
            const auto& ts = Optima::trace::state();
            const bool retried = optimaIterTotal != traceFirstIt;
            double dinf = 0., dinfp = 0.;
            long nfar = 0, nfarp = 0, npresent = 0;
            for( long int j = 0; j < L && (long)traceX0.size() == L; j++ )
            {
                const double d = std::abs( std::log( std::max( traceX0[j], dcFloor ) ) - std::log( std::max( pm.Y[j], dcFloor ) ) );
                dinf = std::max( dinf, d );
                if( d > 1.0 ) nfar++;
                if( pm.Y[j] > 1e-12 * ScFact ) { npresent++; dinfp = std::max( dinfp, d ); if( d > 1.0 ) nfarp++; }
            }
            fprintf( tfp, "ANSWER solve=%ld it_first=%ld succ_first=%d it_retry=%ld succ_retry=%d retried=%d it_reported=%ld f=%.12e err=%.6e"
                          " dx0_inf=%.6e dx0_n_far=%ld dx0_inf_present=%.6e dx0_n_far_present=%ld n_present=%ld"
                          " F_updates_first=%ld F_updates_total=%ld fx_evals_first=%ld fx_evals_total=%ld ok=%d\n",
                     ts.solve, traceFirstIt, traceFirstOk ? 1 : 0, retried ? (long)optimaIterTotal - traceFirstIt : 0L,
                     retried ? ( result.succeeded ? 1 : 0 ) : 0, retried ? 1 : 0, (long)optimaIterTotal, *traceLastF * pm.RT / ScFact, result.error,
                     dinf, nfar, dinfp, nfarp, npresent, traceFupFirst, ts.fupdates, traceFxFirst, ts.fxevals, result.succeeded ? 1 : 0 );
            std::string onset, gone;
            std::vector<char> now;
            for( long int k = 0; k < pm.FI; k++ )
            {
                now.push_back( pm.YF[k] > 1e-12 * ScFact );
                const bool before = k < (long)tracePhPrev.size() && tracePhPrev[k];
                const std::string nm = char_array_to_string( pm.SF[k], MAXPHNAME );
                if( now.back() && !before ) onset += ( onset.empty() ? "" : "," ) + nm;
                if( !now.back() && before ) gone += ( gone.empty() ? "" : "," ) + nm;
            }
            traceCumIt += optimaIterTotal;
            fprintf( tfp, "STEP  solve=%ld it=%ld cum_it=%ld succeeded=%d onset=%s gone=%s\n",
                     ts.solve, (long)optimaIterTotal, traceCumIt, result.succeeded ? 1 : 0, onset.c_str(), gone.c_str() );
            fflush( tfp );
            tracePhPrev = now;
            traceXPrev.assign( pm.Y, pm.Y + L );
        }

        // ---- pa_OptimaZeroAbsent: absent species in the answer (see BASE_PARAM::OptimaZeroAbsent).
        // Runs after every trustworthiness check and only on a clean verdict, so it never
        // turns an accepted answer into a rejected one.
        if( ( pa_p->OptimaZeroAbsent == 2 || pa_p->OptimaZeroAbsent == 3 ) && result.succeeded && allTargetsMet
            && massBalanceBadIC < 0 && kktOk && stabilityOk )
        {
            // Value 2 (and the first part of 3): keep the amounts Optima returned, repair only a
            // failing balance. Same acceptance gate as the zeroing below.
            const long int Zc = N - pm.E;
            // worst relative residual over the ordinary ICs (the CERT record's mb_rel), worst absolute
            // residual over the charge row(s), and the charge rows' tolerance DHBM x total charge carried
            auto worstResiduals = [&]( double& relOut, double& chgOut, double& chgTolOut )
            {
                relOut = 0.; chgOut = 0.; chgTolOut = 0.;
                for( long int i = 0; i < N; i++ )
                {
                    double ci = pm.B[i], scale = 0.;
                    for( long int jj = 0; jj < L; jj++ )
                    {
                        ci -= pm.A[ i + jj*N ] * pm.Y[jj];
                        scale += std::fabs( pm.A[ i + jj*N ] ) * pm.Y[jj];
                    }
                    if( i < Zc )
                    {
                        // a default seed (ICIsDefaultSeed()) neither triggers nor vetoes the repair
                        if( ICIsDefaultSeed( i ) ) continue;
                        const double bar = pm.B[i] * pm.DHBM;
                        relOut = std::max( relOut, bar > 0. ? std::fabs( ci ) / bar : ( ci != 0. ? 1e300 : 0. ) );
                    }
                    else
                    {
                        chgOut = std::max( chgOut, std::fabs( ci ) );
                        chgTolOut = std::max( chgTolOut, scale * pm.DHBM );
                    }
                }
            };
            double relBefore = 0., chgBefore = 0., chgTol = 0.;
            worstResiduals( relBefore, chgBefore, chgTol );
            if( relBefore > 1. )
            {
                std::vector<double> Ysave( pm.Y, pm.Y + L );
                const bool moved = MassBalanceReproject( pm.Y );
                double relAfter = relBefore, chgAfter = chgBefore, chgTolAfter = 0.;
                if( moved )
                    worstResiduals( relAfter, chgAfter, chgTolAfter );
                const int kept = ( moved && relAfter <= 1. && chgAfter <= std::max( chgBefore, chgTol ) ) ? 1 : 0;
                if( moved && !kept )
                    for( long int j = 0; j < L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                if( moved )
                {
                    // MassBalanceReproject() re-synchronised pm.X and the phase totals to ITS state; bring
                    // everything derived back to the state actually returned (kept or restored).
                    for( long int j = 0; j < L; j++ ) pm.X[j] = pm.Y[j];
                    TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                    CheckMassBalanceResiduals( pm.Y );   // pm.C[] describes the returned state
                }
                native_trace_decide( "optimarepair relbefore=%.3e relafter=%.3e chgbefore=%.3e chgafter=%.3e kept=%d",
                                     relBefore, relAfter, chgBefore, chgAfter, kept );
            }
        }
        else if( pa_p->OptimaZeroAbsent != 0 && result.succeeded && allTargetsMet
            && massBalanceBadIC < 0 && kktOk && stabilityOk )
        {
            std::vector<double> Ysave( pm.Y, pm.Y + L );
            long int nZeroed = 0, nPhasesZeroed = 0;

            // Zero only absent phases, then rebalance: a phase is zeroed as a whole, and
            // only if every member passes (a)-(d); species of a present phase are never
            // touched (small dissolved amounts are real).
            const long int Zc = N - pm.E;
            auto residuals = [&]( std::vector<double>& c, std::vector<double>& tol )
            {
                for( long int i = 0; i < N; i++ )
                {
                    double ci = pm.B[i], scale = 0.;
                    for( long int jj = 0; jj < L; jj++ )
                    {
                        ci -= pm.A[ i + jj*N ] * pm.Y[jj];
                        scale += std::fabs( pm.A[ i + jj*N ] ) * pm.Y[jj];
                    }
                    c[(size_t)i] = ci;
                    // ordinary IC: B_i*DHBM, the relative test MBR applies; charge row
                    // (B = 0): DHBM times the total charge carried
                    tol[(size_t)i] = ( i < Zc ) ? pm.B[i] * pm.DHBM : scale * pm.DHBM;
                }
            };
            std::vector<double> cBefore( (size_t)N ), tolBefore( (size_t)N );
            residuals( cBefore, tolBefore );

            for( long int k = 0, j0 = 0; k < pm.FI && j0 < L; j0 += pm.L1[k], k++ )
            {
                const long int j1 = std::min( j0 + pm.L1[k], L );
                bool absent = ( j1 > j0 );
                for( long int j = j0; j < j1 && absent; j++ )
                {
                    // (a) not kinetically required to be present: DLL > 0 is a caller's
                    //     retention floor and must not be discarded.
                    if( pm.DLL[j] > 0. ) { absent = false; break; }
                    // (b) the bound it sits on must be the NUMERICAL floor, not a
                    //     real constraint.
                    const double lo = problem.xlower[j];
                    if( lo > dcFloor * ( 1. + 1e-9 ) ) { absent = false; break; }
                    // (c) actually sitting on it - same expression the KKT check
                    //     above uses to classify a variable as bound-active.
                    if( pm.Y[j] > lo + std::max( dcFloor, lo * 1e-6 ) ) { absent = false; break; }
                    // (d) correctly there: a reduced gradient not below -kktTol (gradSaved is
                    //     already shifted along a free dual direction when the determinacy
                    //     check above resolved one), or a degenerate box (DUL of exactly 0,
                    //     a hard kinetic exclusion).
                    const bool degenerateBox =
                        ( problem.xupper[j] <= lo + std::max( dcFloor, lo * 1e-6 ) );
                    const double g = gradSaved[(size_t)j];
                    if( !degenerateBox && g < -kktTol ) { absent = false; break; }
                }
                if( !absent ) continue;
                long int nThis = 0;
                for( long int j = j0; j < j1; j++ )
                    if( pm.Y[j] != 0. ) { pm.Y[j] = 0.; nThis++; }
                if( nThis > 0 ) { nZeroed += nThis; nPhasesZeroed++; }
            }

            // Rebalance: zeroing removes ~dcFloor per species, which matters for a trace IC.
            // If any IC or the charge row ends past both its residual before zeroing and its
            // tolerance, put the removed amount back onto the present carriers with
            // MassBalanceReproject(). If an IC is still past its limit, the zeroing is undone.
            // A default seed (ICIsDefaultSeed()) neither triggers the rebalance nor the undo.
            int rebalanced = 0, reverted = 0;
            long int nSeeds = 0;
            for( long int i = 0; i < Zc; i++ )
                if( ICIsDefaultSeed( i ) ) nSeeds++;
            if( nZeroed > 0 )
            {
                std::vector<double> cAfter( (size_t)N ), tolAfter( (size_t)N );
                auto brokenIC = [&]() -> long int
                {
                    residuals( cAfter, tolAfter );
                    for( long int i = 0; i < N; i++ )
                        if( std::fabs( cAfter[(size_t)i] ) >
                            std::max( std::fabs( cBefore[(size_t)i] ), tolBefore[(size_t)i] )
                            && !( nSeeds > 0 && ICIsDefaultSeed( i ) ) )
                            return i;
                    return -1;
                };
                if( brokenIC() >= 0 )
                {
                    rebalanced = MassBalanceReproject( pm.Y ) ? 1 : 0;
                    if( brokenIC() >= 0 )
                    {
                        for( long int j = 0; j < L; j++ ) { pm.Y[j] = Ysave[(size_t)j]; pm.X[j] = Ysave[(size_t)j]; }
                        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                        reverted = 1;
                    }
                }
            }
            native_trace_decide( "zeroabsent phases=%ld species=%ld rebalanced=%d reverted=%d of=%ld seeds=%ld",
                                 (long)nPhasesZeroed, (long)nZeroed, rebalanced, reverted, (long)L, (long)nSeeds );
            if( reverted ) nZeroed = 0;

            if( nZeroed > 0 )
            {
                // Self-gate: keep the zeroed state only if it passes the same per-IC test as
                // the accepted state. This test uses an absolute cutoff and does not see a
                // trace-IC loss of ~1e-13 mol; the CERT trace record does.
                if( CheckMassBalanceResiduals( pm.Y ) >= 0 )
                {
                    for( long int j = 0; j < L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                    CheckMassBalanceResiduals( pm.Y );   // restore pm.C[] to the accepted state
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: pa_OptimaZeroAbsent - "
                                       "reverted, zeroing {} absent species would break the mass "
                                       "balance", nZeroed );
                }
                else
                {
                    // Commit, and recompute everything derived from the primal.
                    for( long int j = 0; j < L; j++ ) pm.X[j] = pm.Y[j];
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: pa_OptimaZeroAbsent - "
                                       "{} of {} species reported as exactly zero", nZeroed, L );
                }
            }
        }

        // pa_OptimaTpdAccept (value = TPD tolerance in RT; 0 = off). Optima's own error test
        // demands a non-negative reduced gradient of each end-member of an absent solution
        // phase, whose composition is undetermined at the floor, so a correct answer can
        // keep cycling. Accept the state when the mass balance holds, every species not at a
        // bound is stationary (|gradJ| <= kktTol), every absent species of a single-species,
        // aqueous or gas phase has gradJ >= -kktTol, and every absent non-ideal condensed
        // phase passes the tangent-plane search (NativeTpdPhase(), tpd_min >= -tol).
        // In plain words: accept an answer that is correct by a direct phase-stability check
        // even if the solver's own stopping test was not met.
        {
            const double tpdAccTol = pa_p->OptimaTpdAccept;   // pa_OptimaTpdAccept; 0 = off
            if( tpdAccTol > 0. && ( !result.succeeded || !kktOk || !stabilityOk ) && massBalanceBadIC < 0 && allTargetsMet )
            {
                auto atLo = [&]( long int j ) {
                    return pm.X[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ); };
                auto atHi = [&]( long int j ) { return pm.X[j] >= problem.xupper[j] * ( 1. - 1e-6 ); };
                bool ok = true; double worstInt = 0., worstAbs = 0., worstTpd = 1e300;
                long int nTpd = 0, nVeto = 0, vetoPh = -1, jb = 0;
                for( long int k = 0; k < pm.FI && ok; k++ )
                {
                    const long int n = pm.L1[k];
                    const char ph = pm.PHC[k];
                    bool allLo = true;
                    for( long int j = jb; j < jb + n; j++ ) if( !atLo( j ) ) { allLo = false; break; }
                    const bool tpdPhase = allLo && n > 1 && k < pm.FIs && ph != PH_AQUEL && ph != PH_GASMIX
                                          && ph != PH_PLASMA && ph != PH_FLUID && ph != PH_SORPTION && ph != PH_POLYEL
                                          && ph != PH_ADSORPT && ph != PH_IONEX;
                    if( tpdPhase )
                    {
                        std::vector<double> yb;
                        const double tpd = NativeTpdPhase( k, jb, yb );
                        nTpd++;
                        worstTpd = std::min( worstTpd, tpd );
                        if( !( tpd >= -tpdAccTol ) ) ok = false;   // a search that threw returns 1e300 (rejected below); NaN fails here
                        if( tpd >= 1e299 ) ok = false;
                    }
                    else
                    {
                        // A vanished solution phase has no chemical potentials of its own, so its
                        // species' signs cannot clear it; the phase stability check decides.
                        if( n > 1 && pa_p->OptimaTpdVeto )
                            for( const auto& v : stabRanked )
                                if( v.k == k && v.wasAbsent && !v.clamped ) { nVeto++; vetoPh = k; }
                        for( long int j = jb; j < jb + n; j++ )
                        {
                            const double gj = gradSaved[(size_t)j];
                            if( problem.xupper[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) ) continue;
                            if( atLo( j ) ) { worstAbs = std::max( worstAbs, -gj ); if( gj < -kktTol ) ok = false; }
                            else if( atHi( j ) ) { if( gj > kktTol ) ok = false; }
                            else { worstInt = std::max( worstInt, std::fabs( gj ) ); if( std::fabs( gj ) > kktTol ) ok = false; }
                        }
                    }
                    jb += n;
                }
                native_trace_decide( "tpdaccept ok=%d succeeded=%d kktOk=%d stabilityOk=%d ntpd=%ld worst_tpd=%.4e "
                                     "worst_interior=%.3e worst_absent=%.3e tol=%.1e veto=%ld",
                                     ok ? 1 : 0, result.succeeded ? 1 : 0, kktOk ? 1 : 0, stabilityOk ? 1 : 0,
                                     (long)nTpd, worstTpd, worstInt, worstAbs, tpdAccTol, (long)nVeto );
                // massBalanceBadIC uses an absolute cutoff (min(DHBM*1e10, 1e-2) mol), so the
                // per-IC relative test (|C_i| <= B_i*DHBM, as in MBR) is required as well.
                double mbRelWorst = 0.;
                {
                    const long int Z = pm.N - pm.E;
                    for( long int i = 0; i < Z; i++ )
                    {
                        double c = pm.B[i];
                        for( long int j = 0; j < pm.L; j++ ) c -= pm.A[i + j*pm.N] * pm.X[j];
                        const double bar = pm.B[i] * pm.DHBM;
                        const double r = bar > 0. ? std::fabs( c ) / bar : ( c != 0. ? 1e300 : 0. );
                        mbRelWorst = std::max( mbRelWorst, r );
                    }
                }
                // pa_OptimaAcceptRepair: when the test fails only on the per-IC mass balance,
                // apply MassBalanceReproject() and re-test; restore if mb_rel stays > 1.
                // A certified leveling start may also repair a rounding-size balance error, and so may a tail
                // stop or the solve after it, as an answer that converged would be.
                const bool accRepair = pa_p->OptimaAcceptRepair > 0
                    || g_lvlShortcut || stallWatch->tailStop || stallWatch->vanishedOff;
                if( ok && mbRelWorst > 1. && accRepair )
                {
                    std::vector<double> Ysave( pm.Y, pm.Y + pm.L );
                    const bool moved = MassBalanceReproject( pm.Y );
                    double after = 0.;
                    const long int Z2 = pm.N - pm.E;
                    for( long int i = 0; i < Z2; i++ )
                    {
                        if( ICIsDefaultSeed( i ) ) continue;
                        double c = pm.B[i];
                        for( long int j = 0; j < pm.L; j++ ) c -= pm.A[i + j*pm.N] * pm.Y[j];
                        const double bar = pm.B[i] * pm.DHBM;
                        after = std::max( after, bar > 0. ? std::fabs( c ) / bar : ( c != 0. ? 1e300 : 0. ) );
                    }
                    const bool kept = moved && after <= 1.;
                    if( !kept ) for( long int j = 0; j < pm.L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                    if( moved )
                    {
                        for( long int j = 0; j < pm.L; j++ ) pm.X[j] = pm.Y[j];
                        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                        CheckMassBalanceResiduals( pm.Y );
                    }
                    native_trace_decide( "tpdaccept-repair mb_rel_before=%.3e after=%.3e kept=%d", mbRelWorst, after, kept ? 1 : 0 );
                    if( kept ) mbRelWorst = after;
                }
                if( ok && mbRelWorst > 1. )
                {
                    ok = false;
                    native_trace_decide( "tpdaccept-rejected mb_rel=%.3e", mbRelWorst );
                }
                if( ok && nVeto > 0 )
                {
                    ok = false;
                    native_trace_decide( "tpdveto phase=%s", char_array_to_string( pm.SF[vetoPh], MAXPHNAME ).c_str() );
                }
                // Water driven to its floor means the solve stopped short, not a dry answer: refuse.
                if( ok && pm.FIs > 0 && pm.PHC[0] == PH_AQUEL && pm.LO > 0
                    && pm.X[pm.LO] <= std::max( pm.XwMinM, problem.xlower[pm.LO] * ( 1. + 1e-9 ) ) )
                {
                    ok = false;
                    native_trace_decide( "tpdaccept-dry water=%.3e floor=%.3e", pm.X[pm.LO],
                                         std::max( pm.XwMinM, problem.xlower[pm.LO] ) );
                }
                if( ok )
                {
                    result.succeeded = true; kktOk = true; stabilityOk = true;
                    ipm_logger->debug( "CalculateEquilibriumStateOptima: TPD acceptance - absent non-ideal "
                                       "phases certified by composition search ({} scanned, min TPD {:.3e})", nTpd, worstTpd );
                }
            }
        }

        // An answer whose water phase is reported must look like water: a pH far off the water scale or an
        // absurd ionic strength means the solve stopped on a wrong state, so it is not accepted. Not checked
        // when a pH or Eh target is set: the target sets the state then.
        if( result.succeeded && R == 0 && pm.LO >= 0 && HasAqueousPhase() )
        {
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            const bool insane = !g_aqueousAbsent && ( pm.pH < -2. || pm.pH > 16. || pm.IC > 100. );
            if( insane )
            {
                result.succeeded = false;
                native_trace_decide( "aqsanity pH=%.3f IS=%.3e", pm.pH, pm.IC );
                ipm_logger->info( "Optima: the answer was refused: its water phase has pH {:.2f} and ionic strength {:.3g}",
                                  pm.pH, pm.IC );
            }
        }
        // Absent phases of an accepted answer (phase_floor_parked), and pa_OptimaZeroAbsent = 3.
        phase_floor_parked.assign( (size_t)pm.FI, 0 );
        if( result.succeeded && allTargetsMet && massBalanceBadIC < 0 && kktOk && stabilityOk )
        {
            // The GEMS3K IPM solver's phase-selection rules (PhaseSelectionSpeciationCleanup), with
            // the stability index logSI = pm.Falp[k] from StabilityIndexes() on the converged answer.
            // A phase whose species all sit on the numerical floor is Optima's form of a zero phase:
            // it stays zero unless logSI >= pa_DF, where it would be inserted. A species with an upper
            // limit of 0 (DUL = 0) is zero in the IPM solver, so its floor amount also counts as zero.
            // A phase under any other kinetic limit is not touched. A phase above the floor with
            // logSI <= -pa_DFM is removed by the IPM solver; here only if the log barrier alone holds
            // it above the floor, because no solve follows. Larger ones are only counted.
            long int nParked = 0, nUnstableOffFloor = 0, nSupersatFloor = 0;
            for( long int k = 0, j0 = 0; k < pm.FI && j0 < L; j0 += pm.L1[k], k++ )
            {
                const long int j1 = std::min( j0 + pm.L1[k], L );
                if( j1 <= j0 )
                    continue;
                bool onFloor = true, nearFloor = true, kinConstr = false, allHeldAtZero = true;
                for( long int j = j0; j < j1; j++ )
                {
                    const bool heldAtZero = ( pm.DUL[j] <= 0. && pm.DLL[j] <= 0. );
                    if( !heldAtZero ) allHeldAtZero = false;
                    if( !heldAtZero && ( pm.DUL[j] < 1e6 || pm.DLL[j] > 0. ) ) { kinConstr = true; break; }
                    const double lo = problem.xlower[j];
                    const bool atLo = !( lo > dcFloor * ( 1. + 1e-9 ) ) && pm.Y[j] <= lo + std::max( dcFloor, lo * 1e-6 );
                    // a pure-phase species held just above the floor by the log barrier (as judged in the
                    // stationarity check above)
                    const bool barrierHeld = ( j >= pm.Ls ) && !( lo > dcFloor * ( 1. + 1e-9 ) )
                        && pa_p->LogBarrierTau / std::max( pm.Y[j], dcFloor ) > kktTol;
                    if( !atLo ) onFloor = false;
                    if( !atLo && !barrierHeld ) nearFloor = false;
                }
                if( kinConstr ) continue;
                if( onFloor )
                {   // a zero phase: stays zero unless it would be inserted
                    if( allHeldAtZero || pm.Falp[k] < pa_p->DF ) { phase_floor_parked[(size_t)k] = 1; nParked++; }
                    else nSupersatFloor++;
                }
                else if( pm.Falp[k] <= -pa_p->DFM )
                {   // a present phase the IPM solver would remove; done here only for barrier-held amounts
                    if( nearFloor ) { phase_floor_parked[(size_t)k] = 1; nParked++; }
                    else nUnstableOffFloor++;
                }
            }
            if( pa_p->OptimaZeroAbsent == 3 && nParked > 0 )
            {
                // Set the absent phases to exactly 0, as the GEMS3K IPM solver reports them, and
                // rebalance. If an element still fails its balance, the absent phases that carry
                // it keep their floor amounts and the rest is tried again.
                const long int Zc = N - pm.E;
                std::vector<double> Ysave( pm.Y, pm.Y + L );
                auto residuals = [&]( std::vector<double>& c, std::vector<double>& tol )
                {
                    for( long int i = 0; i < N; i++ )
                    {
                        double ci = pm.B[i], scale = 0.;
                        for( long int jj = 0; jj < L; jj++ )
                        {
                            ci -= pm.A[ i + jj*N ] * pm.Y[jj];
                            scale += std::fabs( pm.A[ i + jj*N ] ) * pm.Y[jj];
                        }
                        c[(size_t)i] = ci;
                        tol[(size_t)i] = ( i < Zc ) ? pm.B[i] * pm.DHBM : scale * pm.DHBM;
                    }
                };
                std::vector<double> cBefore( (size_t)N ), tolBefore( (size_t)N ), cAfter( (size_t)N ), tolAfter( (size_t)N );
                residuals( cBefore, tolBefore );
                // An element is broken when its residual grew past the residual before and past frac x
                // its tolerance. The rebalance starts at 0.1, the share of the tolerance the GEMS3K IPM
                // solver allows a removed phase to use (CutoffDistortionMBR).
                auto brokenIC = [&]( double frac ) -> long int
                {
                    residuals( cAfter, tolAfter );
                    for( long int i = 0; i < N; i++ )
                        if( std::fabs( cAfter[(size_t)i] ) >
                            std::max( std::fabs( cBefore[(size_t)i] ), frac * tolBefore[(size_t)i] )
                            && !ICIsDefaultSeed( i ) )
                            return i;
                    return -1;
                };
                std::vector<char> zeroIt( phase_floor_parked );
                long int nZeroPh = 0, nZeroDC = 0, nKeptPh = 0, rounds = 0;
                int rebalanced = 0, ok = 0;
                for( ; rounds < 16; rounds++ )
                {
                    for( long int j = 0; j < L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                    nZeroPh = 0; nZeroDC = 0;
                    for( long int k = 0, j0 = 0; k < pm.FI && j0 < L; j0 += pm.L1[k], k++ )
                    {
                        if( !zeroIt[(size_t)k] ) continue;
                        const long int j1 = std::min( j0 + pm.L1[k], L );
                        for( long int j = j0; j < j1; j++ )
                            if( pm.Y[j] != 0. ) { pm.Y[j] = 0.; nZeroDC++; }
                        nZeroPh++;
                    }
                    if( nZeroPh == 0 ) break;
                    rebalanced = 0;
                    long int bad = brokenIC( 0.1 );
                    if( bad >= 0 )
                    {
                        rebalanced = MassBalanceReproject( pm.Y ) ? 1 : 0;
                        bad = brokenIC( 1. );
                    }
                    if( bad < 0 ) { ok = 1; break; }
                    // keep at the floor every zeroed phase that carries the failing element
                    long int nKeep = 0;
                    for( long int k = 0, j0 = 0; k < pm.FI && j0 < L; j0 += pm.L1[k], k++ )
                    {
                        if( !zeroIt[(size_t)k] ) continue;
                        const long int j1 = std::min( j0 + pm.L1[k], L );
                        for( long int j = j0; j < j1; j++ )
                            if( pm.A[ bad + j*N ] != 0. ) { zeroIt[(size_t)k] = 0; nKeep++; break; }
                    }
                    nKeptPh += nKeep;
                    if( nKeep == 0 ) break;   // no absent phase carries it: zeroing is not the cause
                }
                // The same per-IC test as the accepted state must still pass.
                if( ok && CheckMassBalanceResiduals( pm.Y ) >= 0 )
                    ok = 0;
                if( !ok )
                    for( long int j = 0; j < L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                else
                {
                    zeroed_floor_amount.assign( (size_t)L, 0. );
                    for( long int j = 0; j < L; j++ )
                        if( pm.Y[j] == 0. && Ysave[(size_t)j] > 0. ) zeroed_floor_amount[(size_t)j] = Ysave[(size_t)j];
                }
                for( long int j = 0; j < L; j++ ) pm.X[j] = pm.Y[j];
                TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                CheckMassBalanceResiduals( pm.Y );   // pm.C[] describes the returned state
                native_trace_decide( "zeroabsent3 parked=%ld zeroed=%ld species=%ld keptfloor=%ld unstable_offfloor=%ld "
                                     "supersat_floor=%ld rounds=%ld rebalanced=%d reverted=%d",
                                     nParked, ok ? nZeroPh : 0L, ok ? nZeroDC : 0L, ok ? nKeptPh : nParked,
                                     nUnstableOffFloor, nSupersatFloor, rounds, rebalanced, ok ? 0 : 1 );
                ipm_logger->debug( "CalculateEquilibriumStateOptima: pa_OptimaZeroAbsent=3 - {} of {} absent "
                                   "phases set to 0{}", ok ? nZeroPh : 0L, nParked,
                                   ok ? "" : " (reverted: the element balance would break)" );
            }
        }

        if( !result.succeeded && pa_p->DW )
        {
            // pa_DW gates the same decision for native's "MBR iterations exceeded" case.
            Error( "E90IPM: Optima solver: ", "Optima::Solver::solve() did not converge (pa_DW, GEMS: Pa_DPV[1], forces this to a hard error)" );
        }
        else if( !result.succeeded || !allTargetsMet || massBalanceBadIC >= 0 || !kktOk || !stabilityOk )
        {
            // Soft failure: a state was produced, but it is not fully trustworthy - the
            // BAD_GEM_* convention (testMulti() reads pm.MK via TNode::GEM_run()).
            pm.MK = 2;
            std::string buf = "Optima solve ";
            buf += result.succeeded ? "reported success" : "did not converge";
            if( !allTargetsMet )
                buf += std::string("; control condition(s) not met: ") + targetMissBuf;
            if( massBalanceBadIC >= 0 )
                buf += "; mass balance residual exceeds tolerance for IC " + name_for_message(pm.SB[massBalanceBadIC],3);
            if( !kktOk )
                buf += "; KKT stationarity residual " + std::to_string(maxKKTResidual)
                     + " exceeds tolerance " + std::to_string(kktTol) + " at species "
                     + ( worstKKTSpecies >= 0 ? char_array_to_string(pm.SM[worstKKTSpecies],MAXDCNAME) : std::string("?") );
            if( !stabilityOk )
            {
                std::ostringstream oss;
                oss.precision(6);
                oss << std::scientific << "; phase-assemblage stability violated (logSI-vs-threshold gap "
                    << worstStabilityViol << ", YF=" << pm.YF[worstStabilityPhase]
                    << ", DSM=" << pm.DSM << ") for phase "
                    << char_array_to_string(pm.SF[worstStabilityPhase],20)
                    << ( worstStabilityWasAbsent ? " (absent but stable - should be present)"
                                                  : " (present but unstable - should be absent)" );
                buf += oss.str();
            }
            setErrorMessage( 21, "W21IPM: Optima solver: ", buf.c_str() );
        }
        // Same read-only system-definition check as native's converged answer. Runs while
        // still in pa_DG's internal scale; the check reports real moles itself.
        if( result.succeeded && pm.MK == 0 )
            StrandedElementCheck();
    }
    catch( TError& xcpt )
    {
        if( pa_p->DG > 1e-5 )
            RescaleSystemFromInternal( ScFact );
        NumIterFIA = pm.ITF;
        NumIterIPM = pm.ITG;
        pm.t_end = clock();
        pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
        throw;      // same error, already logged once
    }
    catch( std::exception& e )
    {
        // Optima reports its own failures by throwing std::runtime_error, not TError. Rescale
        // back to external units and rethrow as TError, so every caller (including HOP's
        // restore path) sees one failure type.
        if( pa_p->DG > 1e-5 )
            RescaleSystemFromInternal( ScFact );
        NumIterFIA = pm.ITF;
        NumIterIPM = pm.ITG;
        pm.t_end = clock();
        pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
        Error( "E91IPM: Optima exception: ", e.what() );
    }

    if( pa_p->DG > 1e-5 )
        RescaleSystemFromInternal( ScFact );

    // Store the total Gibbs energy in pm.FX, which packDataBr() publishes as DATABR.Gs.
    // Done after RescaleSystemFromInternal() (which divides pm.FX by ScFact) and after every
    // post-solve step has committed its amounts. TotalGibbsEnergy() is what
    // TNode::Get_GibbsEnergy() also calls, so the two agree. It writes pm.X/XF/XFA, so these
    // are saved and restored by copy.
    {
        std::vector<double> Xsave( pm.X,   pm.X   + pm.L );
        std::vector<double> XFsave( pm.XF, pm.XF  + pm.FI );
        std::vector<double> XFAsave( pm.XFA, pm.XFA + pm.FIs );
        pm.FX = TotalGibbsEnergy();
        for( long int j = 0; j < pm.L; j++ )    pm.X[j]   = Xsave[(size_t)j];
        for( long int k = 0; k < pm.FI; k++ )   pm.XF[k]  = XFsave[(size_t)k];
        for( long int k = 0; k < pm.FIs; k++ )  pm.XFA[k] = XFAsave[(size_t)k];
    }

    // Count the present phases as the GEMS3K IPM solver does; GEMSGUI sizes its saved record by this count
    pm.FI1 = 0;
    pm.FI1s = 0;
    for( long int k = 0; k < pm.FI; k++ )
        if( pm.YF[k] >= std::min( pm.PhMinM, 1e-22 ) )
        {
            pm.FI1++;
            if( k < pm.FIs )
                pm.FI1s++;
        }

    NumIterFIA = pm.ITF;
    NumIterIPM = pm.ITG;
    pm.t_end = clock();
    pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
    return pm.t_elap_sec;
}

// HOP: native selects the phases (its own cold IPM/MBR/PSSC), then Optima finishes,
// warm-started from native's converged primal and dual. Dispatched from
// TNode::GEM_run() for NEED_GEM_HOP. The Optima leg runs in its own try/catch: if it
// fails after a successful native solve, native's answer is restored and reported as
// BAD_GEM_HOP, so HOP is never worse than a plain native solve.
// In plain words: the original solver finds the answer, then the Optima solver polishes it.
double TMultiBase::CalculateEquilibriumStateHOP( long int& NumIterFIA, long int& NumIterIPM,
                                                 bool warmNative )
{
    long int fiaN = 0, ipmN = 0;
    bool nativeOk = true;
    bool noProgWarm = false;   // pa_IpmNoProgress stopped the native leg: Optima starts warm from its last state
    double calcTime = 0.;
    std::vector<double> Ysave, Usave;
    double FXsave = kTotalGibbsEnergyUnset;   // native's total G, external units (restored below)

    // SHP (warmNative) starts the native leg warm from the state already on this node,
    // which saves most of its cost in a sweep. A node whose pm.U[] is all zero (or
    // non-finite) has never been solved, so its stored speciation is not a warm start
    // and the leg runs cold.
    bool wantWarmNative = warmNative;
    if( wantWarmNative )
    {
        bool usable = false;
        for( long int i = 0; i < pm.N; i++ )
        {
            if( !std::isfinite( pm.U[i] ) ) { usable = false; break; }
            if( pm.U[i] != 0. ) usable = true;
        }
        if( !usable )
        {
            wantWarmNative = false;
            ipm_logger->warn( "Warm start (SHP) requested, but this node has no earlier solution. "
                              "Native leg runs cold. Try: use HOP for a first solve, and SHP only "
                              "after a solve." );
        }
    }

    try
    {
        // pm.pNP = 1 is native's ordinary warm (SIA) start.
        pm.pNP = wantWarmNative ? 1 : 0;
        calcTime = CalculateEquilibriumState( fiaN, ipmN );
        // Snapshot native's converged state before the Optima leg reads and modifies it.
        Ysave.assign( pm.Y, pm.Y + pm.L );
        Usave.assign( pm.U, pm.U + pm.N );
        FXsave = pm.FX;
    }
    catch( TError& werr )
    {
        if( g_ipmNoProgFired && werr.title.rfind( "E06IPM", 0 ) == 0 )
        {
            noProgWarm = true;
            nativeOk = false;
            ipm_logger->info( "HOP: native leg stopped without progress; the Optima leg starts from its last state" );
        }
        // Cold fallback for the warm native leg: native SIA can refuse a state its own
        // cold path returned. Retry the leg cold, so SHP degrades to HOP at the cost of
        // one wasted attempt.
        else if( wantWarmNative )
        {
            ipm_logger->info( "HOP: warm native leg failed ({}: {}); retrying it cold", werr.title, werr.mess );
            wantWarmNative = false;
            fiaN = ipmN = 0;
            try
            {
                pm.pNP = 0;
                calcTime = CalculateEquilibriumState( fiaN, ipmN );
                Ysave.assign( pm.Y, pm.Y + pm.L );
                Usave.assign( pm.U, pm.U + pm.N );
                FXsave = pm.FX;
            }
            catch( TError& nerr2 )
            {
                nativeOk = false;
                ipm_logger->info( "HOP: native leg failed cold too ({}: {}); solving cold with Optima",
                                   nerr2.title, nerr2.mess );
                fiaN = ipmN = 0;
                calcTime = 0.;
            }
        }
        else
        {
            // No assemblage to hand over: degrade to a plain cold Optima solve (an AOP
            // result under a HOP status) and say so.
            nativeOk = false;
            ipm_logger->info( "HOP: native leg failed ({}: {}); solving cold with Optima", werr.title, werr.mess );
            fiaN = ipmN = 0;
            calcTime = 0.;
        }
    }

    // Warm-start Optima from native's converged primal and dual.
    pm.pNP = nativeOk ? 1 : 0;
    long int fiaO = 0, ipmO = 0;
    // Everything in calcTime so far belongs to the native leg.
    const double timeNativeLeg = calcTime;

    // Tell the Optima leg it runs on top of native's assemblage, which lets the
    // dimension reduction run in front of it when pa_OptimaDimReduce > 0 (AUTO does not
    // reach it). Set only when native produced an answer. RAII, so the flag never leaks
    // into a later call.
    struct HopLegGuard {
        TMultiBase* m;
        bool prev;
        HopLegGuard( TMultiBase* mm, bool on ) : m(mm), prev(mm->optima_hop_leg)
            { m->optima_hop_leg = on; }
        ~HopLegGuard() { m->optima_hop_leg = prev; }
    } hopLegGuard( this, nativeOk );

    // The per-leg record (hop_split) is written by a destructor, so a call that throws
    // is recorded too. fiaO/ipmO are real on a throw (set by the Optima path's own catch);
    // timeOptima is not recoverable and stays 0, with `failed` set.
    struct HopLegRecorder {
        TMultiBase* m; bool warm, nok; const long int &fN, &iN, &fO, &iO;
        const double &ct, tN; bool completed = false;
        ~HopLegRecorder()
        {
            m->hop_split.valid      = true;
            m->hop_split.failed     = !completed;
            m->hop_split.warmNative = warm;
            m->hop_split.nativeOk   = nok;
            m->hop_split.fiaNative  = fN;
            m->hop_split.ipmNative  = iN;
            m->hop_split.fiaOptima  = fO;
            m->hop_split.ipmOptima  = iO;
            m->hop_split.timeNative = tN;
            m->hop_split.timeOptima = completed ? ct - tN : 0.;
            // DECIDE record with the iteration split; zero cost when the trace is off.
            // Wall times are kept out of it so traces of identical runs stay identical.
            native_trace_decide( "hop-legsplit warm=%d nativeok=%d failed=%d "
                                 "fian=%ld ipmn=%ld fiao=%ld ipmo=%ld",
                                 warm ? 1 : 0, nok ? 1 : 0, completed ? 0 : 1,
                                 (long)fN, (long)iN, (long)fO, (long)iO );
        }
    } legRecorder{ this, warmNative, nativeOk, fiaN, ipmN, fiaO, ipmO, calcTime, timeNativeLeg };

    if( noProgWarm )
    {
        // Warm from the stopped native state, as the native rescue does; cold if that fails.
        bool ok = false;
        try
        {
            pm.pNP = 1;
            optima_hop_leg = true;
            calcTime += CalculateEquilibriumStateOptima( fiaO, ipmO, false, /* runKinetics */ false );
            optima_hop_leg = false;
            ok = !( pm.MK || pm.PZ );
        }
        catch( TError& ) { optima_hop_leg = false; }
        native_trace_decide( "ipmnoprog-hop warm=%d", ok ? 1 : 0 );
        if( !ok )
        {
            long int fiaC = 0, ipmC = 0;
            pm.pNP = 0;
            try { calcTime += CalculateEquilibriumStateOptima( fiaC, ipmC, false, false ); }
            catch( TError& ) { fiaO += fiaC; ipmO += ipmC; NumIterFIA = fiaN + fiaO; NumIterIPM = ipmN + ipmO; throw; }
            fiaO += fiaC; ipmO += ipmC;
        }
    }
    else if( !nativeOk )
    {
        // No native answer to fall back to: a failure propagates as a plain cold AOP
        // call's would.
        calcTime += CalculateEquilibriumStateOptima( fiaO, ipmO, false,
                              /* runKinetics */ false );  // the native leg already advanced it
    }
    else
    {
        // The Optima leg changes the bulk composition and the control results; keep them.
        const std::vector<double> Bsave( pm.B, pm.B + pm.N );
        const std::vector<EqControlCondition> controlsSave = optima_control_conditions;
        try
        {
            calcTime += CalculateEquilibriumStateOptima( fiaO, ipmO, false,
                              /* runKinetics */ false );  // the native leg already advanced it
        }
        catch( TError& oerr )
        {
            // Restore native's converged primal and dual, and recompute every quantity
            // derived from the primal, which the Optima leg's objective may have changed.
            for( long int j = 0; j < pm.L; j++ ) pm.Y[j] = pm.X[j] = Ysave[(size_t)j];
            for( long int i = 0; i < pm.N; i++ ) pm.U[i] = Usave[(size_t)i];
            // pm.FX too: the Optima leg resets it to kTotalGibbsEnergyUnset. Native's value
            // was saved after its leg rescaled, so it is in external units.
            pm.FX = FXsave;
            for( long int i = 0; i < pm.N; i++ ) pm.B[i] = Bsave[(size_t)i];
            optima_control_conditions = controlsSave;
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            // Soft failure (pm.MK = 2, read by testMulti()): native's answer is reported.
            pm.MK = 2;
            std::string buf = std::string("Optima leg failed after a successful native solve (")
                             + oerr.title + oerr.mess
                             + "); reporting native's own converged state instead";
            setErrorMessage( 21, "W21IPM: HOP: ", buf.c_str() );
            // fiaO/ipmO were set by CalculateEquilibriumStateOptima()'s own catch before it
            // re-threw; the iterations spent are kept, only the state is discarded.
            ipm_logger->warn( "HOP: Optima leg failed after a good native solve ({}: {}). Native's own "
                               "answer is returned as BAD_GEM_HOP. Try: AOP, or check the bulk "
                               "composition for very small amounts", oerr.title, oerr.mess );
        }
    }

    // Report the total cost of both legs; a discarded Optima attempt's iterations are
    // still counted.
    NumIterFIA = fiaN + fiaO;
    NumIterIPM = ipmN + ipmO;

    legRecorder.completed = true;       // see HopLegRecorder above - the record is written
    return calcTime;                    // by its DESTRUCTOR, on this path and on a throw alike
}

// Native left a trace solid phase that still tends to grow (TraceGrowthCheck). Optima, started from
// native's answer, can move it; its answer is kept only if it is better, otherwise native's is put back.
double TMultiBase::NativeOptimaFinish( long int& NumIterFIA, long int& NumIterIPM )
{
    const long int cap = base_param()->NativeOptimaFinish;
    if( cap <= 0 || !trace_growth_stuck || pm.MK || pm.PZ || !optima_control_conditions.empty() )
        return 0.;
    const long int savedPNP = pm.pNP, savedMK = pm.MK, savedPZ = pm.PZ;
    auto snap = [&]( const double* p, long int n ) { return std::vector<double>( p, p + n ); };
    // Native's total G, without disturbing the phase totals it recomputes.
    const std::vector<double> X0 = snap( pm.X, pm.L ), XF0 = snap( pm.XF, pm.FI ), XFA0 = snap( pm.XFA, pm.FIs );
    const double Gnat = TotalGibbsEnergy();
    for( long int j = 0; j < pm.L; j++ )   pm.X[j]   = X0[(size_t)j];
    for( long int k = 0; k < pm.FI; k++ )  pm.XF[k]  = XF0[(size_t)k];
    for( long int k = 0; k < pm.FIs; k++ ) pm.XFA[k] = XFA0[(size_t)k];
    // Everything the Optima leg may change, to put native's answer back exactly.
    const std::vector<double> Ys = snap( pm.Y, pm.L ), Xs = snap( pm.X, pm.L ), Us = snap( pm.U, pm.N ),
        Bs = snap( pm.B, pm.N ), G0s = snap( pm.G0, pm.L ), fDQFs = snap( pm.fDQF, pm.L ), F0s = snap( pm.F0, pm.L ),
        Gs = snap( pm.G, pm.L ), Fs = snap( pm.F, pm.L ), lnGams = snap( pm.lnGam, pm.L ), Gammas = snap( pm.Gamma, pm.L ),
        DULs = snap( pm.DUL, pm.L ), DLLs = snap( pm.DLL, pm.L );
    const double FXs = pm.FX;
    long int fiaO = 0, ipmO = 0;
    double dt = 0., Gopt = 0.;
    bool accept = false;
    std::string why;
    struct Budget {   // the Optima leg's iteration cap; always lifted again, also on a throw
        explicit Budget( long int n ) { g_optimaIterBudget = n; }
        ~Budget() { g_optimaIterBudget = -1; }
    } budget( cap );
    pm.pNP = 1;
    try
    {
        optima_hop_leg = true;
        dt = CalculateEquilibriumStateOptima( fiaO, ipmO, false, /* runKinetics */ false );
        optima_hop_leg = false;
        Gopt = pm.FX;
        if( pm.MK || pm.PZ ) why = "notclean";
        else if( !( Gopt < Gnat - 1e-8 ) ) why = "notlower";
        else accept = true;
    }
    catch( TError& )
    {
        optima_hop_leg = false;
        why = "failed";
    }
    if( !accept )
    {
        for( long int j = 0; j < pm.L; j++ )
        {
            pm.Y[j] = Ys[(size_t)j]; pm.X[j] = Xs[(size_t)j];
            pm.G0[j] = G0s[(size_t)j]; pm.fDQF[j] = fDQFs[(size_t)j]; pm.F0[j] = F0s[(size_t)j];
            pm.G[j] = Gs[(size_t)j]; pm.F[j] = Fs[(size_t)j]; pm.lnGam[j] = lnGams[(size_t)j];
            pm.Gamma[j] = Gammas[(size_t)j]; pm.DUL[j] = DULs[(size_t)j]; pm.DLL[j] = DLLs[(size_t)j];
        }
        for( long int i = 0; i < pm.N; i++ ) { pm.U[i] = Us[(size_t)i]; pm.B[i] = Bs[(size_t)i]; }
        pm.FX = FXs;
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );   // activity coefficients were put back above
        pm.MK = savedMK; pm.PZ = savedPZ;
        ipm_logger->warn( "{}", trace_growth_note );
    }
    else
        ipm_logger->info( "Native solver left a trace phase that still tends to grow; finished with Optima, "
                          "which lowered the energy by {:.3g} RT", Gnat - Gopt );
    pm.pNP = savedPNP;
    native_trace_decide( "nativefinish ok=%d why=%s gnat=%.10f gopt=%.10f fiao=%ld ipmo=%ld cap=%ld",
                         accept ? 1 : 0, accept ? "-" : why.c_str(), Gnat, Gopt, (long)fiaO, (long)ipmO, (long)cap );
    NumIterFIA += fiaO; NumIterIPM += ipmO;
    return dt;
}

// Native stopped at its iteration limit, or (broken) its linear system broke down. Optima solves from native's
// last state, then once from scratch; a clean answer is kept, otherwise the caller reports native's failure.
bool TMultiBase::NativeOptimaRescue( long int& NumIterFIA, long int& NumIterIPM, double& seconds, bool broken )
{
    const long int cap = base_param()->NativeOptimaFinish;
    if( cap <= 0 || !optima_control_conditions.empty() )
        return false;
    const long int savedPNP = pm.pNP;
    long int fiaO = 0, ipmO = 0;
    bool ok = false;
    int attempt = 0;
    for( ; attempt < 2 && !ok; attempt++ )
    {
        struct Budget {   // each attempt's iteration cap; always lifted again, also on a throw
            explicit Budget( long int n ) { g_optimaIterBudget = n; }
            ~Budget() { g_optimaIterBudget = -1; }
        } budget( cap );
        long int fia = 0, ipm = 0;
        pm.pNP = ( attempt == 0 ) ? 1 : 0;
        try
        {
            optima_hop_leg = ( attempt == 0 );
            seconds += CalculateEquilibriumStateOptima( fia, ipm, false, /* runKinetics */ false );
            ok = !( pm.MK || pm.PZ );
        }
        catch( TError& ) {}
        optima_hop_leg = false;
        fiaO += fia; ipmO += ipm;
    }
    pm.pNP = savedPNP;
    native_trace_decide( "nativerescue ok=%d why=%s attempts=%d fiao=%ld ipmo=%ld cap=%ld", ok ? 1 : 0,
                         broken ? "E07" : "E06", attempt, (long)fiaO, (long)ipmO, (long)cap );
    ipm_logger->info( "Native solver {}; Optima {} the answer",
                      broken ? "could not solve its linear system" : "stopped at its iteration limit",
                      ok ? ( attempt == 1 ? "finished" : "found" ) : "did not find" );
    NumIterFIA += fiaO; NumIterIPM += ipmO;
    return ok;
}

#endif // USE_OPTIMA_SOLVER
