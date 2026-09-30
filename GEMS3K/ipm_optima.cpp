//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_optima.cpp
/// Chemical equilibrium via the Optima library's general primal-dual
/// interior-point NLP solver, as an alternative to GEMS3K's own
/// IPM/MBR loop. Reuses GEMS3K's own G0[]/activity-coefficient
/// machinery as the objective function fed to Optima; Optima owns the
/// Newton iteration and box/linear-equality handling entirely.
///
/// Two new solver modes are dispatched here from TNode::GEM_run()
/// (node.cpp): AOP (cold/AIA-equivalent start) and SOP (warm/SIA-
/// equivalent start) - see NODECODECH in databr.h. Both call the one
/// CalculateEquilibriumStateOptima() below, which reads pm.pNP (set by
/// GEM_run() before the call, exactly like the native AIA/SIA path)
/// to choose the starting point.
///
/// Optional equilibrium "control conditions" (pH, Eh, and - by design -
/// any future closed-form target such as fixed fugacity or fixed
/// species activity) are folded into the same joint Newton solve as
/// extra unknowns with a fixed objective-gradient, matching Reaktoro's
/// own EquilibriumSpecs mechanism (see ipm_optima.h). With none
/// registered this is a plain equilibrium solve, architecturally
/// equivalent to AIA/SIA but solved via Optima instead of IPM/MBR.
///
/// Scope, deliberately: no phase-stability (PSSC) pass, no volume/T,P-
/// as-unknown control conditions.
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

#include <chrono>
#include "ms_multi.h"

#ifdef USE_OPTIMA_SOLVER

#include "node.h"
#include <Optima/Optima.hpp>
#include <Optima/Sensitivity.hpp>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>

// pa_OptimaFDDiagFloor hit counter (plan v5 s139.5): FD-Hessian columns whose diagonal was restored to the
// analytic value, accumulated across the reduced pre-solve and the full solve and reported by a DECIDE
// fddiagfloor record after each full solve, then reset.
static long g_fdDiagFloorHits = 0, g_fdDiagFloorCols = 0;

// pa_OptimaLineSearch - put Optima's merit line search on the unmasked error at the given trigger factor.
// Only meaningful with the local Optima fix in ErrorControl::execute (plan v5 s139.6).
static void apply_optima_linesearch( Optima::Options& o, double factor )
{
    if( !( factor > 0. ) ) return;
    o.linesearch.enabled = true;
    o.linesearch.use_unmasked_error = true;
    o.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = factor;
}

namespace {
// Symmetric eigenvalue floor for a small dense block, by cyclic Jacobi
// rotations (n is the number of end-members of one solution phase - 3 for a
// ternary feldspar, a few dozen at most).
//
// Purpose: the exact (finite-differenced) curvature of a solution phase inside
// a miscibility gap is INDEFINITE by construction - that negative eigenvalue
// along the unmixing direction is what a spinodal IS. Newton needs a positive-
// definite model to produce a descent direction at all, and the magnitude of
// the floor sets how far it steps along that soft direction: too large and the
// step is as short as the ideal-mixing model's (the near-critical crawl at rate
// 0.9994); too small and the step is unbounded and annihilates the phase.
// Flooring the eigenvalues at `ratio` times the block's own largest |eigenvalue|
// is therefore a well-scaled trust region expressed in the only units that mean
// anything here - the phase's own curvature - rather than an absolute delta.
//
// Returns false and leaves `a` untouched if the rotations do not settle.
bool SymEigFloorInPlace( std::vector<double>& a, int n, double ratio )
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
    const double floorVal = ratio * lmax;
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

// The system size Optima's DEFAULT boxes are scaled by (species upper bounds with no DUL,
// the solvent-reseed ceiling, the control-condition titrant bound). On the rescaled path
// (pa_DG > 1e-5, ScaleSystemToInternal()) that is pm.SMols = pa_DG, the internal total IC
// moles. With rescaling DISABLED pm.SMols is never set - it stays 0 from ms_multi_file.cpp -
// so every default upper bound collapsed to max(0,1)*10 = 10 mol whatever the system held.
// Measured 2026-09-14 on a scratch copy of Cu-Pourbaix_G_pHtitr with pa_DG = 0 (55.5 mol
// H2O): plain AOP returned ERR_GEM_AOP after 13000 iterations at G = -2375.79 against
// -5327.89 (native unaffected, same G either way), and pH/Eh targeting returned
// ERR_GEM_AOP at every target with both titrants pinned at +-2. Unscaled, the equivalent of
// pm.SMols is the actual total IC moles - the same sum SystemTotalMolesIC() takes - so the
// rescaled path is bit-identical to before and the unscaled one gets the same 10x headroom.
static double optima_default_box_moles( const MULTI& pm, double DG )
{
    if( DG > 1e-5 )
        return pm.SMols;
    double tot = 0.;
    for( long int i = 0; i < pm.N - pm.E; i++ )
        tot += pm.B[i];
    return tot;
}

// Small, self-contained two-phase primal simplex (dense tableau, Bland's
// rule throughout for BOTH entering- and leaving-variable selection - the
// textbook anti-cycling guarantee, needed here since this runs unattended
// with no opportunity for a human to notice a hang): solves
// minimize sum_j x_j  s.t.  A*x = b (A: nRows x nCols, row-major via the
// `a` accessor), x >= 0. Used once per ROP solve (not per Newton
// iteration, see TMultiBase::LPFeasibilitySeed() below), so a plain dense
// O(rows*cols) tableau - no sparsity exploitation - is adequate: nRows/
// nCols are at most a few dozen/roughly a hundred for every GEMS3K project
// checked so far.
//
// Returns false without modifying xOut - meaning "don't trust this",
// covering both a genuinely infeasible system and an iteration-cap/
// numerical failure (the two are not distinguished; either way the caller
// falls back to a simpler seed rather than trust a partial result) - and
// true with xOut sized nCols on success.
//
// `cost` (optional, null = all ones) supplies the objective; `yOut` (optional,
// null = not wanted) receives the optimal DUAL of the equality rows, i.e. the
// y for which c_j - sum_i y_i a_ij >= 0 holds for every column at the optimum.
// With cost = pm.G0 that dual is a genuine first-order estimate of pm.U[], and
// pricing every species against it is what OptimaReducedPreSolve() uses to pick
// a generous-but-solve-independent initial active set (see there).
bool TwoPhaseSimplexMinSum( long int nRows, long int nCols,
                             const std::function<double(long int,long int)>& a,
                             const double* b, std::vector<double>& xOut,
                             const double* cost = nullptr, double* yOut = nullptr )
{
    if( nRows <= 0 || nCols <= 0 )
        return false;
    const double eps = 1e-9;
    const long int nVar = nCols + nRows; // real (0..nCols) + artificial (nCols..nVar) columns
    const long int nCon = nRows;
    // Dense row-major tableau: T[i*(nVar+1)+j], i in [0,nCon), j in [0,nVar]
    // (column nVar is the RHS).
    std::vector<double> T( (size_t)nCon * (size_t)(nVar+1), 0. );
    auto Tref = [&]( long int i, long int j ) -> double& { return T[ (size_t)i*(size_t)(nVar+1) + (size_t)j ]; };
    std::vector<long int> basis( nCon );
    std::vector<double> rowSgn( (size_t)nCon, 1. );   // per-row flip, to un-flip the dual

    double bScale = 1.;
    for( long int i = 0; i < nCon; i++ )
        bScale = std::max( bScale, std::fabs(b[i]) );
    const double feasTol = std::max( 1e-7, 1e-9 * bScale );

    // Flip each row's sign so its RHS is non-negative (b_i==0 keeps sign
    // +1, arbitrarily), then attach that row's own artificial column.
    for( long int i = 0; i < nCon; i++ )
    {
        const double sgn = ( b[i] < 0. ) ? -1. : 1.;
        rowSgn[(size_t)i] = sgn;
        for( long int j = 0; j < nCols; j++ )
            Tref(i,j) = sgn * a(i,j);
        Tref(i, nCols + i) = 1.;
        Tref(i, nVar) = sgn * b[i];
        basis[i] = nCols + i;
    }

    // Phase-1 objective row (reduced costs c_j - z_j, minimizing sum of
    // artificials): with every row's basic variable currently costing 1
    // (its own artificial), canonical-form reduced cost for a real column
    // j is -sum_i T[i][j]; every artificial column is already exactly
    // canonical (reduced cost 0, left at its zero-initialized default -
    // c_j=1 for artificial j, z_j=1*T[i][j]=1 in its own row, 0 elsewhere,
    // so c_j-z_j=0 identically).
    std::vector<double> row0( nVar + 1, 0. );
    for( long int j = 0; j < nCols; j++ )
    {
        double s = 0.;
        for( long int i = 0; i < nCon; i++ ) s += Tref(i,j);
        row0[j] = -s;
    }
    {
        double s = 0.;
        for( long int i = 0; i < nCon; i++ ) s += Tref(i,nVar);
        row0[nVar] = -s; // -(current phase-1 objective value)
    }

    auto pivot = [&]( long int prow, long int pcol )
    {
        const double piv = Tref(prow,pcol);
        for( long int j = 0; j <= nVar; j++ )
            Tref(prow,j) /= piv;
        for( long int i = 0; i < nCon; i++ )
        {
            if( i == prow ) continue;
            const double f = Tref(i,pcol);
            if( f == 0. ) continue;
            for( long int j = 0; j <= nVar; j++ )
                Tref(i,j) -= f * Tref(prow,j);
        }
        const double f0 = row0[pcol];
        if( f0 != 0. )
            for( long int j = 0; j <= nVar; j++ )
                row0[j] -= f0 * Tref(prow,j);
        basis[prow] = pcol;
    };

    // limitCol: only columns [0,limitCol) may enter - phase 2 passes
    // nCols to permanently exclude artificials from re-entering.
    auto runSimplex = [&]( long int limitCol, long int maxIter, double entEps ) -> bool
    {
        for( long int iter = 0; iter < maxIter; iter++ )
        {
            long int q = -1;
            for( long int j = 0; j < limitCol; j++ ) // Bland's rule: smallest negative-reduced-cost index
                if( row0[j] < -entEps ) { q = j; break; }
            if( q < 0 )
                return true; // optimal for this phase
            long int prow = -1; double bestRatio = 0.;
            for( long int i = 0; i < nCon; i++ )
            {
                if( Tref(i,q) <= eps ) continue;
                const double ratio = Tref(i,nVar) / Tref(i,q);
                if( prow < 0 || ratio < bestRatio - 1e-12 ||
                    ( std::fabs(ratio-bestRatio) <= 1e-12 && basis[i] < basis[prow] ) ) // Bland's tie-break
                { prow = i; bestRatio = ratio; }
            }
            if( prow < 0 )
                return false; // unbounded - cannot actually happen for a nonnegative-cost minimization
                               // over x>=0 (bounded below by 0); treat as a numerical failure if hit
            pivot( prow, q );
        }
        return false; // iteration cap hit
    };

    const long int maxIter = std::max( (long int)2000, 20*nVar );
    if( !runSimplex( nVar, maxIter, eps ) )
        return false;
    if( -row0[nVar] > feasTol ) // phase-1 optimum not ~0: genuinely infeasible
        return false;

    // Drive out any artificial still basic at (necessarily zero) level,
    // wherever a genuinely nonzero real-column pivot is available - a
    // redundant/rank-deficient row is the only case where none exists,
    // and is harmless to leave as-is (its value is 0 regardless).
    for( long int i = 0; i < nCon; i++ )
    {
        if( basis[i] < nCols ) continue;
        long int q = -1;
        for( long int j = 0; j < nCols; j++ )
            if( std::fabs(Tref(i,j)) > 1e-7 ) { q = j; break; }
        if( q >= 0 )
            pivot( i, q );
    }

    // Phase 2: minimize sum of REAL variables only; artificial columns
    // permanently excluded from entering (limitCol=nCols below), so any
    // still basic-at-zero from the step above just sit inert.
    double cMax = 1.;
    if( cost )
        for( long int j = 0; j < nCols; j++ ) cMax = std::max( cMax, std::fabs( cost[j] ) );
    std::fill( row0.begin(), row0.end(), 0. );
    for( long int j = 0; j < nCols; j++ ) row0[j] = cost ? cost[j] : 1.;
    for( long int i = 0; i < nCon; i++ )
    {
        if( basis[i] >= nCols ) continue; // degenerate artificial: cost 0, no canonicalization needed
        const double cB = cost ? cost[ basis[i] ] : 1.;
        if( cB == 0. ) continue;
        for( long int j = 0; j <= nVar; j++ )
            row0[j] -= cB * Tref(i,j);
    }
    if( !runSimplex( nCols, maxIter, eps * cMax ) )
        return false;

    xOut.assign( (size_t)nCols, 0. );
    for( long int i = 0; i < nCon; i++ )
        if( basis[i] < nCols )
            xOut[ (size_t)basis[i] ] = std::max( 0., Tref(i,nVar) );

    // Dual of the equality rows. Artificial column i is e_i in the FLIPPED
    // system and costs 0, so its reduced cost is c-z = -y'_i; un-flip with the
    // row's own sign. A row whose artificial is still basic (a redundant row)
    // keeps reduced cost 0 and so gets y_i = 0, which is correct for it.
    // The sign convention is then VERIFIED rather than assumed: the LP optimum
    // must satisfy c_j - y.A_j >= 0 for every column, so if it does not, try the
    // opposite sign, and if that fails too report failure rather than hand back
    // a dual nothing downstream could trust.
    if( yOut )
    {
        std::vector<double> y( (size_t)nCon );
        for( long int i = 0; i < nCon; i++ )
            y[(size_t)i] = -rowSgn[(size_t)i] * row0[nCols + i];
        const double dualTol = 1e-6 * cMax;
        bool ok = false;
        for( int attempt = 0; attempt < 2 && !ok; attempt++ )
        {
            if( attempt == 1 )
                for( long int i = 0; i < nCon; i++ ) y[(size_t)i] = -y[(size_t)i];
            ok = true;
            for( long int j = 0; j < nCols && ok; j++ )
            {
                double z = 0.;
                for( long int i = 0; i < nCon; i++ ) z += y[(size_t)i] * a(i,j);
                if( ( cost ? cost[j] : 1. ) - z < -dualTol ) ok = false;
            }
        }
        if( !ok )
            return false;
        for( long int i = 0; i < nCon; i++ ) yOut[i] = y[(size_t)i];
    }
    return true;
}
} // namespace

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
    // Direct electron/charge (Zz-row) titrant, matching Reaktoro's own
    // qvar.substance="e-" choice. An O2 titrant (mass-balance coupling on
    // the O row only) would have no direct algebraic relationship to Zz's
    // dual/Eh, so it wouldn't admit a closed-form fixed gradient here -
    // the direct charge titrant is what makes the one-shot joint solve
    // possible.
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

/// Does this system actually have an aqueous phase?
///
/// `pm.LO` (index of the water solvent) is initialised to **0** in
/// ms_multi_file.cpp and only reassigned when some phase carries `ccPH == 'a'`
/// (ms_multi_format.cpp). On an aqueous-free system it therefore STAYS 0 - which
/// is a perfectly valid DC index - so a range test like `pm.LO >= 0 && pm.LO < pm.L`
/// or `pm.LO >= j0 && pm.LO < j1` **cannot distinguish "there is no solvent" from
/// "the solvent is DC 0"**. Every guard here used to make exactly that mistake.
///
/// Measured consequence (2026-08-29, `LBE-6`, the one aqueous-free project we have,
/// `gssss` = a 12-species gas phase plus four single-species solids): the adapter
/// treated **gas species 0 as the aqueous solvent**, applying the solvent-coupling
/// Hessian block to the gas phase and pointing the solvent-collapse detector at the
/// wrong species. It still converged, which is why one project hid this.
///
/// The phase classifier is the correct test. Aqueous systems are unaffected: for them
/// `PH_AQUEL` is present and `pm.LO` is genuinely the solvent, so every guard behaves
/// exactly as before.
bool TMultiBase::HasAqueousPhase() const
{
    for( long int k = 0; k < pm.FIs; k++ )
        if( pm.PHC[k] == PH_AQUEL )
            return true;
    return false;
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
            if( x[pm.LO] >= otherTotal )
                return false; // solvent already dominates - no correction needed
            // Reseed VALUE: NOT simply "total of the other species already
            // in this phase" - confirmed wrong in general (works for one
            // project purely by coincidence, undershoots badly on
            // another where bIC is essentially pure water - see the
            // detailed trace in GEMS3K's CLAUDE.md, "aqueous solvent
            // collapses to its floor", 2026-08-24). Use a physically-
            // grounded bound instead: the most water any system could
            // possibly have is limited by its own total H and O content.
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
                                                long int excludePhaseIdx )
{
    reseedsOut.clear();
    long int j0 = 0;
    for( long int k = 0; k < pm.FIs; k++ ) // multicomponent phases only
    {
        const long int j1 = j0 + pm.L1[k];
        const long int nEnd = j1 - j0;
        if( k != excludePhaseIdx && nEnd > 1 ) // single-species phases have their own PhMinM guard, not this one
        {
            double phaseTotal = 0.;
            for( long int j = j0; j < j1; j++ )
                phaseTotal += x[j];
            // Aqueous carries an EXTRA, stricter threshold
            // (Y[LO]<=pm.XwMinM) on top of the generic pm.DSM one - both
            // guard the SAME "phase reported absent" cliff in
            // PrimalChemicalPotentials(), aqueous is just a strictly
            // stricter instance of the same trap, not a separate one.
            const double threshold = ( pm.PHC[k] == PH_AQUEL ) ? std::max( pm.DSM, pm.XwMinM ) : pm.DSM;
            // Safety margin - require comfortably clear of the cliff, not
            // merely nonzero (matches the spirit of the aqueous-specific
            // "does the solvent DOMINATE its phase" check this
            // generalizes, not just "is it nonzero").
            if( phaseTotal < threshold * 10. )
            {
                // Bound this phase's plausible total by the bulk
                // composition actually available for the ICs its own
                // end-members touch - the same generalization of the
                // water-specific min(bulk-H/2,bulk-O) formula, applied to
                // any phase's own stoichiometry rather than water's
                // specific H:O=2:1 one. Bounded PER END-MEMBER, not once
                // for the whole phase.
                // Taking a single min over every IC that ANY end-member
                // touches lets one trace IC, present in only one minor
                // end-member, cap the entire phase's seed. Measured
                // instance (2026-08-26, the sanidine-albite solvus
                // benchmark): the two ternary feldspars carry an
                // Anorthite (CaAl2Si2O8) end-member while the bulk holds
                // bIC[Ca] = 1.78e-8 mol, so the phase-wide bound came out
                // at Ca's own 3.05e-6 (internally rescaled) and the whole
                // second feldspar was seeded at 3e-6 mol against a true
                // equilibrium amount near 20 - i.e. it was "reseeded"
                // straight back onto the phase-absence cliff this
                // function exists to clear. Bounding each end-member by
                // ITS OWN limiting IC gives Albite 30.5 / Sanidine 14.9 /
                // Anorthite 3e-6, which is both correct and, unlike a
                // phase-wide bound, cannot be dragged down by an
                // end-member the bulk composition cannot support.
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
                    // Each end-member gets its own bound divided by the
                    // end-member count - i.e. still no guess at which one
                    // the true equilibrium favors beyond what the bulk
                    // composition itself already dictates, and still a
                    // long way from a sparse LP vertex (which would
                    // concentrate mass in as few species as possible -
                    // exactly the failure mode this mechanism exists to
                    // avoid). Newton's own gradient-driven iteration is
                    // expected to correct the actual mix from there.
                    for( long int j = j0; j < j1; j++ )
                    {
                        const double b = emBound[(size_t)(j-j0)];
                        if( b <= 0. ) continue;
                        const double perSpecies = std::max( b, threshold * 10. ) / (double)nEnd;
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

bool TMultiBase::LPFeasibilitySeed( std::vector<double>& nOut )
{
    const long int N = pm.N;
    const long int L = pm.L;
    auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
    if( !TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nOut ) )
        return false;
    if( (long int)nOut.size() != L )
        return false;

    // Self-verification before ever trusting the simplex's own output -
    // this method's whole value proposition (a single deterministic seed,
    // no retry-ordering risk) depends on never silently handing back a bad
    // point. Same relative tolerance discipline as the rest of this file
    // (feasTol-equivalent, scaled by the bulk composition's own magnitude).
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
        ipm_logger->warn( "LPFeasibilitySeed: self-check failed (maxResid={}, minVal={}, tol={}) - discarding",
                           maxResid, minVal, tol );
        return false;
    }
    return true;
}

// Dual of the LINEARISED-Gibbs LP:  min sum_j G0[j]*n_j  s.t.  A n = b, n >= 0.
//
// Same simplex, same rows and same feasibility guarantee as LPFeasibilitySeed()
// - only the objective differs, and that difference is the point. The
// feasibility LP's dual is an artefact of "minimise total moles" and says
// nothing about chemistry; this one is the exact dual of a first-order model of
// the real objective, so pricing a species against it
// (s_j = G0[j] - sum_i y_i A[i,j]) is a meaningful "how far from stable is this
// species". Measured on three projects (2026-09-02): |y - U_converged| is
// within 6-12% of the converged dual's own magnitude.
//
// Returns false without touching yOut if the LP or its dual is not trustworthy.
bool TMultiBase::LPGibbsDual( std::vector<double>& yOut, const double* cost )
{
    // `cost` defaults to pm.G0, which is what every shipped caller uses. Passing pm.G instead
    // prices against the CURRENT chemical potentials (G0 + fDQF + F0, i.e. including the mixing
    // and activity terms) rather than the pure standard state - the difference the 29.3 RT dual
    // gap looked like it was made of, and measurably is NOT (plan v5 137.8c/d: re-solving at the
    // CONVERGED potentials leaves the dual 26.5 RT out, because an LP fixes its objective and not
    // its dual). Used only by LpFillProbeReport()'s LPRELP record.
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
                                                   PhStabCensus* censusOut )
{
    const BASE_PARAM *pa_p = base_param();
    const long int L = pm.L;
    std::vector<PhStabViolation> ranked;
    PhStabCensus census;
    census.phases = pm.FI;

    // Phase-assemblage stability check - a second, structurally
    // different correctness signal from the per-species KKT check run
    // at the call site: that check confirms the CURRENT assemblage (which
    // phases are present/absent) is a local stationary point, but
    // cannot by itself tell that assemblage apart from a DIFFERENT,
    // wrong local optimum where some phase that should be present was
    // never let in at all (or vice versa) - exactly the class of error
    // a non-convex NLP can converge to, confirmed real on this
    // integration by the j_CASHNK finding documented in GEMS3K's
    // CLAUDE.md, 2026-08-24 ("Why ROP doesn't just try the fast
    // (water-first) option first"). Native GEMS3K guards against this
    // via PhaseSelectionSpeciationCleanup()'s own PhaseSelect() logic
    // (ipm_chemical.cpp), which reads exactly the same
    // TMultiBase::StabilityIndexes() this reuses - a pure diagnostic
    // computed directly from the already-converged dual (no MBR/IPM/
    // PSSC retry machinery), so calling it here does not violate ROP/
    // AOP's "no native solver calls" design (see this method's own
    // doc comment, ms_multi.h).
    //
    // Confirmed safe to call standalone by reading
    // PhaseSelectionSpeciationCleanup() closely first, per instruction
    // - not just declared-reachable: StabilityIndexes() needs
    // pm.Gamma[] (current, from the LINK_UX_MODE refresh above),
    // pm.Y_la[] (populated by CalculateConcentrations() above, using
    // the FINAL committed pm.U[] - confirmed by tracing
    // CalculateConcentrationsInPhase()'s own Muj=DC_DualChemicalPotential(pm.U,...)
    // call), pm.fDQF[]/pm.sMod[] (loaded once at project init, call-
    // invariant), and pm.K2 (never incremented on this Optima-only
    // path, stays 0 for the whole call - the one branch that reads
    // pm.GamFs[] only fires when K2!=0, so it's simply never taken
    // here, not silently wrong). The one prerequisite nothing above
    // already provides: pm.YF[]/pm.YFA[] are NOT the same arrays as
    // pm.XF[]/pm.XFA[] already refreshed above - "approximation for
    // the next IPM iteration", a genuinely separate pair -
    // StabilityIndexes() reads pm.YF[0]/pm.YFA[0] directly, and
    // PhaseSelectionSpeciationCleanup() always refreshes them from
    // pm.Y right before its own StabilityIndexes() call. pm.Y==pm.X
    // exactly at this point (both set to state.x above), so this is a
    // cheap, exact refresh, not an approximation of one.
    TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
    StabilityIndexes();

    double worstStabilityViol = 0.;
    long int worstStabilityPhase = -1;
    bool worstStabilityWasAbsent = false;
    long int j0stab = 0;
    for( long int k = 0; k < pm.FI; k++ )
    {
        const long int j1stab = j0stab + pm.L1[k];
        // A phase with ANY end-member under a genuine, non-trivial
        // kinetic restriction - DUL[j] < 1e6 (an upper ceiling, however
        // small, including exactly 0) OR DLL[j] > 0.0 (a lower floor
        // forcing presence) - is exempt from this check entirely.
        // Mirrors native's own KinConstrDC/KinConstrPh criterion
        // verbatim (ipm_chemical.cpp:1450-1455: "if(pm.DUL[j]<1e6)
        // KinConstrDC=true; if(pm.DLL[j]>0.0) KinConstrDC=true;", which
        // then `goto NextPhase` - skips PhaseSelect()'s own stability
        // decision for the whole phase entirely). A deliberately,
        // kinetically excluded phase (e.g. Quartz, DUL=DLL=0) is
        // EXPECTED to look thermodynamically "stable but absent"; a
        // phase kinetically forced to REMAIN present via a DLL floor
        // (e.g. a relict mineral the fluid is undersaturated with
        // respect to but hasn't kinetically finished dissolving) is
        // EXPECTED to look "present but unstable" - both are the whole
        // POINT of the restriction, not evidence of a wrong assemblage.
        // The first version of this fix only caught the FULLY fixed
        // (DUL==DLL) case, via a box-degeneracy test - too narrow: a
        // one-sided DLL floor with DUL left unrestricted (a genuine,
        // very common metastability pattern, e.g. Gibbsite in
        // o_Kaolinite: DUL=1e-7, and a synthetic DLL=1e-8 test) is NOT
        // a degenerate box, sits cleanly at its lower bound, and was
        // still misflagged "present but unstable - should be absent"
        // until this broader test replaced it. Found by testing the
        // DLL side directly, GEMS3K/CLAUDE.md 2026-08-26.
        bool kinConstrPh = false;
        for( long int j = j0stab; j < j1stab; j++ )
            if( pm.DUL[j] < 1e6 || pm.DLL[j] > 0.0 )
            { kinConstrPh = true; break; }
        // Same exemption, same reason, for a phase THIS call's
        // phase-extinction retry deactivated: it was removed only after
        // being proved interchangeable with a present twin (identical
        // stoichiometry and identical G0 for every end-member), so it lies
        // on the SAME tangent plane as that twin and therefore ALWAYS
        // reports logSI ~ 0, i.e. "absent but stable". That is a tautology
        // for a redundant duplicate, not evidence of a wrong assemblage -
        // re-admitting it could not lower G, since its twin already offers
        // exactly the same end-members at exactly the same potentials.
        // Confirmed on j_CASHNK, where this check otherwise rejected a
        // result whose composition, volume and G match native's.
        if( !kinConstrPh && L > 0 && exemptSpecies != nullptr )
            for( long int j = j0stab; j < j1stab && j < L; j++ )
                if( exemptSpecies[(size_t)j] ) { kinConstrPh = true; break; }
        j0stab = j1stab;
        if( kinConstrPh ) { census.exempt++; continue; }
        census.scanned++;

        const double logSI = pm.Falp[k];
        // Is this phase's stability index a MEASUREMENT or a guard value?
        // StabilityIndexes() clamps each species' dual activity to
        // [-608, +609] before exponentiating (ipm_chemical.cpp, and it warns
        // on the gems3k channel each time), so a phase whose sum is dominated
        // by a clamped species carries a logSI of ~609/ln10 that says only
        // "overflowed", not "this far from stable". Reconstructed here rather
        // than counted at the source, because a MULTI member would change the
        // ABI: pm.NMU[j] is log(exp(ln_ax_dual)/gamma), and on this path
        // pm.K2 is always 0, so the effective gamma is pm.Gamma[j] (replaced
        // by 1 outside [1e-33, 1e33], exactly as StabilityIndexes() does).
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
        // "Present" is decided STRUCTURALLY as well as by magnitude: a phase
        // every one of whose species sits at its own lower bound was pinned out
        // by the solver and is absent, whatever its total adds up to. The
        // magnitude rule alone scales with the phase's species COUNT while the
        // threshold does not, so a wide phase held entirely at the floor reads
        // as present. Measured on 07PSIna_G_complex_1 (pa_DS = 1e-20, so
        // presenceThreshold collapses onto dcFloor*10 = 1e-12): gas_gen's ten
        // species at exactly dcFloor = 1e-13 sum to exactly 1e-12 and were
        // reported "present but unstable" - the single reason the project
        // returned BAD rather than OK once the dimension reduction let it
        // converge at all. The same signature is on record for
        // 07PSIna_G_simple_0_0_0_150.
        //
        // The test is the one OptimaReducedPreSolve() already uses to read a
        // seed's own support, and it can only move a phase from present to
        // absent, never the other way - so it cannot manufacture a
        // "present but unstable" verdict, only retire one that the magnitude
        // rule invented.
        bool anyOffFloor = false;
        for( long int j = j1stab - pm.L1[k]; j < j1stab && j < L; j++ )
        {
            const double lo = std::max( pm.DLL[j], dcFloor );
            if( pm.Y[j] > lo * ( 1. + 1e-6 ) ) { anyOffFloor = true; break; }
        }
        // ... and a third, PHYSICAL rule, OR-ed with the two above: a phase
        // every one of whose species holds a negligible FRACTION of the most
        // of that species the bulk composition could possibly support is
        // absent, whatever its absolute amount.
        //
        // Why this is needed even with the two rules above. presenceThreshold
        // is max(pm.DSM, dcFloor*10), i.e. an ABSOLUTE amount, and on a project
        // that sets pa_DS = 1e-20 it collapses onto the numerical floor - so a
        // species three orders of magnitude OFF its floor (hence structurally
        // present) but at 1e-16 mol still reads present. Measured case,
        // plan v5 section 44.1: 07PSIna_G_simple_0_0_0_150 at
        // pa_OptimaDcFloor = 1e-20 reports AlOOH(cr) "present but unstable" at
        // YF = 9.64e-17 against a threshold of 1e-19. Its ceiling is set by
        // aluminium (bIC[Al] = 0.098 mol), so x/n_max = 9.8e-16 - negligible by
        // fifteen orders of magnitude, and absent under any eps in the whole
        // usable band.
        //
        // The ceiling: n_max(j) = min over ICs i with a(j,i) > 0 of B[i]/a(j,i)
        // - a pure stoichiometric BOUND (it ignores competition between species
        // for the same element), which is why the threshold sits so far below 1.
        // It is scale-invariant: pm.B[] and pm.Y[] are both in the same
        // internally rescaled frame, so the ratio does not depend on pa_DG.
        //
        // Two traps, both real, both handled (plan v5 section 44.2):
        //  (a) the CHARGE IC must be excluded. On a charge-balanced system
        //      B[Zz] = 0, so every cation would get a ceiling of 0 and
        //      "x < eps*0" would never be true - nothing would ever be absent.
        //      Note DetectPhaseCollapseAndReseed()'s own ceiling loop above does
        //      NOT exclude it; that is a different consumer with a different
        //      guard (b <= 0 skips the end-member), deliberately left alone.
        //  (b) a species whose only positive-coefficient ICs genuinely have
        //      B[i] = 0 has a real ceiling of 0 and must SHORT-CIRCUIT to
        //      absent, since "x < eps*0" is false for any positive x. Not
        //      exercised by any project in the corpus - written anyway.
        //
        // MEASURED on six projects spanning 23-265 species (section 44.3): the
        // smallest ratio among genuinely present species is 1.2e-6 to 3.6e-6 -
        // remarkably constant across utterly different chemistry, which is the
        // scale invariance this rule was proposed for showing up in data - and
        // the largest among numerically negligible ones is <= 7.3e-19. Any eps
        // between about 1e-18 and 1e-7 classifies all six correctly. 1e-12 is
        // the log-centre of that band, roughly six decades clear on each side.
        //
        // OR, not AND, and that is a fact about this consumer rather than a
        // preference (section 44.4): in the AlOOH(cr) case the absolute rule
        // says "present" and only the relative one can overrule it, so an AND
        // would leave the false violation exactly as it is. The DIAGNOSTIC
        // consumer (tools/dimension_ceiling) wants the opposite; the two must
        // not share one predicate.
        //
        // Gateable for the same reason section 33.5's structural rule was: this
        // can only move a phase present -> absent, so it can only RETIRE a
        // "present but unstable" violation, never manufacture one.
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
        const bool present = anyOffFloor && ( pm.YF[k] >= presenceThreshold ) && anySignificant;
        // GEMS3K_PHSTAB_PROBE=<path>: one line per phase per evaluation - amount,
        // stability index, the three presence clauses and the verdict this check
        // would return. Zero-cost when unset (one getenv in a static initialiser,
        // then a null test), same shape and same footing as OPTIMA_BETA_PROBE.
        //
        // KEPT because it answered a question nothing else could (plan v5 section
        // 91): this function returns only the WORST violation, so from outside it
        // is impossible to tell "no violation" from "a violation the repair loop
        // is not allowed to act on, with actionable ones ranked behind it". On
        // 07PSIna_G_complex_1 @ 80 C it showed 289 of 314 phases saturated at
        // StabilityIndexes()' own overflow clamp (609/ln10 - DF = 264.4753),
        // ZERO present-but-unstable, and the gas phase the project's whole
        // failure is about reading present-and-stable at logSI = +60 where
        // native's converged value is -0.332. Reach for it whenever a
        // phase-assemblage verdict needs explaining rather than just observing.
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
        // Same thresholds/logic as native's own PhaseSelect()
        // (ipm_chemical.cpp): a stable phase (logSI > DF) not
        // currently in the assemblage, or an unstable one
        // (logSI < -DFM) that IS, both indicate the wrong assemblage
        // was reached.
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
    // Descending by size, ties broken by phase index - a stable, reproducible
    // order, which matters more than usual here because clamped entries all
    // carry the SAME viol and would otherwise be ordered by nothing.
    std::sort( ranked.begin(), ranked.end(),
               []( const PhStabViolation& a, const PhStabViolation& b )
               { return a.viol != b.viol ? a.viol > b.viol : a.k < b.k; } );
    {
        static FILE* dbg = []() -> FILE* {
            const char* fn = std::getenv( "GEMS3K_PHSTAB_PROBE" );
            return fn ? fopen( fn, "a" ) : nullptr; }();
        if( dbg )
            fprintf( dbg, "CENSUS phases=%ld exempt=%ld scanned=%ld clamped=%ld "
                          "absent-but-stable=%ld present-but-unstable=%ld worst=%ld\n",
                     (long)census.phases, (long)census.exempt, (long)census.scanned,
                     (long)census.clamped, (long)census.absentStable,
                     (long)census.presentUnstable, (long)worstStabilityPhase );
    }
    if( rankedOut != nullptr ) rankedOut->swap( ranked );
    if( censusOut != nullptr ) *censusOut = census;
    violOut = worstStabilityViol;
    wasAbsentOut = worstStabilityWasAbsent;
    return worstStabilityPhase;
}


// ---------------------------------------------------------------------------
// Species-level dimension-reduction pre-solve  (pa_OptimaDimReduce, default 0)
// ---------------------------------------------------------------------------
//
// See BASE_PARAM::OptimaDimReduce (ms_multi.h) for the motivation, the measured
// prize, and why this must omit species rather than pin them. The short form:
// native's linear system is N x N over independent components while the Optima
// path's is (L + R) x (L + R) over species PLUS components, and on the corpus's
// large systems 68-92% of those species are absent at the answer. Pinning them
// (pa_OptimaPhaseCompaction) leaves the factorisation exactly as large - which
// is why that field, measured, does nothing for the size wall. Omitting them is
// the only thing that can.
//
// The mathematics is an ordinary active-set method and is exact at its fixed
// point, not an approximation:
//
//   Hold every omitted species j at its own lower bound l_j. The remaining
//   mass-balance rows are then
//       sum_{s in S} A[i, s] * x_s  =  b_i - sum_{j not in S} A[i, j] * l_j
//   and the reduced problem is minimise G over x_S subject to those rows and
//   the surviving boxes. Its dual y is a full N-vector, because every IC row
//   survives - only columns are dropped - so the omitted columns can be PRICED
//   on it exactly as a simplex prices non-basic columns:
//       s_j = dG/dn_j - sum_i y_i A[i,j] = F[j] - sum_i U[i] A[i,j].
//   s_j >= 0 means "at its lower bound and wanting to stay there", which is the
//   full problem's own KKT condition for a bound-active variable, and is
//   literally the test Optima itself uses (Stability.cpp's is_lower_unstable).
//   So when a pass prices nothing back in, the reduced answer satisfies the
//   FULL problem's KKT conditions - it is a solution, not an approximation.
//
// Two deliberate design choices, both from measured failures on this branch:
//
//  1. THE INITIAL SET COMES FROM TWO LPs, NOT FROM A SOLVER STATE. The tempting
//     alternative - run a short probe and drop whatever sits at the floor - is
//     measured-unsafe: on Resources/gems3k/j_10TH_G_seawater, 4 of the 60
//     species pinned at the floor in the failing run are genuinely present in
//     the answer, dolomite among them at 1.0e-2 mol. A set derived that way
//     omits dolomite and converts an honest failure into a silent wrong answer.
//     So the set is the union of (a) the LP-FEASIBILITY seed's own support,
//     which is feasible by construction (A*Y = b exactly on it), and (b) every
//     species priced below OptimaDimReduceTol against the LP-GIBBS dual (see
//     BASE_PARAM::OptimaDimReduceTol and the call site below). Both are LPs
//     over pm.A/pm.B alone, so both are wholly independent of any solve.
//     (b) exists because (a) on its own is a VERTEX - N species - and that
//     sparsity, measured 2026-09-02, is what broke this feature on the two
//     largest projects in the corpus.
//
//  2. READMISSION IS MONOTONE - a species that has been activated is never
//     dropped again. That is what bounds the loop, and it keeps this free of
//     the active-set cycling that has repeatedly bitten this solver whenever a
//     set was allowed to shrink again (see CLAUDE.md's retry-ordering history).
//     Removal is already covered, after the fact, by the phase-selection repair
//     loop and the phase-extinction tier in CalculateEquilibriumStateOptima().
//
// This never returns an answer on its own. It leaves its primal in pm.Y[] and
// its dual in pm.U[], and the ordinary full-dimension solve runs warm-started
// from both - which costs O(1) Optima iterations for a correct (x,y) pair
// (measured: 3 on f_TestPNTDB, plan-v5 section 27) and keeps every existing
// retry, KKT, mass-balance and phase-stability check running at full dimension
// and full index. A pre-solve that fails, or that cannot omit anything, is
// simply discarded.
bool TMultiBase::OptimaReducedPreSolve( long int maxPasses, double dcFloor, double dimTol,
                                        long int passBudget,
                                        long int& iterationsOut, long int& activeOut )
{
    iterationsOut = 0;
    activeOut = 0;
    const long int N = pm.N;
    const long int L = pm.L;
    if( L <= 0 || N <= 0 || maxPasses <= 0 )
        return false;

    const BASE_PARAM* pa_p = base_param();
    const bool kMoleFracHessian = ( pa_p->OptimaMoleFracHessian != 0 );
    const bool kFDHessian       = ( pa_p->OptimaFDHessian != 0 );
    const bool kFDDiagFloor     = ( pa_p->OptimaFDDiagFloor != 0 );
    const double kLogBarrierTau = pa_p->LogBarrierTau;
    const double kPhaseHessianFloor = pa_p->PhaseHessianFloor;
    const bool hasAq = HasAqueousPhase();
    // pa_OptimaReadmitSeed: 0 = readmit at the floor (the behaviour this has
    // always had); > 0 = seed at White 1958's Eq. 12a predicted amount, with
    // this value capping the growth exponent in RT units. See
    // BASE_PARAM::OptimaReadmitSeed (ms_multi.h).
    const double kReadmitSeed = pa_p->OptimaReadmitSeed;

    // Box bounds, built exactly as the full path builds them - including the
    // "pm.DUL[j] < 1e6 without a `> 0.` guard" convention, so a DUL of exactly
    // zero stays the hard exclusion GEMS3K means it to be.
    std::vector<double> xlo( (size_t)L ), xhi( (size_t)L );
    for( long int j = 0; j < L; j++ )
    {
        xlo[(size_t)j] = std::max( pm.DLL[j], dcFloor );
        xhi[(size_t)j] = ( pm.DUL[j] < 1e6 )
                          ? std::max( pm.DUL[j], dcFloor ) : std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 10.;
    }

    // Initial active set: the seed's own support. A species whose box is
    // degenerate (xhi <= xlo - a hard kinetic exclusion) is never active:
    // omitting it and holding it at its bound is exactly equivalent, and it
    // cannot be priced back in either.
    std::vector<char> act( (size_t)L, 0 );
    for( long int j = 0; j < L; j++ )
    {
        if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;
        if( pm.Y[j] > xlo[(size_t)j] * ( 1. + 1e-6 ) ) act[(size_t)j] = 1;
    }
    // The solvent is never omitted, whatever the seed says about it - the
    // solvent-collapse trap this file documents at length is exactly a
    // transient state reporting water as absent.
    if( hasAq && pm.LO >= 0 && pm.LO < L && xhi[(size_t)pm.LO] > xlo[(size_t)pm.LO] )
        act[(size_t)pm.LO] = 1;

    // ... and then widen that set by PRICING every species against the dual of
    // the linearised-Gibbs LP (see BASE_PARAM::OptimaDimReduceTol for the
    // measurements and the physical reading of the threshold). The seed's own
    // support is a VERTEX - N species - and on 2026-09-02 that sparsity, not
    // the pricing loop, was measured to be what broke this feature on the two
    // largest projects: pass 0 was either unsolvable (07PSIna_G_complex_1, 69
    // species of 1392) or produced a dual so poor that 427 of 641 omitted
    // columns priced back in one step (f_TestPNTDB, 49 of 690). The LP-Gibbs
    // dual is just as independent of any solver state - same simplex, same
    // rows, only the objective differs - and is within 6-12% of the converged
    // dual's own magnitude on all three projects measured.
    //
    // NOT ADOPTED, measured the same day: additionally seeding pass 0's dual
    // from this same LP dual. It helped at one threshold (261 -> 198 iterations
    // on 07PSIna_G_mid_1 at 5) and hurt at the next (269 -> 410 at 10) - not
    // monotone in its own parameter, this solver's familiar signature for a
    // knob that should not be tuned. Pass 0 therefore still starts from a zero
    // dual, and only later passes inherit their predecessor's.
    {
        std::vector<double> lpDual;
        if( dimTol != 0. && pm.G0 != nullptr && LPGibbsDual( lpDual ) )
        {
            // Price every candidate once; the sign of OptimaDimReduceTol then
            // selects how the prices are turned into a set. A POSITIVE value is
            // an absolute threshold in RT (the physical reading documented at
            // that field). A NEGATIVE value is the dimensionless RANK rule
            // -|tol| x N: admit the cheapest-priced species until the active set
            // reaches that multiple of the IC count. The rank rule exists
            // because the three projects the threshold was calibrated on landed
            // at 3.4-4.6 x N regardless of their size, which would make the set
            // size, not the RT cut, the quantity that actually matters.
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
                ipm_logger->info( "OptimaReducedPreSolve: LP-Gibbs pricing admitted {} extra species "
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
                ipm_logger->info( "OptimaReducedPreSolve: LP-Gibbs rank rule admitted {} extra species "
                                   "(target {} = {} x N, seed support {})",
                                   added, target, -dimTol, already );
            }
        }
    }

    std::vector<long int> nxToJ, jToNx( (size_t)L, -1 );
    std::vector<double>   Fbase( (size_t)L, 0. );

    // ONLY A FIXED POINT IS EVER HANDED OVER. An intermediate pass's answer is
    // the solution of a DIFFERENT (smaller) problem, so its dual is confidently
    // wrong about every species not yet active - and this branch has measured
    // twice that a confidently wrong or internally inconsistent dual costs more
    // than no dual at all (the cold-start dual estimate: 4184 it -> 60000; the
    // partial warm dual: 3 it -> 837). Measured here too, before this guard
    // existed: on f_TestPNTDB pass 1 failed to converge and pass 0's 49-species
    // answer was handed over anyway, taking the full solve from 2736 iterations
    // to 4048 and the project from 330 s to 498 s. So on any exit other than
    // "priced nothing back in", restore what the caller had.
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
    apply_optima_linesearch( options, pa_p->OptimaLineSearch );

    // ---- Stall / wall-clock guard for the pre-solve itself ----
    // The full solve in CalculateEquilibriumStateOptima() has carried a stall
    // watch and a pa_OptimaMaxSeconds budget since 2026-08; this pre-solve had
    // neither, and on 2026-09-02 that was measured to be the whole of its worst
    // regression. On f_CASHNK the reduction does not settle: pass 1 runs a
    // 28-of-35-species problem to the full pa_IIM budget and is then discarded,
    // so the project pays ~10060 wasted iterations on top of its ordinary 643
    // (which the full path reaches only via its own stall rescue). 643 -> 10703
    // iterations, 343 ms -> 7464 ms, for an answer that is discarded anyway.
    // Capping a pass that has stopped improving costs nothing - the pass is
    // thrown away in either case - and it bounds what a non-settling reduction
    // can cost. The budget is deliberately SHARED across all passes: a project
    // whose passes each stall should pay one window, not maxPasses windows.
    //
    // Same two non-obvious points as the full path's watch, for the same
    // reasons: `check` returning true means CONVERGED to Optima, so the flag is
    // folded back by forcing succeeded=false (here that lands on the existing
    // `if( !result.succeeded ) return discard()`, which is the correct
    // response); and the test is deliberately conservative about what counts as
    // a stall - it asks whether the pass is DEADLOCKED (best-so-far has not
    // moved AND the live iterate has not moved either, over a whole window),
    // not whether it is converging fast enough. See the test body for the
    // measurement that forced that distinction: a per-step "did best-so-far
    // improve" question kills f_TestSUP98's own converging pre-solve pass,
    // which reaches its answer through a 661-iteration excursion.
    struct PreStallWatch {
        long int window = 0;
        double bestErr = 0., bestComp = 0.;
        // Best-so-far as it stood one whole WINDOW ago, plus the live error's
        // range within the current window - see the test itself for why both
        // are needed and what each one protects against.
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
    if( preWatch->window > 0 || preWatch->maxSeconds > 0. )
        options.convergence.check =
            [preWatch]( Optima::ConvergenceCheckArgs const& args ) -> bool
            {
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
                // The question this asks is "is this pass DEADLOCKED", and it is
                // deliberately not "is it converging fast enough". Over one whole
                // window it declares a stall only when BOTH
                //   (a) best-so-far has not fallen by a meaningful relative amount
                //       against its value at the window's start, AND
                //   (b) the live error has not even MOVED - its range within the
                //       window is below the same relative amount.
                // Either one alone is measurably wrong; the pair is not.
                //
                // WHY (a) ALONE IS WRONG - this is a real false positive, found
                // 2026-09-03 and the reason this test was rewritten. The previous
                // version asked only "did best-so-far improve meaningfully on this
                // STEP", with the same window. On f_TestSUP98 (923 species) the
                // reduced pre-solve's pass 0 CONVERGES, in 1868 iterations, from
                // e=2.60e+03 to 1.06e-08 - but it does so through a 661-iteration
                // excursion (iterations 918-1579) in which best-so-far is EXACTLY
                // frozen at 9.671475 while the live error wanders up to 4.5e+04 and
                // back. A step-based rule at a 500 window fires at iteration 1418,
                // i.e. kills a converging pass three-quarters of the way through.
                // The cumulative-only variant fires at 1500 - same fate. That
                // project does not set pa_OptimaStallWindow today, so nothing
                // regressed; what it blocked was arming the field by default.
                //
                // WHY (b) ALONE IS WRONG: a genuinely slow but SMOOTH converger
                // descending even 0.5 % per iteration traverses a factor of ~12
                // over 500 iterations, so a per-step movement test would call it
                // static and kill it. Measuring the RANGE over the whole window
                // instead of the per-step change is what removes that risk: a run
                // that is going anywhere at all cannot look static over 500 steps.
                //
                // WHAT IT CATCHES, and it is exactly one thing: a pass whose
                // ITERATE HAS STOPPED MOVING. On j_TestPNTDB at
                // pa_OptimaDimReduceTol = -2 - whose pass 0 is discarded after 7001
                // iterations - the live error is frozen to 1.2e-13 relative for
                // 6500 consecutive iterations, so both clauses hold from the second
                // window on and it fires at 1000 (against 7001 unbounded).
                //
                // THE THRESHOLD IS NOT A TUNING PARAMETER. Measured per window
                // across nine passes of five projects, the two quantities are
                // bimodal with a NINE-order gap: a deadlocked window scores 1.2e-13
                // on both clauses, and the smallest score any live window produces
                // is 7.3e-4. Every value from 1e-12 to 1e-4 gives the identical
                // answer on every pass. 1e-8 is used because pa_OptimaTol's own
                // default already is - no new number enters the file.
                //
                // WHAT IT DELIBERATELY DOES NOT CATCH, and why that is correct:
                // 07PSIna_G_vcomplex at 80 C was on record as a third "crawling"
                // regime - real gains, uselessly small, bounded only by pa_IIM. The
                // per-window measurement does not support that reading: its
                // best-so-far falls by 7.3e-4, then 5.1e-2, then 3.7e-1 over three
                // consecutive windows, i.e. the rate is ACCELERATING and the pass
                // is converging - it simply needs more iterations than its budget
                // allows. Firing on it would be a false positive, not a catch. A
                // pass that is genuinely converging too slowly for its budget is a
                // BUDGET question (pa_IIM, pa_OptimaMaxSeconds), not a stall one,
                // and no stall watch should be asked to answer it.
                // 07PSIna_G_edt_2's discarded pass 1 is likewise not caught: its
                // best-so-far is frozen at 1.2e-13 per window, but its live error
                // ranges over a factor of 27, i.e. it carries the same signature as
                // f_TestSUP98's excursion, which recovered. Whether it would
                // recover given budget is unknown, so it is treated as unproven
                // rather than as dead. That costs the bound this watch used to put
                // on that project (1241 iterations); the safety is worth more.
                //
                // `best` is updated on every genuine fall, even a fall too small to
                // count, so the reference the next window is measured against never
                // drifts upward on accumulated noise.
                //
                // Deliberately NOT applied to the full solve's own stall watch
                // below, which has the same shape: that watch is what rescues
                // f_/j_CASHNK (plan v5 section 30.5 - they converge only via it, at
                // ~642 iterations), so changing when it fires changes those
                // projects' results. Different question, bigger blast radius.
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
            // Nothing left to omit - the reduced problem IS the full one, so
            // there is no point paying for it twice. Whatever the previous
            // pass produced (if any) still stands as a warm start.
            ipm_logger->info( "OptimaReducedPreSolve: active set reached the full species "
                               "count ({}) at pass {} - no reduction available", L, pass );
            // The previous pass CONVERGED and re-admission has brought every species back: that answer
            // (re-admitted species at their lower bound) is the warm start for the full solve, as the
            // comment above says. Until 2026-09-14 this branch called discard(), which restores the SEED,
            // so a converged pass was thrown away and the full solve started from scratch. Measured on the
            // one project where it happens, T8ax2_nIC45 (130 species, pass 0 converges on 94 and re-admits
            // 36), AOP over 9 draws at 1e-15: 330/348 median/max iterations -> 161/169, against 183/202
            // with the pre-solve off; unchanged on the 12 other projects of that trial
            // (gems-benchmark Docs/presolve-trial-2026-09-14c.txt).
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
            // Omitted species are held at their lower bound; their (constant)
            // contribution to each IC row moves to the right-hand side.
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

        // Objective. Same physics as the full path's - the gradient of the
        // FULL Gibbs energy with respect to an active species is still just
        // that species' chemical potential (the mu-dependence terms cancel by
        // Gibbs-Duhem), so no reduction-specific correction is needed. The
        // only differences are index mapping and that the phase-decay counters
        // (pa_MbTrendPhaseDecay, default off) are not maintained here; the
        // full solve that follows maintains its own.
        problem.f = [this, L, nS, dcFloor, kLogBarrierTau, kPhaseHessianFloor,
                     kFDHessian, kFDDiagFloor, kMoleFracHessian, hasAq, &nxToJ, &jToNx, &xlo, &act, &Fbase]
                    ( Optima::ObjectiveResultRef res, Optima::VectorView x,
                      Optima::VectorView /*p*/, Optima::VectorView /*c*/,
                      Optima::ObjectiveOptions opts )
        {
            for( long int j = 0; j < L; j++ )
                pm.X[j] = act[(size_t)j] ? x[ jToNx[(size_t)j] ] : xlo[(size_t)j];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );

            // f is the FULL Gibbs energy (omitted species included) so the
            // reported value stays comparable with the full path's; the
            // omitted terms are a constant shift as far as the optimiser is
            // concerned.
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
                // Per-multicomponent-phase ideal-mixing curvature, gathered
                // into the reduced slots. Identical expressions to the full
                // path's - see the long derivation comments there, especially
                // the aqueous solvent ROW, which uses the ideal water-activity
                // convention and is NOT d ln x_w/d n (do not "fix" it).
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
                        // Same present-only restriction as the main objective -
                        // see the long comment there for the measurement.
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
                    const double dIdeal = res.fxx(si,si);
                    pm.X[i] = Xi + h;
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                    for( long int sr = 0; sr < nS; sr++ )
                    {
                        const long int r = nxToJ[(size_t)sr];
                        res.fxx(sr,si) = ( pm.F[r] - Fbase[(size_t)r] ) / h;
                    }
                    // pa_OptimaFDDiagFloor: an FD diagonal that is not positive (exactly 0 when X_i is below
                    // PrimalChemicalPotentials()'s recompute threshold) replaces a positive analytic one and
                    // makes the reduced Hessian indefinite; put the analytic value back (plan v5 s139.5).
                    if( kFDDiagFloor ) { g_fdDiagFloorCols++;
                        if( !( res.fxx(si,si) > 0. ) && dIdeal > 0. ) { res.fxx(si,si) = dIdeal; g_fdDiagFloorHits++; } }
                    pm.X[i] = Xi;
                }

                // pa_PhaseHessianFloor: exact, eigenvalue-floored curvature for
                // the non-aqueous multicomponent phases. Only end-members that
                // are both PRESENT and ACTIVE take part.
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
                                if( SymEigFloorInPlace( blk, nP, regRatio ) )
                                    for( int a = 0; a < nP; a++ )
                                        for( int b = 0; b < nP; b++ )
                                            res.fxx(pres[(size_t)a],pres[(size_t)b]) = blk[(size_t)a*nP+b];
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

        Optima::State state( dims );
        for( long int s = 0; s < nS; s++ )
        {
            const long int j = nxToJ[(size_t)s];
            state.x[s] = std::min( std::max( pm.Y[j], problem.xlower[s] ), problem.xupper[s] );
        }
        // Seed the dual from the previous pass's own answer. ALL OR NOTHING,
        // per the measurement recorded at the full path's own dual seed: a
        // partial dual is worse than no dual at all. Every IC row survives the
        // reduction, so this dual is complete by construction.
        if( haveAnswer )
        {
            bool duals_usable = true;
            for( long int i = 0; i < N; i++ )
                if( !std::isfinite( pm.U[i] ) ) { duals_usable = false; break; }
            if( duals_usable )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -pm.U[i];
        }


        // Per-pass Optima trace for the REDUCED problem, same env var and the
        // same append-mode file as the full path's block in
        // CalculateEquilibriumStateOptima(). Off by default; one getenv per
        // pass on a path that already builds an N x nS matrix.
        //
        // WHY THIS EXISTS (plan v5 section 109.2, work item 21). The full path
        // carried options.output and this one did not, so every trace taken on
        // a project where pa_OptimaDimReduce engages recorded only the tail:
        // ew_exact_fraction.py read `traced = 2` against an ITG of 112-11215 on
        // 13 of 28 gems3k-psina projects, and the correlation with dimension
        // reduction was 28 of 28 with no exception in either direction. AUTO
        // engages at >= 200 species, which is precisely the set the size wall
        // lives on - so the instrument was blind on exactly the projects it was
        // built to measure. A row whose `traced` is far below its `ITG` is that
        // blindness, not a short solve.
        //
        // The names are the ACTIVE SET's, so they change from pass to pass and
        // must be rebuilt here rather than once beside options.maxiters. The
        // Outputter opens with std::ios_base::app (Optima/Outputter.cpp), so
        // each pass appends its own '=' rule and 'Iteration' header and
        // ew_exact_fraction.py splits on those - a reduced block and a full
        // block differ only in their x[..] columns, and the six columns that
        // tool reads (Iteration, f, Error, ||ex||max, ||ep||max, ||ew||max) are
        // in the same places in both.
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
        Optima::Result result = solver.solve( problem, state );
        iterationsOut += (long int)result.iterations;
        if( preWatch->stalled || preWatch->timedOut ) result.succeeded = false;

        if( !result.succeeded )
        {
            // GEMS3K_PRESOLVE_RESID_PROBE=<file>: on a failing pass, append the
            // reduced-gradient ranking of the state Optima gave up at. Kept as
            // zero-cost-when-off infrastructure (one getenv on a path already
            // discarding thousands of iterations), on the same footing as the
            // local Optima checkout's betamin probe, because it answers "WHICH
            // species is this pass stuck on" in a single run and nothing else
            // does. Recomputes exactly the residual the post-solve KKT check
            // computes, so the ranking is the solver's own, not an approximation.
            //
            // What it found (plan v5 section 57): on
            // Resources/gems3k-psina/07PSIna_G_complex_1_0_1_80_0 ONE species -
            // H2O(g), INTERIOR at 18.93 internal units, five orders from either
            // bound - holds 99.0% of the residual, and its value 7.647793e-01 is
            // that species' own driving force out of existence at native's
            // converged dual (7.645308e-01, three significant figures). Not the
            // absent-species/inaccurate-dual shape that had been assumed: the
            // dual is essentially right and the solver simply cannot shrink a
            // phase it should never have carried.
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
            ipm_logger->info( "OptimaReducedPreSolve: pass {} did not converge on {} of {} "
                               "species{} - discarding the reduced pre-solve", pass, nS, L,
                               preWatch->timedOut ? " (wall-clock budget)"
                                                  : ( preWatch->stalled ? " (stalled)" : "" ) );
            // PER-PASS cost, the decomposition the caller's `dimreduce` record
            // cannot carry. `iters` there is a TOTAL over all passes of both
            // attempts, so plan v5 section 120.7's proposed first-attempt budget
            // could not be sized from it: "111-2101 iterations" over 8 passes is
            // compatible with a flat 14-262 per pass and with one 2000-iteration
            // pass among seven cheap ones, and the two imply completely different
            // budgets. `tol` identifies WHICH attempt this pass belongs to
            // (configured rule vs the fallback the caller retries at) without
            // relying on position in the trace.
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

        // Price the omitted columns on this dual and readmit every one whose
        // reduced gradient is negative - the full problem's own KKT test for a
        // variable sitting at its lower bound, and literally Optima's own
        // is_lower_unstable condition.
        //
        // Readmitted species stay AT their lower bound rather than being seeded
        // strictly inside the box. Both alternatives to this plain rule were
        // implemented and measured on 2026-09-02, and both were WORSE - see
        // plan-v5 for the numbers:
        //   - partial pricing (readmit only the max(nS,N) most negative per
        //     pass, so an early inaccurate dual cannot overshoot): the smaller
        //     intermediate sets are themselves harder to solve - on
        //     07PSIna_G_mid_1 a 52-species pass failed where the 93-species one
        //     it replaced converged in 303 iterations;
        //   - interior seeding of readmitted species at their own
        //     bulk-composition bound x 1e-6, which is what ORCHESTRA's
        //     insertion pass does (activate at 1e-3, never at a floor): on the
        //     same project pass 1 went 303 -> 2263 iterations and pass 2 then
        //     failed.
        // Neither response was monotone in its own parameter, which is this
        // solver's now-familiar signature for a knob that should not be tuned.
        // Do not re-try either without a genuinely new hypothesis.
        //
        // THE ONE EXCEPTION, and it is the hypothesis that comment asked for:
        // pa_OptimaReadmitSeed (default 0 = off, so the paragraph above still
        // describes the shipped behaviour) seeds at White 1958 Eq. 12a's
        // PREDICTED amount rather than at an arbitrary fraction of the bulk
        // bound. The 1e-6 measurement above is evidence against an arbitrary
        // seed; it says nothing about the thermodynamically predicted one.
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

            // pa_OptimaReadmitSeed (default 0 = off): start the readmitted
            // species at the amount its own chemical potential predicts,
            // instead of at the floor. White, Johnson & Dantzig (1958) Eq. 12a;
            // the derivation, the species-class restriction and the reason the
            // knob is the exponent are all at BASE_PARAM::OptimaReadmitSeed.
            //
            // pm.Y[j] is what the next pass's Optima::State is built from (the
            // clamp into [xlower, xupper] happens there), so writing it here is
            // the whole mechanism - nothing else in this function changes.
            if( kReadmitSeed > 0.
                && ( pm.DCCW[j] == DC_SYMMETRIC || pm.DCCW[j] == DC_ASYM_SPECIES ) )
            {
                const double xcur = pm.X[j];             // == xlo[j] for an omitted species
                if( xcur > 0. )
                {
                    // Stoichiometric ceiling: the most of this species the bulk
                    // composition could support. Same bound
                    // DetectPhaseCollapseAndReseed() computes, and the only one
                    // of the three clamps that is a statement about the system
                    // rather than about arithmetic.
                    double nmax = std::numeric_limits<double>::infinity();
                    for( long int i = 0; i < N; i++ )
                    {
                        const double aij = pm.A[ i + j*N ];
                        if( aij > 0. )
                            nmax = std::min( nmax, pm.B[i] / aij );
                    }
                    const double grow = std::min( -sj, kReadmitSeed );   // sj < 0 here
                    double xpred = xcur * std::exp( grow );
                    if( nmax > 0. && xpred > nmax ) xpred = nmax;
                    xpred = std::min( std::max( xpred, xlo[(size_t)j] ), xhi[(size_t)j] );
                    if( xpred > pm.Y[j] ) { pm.Y[j] = xpred; seeded++; }
                }
            }
        }

        ipm_logger->info( "OptimaReducedPreSolve: pass {} - {} of {} species active, "
                           "{} Optima iterations, {} readmitted ({} seeded above the floor)",
                           pass, nS, L, result.iterations, readmitted, seeded );
        // The converged half of the same decomposition. `readmit=0` marks the
        // fixed point that ends the loop, so a reader can tell a pre-solve that
        // SETTLED from one that merely ran out of passes without consulting the
        // caller's record.
        native_trace_decide( "dimreducepass pass=%ld tol=%.6g ns=%ld of=%ld iters=%ld ok=1 "
                             "readmit=%ld seeded=%ld budget=%ld",
                             (long)pass, dimTol, (long)nS, (long)L, (long)result.iterations,
                             (long)readmitted, (long)seeded, (long)passBudget );

        if( readmitted == 0 )
            return true;   // fixed point: the reduced answer satisfies the full KKT conditions
    }

    ipm_logger->warn( "OptimaReducedPreSolve: readmission did not settle within {} passes - "
                       "discarding (an unsettled active set is not a fixed point, and its dual is "
                       "wrong about everything still omitted)", maxPasses );
    return discard();
}

double TMultiBase::CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool reaktoroMode,
                                                    bool runKinetics )
{
    // Disable the IPM-2 chemical-potential smoothing for the whole of this
    // call. That blend (ipm_chemical.cpp, DC_PrimalChemicalPotentialUpdate())
    // mixes the current potential with the value pm.F0[j] held from the
    // PREVIOUS evaluation, and Optima's objective callback re-enters
    // CalculateActivityCoefficients(LINK_UX_MODE) - hence SetSmoothingFactor()
    // and that blend - on every Newton iteration. With a smoothing factor
    // s < 1 the objective is then a moving, history-dependent target rather
    // than a function of x alone, which no Newton method can converge
    // against. See TMultiBase::optima_disable_smoothing (ms_multi.h) and
    // CLAUDE.md 2026-08-25 (plan-v5 Phase A / A.1) for the full trace.
    //
    // RAII rather than a plain assignment because this function has many
    // exits, including thrown Error(...) on non-convergence - the flag must
    // not leak into a subsequent native AIA/SIA call on the same TMultiBase.
    struct SmoothingGuard {
        TMultiBase* m;
        bool prev;
        explicit SmoothingGuard( TMultiBase* mm )
            : m(mm), prev(mm->optima_disable_smoothing)
{ m->optima_disable_smoothing = true; }
        ~SmoothingGuard() { m->optima_disable_smoothing = prev; }
    } smoothingGuard( this );

    double ScFact = 1.;
    const BASE_PARAM* pa_p = base_param();
    // Reuse GEMS3K's own tolerances rather than adding parallel
    // Optima-specific BASE_PARAM fields: DHB is the numerical DC-amount
    // floor, IIM/DK are the "max iterations"/"convergence tolerance"
    // knobs (below), DW gates hard-error vs. soft-BAD on non-convergence.
    // pa_OptimaDcFloor > 0 decouples this from pa_DHB - see that field for why
    // reusing pa_DHB made a low floor untestable (it is also native's
    // mass-balance tolerance, so lowering it breaks the reference first).
    const double dcFloor = pa_p->OptimaDcFloor > 0. ? pa_p->OptimaDcFloor
                                                    : std::max( pa_p->DHB, 1e-300 );

    InitalizeGEM_IPM_Data();

    // Redundant species are held at zero for this call only (ExcludeRedundantDCs, ipm_main.cpp);
    // the guard restores their metastability settings on every exit, exceptions included.
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

    // One kinetics/metastability time step, at the SAME position native runs it
    // (ipm_simplex.cpp): after InitalizeGEM_IPM_Data() and ExcludeRedundantDCs(), and BEFORE the
    // internal rescaling, so TKinMet sees the caller's real units. Added 2026-09-20 - until then
    // this path never called it, so the Additional Metastability Restrictions were never updated
    // and a kinetically controlled phase could not change on AOP/SOP at all.
    // runKinetics is false only for CalculateEquilibriumStateHOP()'s Optima leg, whose native leg
    // has already advanced the step.
    if( runKinetics )
        RunKineticsStep();

    if( pa_p->DG > 1e-5 )
    {
        ScFact = SystemTotalMolesIC();
        ScaleSystemToInternal( ScFact );
    }
    // An element with less material than its species' floor amounts can hold (ms_multi.h). Once per call, after
    // rescaling, so pm.B, pm.DLL and dcFloor are in the same internal units.
    SubFloorElementCheck( dcFloor );

    try
    {
        // Allocates/parametrizes each multicomponent phase's TSolMod
        // instance (phSolMod[k]) at the current T,P - a prerequisite for
        // any later LINK_UX_MODE call. The native path gets this from its
        // own GEM_IPM_Init() (ipm_simplex.cpp); this path has no such
        // step, so it must be called explicitly here.
        CalculateActivityCoefficients( LINK_TP_MODE );

        // pm.pNP is set by TNode::GEM_run() before this call, exactly as it
        // is for the native AIA/SIA path: 0 = cold start (AOP, like AIA),
        // 1 = warm start reusing the existing pm.Y[] (SOP, like SIA).
        // Deliberately Optima-only: no native MBR/IPM/PSSC refinement is
        // used to seed or correct this solve (see GEMS3K's CLAUDE.md,
        // 2026-08-23, for why a native-seeded two-pass variant was tried
        // and then rejected in favor of keeping this path a genuine,
        // self-contained alternative solver rather than a wrapper around
        // native GEMS3K).
        // Cold start (pNP==0), for BOTH AOP and ROP: the LP-feasibility seed.
        // Adopted 2026-08-25 after a one-variable A/B measured on the full
        // 25-project sweep (GEMS3K/CLAUDE.md) - it replaced
        // AutoInitialApproximation()+DC_RaiseZeroedOff() on the AOP path, which
        // produced a degenerate LP vertex (all of a solid solution's mass in one
        // end-member, every other end-member pinned at the insertion floor).
        // Paired with the now-unconditional FD Hessian below it removed all four
        // "OK at a worse G" Solvus rows AND cut total iterations 47862 -> 43389.
        // SOP's warm start (pNP!=0, reuse the incoming pm.Y) is untouched.
        //
        // ...EXCEPT when there is nothing to warm-start FROM. Added 2026-08-28
        // after NEED_GEM_SOP on a FRESH node was found to be a silent foot-gun,
        // strictly worse than either real mode. TNode::GEM_run() sets pm.pNP=1
        // for NEED_GEM_SOP unconditionally (node.cpp), so on a node that has
        // never solved anything this branch was skipped and the solve started
        // from whatever speciation the project's .dbr file happened to ship -
        // frozen at that file's own state point, with an all-zero dual, and no
        // LP-feasibility seed. That is not a warm start; it is a cold start
        // from stale data. It is what made solvus_sweep_test --modes sop look
        // catastrophic (104 of 301 temperatures failing, Tc reported 105 C too
        // low) - a harness that builds a fresh TNode per temperature BY DESIGN,
        // so no state can carry over and the .dbr speciation gets progressively
        // more wrong as T walks away from where it was saved. That result was
        // originally misread as warm-start hysteresis; there was no carry-over
        // for hysteresis to happen in. See GEMS3K/CLAUDE.md, 2026-08-27/28.
        //
        // Detector: pm.U[] identically zero. It is the dual - IC chemical
        // potentials over RT - so after ANY real solve by ANY solver on this
        // instance at least one entry is nonzero and large; all-zero means no
        // solver has ever run here. Deliberately NOT "has this Optima path run
        // before": native AIA/SIA writes pm.U[] too (ipm_main.cpp), and handing
        // a converged NATIVE state to SOP is a legitimate, useful warm start
        // that must keep working (it is the cheapest leg of the CTest guard in
        // gems-benchmark's optima_regression.h). Non-finite is treated the same
        // way, for the same reason the dual seed below refuses it.
        //
        // On detection: log loudly AND fall back to the cold path, which is
        // strictly better than what happened before (an LP-feasible seed with a
        // self-consistent zero dual, rather than stale data with a zero dual).
        // Note the two halves must agree - the dual seed below is gated on the
        // same flag, because seeding a dual for a primal we just replaced would
        // be exactly the inconsistent-start case measured as worse than either
        // consistent one (see that block's own comment).
        bool warmStateUsable = false;
        if( pm.pNP != 0 )
        {
            for( long int i = 0; i < pm.N; i++ )
            {
                if( !std::isfinite( pm.U[i] ) ) { warmStateUsable = false; break; }
                if( pm.U[i] != 0. ) warmStateUsable = true;
            }
            if( !warmStateUsable )
                ipm_logger->warn( "CalculateEquilibriumStateOptima: warm start (SOP) requested but "
                                  "this node carries no previous solution (pm.U[] is all zero) - "
                                  "the .dbr file's stored speciation is NOT a warm start. Falling "
                                  "back to the cold (LP-feasibility) seed. Use AOP for a first "
                                  "solve, or SOP only on a node that has already solved." );
        }

        if( pm.pNP == 0 || !warmStateUsable )
        {
            // ROP's own initial guess - superseded, same day, 2026-08-24:
            // a plain uniform-tiny seed (Reaktoro's own ChemicalState
            // default, "n.setConstant(nspecies, 1e-16)") is no longer used
            // here directly. It needed two retries (toggle, then a
            // GEMS3K-specific water reseed) to reach a trustworthy answer
            // on most tested systems, and reordering those retries turned
            // out to be UNSAFE in general (see GEMS3K's CLAUDE.md,
            // 2026-08-24, "Why ROP doesn't just try the fast (water-first)
            // option first" - two independent reordering attempts both
            // broke j_CASHNK, because different retry orders are different
            // Newton starting points that can converge to different local
            // optima on this non-convex NLP). The real fix identified
            // there: replace the SEED itself with something already close
            // to a physically sensible point, removing the need for most
            // retries (and the ordering risk they carry) rather than
            // trying a third retry ordering.
            //
            // LPFeasibilitySeed() (ms_multi.h/above) solves the genuine LP
            // min sum(n) s.t. A*n=b, n>=0 - the same class of computation
            // Reaktoro's own compared-against iteration counts (38/77 on
            // f_CalcDolo/j_Flowline) actually started from
            // (reaktoro_bench.py's build_seed(), scipy.optimize.linprog),
            // not a bare uniform value. A plain min-sum LP vertex is known
            // to concentrate mass in as few species as possible - the same
            // mechanism that made AutoInitialApproximation()'s own
            // LP-simplex seed starve the aqueous phase in the original
            // solvent-collapse trap - so the LP output is run through the
            // SAME solvent-/phase-collapse detection-and-redistribution
            // already validated for the retries below, applied ONCE, up
            // front, deterministically (not as a retry: this seed
            // construction has no ordering choice to get wrong, unlike two
            // separate solver calls racing different starting points).
            //
            // Falls back to the old uniform-tiny seed (unchanged) if the
            // LP itself fails (infeasible/numerical issue, self-checked
            // inside LPFeasibilitySeed()) - in that case the retry logic
            // below is still there as a safety net, exactly as before this
            // change.
            std::vector<double> lpSeed;
            if( LPFeasibilitySeed( lpSeed ) )
            {
                for( long int j = 0; j < pm.L; j++ )
                    pm.Y[j] = lpSeed[j];
                if( HasAqueousPhase() && pm.LO >= 0 && pm.LO < pm.L )
                {
                    // pm.DUL[j] < 1e6 is GEMS3K's own established convention for
                    // "this DC carries a real (possibly zero) kinetic/metastability
                    // upper restriction" - matches ipm_chemical.cpp's own
                    // "if( pm.DUL[j] < 1e6 ) KinConstrDC = true;" verbatim. A bare
                    // `> 0.` guard (as this read before) wrongly treats a DUL of
                    // EXACTLY 0 - a hard exclusion, e.g. a kinetically suppressed
                    // mineral - as "unset", falling through to the large default and
                    // silently discarding the restriction. See the matching fix at
                    // problem.xupper[j]'s own construction below for the case this
                    // was actually caught on (o_/t_Kaolinite's Quartz, DUL=0).
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
                if( DetectPhaseCollapseAndReseed( pm.Y, seedReseeds, aqueousPhaseIdx ) )
                    for( const auto& rs : seedReseeds )
                        pm.Y[ rs.first ] = rs.second;
            }
            else
            {
                ipm_logger->warn( "CalculateEquilibriumStateOptima: LP-feasibility seed unavailable "
                                   "(infeasible or numerical issue) - falling back to the uniform seed" );
                const double reaktoroSeed = std::max( dcFloor, 1e-16 * ScFact );
                for( long int j = 0; j < pm.L; j++ )
                    pm.Y[j] = reaktoroSeed;
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
        // Off by default. When on, this solves the equilibrium over a reduced
        // set of species first and leaves its primal in pm.Y[] and its dual in
        // pm.U[]; the ordinary full-dimension solve below then warm-starts from
        // both and, because a correct (x,y) pair costs O(1) Optima iterations,
        // becomes a cheap verification that keeps every existing retry and
        // post-solve check running at full dimension and full index. Skipped
        // with control conditions active (their virtual slots are a separate
        // mechanism) and in ROP, which is a faithful port of Reaktoro's own.
        // See BASE_PARAM::OptimaDimReduce and OptimaReducedPreSolve() above.
        // COLD STARTS ONLY (pm.pNP == 0). A warm start already carries the very
        // thing this pre-solve exists to manufacture - a consistent (primal,
        // dual) pair - and costs O(1) Optima iterations because of it (3 on
        // f_TestPNTDB, section 27). Running a reduction in front of that is
        // pure waste when it discards, and actively destructive when it
        // settles, since it would overwrite a perfect warm start with the
        // reduced problem's own answer. Without this gate, setting
        // pa_OptimaDimReduce on a project would silently degrade every SOP call
        // on it - which is exactly what the suite's own warm-restart case for
        // f_TestPNTDB guards, and that project is the first one the field is
        // enabled on.
        //
        // ...WITH ONE EXCEPTION: the HOP leg (optima_hop_leg, ms_multi.h,
        // where the reasoning is written out in full). There the incoming
        // warm pair is NATIVE's, not an Optima fixed point, and native decides
        // absence twenty orders of magnitude below Optima's box floor - so
        // every species native calls absent arrives sitting AT that floor with
        // a reduced gradient Optima has never tested, and on a large system
        // those are exactly what its max-norm residual is made of. The
        // pre-solve there is not manufacturing a (primal, dual) pair that
        // already exists; it is re-expressing native's own assemblage at a
        // dimension where the residual is not dominated by species that are
        // not present. It needs no new selection rule to do it:
        // OptimaReducedPreSolve() builds its initial active set from pm.Y[],
        // which on this path is native's converged answer.
        //
        // pa_OptimaDimReduce is THREE-VALUED: > 0 is an explicit pass count,
        // < 0 is explicitly off, and 0 - the compiled default, and what every
        // project that has never heard of this field carries - is AUTO: on
        // above a species-count gate, off below it. The gate exists because
        // section 33.1's corpus sweep found the payoff to be entirely a
        // function of how much the reduction can OMIT, which tracks size:
        // 85-96 % of species stay active below ~35 species, so the pre-solve
        // there prices and re-solves almost the whole list and the ordinary
        // solve then repeats the work; from 122 species up only 49-56 % stay
        // active. kDimReduceAutoMinDC sits above every measured loss
        // (j_GEOTHERM, 154 species, 2.7x) and below every large win (mid_1,
        // 265 species, 60x; f_TestPNTDB, 690, 140x; f_/j_TestSUP98, 923, from
        // never-finished to under a minute). Nothing in the corpus lies
        // between 154 and 265 species, so its exact value is interpolated
        // rather than measured - which is also why it is a named constant with
        // an explicit off switch rather than a silent hardcode.
        // The gate itself now lives in ms_multi.h (optima_dimreduce_passes) so
        // that native_trace_run_header()'s EFF line resolves it exactly as this
        // call site does. Before that they could not disagree only by luck - and
        // nothing in the trace said which value had run. Plan v5 section 95.4.

        // Fallback initial-set rule, used ONCE if the configured rule's
        // pre-solve is discarded. See BASE_PARAM::OptimaDimReduceTol for the
        // measurement this implements: the threshold rule and the rank rule
        // fail on DISJOINT projects, and no single value of either is safe
        // across the corpus, so the way to make them both usable is to retry
        // under the other rule rather than to keep tuning one of them.
        // A discarded pre-solve has already restored pm.Y/pm.U, so a second
        // attempt costs only the wasted pass - the same probe-then-commit
        // pattern pa_OptimaFDHessianDelay already uses.
        //
        // WHAT THE WASTED PASS ACTUALLY COSTS, measured rather than assumed:
        // one FULL pre-solve budget, max(2000, pa_IIM). It is NOT bounded by
        // the stall watch in general - on j_TestPNTDB, which does set
        // pa_OptimaStallWindow = 500, the discarded pass still ran to 7001
        // iterations because it was still improving by the watch's two-signal
        // test while not converging. The trade is still strongly favourable
        // there (7496 iterations / 18.3 s with the fallback, against 9798 /
        // ~400 s without it, same G), but size the expectation from the budget,
        // not from the watch.
        static const double kDimReduceDefaultTol  =  10.;   // the shipped threshold
        static const double kDimReduceFallbackRank = -3.;   // 3 x N: of the rank values
                                                            // measured, the only one that never
                                                            // fails where another rank value
                                                            // succeeds (2 N fails j_TestPNTDB,
                                                            // 2.5 N fails f_TestSUP98). It does
                                                            // not rescue 07PSIna_G_complex_1 at
                                                            // 80 C, but nothing does.

        const long int dimReducePasses = optima_dimreduce_passes( pa_p->OptimaDimReduce, L );

        // On the HOP leg the reduction is EXPLICIT OPT-IN ONLY - AUTO does not
        // reach it. Measured 2026-09-04 on 07PSIna_G_mid_1 (265 species, so
        // above the AUTO gate), same binary, native's own answer handed over
        // either way:
        //     warm verification at full dimension :  2 Optima iterations,  6 ms
        //     reduction first, then verification  : 35 Optima iterations, 85 ms
        // i.e. exactly the "waste when it discards, destructive when it
        // settles" the cold-start-only rule above predicts - HOP's warm
        // verification there is ALREADY O(1), so there is nothing for a
        // reduction to buy and it charges 14x on the Optima leg for the same
        // answer. The reduction pays on the HOP leg only where that
        // verification does NOT work, and section 39.4 established that no
        // static property of a project predicts which of those it is (set size
        // least of all: 164 species solvable, 205 not, 246 solvable). So this
        // is per-project, like every other knob on this path whose response is
        // not uniformly favourable - and it reuses pa_OptimaDimReduce's own
        // three-valued convention rather than adding a field: > 0 means "yes,
        // on this project, including its HOP leg", 0 (AUTO) stays cold-start
        // only. Default HOP behaviour is therefore byte-identical to before.
        const bool hopReduce = optima_hop_leg && pa_p->OptimaDimReduce > 0;

        long int dimReduceIters = 0;
        bool dimReduceDone = false;
        if( !reaktoroMode && R == 0 && ( pm.pNP == 0 || hopReduce ) && dimReducePasses > 0 )
        {
            // Re-establish the same consistent (Y, X, XF/XFA, activity
            // coefficients) state the seed block above leaves behind. The
            // pre-solve's own objective callback mutates all of those while it
            // runs, and on the discard path pm.Y[] is still the seed - so this
            // is correct whether it succeeded or not, and makes "discarded"
            // mean genuinely discarded. It also has to run BETWEEN the two
            // attempts, so the fallback starts from exactly the state the first
            // attempt did.
            auto reestablish = [&]() {
                TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
                for( long int j = 0; j < L; j++ )
                    pm.X[j] = pm.Y[j];
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
            };

            const double configuredTol = pa_p->OptimaDimReduceTol;
            long int nActive = 0;
            // First attempt at the lowered per-pass budget (pa_OptimaPreSolveFirstIters,
            // plan v5 s122.7); the fallback below keeps max(2000, pa_IIM).
            dimReduceDone = OptimaReducedPreSolve( dimReducePasses, dcFloor, configuredTol,
                                                   optima_presolve_pass_budget(
                                                       pa_p->OptimaPreSolveFirstIters,
                                                       (long int)pa_p->IIM, true ),
                                                   dimReduceIters, nActive );
            reestablish();

            // WASTED work, tracked separately from total. A discarded attempt's
            // iterations are still paid and still land in `iters`, so a total
            // alone cannot distinguish "this problem is hard" from "the first
            // admission rule did not converge and everything it spent was thrown
            // away". See the DECIDE records below - plan v5 section 115.
            long int dimReduceWasted = 0;
            long int dimReduceAttempts = 1;

            if( !dimReduceDone )
            {
                // The configured rule was discarded - try the other one once.
                // A negative configured value is the rank rule, so its fallback
                // is the shipped threshold; anything else (the threshold rule,
                // or 0 = seed support only) falls back to the rank rule.
                const double fallbackTol = ( configuredTol < 0. )
                                            ? kDimReduceDefaultTol : kDimReduceFallbackRank;
                long int fallbackIters = 0;
                ipm_logger->info( "CalculateEquilibriumStateOptima: dimension-reduction pre-solve "
                                   "discarded at tol={} after {} iterations - retrying once at "
                                   "tol={}", configuredTol, dimReduceIters, fallbackTol );
                // THE RECORD THIS SITE WAS MISSING, and it is the single most
                // expensive decision the solver makes on a large system.
                // Measured on T14_ball000 (plan v5 section 115): the project's AOP
                // cost is BIMODAL under a 1e-15 bIC nudge - 405 iterations on six
                // draws of nine and ~10748 on two - and the entire difference is
                // whether THIS branch is taken. On the expensive draw pass 1 fails
                // to converge on 299 of 1167 species, the whole attempt is
                // discarded after 10353 iterations, and the retry then produces
                // the SAME answer in ~394. So ~96 % of that run is work thrown
                // away, `G` is bit-identical either way, and until now the only
                // trace of it was an spdlog INFO line - absent from the trace,
                // absent from the freeze's `# dec` column, and therefore invisible
                // to anything scoring a benchmark. Both draws printed an
                // IDENTICAL-looking `dimreduce done=1 active=299 of=1167 passes=8`
                // and differed only in `iters`, which reads as "this run was
                // harder" rather than "this run wasted 10353 iterations".
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
                ipm_logger->info( "CalculateEquilibriumStateOptima: dimension-reduction pre-solve "
                                   "produced a warm start over {} of {} species in {} iterations",
                                   nActive, L, dimReduceIters );
            // This is the record whose absence cost plan v5 section 95.4: the whole
            // T14 ladder was measured against pa_OptimaPhaseCompaction while this ran
            // on every rung, and nothing anyone diffs said so.
            // `attempts` and `wasted` added 2026-09-11 (section 115) so the total is
            // DECOMPOSABLE: iters = wasted + productive. A row whose `wasted` is most
            // of its `iters` is not a hard problem, it is a discarded attempt, and the
            // two were indistinguishable in this record for as long as it has existed.
            //
            // AND IF THE RETRY FAILED TOO, ALL OF IT IS WASTE. Measured on T-cement
            // the next day (section 119): both attempts run to the 10000-iteration cap
            // and end at `active=0 of=207 done=0`, so the pre-solve hands over nothing
            // and 20000 of that call's 33000 iterations buy nothing at all - while the
            // field as first written reported `wasted=10000`, counting attempt 1 only.
            // A decomposition that under-reports waste on exactly the runs where the
            // waste is total is the same defect this record was added to fix, one level
            // down: `done=0` is what makes the remainder unproductive, so read it here
            // rather than leaving the reader to multiply.
            if( !dimReduceDone )
                dimReduceWasted = dimReduceIters;
            native_trace_decide( "dimreduce done=%d active=%ld of=%ld iters=%ld passes=%ld "
                                 "attempts=%ld wasted=%ld",
                                 dimReduceDone ? 1 : 0, (long)nActive, (long)L,
                                 (long)dimReduceIters, (long)dimReducePasses,
                                 (long)dimReduceAttempts, (long)dimReduceWasted );
        }

        // GEMS3K_PRESOLVE_HANDOVER_PROBE=<path>: dump EVERY piece of state the main
        // solve is about to start from, so the pre-solve-on and pre-solve-off arms
        // can be diffed directly.
        //
        // THE QUESTION IT EXISTS FOR (plan v5 s120.5, handoff 2026-09-11c item 3):
        // a DISCARDED pre-solve still changes the primary solve - on T-cement at
        // k=+2 the run converges in the AUTO arm and not with the reduction off,
        // although discard() restores pm.Y and pm.U element for element. So
        // something OTHER than the primal and dual is carried across, and the
        // available evidence could not say what. Reasoning about it from the source
        // has already produced two wrong candidates (the IPM-2 smoothing blend,
        // which SmoothingFactor() makes an algebraic no-op on this path; and pm.X,
        // which reestablish() rewrites from pm.Y) - so dump the state and diff it
        // rather than arguing about it. Placed OUTSIDE the pre-solve block on
        // purpose: in the off arm the block does not run, and the whole point is to
        // compare the two arms at the same instant.
        //
        // Deliberately NOT a checksum. A checksum answers "did anything move",
        // which is already known; the open question is WHICH array, and a
        // per-element dump lets an ordinary diff localise it and name the species.
        if( const char* hp = std::getenv( "GEMS3K_PRESOLVE_HANDOVER_PROBE" ) )
        {
            if( FILE* fh = fopen( hp, "a" ) )
            {
                // The header has to say WHICH CALL this is. The probe fires on
                // every Optima call that reaches this point - AOP, SOP, HOP and
                // SHP all do - so a run of rop_compare appends four blocks, and
                // without pNP/hop/reaktoro to separate them a diff would be
                // comparing a cold arm against a warm one and reporting the mode
                // difference as the finding.
                fprintf( fh, "# HANDOVER L=%ld N=%ld FI=%ld pNP=%ld hop=%d reaktoro=%d "
                             "IT=%ld ITG=%ld ITF=%ld K2=%ld FitVar3=%.17g FitVar4=%.17g\n",
                         (long)L, (long)N, (long)pm.FI, (long)pm.pNP,
                         optima_hop_leg ? 1 : 0, reaktoroMode ? 1 : 0, (long)pm.IT,
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

        // Resolve each active condition's fixed objective-gradient value
        // once, up front (needs this call's G0[]/T, only valid now that
        // InitalizeGEM_IPM_Data() has run) - assign its unknown slot.
        std::vector<double> fixedGrad( R, 0. );
        for( long int k = 0; k < R; k++ )
        {
            conditions[k].slot = L + k;
            fixedGrad[k] = conditions[k].fixedGradientFn( conditions[k].target );
        }

        Optima::Dims dims;
        dims.x  = L + R;
        dims.be = N;
        // Sensitivity parameters c := the bulk composition b, so that
        // Sensitivity::xc is dn/db directly (see optima_want_sensitivity).
        // Declared only on request - dims.c = 0 keeps the problem exactly as
        // it was, and Optima skips the sensitivity solve entirely.
        const bool wantSens = optima_want_sensitivity;
        if( wantSens ) dims.c = N;

        Optima::Problem problem( dims );

        for( long int i = 0; i < N; i++ )
            for( long int j = 0; j < L; j++ )
                problem.Aex(i, j) = pm.A[ i + j*N ];

        if( wantSens )
        {
            // be[i] IS b[i], so d(be)/dc is the identity and c carries the
            // current bulk composition. Both are set before the solve; Optima
            // needs bec to form the sensitivity system, and c only so the
            // parameters have a meaningful value.
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

        // Box constraints: reuse GEMS3K's own DC kinetic-restriction
        // bounds directly, floored by dcFloor (pa_DHB) for well-posedness
        // of the log-based Hessian approximation.
        //
        // pm.DUL[j] < 1e6 (no `> 0.` guard) is GEMS3K's own established
        // convention for "this DC carries a real, possibly-zero kinetic/
        // metastability upper restriction" - see ipm_chemical.cpp's own
        // "if( pm.DUL[j] < 1e6 ) KinConstrDC = true;". The previous `> 0.`
        // guard here silently discarded a DUL of EXACTLY 0 - a hard, common
        // exclusion (e.g. a kinetically suppressed mineral, DUL=DLL=0) -
        // falling through to the large "unrestricted" default and letting
        // Optima grow a formally-excluded species freely. Confirmed on
        // o_/t_Kaolinite: Quartz ships DUL=DLL=0 (kinetically forbidden in
        // favor of metastable SiO2(am), which native respects), and AOP was
        // growing it to ~0.033 mol regardless - not a "different, better
        // assemblage" as earlier sessions assumed, but this bug. Floored by
        // dcFloor, matching xlower's own floor immediately below, so a
        // DUL=DLL=0 species gets the well-posed degenerate box
        // [dcFloor,dcFloor] rather than an inverted, infeasible one
        // (xlower=dcFloor > xupper=0).
        for( long int j = 0; j < L; j++ )
        {
            problem.xlower[j] = std::max( pm.DLL[j], dcFloor );
            problem.xupper[j] = ( pm.DUL[j] < 1e6 )
                                 ? std::max( pm.DUL[j], dcFloor ) : std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 10.;
        }
        // NOTE on an experiment that was tried here and reverted: raising
        // the solvent's (pm.LO) own lower box bound above pm.XwMinM, to
        // stop Optima's undamped Newton step from ever driving it into
        // PrimalChemicalPotentials()'s "this aqueous phase is absent"
        // regime (see the still-open item on this in GEMS3K's CLAUDE.md,
        // 2026-08-23, "water solvent collapsing to its floor..."). Tried
        // pm.XwMinM*1e6 on Resources/gems3k/j_Flowline_G_series1_... -
        // made things WORSE, not better: the solve now hits its iteration
        // cap and diverges outright (VXc going negative, pH~5.7e9) rather
        // than landing on the previous wrong-but-plausible-looking answer,
        // and Albite - which the Hessian fixes below get right on their
        // own - collapsed back to its own floor too. Reverted rather than
        // kept and tuned further, consistent with this file's own
        // established caution about single-species bound/step-size
        // "magic number" patches (see the old MBR-based Tier A session
        // history above) - a real fix needs a properly scoped globalization
        // strategy, not a guessed threshold on one species.
        // Titrant unknowns are free (can be positive or negative), unlike
        // ordinary species amounts - bounded only generously, as a
        // numerical safety net, not a physical constraint.
        const double titrantBound = std::max( optima_default_box_moles( pm, pa_p->DG ), 1.0 ) * 2.0;
        for( long int k = 0; k < R; k++ )
        {
            problem.xlower[L+k] = -titrantBound;
            problem.xupper[L+k] =  titrantBound;
        }

        // Objective: Gibbs energy, in the same RT-normalized units GEMS3K's
        // own G0[]/G[] already use. fx = chemical potentials, recomputed
        // via GEMS3K's real activity-coefficient machinery every Optima
        // Newton iteration (not frozen/linearized the way MBR's Schur-
        // complement weights are). Virtual control-condition slots get a
        // FIXED (not composition-recomputed) gradient equal to their
        // target's own closed-form chemical potential - the whole
        // mechanism (see ipm_optima.h).
        //
        // Pure single-species phases (pm.Ls..pm.L - GEMS3K's own convention
        // of listing multicomponent-phase species first, single-component
        // ones after, per pm.Ls's own doc comment "total number of DC in
        // multi-component phases") get a small logarithmic-barrier term
        // added on top, ported directly from Reaktoro's own
        // EquilibriumSetup.cpp (updateGibbsEnergy()/updateGradX()/
        // updateHessX(), read from source - not guessed): Reaktoro adds
        // exactly this term for what it calls "pure phase species...whose
        // chemical potentials do not depend on composition" - precisely
        // our single-species mineral phases. tau (pa_p->LogBarrierTau)
        // defaults to Reaktoro's own default (EquilibriumOptions::epsilon *
        // logarithm_barrier_factor = 1e-16 * 1.0). Motivated by, but not
        // yet independently re-validated against, a real failure on
        // Resources/gems3k/j_Flowline_G_series1_...: without this term,
        // Optima's cold start converges to a genuine but wrong KKT point
        // where Albite (among other minerals that should precipitate)
        // never leaves its near-zero starting value - the plain diagonal
        // Hessian's curvature at small X[j] discourages growth rather than
        // encouraging it, and this barrier's negative gradient contribution
        // (-tau/X[j]) is intended to counteract that, matching Reaktoro's
        // own mechanism for the identical class of species (see GEMS3K's
        // CLAUDE.md, 2026-08-23, "should users switch between native and
        // Optima solvers, not use them together") - re-test on that project
        // before trusting this comment's own "fixes it" framing.
        // NOTE (plan v5 A.0, measured 2026-08-25): Reaktoro couples tau to the
        // species lower bound - `tau = epsilon * logarithm_barrier_factor`
        // where `epsilon` IS that bound - so its barrier pushes off the floor
        // with an O(1) gradient. GEMS3K's floor is dcFloor (pa_p->DHB,
        // 1e-13..1e-15), not this 1e-16 default, so as shipped the push is
        // ~1e-3 against mu0/RT of O(1..1000) - i.e. effectively inert, which
        // is what "tau does nothing" measurements have been seeing.
        // Restoring the coupling (tau = dcFloor) WAS A/B'd over the full
        // suite and is NOT adopted: it converts j_Solvus ROP and j_TiQ_PRSV
        // ROP from FAIL to OK (the latter in 20 iterations) and cuts AOP
        // iteration counts 25-45% on the Flowline/GEOTHERM tier, but breaks
        // f_Flowline ROP and o_/t_Solvus AOP - net AOP 17->15 OK,
        // ROP 12->13 OK. That fails the "no currently-converging project
        // regresses" gate, so the default stays as shipped. No code change is
        // needed to explore this: pa_LogBarrierTau is already a per-project
        // ipm-dat field, so setting it to the project's own pa_DHB reproduces
        // the tested configuration exactly. See CLAUDE.md 2026-08-25.
        const double kLogBarrierTau = pa_p->LogBarrierTau;
        // No reaktoroMode capture: since the FD PartiallyExact Hessian was made
        // unconditional (2026-08-25) the objective, gradient and Hessian are
        // byte-identical for AOP/SOP and ROP. The two modes now differ ONLY in
        // Optima::Options and in their retry chains - see below.
        const double kPhaseHessianFloor = pa_p->PhaseHessianFloor;
        // Captured by value like the constants above - pa_p is not in scope inside the lambda.
        const bool kFDHessian = ( pa_p->OptimaFDHessian != 0 );
        const bool kFDDiagFloor = ( pa_p->OptimaFDDiagFloor != 0 );
        const bool hasAq = HasAqueousPhase();

        // Leal 2014 §2.3.3's
        // two-clause unstable-phase test - below threshold AND DECREASING since
        // the last iterate. The phase-extinction tier below uses an exactness
        // argument (interchangeable twins) precisely because both magnitude-only
        // detectors fail on f_CASHNK; a trend test is the third option neither
        // considered. Tracked here because only the objective sees every iterate.
        auto phLast = std::make_shared<std::vector<double>>( pm.FIs, -1. );
        auto phDec  = std::make_shared<std::vector<long int>>( pm.FIs, 0 );
        auto phMax  = std::make_shared<std::vector<double>>( pm.FIs, 0. );
        // Tiered Hessian (pa_OptimaFDHessianDelay, default 0 = off). Set while the
        // short cheap-Hessian ATTEMPT below is running; the objective then behaves
        // exactly as pa_OptimaFDHessian = 0. shared_ptr because the objective lambda
        // outlives this scope inside Optima.
        //
        // It is an attempt-and-restart, NOT an in-place switch, and that is a
        // measured requirement rather than a style choice. An in-place version -
        // count iterations in the convergence hook, flip the FD loop on at N and
        // carry straight on - was implemented first and FAILS on the 301-point
        // solvus sweep at T = 560 C: with FD from the start that point converges in
        // 556 iterations, with the cheap Hessian alone it runs to the 7001 cap, and
        // switching FD on mid-trajectory does not rescue it at ANY delay tried
        // (300, 1000, 1500, 3000 all hit the cap). The first cheap iterations move
        // the iterate somewhere the exact columns cannot recover from, so the
        // trajectory has to be discarded rather than corrected.
        auto fdSuppress = std::make_shared<bool>(false);
        // AOP/SOP only. ROP is a deliberately faithful port of Reaktoro's own
        // arrangement and leaves Optima::Options untouched, so its budget is the
        // library default 200 per attempt; any delay comparable to the ~1100
        // iterations the cheap Hessian needs on f_/j_GEOTHERM would simply prevent
        // ROP from ever reaching the exact columns. Measured: with the delay applied
        // to ROP as well, j_TiQ_PRSV ROP goes OK 20 it -> FAIL 402.
        // COLD START ONLY, for the same reason OptimaReducedPreSolve() is (section 65.4):
        // a warm call already carries a consistent (primal, dual) pair, which is what the
        // cheap-Hessian probe exists to manufacture, so on a warm call the probe is pure
        // waste. Measured before gating it: a corpus freeze with the delay on taxed the
        // warm modes systematically - about twenty SOP/SHP rows went from 1 iteration to
        // 2, taking SHP +42 % and SOP +2 % overall - while AOP, the cold mode it is
        // actually for, went -5 %. HOP's Optima leg is warm (pm.pNP = 1) and was +7 %.
        const long int kFDDelay = ( !reaktoroMode && pm.pNP == 0 && pa_p->OptimaFDHessianDelay > 0 )
                                  ? pa_p->OptimaFDHessianDelay : 0;
        const bool kMoleFracHessian = ( pa_p->OptimaMoleFracHessian != 0 );
        problem.f = [this, L, R, dcFloor, &fixedGrad, kLogBarrierTau, kPhaseHessianFloor, kFDHessian, kFDDiagFloor, kMoleFracHessian, fdSuppress, hasAq, phLast, phDec, phMax]
                    ( Optima::ObjectiveResultRef res, Optima::VectorView x,
                      Optima::VectorView /*p*/, Optima::VectorView /*c*/,
                      Optima::ObjectiveOptions opts )
        {
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
            if( opts.eval.fxx )
            {
                res.fxx.setZero();
                // 1/X[j] ("ideal mixing") curvature applies only to species in
                // TRUE multicomponent phases (j < pm.Ls), where mu_j actually
                // depends on the phase's own composition via ln(x_j) (see
                // DC_PrimalChemicalPotential()'s DC_SYMMETRIC/DC_ASYM_* cases
                // above, ipm_chemical.cpp). For j in [pm.Ls, L) - single-
                // component ("pure") phases - PrimalChemicalPotentials()
                // itself already returns F[j]=G (DC_SINGLE case): a pure
                // phase's mu does NOT depend on its own amount at all, so the
                // TRUE second derivative there is exactly zero, not 1/X[j].
                // Applying 1/X[j] here anyway (as an earlier version of this
                // code did) fabricates curvature the gradient formula above
                // doesn't have: at a tiny cold-start X[j] (e.g. 1e-5 mol) it
                // manufactures an enormous artificial resistance
                // (1/1e-5=1e5) to that species growing, even when its real
                // driving force (F[j]-dual) strongly favors precipitation -
                // exactly the mechanism that left Albite pinned at its
                // near-zero starting value on Resources/gems3k/
                // j_Flowline_G_series1_... (GEMS3K's CLAUDE.md, 2026-08-23).
                // Confirmed against Reaktoro's own EquilibriumSetup.cpp/
                // EquilibriumHessian.cpp (read from source): its equivalent
                // "approximate"/"diagonal" Hessian computes an ideal-mixing
                // block PER PHASE (lnMoleFractionsJacobian et al.) - for a
                // single-species phase this is exactly a 1x1 block of zero
                // (ln(x)=ln(1)=0 identically, so its derivative is zero too),
                // and the ONLY curvature Reaktoro's own default (Partially-
                // Exact) Hessian gives such a species is its own tiny
                // log-barrier term (tau/n^2) - matched here by keeping that
                // term (below) while dropping the erroneous 1/X[j] one.
                // Per-multicomponent-phase curvature. Default: same ideal-
                // mixing diagonal as before (1/X[j]). For the ONE phase that
                // contains the water-solvent species (pm.LO) - i.e. the
                // aqueous phase, identified structurally rather than by name
                // so this isn't tied to any particular project's phase
                // ordering - also fill the off-diagonal row/column coupling
                // every solute's chemical potential to the solvent's own
                // amount, matching Reaktoro's own EquilibriumHessian.cpp
                // "approximate()" aqueous block (lnMolalitiesJacobian, read
                // from source): ln(molality_i) = ln(X[i]) - ln(X[w]) + const
                // for solutes i!=w, so d F[i]/dX[w] = -1/X[w] (and
                // symmetrically d F[w]/dX[i] = -1/X[w], matching
                // DC_ASYM_CARRIER's own F[w] ~ ln(X[w]) - ln(phase total)
                // form in DC_PrimalChemicalPotential() above); the diagonal
                // for w itself is (Xf-X[w])/X[w]^2, not bare 1/X[w].
                // A pure diag(1/X[j]) Hessian - what this code used
                // unconditionally before this fix - has NO off-diagonal
                // coupling at all between a trace redox pair (e.g.
                // H2(aq)/O2(aq)) and the solvent, even though a mass-
                // balance-neutral "swap water for an equivalent H2+O2
                // mixture" direction exists (H2's H:O=2:0 and O2's H:O=0:2
                // sum to water's own 2:1 ratio) - confirmed by direct
                // tracing on Resources/gems3k/j_Flowline_G_series1_... that
                // this is exactly the direction Optima's cold-start
                // iteration runs away along once the pure-phase Hessian fix
                // above lets the rest of the system actually move (H2(aq)/
                // O2(aq) blowing up to tens of mol, in the exact 2:1 ratio,
                // while mass balance and the plain diagonal-Hessian KKT
                // check both still report a clean solve) - the same near-
                // rank-1-singular "water ties H and O at a fixed 2:1 ratio"
                // mechanism already diagnosed for native MBR's Schur-
                // complement matrix (GEMS3K's CLAUDE.md, "Non-diagonal
                // ill-conditioning root cause", 2026-08-21). This off-
                // diagonal term is exactly what a pure diagonal Hessian is
                // structurally blind to.
                // Ideal-mixing curvature per multicomponent phase.
                //
                // AQUEOUS phase (the one containing pm.LO, identified
                // structurally rather than by name): solute rows are the
                // molality Jacobian, d ln m_i/d n = delta_ij/n_j - delta_jw/n_w;
                // the SOLVENT row uses the ideal water-activity convention
                //     ln a_w = -(1 - x_w)/x_w = -(nSum - n_w)/n_w
                //  => d ln a_w/d n_i = -1/n_w (i != w), (nSum-n_w)/n_w^2 (i = w)
                // This is NOT d ln x_w/d n_i, and the difference is deliberate.
                // Verified 2026-08-27 against Reaktoro's own
                // EquilibriumHessian.cpp: its aqueous approxfuncs calls
                // lnMolalitiesJacobian() (which zeroes the solvent row) and then
                // OVERWRITES that row with exactly these two expressions, under
                // exactly this comment. An earlier reading this session looked
                // only at the lnMolalitiesJacobian helper, concluded the solvent
                // row was a factor 1/x_w too large, "fixed" it to the
                // mole-fraction form, and measured a clear regression (10 new
                // failures on the 301-point solvus sweep, iterations 49k -> 109k,
                // f_CASHNK unaffected). Reverted. DO NOT "fix" this again -
                // read EquilibriumHessian.cpp's aqueous branch itself, not just
                // the helper it calls.
                //
                // NON-AQUEOUS solution phases: pa_p->OptimaMoleFracHessian
                // selects the form.
                //   0 (default) - diag(1/X[j]) only, no off-diagonal.
                //   1           - the full ideal mole-fraction Jacobian
                //                 d ln x_j/d n_i = delta_ij/X[j] - 1/Xf,
                //                 which is Leal/Kulik/Smith/Saar (2017) Eq. 80
                //                 (Docs/literature/), is what Reaktoro's own
                //                 non-aqueous approxfuncs assembles
                //                 (lnMoleFractionsJacobian, NOT overwritten,
                //                 unlike the aqueous branch above), and is the
                //                 exact derivative of the F[] this code already
                //                 computes for DC_SYMMETRIC species
                //                 (F = G + ln n_j - ln nSum, ipm_chemical.cpp).
                //
                // The default is 0 - i.e. formally the wrong derivative -
                // because it measures better on this corpus, not because it is
                // right. Setting 1 subtracts a rank-1 (1/Xf)*ones*ones^T from
                // each block. Derived exactly, that term is IDENTICALLY ZERO
                // along the unmixing direction (v=(1,-1) gives 1/n1+1/n2 either
                // way) and is the whole curvature along the phase-SCALING
                // direction (v=(1..1)), where the true block is singular by
                // construction: scaling a phase at fixed composition changes no
                // ln x. So it only affects how freely a phase's TOTAL amount
                // can move.
                //
                // MEASURED 2026-08-27, and it is a trade, not a win:
                //   f_CASHNK AOP 10003 -> 643 it (15.6x, same pH/Eh/Vs/Ms) -
                //     the interchangeable-Berman-twin case, whose vestigial twin
                //     decays toward extinction one order of magnitude per ~2000
                //     iterations. Freeing the scaling direction is what that
                //     case needs.
                //   f_Solvus AOP 145 -> 101 it ; j_Solvus 108 -> 100
                //   BUT the 301-point solvus sweep goes 0 -> 13 convergence
                //   failures, worst limb error 2.7e-3 -> 2.7e-1, iterations
                //   49k -> 187k; and j_CASHNK 563 -> 1600, f_/j_CalcDolo ~+10%.
                // Hence per-project opt-in, default off, same pattern as
                // pa_OptimaFDHessian / pa_OptimaStallWindow.
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
                    // Consecutive-decrease counter per phase (pa_MbTrendPhaseDecay)
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
                        // The rank-1 -(1/Xf) coupling is applied ONLY over the
                        // end-members that are PRESENT, by exactly the criterion
                        // pa_PhaseHessianFloor's block below uses to decide which
                        // columns deserve an exact (FD) curvature. An absent
                        // end-member keeps the plain diag(1/X_j).
                        //
                        // MEASURED 2026-09-06, and this is a defect fix, not a
                        // tuning choice. Applying the coupling to an end-member
                        // sitting at the numerical floor asserts curvature the
                        // physical model does not have (its X_j is a clamp, not
                        // an amount), and it opens a self-reinforcing trap:
                        //   - the coupling lets a whole phase be driven to the
                        //     floor, so nPresent falls to 0 or 1;
                        //   - the eigenvalue regularisation below requires
                        //     nP > 1 and is therefore SKIPPED exactly there;
                        //   - what Optima then sees is the raw analytic block,
                        //     which is exactly singular by construction (its
                        //     null vector is the phase's own amounts: scaling a
                        //     phase at fixed composition changes no ln x), mixed
                        //     with whatever FD columns landed on it - i.e.
                        //     indefinite - so the phase cannot recover.
                        // On the j_Solvus 61-point sweep at pa_OptimaMoleFracHessian
                        // = 1, T = 540 C: 2332 objective evaluations in the
                        // unregularised nPresent<=1 regime against 43 with the
                        // field off (54x), and lmin < 0 in 1015 of them (44%).
                        // Where the regularisation does run (nPresent >= 2) lmin
                        // was never negative in either arm, 13847 evaluations.
                        // Restricting the coupling to the present set closes the
                        // trap by construction: with nPresent <= 1 the coupled
                        // sub-block is at most 1x1, so the phase block is
                        // diag(1/X_j) and positive definite.
                        double phTot = 0.;
                        for( long int j = j0; j < j1; j++ ) phTot += std::max( pm.X[j], 0. );
                        const double presThr = std::max( dcFloor * 1e3, phTot * 1e-6 );
                        for( long int j = j0; j < j1; j++ )
                        {
                            const double Xj = std::max( pm.X[j], dcFloor );
                            if( kMoleFracHessian && pm.X[j] > presThr )
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
                // Virtual slots get zero direct curvature - a fixed
                // gradient has no second derivative in its own value; all
                // coupling to the rest of the system comes through the
                // shared Aex mass-balance rows.

                // ROP only: PartiallyExact (Reaktoro's own DEFAULT
                // GibbsHessian mode - EquilibriumSetup.cpp's
                // updateGradX(), read from source, not guessed). Reaktoro
                // starts from the exact same per-phase "approximate"
                // block as above (its own EquilibriumHessian::
                // approximate(), plus the identical log-barrier addition
                // for pure phases), THEN overwrites the FULL COLUMN
                // (every row i, not just the diagonal) for each Optima-
                // reported BASIC variable j (opts.ibasicvars - Optima
                // itself tells us which, so this doesn't need to be
                // re-derived from box bounds) with the true derivative
                // d(chem.potential_i)/dX_j:
                //   for(auto i : ibasicvars) { updateFx(i); Hxx.col(i) = grad(F); }
                // GEMS3K has no autodiff, so this is a forward-finite-
                // difference port of that exact mechanism - an
                // approximation of Reaktoro's own exact analytic
                // derivative, not a claim of bit-identical numerics, but
                // structurally the same "exact column for basic
                // variables, approximate elsewhere" strategy, which is
                // ROP's whole point (matching Reaktoro's MECHANISM, not
                // just tuning AOP's existing approximate Hessian).
                // TEMPORARY diagnostic ablation, added 2026-08-24 - remove
                // before this is considered done. Isolates the FD
                // PartiallyExact overwrite's effect on convergence rate,
                // holding ROP's own uniform seed / untouched Options /
                // Reaktoro-style retry all fixed - i.e. this env var skips
                // ONLY the FD column overwrite below, falling back to the
                // same base analytic (ideal-mixing + log-barrier) Hessian
                // AOP already uses.
                // Unconditional since 2026-08-25 (was ROP-only). This is not a
                // convergence-rate optimisation, it is a correctness requirement
                // for any phase with a miscibility gap: the base block above is
                // ideal-mixing only, hence positive-definite in composition
                // everywhere, i.e. it describes a SINGLE well. The real objective
                // is a double well - the gradient pm.F does carry the true
                // activity coefficients, and it is precisely the excess term's
                // negative (spinodal) curvature that opens the gap. Handing Newton
                // a convex Hessian for a non-convex objective leaves no barrier
                // between the two limbs, so an unmixed state slides back to the
                // homogeneous one and the gap is simply absent from the solver's
                // model of the problem. Measured on the Solvus family, where two
                // feldspar phases are deliberately defined with IDENTICAL
                // end-members and distinguished only by their major/junior
                // ('M'/'J') marks: without this block AOP returned one homogeneous
                // composition duplicated across both phase slots (66.3% albite in
                // each) where the true answer is two limbs at 79.1% and 49.9%, and
                // reported it as OK because the energetic penalty is only 2e-4 RT
                // near the critical point. See GEMS3K/CLAUDE.md, 2026-08-25.
                {
                    // pa_OptimaFDHessian = 0 skips this loop entirely - see that
                    // field's declaration comment (ms_multi.h) for what it costs
                    // and the measurement suggesting it may be redundant.
                    // pa_OptimaFDHessianDelay suppresses this loop for the
                    // duration of the short cheap-Hessian attempt made just
                    // before the primary solve - see fdSuppress's declaration.
                    const bool skipFD = !kFDHessian || *fdSuppress;
                    std::vector<double> Fbase( L );
                    for( long int j = 0; j < L; j++ ) Fbase[j] = pm.F[j];
                    for( Optima::Index bk = 0; !skipFD && bk < opts.ibasicvars.size(); bk++ )
                    {
                        const long int i = (long int)opts.ibasicvars[bk];
                        if( i >= L ) continue; // a control-condition virtual slot, not a real species - matches Reaktoro's own "i>=Nn: continue" guard
                        const double Xi = pm.X[i];
                        const double h = std::max( std::fabs(Xi) * 1e-7, dcFloor * 10. );
                        const double dIdeal = res.fxx(i,i);
                        pm.X[i] = Xi + h;
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                        for( long int j = 0; j < L; j++ )
                            res.fxx(j,i) = ( pm.F[j] - Fbase[j] ) / h;
                        // pa_OptimaFDDiagFloor - same rule as the reduced path (plan v5 s139.5).
                        if( kFDDiagFloor ) { g_fdDiagFloorCols++;
                            if( !( res.fxx(i,i) > 0. ) && dIdeal > 0. ) { res.fxx(i,i) = dIdeal; g_fdDiagFloorHits++; } }
                        pm.X[i] = Xi;
                    }
                    // EXACT, REGULARISED curvature for the non-aqueous
                    // multicomponent (solution) phases. Optima's basic set is
                    // exactly pm.N variables, so on a system with two ternary
                    // feldspars at most ONE of the two phases ever gets exact
                    // columns; the other keeps the ideal-mixing 1/X_j, which is
                    // positive and O(1/20) where the true curvature along the
                    // unmixing direction goes to zero at criticality. That
                    // mismatch is what makes Newton contract at 1 - H/B -> 1
                    // near Tc (measured: 0.9424 at 600 C, 0.9946 at 640,
                    // 1.0000 at 650, with FULL Newton steps throughout).
                    //
                    // Taking every column of such a phase exactly fixes the step
                    // magnitude (651 C: 288 iterations against >7001) but the
                    // exact block is indefinite inside the gap and an unmodified
                    // Newton step on it is unbounded - it proposes x = -2e11 and
                    // annihilates the phase. So symmetrise the block and floor
                    // its eigenvalues at a fraction of its own largest, which is
                    // the standard modified-Newton remedy and is well scaled in
                    // the phase's own curvature units.
                    //
                    // The AQUEOUS phase is deliberately excluded: its analytic
                    // ideal-molality block is already the right shape, and FD
                    // columns for its many near-floor trace species are noise
                    // (measured: extending this to the aqueous phase fails 10 of
                    // 11 near-critical points and returns a collapsed answer on
                    // the one that "converges").
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
                                // Only end-members that are actually PRESENT.
                                // An end-member at (or near) the numerical floor
                                // has 1/X_j up to 1e13, which would set the
                                // block's largest eigenvalue and so set the
                                // floor - flattening the physically meaningful
                                // O(1/X_phase) curvatures by ten orders of
                                // magnitude. Measured directly: with Anorthite
                                // at 1.4e-13 the feldspar block's max |entry|
                                // was 2.4e12 and the floor destroyed the model.
                                // Its FD column is meaningless anyway - the step
                                // h = max(|X|*1e-7, dcFloor*10) is larger than
                                // the species itself there, so the "derivative"
                                // is a 10x-perturbation secant.
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
                                    const bool okEig = SymEigFloorInPlace( blk, nP, regRatio );
                                    if( okEig )
                                        for( int a = 0; a < nP; a++ )
                                            for( int b = 0; b < nP; b++ )
                                                res.fxx(pres[(size_t)a],pres[(size_t)b]) = blk[(size_t)a*nP+b];
                                }
                            }
                            p0 = p1;
                        }
                    }

                    // Restore the base state the gradient/fval above were
                    // computed from - the FD loop above always restores
                    // pm.X[i] itself right after each perturbation, but
                    // XF/XFA/activity coefficients/pm.F were left at the
                    // LAST perturbation's values and must be refreshed.
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    for( long int j = 0; j < L; j++ )
                        pm.F[j] = Fbase[j];
                }
            }
            res.succeeded = true;
        };

        Optima::State state( dims );
        for( long int j = 0; j < L; j++ )
            state.x[j] = std::max( pm.Y[j], problem.xlower[j] );
        for( long int k = 0; k < R; k++ )
            state.x[L+k] = 0.0;

        // WARM START MUST SEED THE DUAL TOO, not just the primal. Optima is a
        // primal-DUAL interior-point solver: at (x, y) the reduced gradient is
        // s = g - Aex^T*y, so handing it a perfect x with y = 0 (which is what
        // Optima::State's constructor leaves behind) presents an iterate whose
        // stationarity is violated by the full size of g. The solver then has
        // to rebuild the entire dual from scratch, and the primal - already
        // correct - simply waits.
        //
        // Measured on Resources/gems3k/f_TestPNTDB (2026-08-27): re-solving via
        // SOP from native's own exact converged answer took **493 iterations /
        // 70 s**, against native's 224 / 70 ms. The per-iteration trace
        // (GEMS3K_OPTIMA_TRACE_FILE) shows exactly why, and it is not slow
        // convergence: ||ex||max sits at **3576, flat to four significant
        // figures, for 492 consecutive iterations**, then collapses to 2.9e-9
        // in the final one. x[compreignacite] never moves off its 1e-13 floor
        // the entire time and ||ep||max is 0 throughout - the primal was right
        // from iteration 0. What moves is y: every dual starts at EXACTLY 0,
        // is still ~5e-6 at iteration 250, and only reaches its true O(250)
        // magnitude at the end. The whole 492-iteration plateau is the dual
        // climbing geometrically from zero, and the 3576 is just
        // s = g - Aex^T*0 = g for the one species furthest from stationarity.
        //
        // GEMS3K already has that dual in pm.U[] - native's IPM writes it, this
        // function writes it (below), and unpackDataBr() does not touch it
        // (grepped: node.cpp never assigns pm.U[]). So it costs nothing to hand
        // it over. Sign convention is the inverse of the unpack further down
        // ("Optima's Lagrange multipliers ye use the opposite sign convention
        // from GEMS3K's own dual U[]"), hence the negation here.
        //
        // WARM PATH ONLY (pm.pNP != 0). On a cold AOP call pm.U[] is whatever
        // some earlier, unrelated calculation left behind - seeding from it
        // would be worse than zero, which is a defensible neutral start.
        // `warmStateUsable` (above) is false when SOP was requested on a node
        // that has never solved - the primal was then replaced by the cold LP
        // seed, so seeding a dual here would pair a fresh primal with a stale
        // (all-zero, hence meaningless) dual. Both halves of the start must
        // come from the same place; see that block's own comment.
        // dimReduceDone: the reduced pre-solve just wrote BOTH pm.Y[] and
        // pm.U[], so this is a warm start in every sense even on a cold call -
        // and seeding the primal without the dual is exactly the 493-iteration
        // failure the paragraph above documents.
        if( ( pm.pNP != 0 && warmStateUsable ) || dimReduceDone )
        {
            bool duals_usable = true;
            for( long int i = 0; i < N; i++ )
                if( !std::isfinite( pm.U[i] ) ) { duals_usable = false; break; }
            // ALL OR NOTHING - a PARTIAL dual is worse than no dual at all.
            // Measured 2026-08-27 on f_TestPNTDB by withholding the dual for the
            // trace ICs only (the exact thing a staged major/trace decomposition
            // could hand a second stage), SOP iterations:
            //     0 withheld   ->    3        1 of 42 withheld ->    3
            //    21 of 42      ->  837       22 of 42          ->  789
            //    28 of 42      ->  670       all (ye = 0)      ->  493
            // A dual that is right for the majors and zero for the traces is
            // INCONSISTENT - further from the central path than the uniformly
            // zero start, which is at least self-consistent - and costs more
            // than starting from scratch. Same lesson as the rejected
            // cold-start dual estimate below, one step sharper: it is not
            // "wrong is bad", it is "internally inconsistent is worst".
            // So never seed a subset of ye; seed all of it or none.
            if( duals_usable )
                for( long int i = 0; i < N; i++ )
                    state.ye[i] = -pm.U[i];
        }

        // A COLD-start dual ESTIMATE was tried here and REMOVED, 2026-08-27 - do
        // not re-add without reading this. The idea: with no previous dual to
        // reuse, fit one to the seed's own gradient by weighted least squares
        // (g_j = sum_i U_i*a(j,i), normal equations weighted by X[j]), the same
        // form native's MBR assembles. Implemented and A/B'd on
        // Resources/gems3k-psina/07PSIna_G_mid_1: **OK 4184 it / 42.0 s with a
        // zero dual, FAIL 60000 it / 254 s with the estimate.** A confidently
        // WRONG dual is far worse than a neutral zero one - the LP seed is a
        // sparse vertex with most species at the floor, so the fit is badly
        // determined, and a wrong dual misclassifies which species are
        // stable/unstable from iteration 0. Zero is not merely a lazy default
        // here; it is a defensible neutral start. (The warm path above is a
        // different case entirely: there the dual is not estimated but KNOWN.)
        Optima::Options options;
        if( !reaktoroMode )
        {
            // Reuse GEMS3K's own equivalents rather than adding parallel
            // Optima-only fields wherever one exists (per the AOP/SOP design
            // note in ms_multi.h): IIM is the native IPM loop's own iteration
            // cap, OptimaTol its convergence tolerance (a real trailing
            // BASE_PARAM field, default 1e-8 = Optima's own library default -
            // see its declaration for why DK is NOT reused here).
            //
            // maxiters is FLOORED AT 2000 rather than taken straight from IIM:
            // adopting the LP seed + unconditional FD Hessian (2026-08-25) left
            // f_GEOTHERM, whose project sets pa_IIM=1000, needing 1046
            // iterations and reporting FAIL at an otherwise correct state -
            // clipped 46 iterations short, confirmed by patching that project's
            // own pa_IIM. Raising a cap can only turn a budget-exhausted FAIL
            // into a converged result, never make a converging system diverge,
            // so it needs no opt-in flag.
            options.maxiters = (unsigned)std::max( 2000L, (long int)pa_p->IIM );
            // pa_OptimaEarlyStabilityAt's POSITIVE (cap) form is deliberately
            // NOT a maxiters clamp any more (2026-09-09, plan v5 section 102).
            // It was one until then, and that is precisely what made it
            // ungateable: a budget cannot be conditional. The solver runs out of
            // iterations inside its own stepping() loop, so the capped form never
            // entered the convergence hook and could not consult the dual-settled
            // test that the NEGATIVE (trend) form has had since section 87 - the
            // one thing that tells a trustworthy early stability verdict from a
            // premature one.
            //
            // The cap is now enforced in that same hook ("stop once iteration >=
            // N AND the dual has settled"), so both forms of this field are one
            // rule evaluated in one place. maxiters therefore stays at the full
            // budget here, and there is nothing to restore after the first
            // attempt - arming is what scopes the cap to it now.
            options.convergence.tolerance = pa_p->OptimaTol;
            // Trust region on per-variable Newton-step growth - a field this
            // branch added to its local Optima checkout. Default 0.0 = off;
            // measured HARMFUL at every nonzero value tried (GEMS3K/CLAUDE.md
            // 2026-08-24), kept only as re-runnable infrastructure.
            options.backtracksearch.max_step_ratio = pa_p->OptimaMaxStepRatio;
            apply_optima_linesearch( options, pa_p->OptimaLineSearch );
        }
        // else (reaktoroMode): leave Optima::Options() entirely at the
        // library's own untouched defaults - matching Reaktoro's own
        // practice exactly (Reaktoro/Equilibrium/EquilibriumOptions.hpp
        // carries a plain `Optima::Options optima` member with no
        // tolerance/iteration/step-ratio override of its own). This is
        // deliberate, not an oversight: applying GEMS3K's own DK/IIM-
        // derived overrides here would defeat the purpose of ROP, which
        // is to reproduce Reaktoro's actual numerics, not AOP's tuned
        // ones.

        // Optional per-iteration Optima trace, off by default (zero cost
        // when the env var is unset - Optima::Outputter itself is a stock
        // feature, on both this checkout and Reaktoro's own, so no
        // rebuild is needed on either side to use it). Set
        // GEMS3K_OPTIMA_TRACE_FILE to dump the same table format Reaktoro
        // exposes via its own Options.optima.output - this is exactly
        // what let a 2026-08-24 session diff GEMS3K+Optima's trajectory
        // against Reaktoro's own on the same chemical system and find the
        // solvent-seed root cause fixed above (see debug-optima-vs-
        // reaktoro/reaktoro_flowline.py for the Reaktoro-side equivalent).
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

        // ---- Stall/freeze limit (pa_OptimaStallWindow; 0 = off, the default) ----
        // TRIAGE, not a fix. Every remaining AOP failure in the benchmark
        // corpus is a FROZEN iterate - f/j_TestPNTDB, f/j_TestSUP98 and the two
        // largest gems3k-psina projects all have f and Error bit-identical from
        // iteration 0 - so each burns its whole budget and returns nothing (40
        // minutes at 1392 species, against a native solve of 1.3 s). This makes
        // that an early, honest failure. See Docs/gems3k-optima-plan-v5.md 15-17.
        //
        // Implemented through Optima::ConvergenceOptions::check, a stock
        // caller-supplied hook (Optima/ConvergenceOptions.hpp) that GEMS3K did
        // not previously use - so this needs NO change to the Optima checkout.
        //
        // CAREFUL, two non-obvious points:
        //  1. `check` returning true means CONVERGED (Convergence.cpp:56 ORs it
        //     into the tolerance test), and every retry tier below is gated on
        //     !result.succeeded. So a naive stall exit would SUPPRESS the
        //     phase-extinction retry and turn f_/j_CASHNK from OK into FAIL.
        //     The flag is therefore folded back by forcing succeeded=false, so
        //     all existing gating and the final verdict work unchanged.
        //  2. The test is on the BEST-SO-FAR error, not on the objective and not
        //     on iterate displacement. Both of those were measured and do not
        //     discriminate: f is frozen to 6 significant figures on f_CASHNK
        //     (which is progressing - its vanishing phase decays four orders of
        //     magnitude) exactly as on complex_1 (which is not), and complex_1's
        //     iterate moves MORE per 50-iteration window than f_CASHNK's.
        //     Nor does ||ex||inf on its own: it INCREASES on f_CASHNK.
        //
        // The state must be reset before EVERY solve() call - the retries reuse
        // this same `options` object, so a carried-over best would make the
        // first retry look instantly stalled.
        struct StallWatch {
            long int window = 0;
            double bestErr = 0.;    // best-so-far ||ex||inf (Optima's own criterion)
            double bestComp = 0.;   // best-so-far max_j |ex_j| * x_j (complementarity)
            // Best-so-far as it stood one whole WINDOW ago. The test is
            // cumulative over a window rather than per-step - see the test body.
            double refErr = 0., refComp = 0.;
            long int run = 0;       // iterations elapsed in the current window
            bool stalled = false;
            // Wall-clock budget (pa_OptimaMaxSeconds). Deliberately NOT cleared by
            // reset(): the budget covers the whole call - primary solve plus every
            // retry - not each attempt separately, so the deadline is set once.
            double maxSeconds = 0.;
            std::chrono::steady_clock::time_point started;
            bool timedOut = false;
            void reset() {
                bestErr = bestComp = std::numeric_limits<double>::infinity();
                refErr  = refComp  = std::numeric_limits<double>::infinity();
                run = 0; stalled = false;
            }
            bool overBudget() const {
                if( maxSeconds <= 0. ) return false;
                return std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - started ).count() > maxSeconds;
            }
        };
        auto stallWatch = std::make_shared<StallWatch>();
        stallWatch->window = pa_p->OptimaStallWindow;
        stallWatch->maxSeconds = pa_p->OptimaMaxSeconds;
        stallWatch->started = std::chrono::steady_clock::now();
        stallWatch->reset();

        // pa_OptimaEarlyStabilityAt < 0 selects a TREND trigger instead of a
        // fixed iteration cap: -N means "end the first attempt as soon as some
        // multicomponent phase has fallen monotonically for N consecutive
        // objective evaluations", so the phase-selection repair loop gets its
        // look while that phase is still dissolving rather than after the
        // primary solve has spent its whole budget reaching an assemblage it
        // then has to correct.
        //
        // WHY A TREND AND NOT A COUNT. A fixed cap is spent on EVERY call,
        // including the great majority whose assemblage is fine and which simply
        // need their budget - measured at ~1.24x on an ordinary cold solve
        // (plan v5 section 58.2), which is exactly why that field stays default
        // 0 despite being worth 33.5x where it acts. A trend trigger is spent
        // only on a run that shows the signature, so an ordinary cold solve -
        // where phases GROW into place from the LP seed and no sustained fall
        // exists - never pays for it at all.
        //
        // Same sign-overload convention as pa_OptimaDimReduceTol, and the same
        // reason: it is a different RULE for one decision, not a second knob.
        //
        // The counters are the ones pa_MbTrendPhaseDecay already maintains in
        // the objective callback (the only place that sees every iterate). This
        // is that field's criterion reaching a gate a CONVERGING run can open -
        // its own consumer sits inside `if( !result.succeeded )`, a failure-only
        // tier, which is why plan v5 section 45.1 found it could never fire on
        // the decay case it was built for.
        //
        // A false trigger costs the abandoned attempt and nothing else: the
        // safety net below re-solves once at the full budget when the early
        // attempt found nothing to repair. That is what makes a liberal trigger
        // safe, and it is why this could not have been built before the net.
        // THREE-VALUED, and 0 is AUTO rather than off since 2026-09-10 - both
        // forms therefore read the RESOLVED value, never pa_p directly, so the
        // solver and native_trace_run_header()'s EFF line cannot disagree about
        // what ran (ms_multi.h, optima_earlystability_at). AUTO can only ever
        // resolve to a positive cap, so the trend form is unaffected by it; it is
        // read through the same call anyway, because a second reader of a
        // three-valued field is exactly how one of them drifts.
        const long int earlyStabilityAt =
            optima_earlystability_at( pa_p->OptimaEarlyStabilityAt,
                                      optima_multisite_phase_count( pm.sMod, pm.FIs ),
                                      pm.pNP != 0 );
        const long int earlyTrendN = earlyStabilityAt < 0 ? -earlyStabilityAt : 0;
        // Only a guard that the phase really is below its own running peak, NOT
        // a magnitude test - measured, and the measurement is the point. On the
        // decay case the dissolving phase is still at 14% of its peak 6638
        // consecutive falls in, so any 10x-or-more clause becomes true only near
        // the very end of the run and the "early" look is not early at all. That
        // is the opposite of the phase-extinction tier's own `decaying`, which
        // wants 100x precisely because it is deciding EXTINCTION rather than
        // asking for a look.
        //
        // The consecutive-decrease count is what discriminates, and it does so
        // cleanly: on the same run the aqueous phase's longest monotone fall is
        // 30 evaluations against the dissolving phase's 6638.
        const double kEarlyTrendDropRatio = 0.999;
        // THE DUAL-SETTLED GATE (plan v5 section 87). The trend trigger asks the
        // phase-selection repair loop to look early, and that loop's verdict is a
        // stability index computed from the DUAL - so it is only as trustworthy as
        // the dual is. Acting on a dual that is still moving is exactly how the
        // trend form went wrong: CSHSnplus OK 1298 -> BAD 197 with G worse by
        // 0.8 RT (section 59.3), and section 58's safety net cannot catch it,
        // because the net only re-solves when the look found NOTHING.
        //
        // Optima's convergence hook receives BOTH iterates of u = (x, p, w), and w
        // is the dual - so the dual's own relative movement is free here and needs
        // no new plumbing. Measured at the moment the trigger fires, on the two
        // projects that decide this field:
        //     j_CASHNK   acting is RIGHT   max|dw|/max|w| = 2.97e-15
        //     CSHSnplus  acting is WRONG   max|dw|/max|w| = 1.33e-04
        // Eleven orders of magnitude apart - "settled to machine precision" against
        // "still moving" - so the threshold is not a tuned knob. It is
        // pa_OptimaTol, which is not an arbitrary pick either: it says the dual has
        // stopped moving at the scale the solve is trying to converge to, and it
        // sits ~7 decades above the RIGHT case and ~4 below the WRONG one.
        const double kEarlyTrendDualSettled = pa_p->OptimaTol;
        // THE RATE CLAUSE (plan v5 section 92) - the discriminator the dual gate
        // cannot be. Section 87.4 left one false positive standing: on a WARM leg
        // the trigger fires on a phase that has fallen 0.09 % over 342
        // evaluations - settling, not dissolving - and the dual gate is a no-op
        // exactly there, because a warm start inherits an accurate dual. Section
        // 87.5 showed the obvious LEVEL clause ("the phase must be below some
        // fraction of its peak") is blocked: the usable window is 0.14 to 0.91, a
        // tuned constant with ~2x margin.
        //
        // A RATE is scale-free in the length of the fall where a level is not.
        // Measured with GEMS3K_EARLYTREND_PROBE on the three cases that decide
        // this field (all three reproduce their recorded figures exactly - the
        // dual movements below are section 87.2's 2.97e-15 and 1.33e-04 to three
        // digits, from a different instrument):
        //
        //   case                 falls  frac of peak   rate/eval   dual     acting
        //   j_CASHNK  AOP cold      50      1.5e-12      2.0e-02   settled  RIGHT
        //   j_CASHNK  HOP warm     342      0.99900      2.9e-06   settled  WRONG
        //   CSHSnplus AOP cold      50      0.92401      1.5e-03   moving   WRONG
        //
        // The rate only has to separate what the dual gate LETS THROUGH, and
        // there it separates right from wrong by 6800x. 1e-4 sits ~200x below the
        // right case and ~34x above the wrong one - three decades of margin,
        // against the level clause's two-fold.
        //
        // Two things to be honest about. (1) rate = (1 - frac)/falls collapses to
        // 1/falls once a phase has essentially vanished, so on the RIGHT case the
        // number is set by how long the fall has run, not by how deep it is -
        // which is the intent ("has this phase lost a fraction of itself
        // comparable to the time it has been falling"), but it means the clause is
        // not independent of the trigger's own N. (2) Section 87.5's fourth case -
        // the f_Solvus decay run, "still at 14 % of peak 6638 consecutive falls
        // in", which is what BLOCKED the level clause - does not reproduce at
        // anything like that LENGTH: surveyed across the 61-point solvus AOP
        // sweep, the longest monotone fall of a non-solvent phase is between 40
        // and 45 evaluations (2 phases fire at N = 40, none at N = 45), not 6638.
        //
        // But the SIGNATURE is there, and it tests this clause rather than being
        // absent: those two reach 40 falls at 0.77 % and 11.7 % of their own peak
        // - the second is essentially 87.5's 14 % case - and their rates are
        // 2.5e-02 and 2.2e-02. Across all 56 non-solvent firings at N = 20 the
        // rate runs 2.1e-03 to 5.0e-02, so the SMALLEST genuine dissolution on
        // this corpus sits 21x above kEarlyTrendMinRate and 730x above the warm
        // false positive's 2.9e-06. On everything now measurable the clause fires
        // where it should and not where it should not.
        //
        // Both clauses separate the cases that remain; the rate's margin is 730x
        // against the level's ~8x, and 87.5's specific objection (a level clause
        // only becomes true near the end of a 6638-evaluation fall) has no case
        // left to stand on at these lengths. Section 92.4, corrected in 94.2.
        const double kEarlyTrendMinRate = 1e-4;
        auto earlyTrend      = std::make_shared<bool>( false );
        auto earlyTrendArmed = std::make_shared<bool>( earlyTrendN > 0 );
        // THE POSITIVE (cap) FORM, evaluated in the same hook and gated on the
        // same dual-settled test (plan v5 section 102). "Stop the first attempt
        // once iteration >= N, PROVIDED the dual has settled" - not "give the
        // first attempt N iterations of budget", which is what it was until
        // 2026-09-09 and why it could not be gated at all.
        //
        // WHY THE GATE BELONGS ON THIS FORM TOO, and it is the same argument
        // section 87 makes for the trend form: the early look's verdict is a
        // stability index computed FROM the dual, so it is worth exactly what the
        // dual is worth. Measured 2026-09-09 on the one row the cap form loses,
        // gems3k-fail/CSHSnplus_G_CSH1_5_bufs AOP, with the cap armed so every
        // sample is necessarily from an evaluation <= 200:
        //
        //   consecutive falls surveyed   1        5        20
        //   smallest max|dw|/max|w|      1.45e-01 5.25e-02 1.18e-03
        //
        // The smallest dual movement anywhere inside the capped window is 1.18e-03
        // against a 1e-8 threshold - five orders above it - and at the earliest
        // sample the dual is still moving by 100 % of its own magnitude. So the
        // gate blocks here by a wide margin rather than a close one.
        //
        // What that costs the row when it is NOT gated, measured the same day and
        // correcting a claim made confidently before it: the cap's loss on
        // CSHSnplus is NOT a budget shortfall that the safety net below should
        // have absorbed. The net fires, the full-budget re-solve converges in 727
        // iterations, and the row still fails - on the ASSEMBLAGE. At iteration
        // 200 the early look sees CASH+Sn at a logSI gap of 0.044 and deactivates
        // it (DECIDE phasesel loop=0 deactivate=2 logsigap=4.378700e-02); the
        // converged state contradicts that by two orders of magnitude (9.01) and
        // the phase is reported "absent but stable - should be present". G lands
        // 0.79 J high. A wrong decision on noise, not a truncated solve.
        //
        // The threshold is pa_OptimaTol for the same reason as the trend form's -
        // it says the dual has stopped moving at the scale the solve is trying to
        // converge to - and the two forms share kEarlyTrendDualSettled rather than
        // each carrying their own constant, because they are one rule.
        const long int earlyCapN = earlyStabilityAt > 0 ? earlyStabilityAt : 0;
        auto earlyCap      = std::make_shared<bool>( false );
        // Armed IMMEDIATELY BEFORE the primary solve, not here, and that is not
        // cosmetic. `options` carries this hook, and every solve() before the
        // primary one copies `options`: the phase-compaction probe and the
        // cheap-Hessian attempt both do. A cap firing inside either would set the
        // flag, and the fold-back after the primary solve would then mark THAT
        // solve failed however well it went - the exact silent failure the trend
        // form's own guard around the cheap attempt was added for
        // (j_10TH_G_seawater, OK 229 -> ERR 1226). Arming late makes the cap
        // unreachable from any earlier solve by construction rather than by
        // remembering to disarm at each one, and it reproduces the old maxiters
        // clamp's scope exactly: both of those probes overrode maxiters, so the
        // clamp never applied to them either.
        auto earlyCapArmed = std::make_shared<bool>( false );
        // The dual movement at the moment either form fired, for the DECIDE
        // record after the primary solve - the trace is what makes this decision
        // recoverable from the data rather than asserted (GEMS3K/CLAUDE.md s4).
        auto earlyStopDwRel = std::make_shared<double>( -1. );
        // GEMS3K_EARLYTREND_PROBE=<path>: one line the FIRST time each phase
        // satisfies the trend trigger's own count-and-drop clause, WHETHER OR
        // NOT the field is armed and whether or not the dual gate would let it
        // fire. It exists because the trigger's log line prints only when the
        // trigger actually fires, so from outside there is no way to see the
        // cases where it would have fired and was gated - which is exactly the
        // population a discriminator has to be measured against (plan v5
        // section 87.5's rate clause is on record as measured on TWO points).
        //
        // Prints the two quantities that separate the four known cases: the
        // phase's fraction of its own peak, and the average fractional fall per
        // evaluation over the monotone run - plus the dual's own movement, so
        // the existing gate and any candidate successor can be scored from one
        // file. Zero-cost when unset. GEMS3K_EARLYTREND_PROBE_N sets the
        // consecutive-fall count to survey at when the field itself is off
        // (default 200); when it is armed, the field's own -N is used, because
        // the question is what the trigger would see.
        static FILE* earlyTrendDbg = []() -> FILE* {
            const char* fn = std::getenv( "GEMS3K_EARLYTREND_PROBE" );
            return fn ? fopen( fn, "a" ) : nullptr; }();
        static const long int kEarlyTrendProbeN = []() -> long int {
            const char* v = std::getenv( "GEMS3K_EARLYTREND_PROBE_N" );
            return v ? std::atol( v ) : 200; }();
        const long int earlyTrendProbeN = earlyTrendN > 0 ? earlyTrendN : kEarlyTrendProbeN;
        auto earlyTrendSeen = std::make_shared<std::vector<char>>( pm.FIs, (char)0 );
        // The SOLVENT phase is excluded, and not as a tuning choice: the
        // phase-selection repair loop refuses to act on it by construction
        // ("never remove or reseed the solvent phase here"), so a trigger on it
        // can only ever spend a probe for nothing. Measured before this
        // exclusion: f_GEOTHERM fired on the aqueous phase settling by 7% -
        // 332.9 -> 308.8 over 50 evaluations, ordinary convergence, not
        // dissolution - and LOST its answer (OK 1765 -> ERR 3039).
        long int aqPhaseIdxTrend = -1;
        {
            long int j0aq = 0;
            for( long int k = 0; k < pm.FIs; k++ )
            {
                if( hasAq && pm.LO >= j0aq && pm.LO < j0aq + pm.L1[k] ) { aqPhaseIdxTrend = k; break; }
                j0aq += pm.L1[k];
            }
        }
        if( stallWatch->window > 0 || stallWatch->maxSeconds > 0. || earlyTrendN > 0
            || earlyCapN > 0 )
            options.convergence.check =
                [stallWatch, earlyTrendN, kEarlyTrendDropRatio, kEarlyTrendDualSettled,
                 kEarlyTrendMinRate,
                 earlyTrend, earlyTrendArmed, earlyTrendProbeN, earlyTrendSeen,
                 earlyCapN, earlyCap, earlyCapArmed, earlyStopDwRel,
                 aqPhaseIdxTrend, phDec, phMax, phLast]
                ( Optima::ConvergenceCheckArgs const& args ) -> bool
                {
                    // The dual's own relative movement over the last step -
                    // max|dw|/max|w| on u = (x, p, w). Optima's convergence hook
                    // receives BOTH iterates, so this is free and needs no new
                    // plumbing. Defined once and used by all three consumers (the
                    // probe below, the trend trigger, the cap): the whole argument
                    // for gating the cap on section 87's measurements only holds
                    // if the cap gates on the SAME quantity those measurements
                    // were taken of, and one definition is how that stays true.
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
                    // The probe (above): record what the trigger WOULD see, once
                    // per phase, independently of arming and of the gate.
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
                    // Trend trigger first: it is cheap, independent of the stall
                    // signals, and armed only for the first attempt.
                    if( *earlyTrendArmed && !*earlyTrend )
                    {
                        // Is the dual settled enough for the repair loop's stability
                        // index to be worth acting on? If not, do not fire - the phase
                        // is still decaying, so the trigger will be re-tested next
                        // evaluation and fires as soon as the dual does settle.
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
                                ipm_logger->info( "CalculateEquilibriumStateOptima: phase {} has "
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
                    // THE CAP FORM (pa_OptimaEarlyStabilityAt > 0), same fold-back
                    // and same safety net. Two conditions, and the second is the
                    // whole point of moving this rule in here: the first attempt
                    // has reached iteration N, AND the dual has settled, so the
                    // stability index the repair loop is about to compute from it
                    // is worth acting on.
                    //
                    // args.result.iterations is the solver's own count and is
                    // current at this point - Optima checks convergence before its
                    // budget test (MasterSolver::stepping), so the hook sees every
                    // iteration.
                    //
                    // ONE SHOT, AT N - the gate is asked once and the cap is then
                    // disarmed whatever it answered. NOT "wait at N until the dual
                    // settles", which is what this was first built as on
                    // 2026-09-09 and which MEASURES WORSE THAN DOING NOTHING:
                    //
                    //   f_Solvus_G_Test1 AOP   baseline 823   deferred 1530
                    //   CSHSnplus        AOP   baseline 1298  deferred 2119
                    //
                    // The reason is structural rather than particular to those two,
                    // which is why it decides the shape here. On a CONVERGING run
                    // the dual settles when the run is nearly over - f_Solvus's
                    // settles at iteration 707 of 823, at 9.65e-09 against the 1e-08
                    // threshold - so a look deferred until then is not an early look
                    // at all, and the probe it spends is very nearly a whole solve.
                    // Deferring converts the cap's honest, bounded cost ("waste N
                    // iterations") into an unbounded one ("waste however long
                    // convergence takes"), and it does so on exactly the ordinary
                    // well-behaved runs that have nothing to repair.
                    //
                    // Asking once at N keeps the bound: the look happens where the
                    // dual is ALREADY settled at N - a warm leg, or a run whose dual
                    // is determined early, which is precisely j_CASHNK HOP's 2.97e-15
                    // and where the cap's wins live - and otherwise the call is left
                    // alone at its baseline cost and baseline answer. Both of the
                    // rows above return to their baselines exactly.
                    //
                    // The latch is a `>=` test that disarms rather than an `==` one
                    // so that a hook the solver happens not to call at exactly
                    // iteration N cannot silently turn the field off.
                    if( *earlyCapArmed
                        && (long int)args.result.iterations >= earlyCapN )
                    {
                        *earlyCapArmed = false;   // asked and answered, either way
                        const double dwRelCap = dualRelMove();
                        // Recorded whichever way the gate answers, so the DECIDE
                        // record after the primary solve can report a BLOCKED cap
                        // as well as a fired one. A gate that declines is the
                        // interesting half of this mechanism and is otherwise
                        // invisible from outside.
                        *earlyStopDwRel = dwRelCap;
                        if( dwRelCap <= kEarlyTrendDualSettled )
                        {
                            *earlyCap = true;
                            ipm_logger->info( "CalculateEquilibriumStateOptima: reached iteration {} "
                                              "with the dual settled - ending the first attempt so "
                                              "the phase-selection loop can look "
                                              "(pa_OptimaEarlyStabilityAt = {}; max|dw|/max|w| = "
                                              "{:.3e} <= {:.3e})",
                                              (long)args.result.iterations, earlyCapN,
                                              dwRelCap, kEarlyTrendDualSettled );
                            return true;   // folded back to a failure below
                        }
                        ipm_logger->info( "CalculateEquilibriumStateOptima: reached iteration {} but "
                                          "the dual is still moving (max|dw|/max|w| = {:.3e} > "
                                          "{:.3e}) - NOT ending the first attempt; a stability "
                                          "verdict taken here would be read off noise "
                                          "(pa_OptimaEarlyStabilityAt = {})",
                                          (long)args.result.iterations, dwRelCap,
                                          kEarlyTrendDualSettled, earlyCapN );
                    }
                    // Two signals, and BOTH must stagnate. Each alone gives a
                    // false stall on a project that genuinely converges:
                    //  - ||ex||inf alone flags f_CASHNK, whose best-so-far error
                    //    stops improving at iteration ~220 and then INCREASES
                    //    (worst 200-window ratio 1.121) while the solve is in
                    //    fact progressing - its redundant phase decays four
                    //    orders of magnitude. Acting on that measured FAIL,
                    //    because the phase-extinction retry keys on that phase
                    //    becoming 1000x smaller than its twin, which needs
                    //    ~4600 iterations to develop.
                    //  - complementarity alone flags f_GEOTHERM, where it rises
                    //    144x over one 200-iteration window yet the run
                    //    converges at 1765.
                    // On a genuinely frozen run (complex_1, f_TestPNTDB) both are
                    // bit-identical from iteration 0, so both stagnate at once.
                    // Wall-clock guard first - it must fire even while the
                    // iterate is still improving, which is exactly the case it
                    // exists for (f_/j_TestSUP98 progress slowly and forever, so
                    // the stall test below never triggers).
                    if( stallWatch->overBudget() )
                    {
                        stallWatch->timedOut = true;
                        return true;   // stop now; folded back to a failure below
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
                    // The test is CUMULATIVE over a whole window, not per-step:
                    // fire only when best-so-far has failed to fall by a
                    // meaningful RELATIVE amount against its value one whole
                    // window ago. Rewritten 2026-09-03; the previous version
                    // counted consecutive iterations in which best-so-far did not
                    // improve AT ALL, and that had a measured false positive that
                    // is the entire reason pa_OptimaStallWindow's compiled default
                    // was reverted from 500 to 0 (see that field in ms_multi.h).
                    //
                    // WHY THE PER-STEP VERSION WAS WRONG. On the 301-point solvus
                    // temperature sweep the point at 588 C converges in 1074
                    // iterations, but reaches its answer through a plateau in which
                    // best-so-far improves only microscopically: over the window
                    // [500,1000) it falls by 0.21 % (0.012713 -> 0.012686), which
                    // is real progress but never a large enough single step to
                    // reset a per-step counter. The per-step rule at a 500 window
                    // fires at iteration 701 and turns a converging point into a
                    // failed one - measured directly, and it reproduces the
                    // reversal exactly. The cumulative test sees the 0.21 % and
                    // survives, with five orders of margin over the threshold.
                    //
                    // WHY THERE IS NO RANGE CLAUSE HERE, unlike the pre-solve
                    // watch in OptimaReducedPreSolve() which uses the same window
                    // machinery with an extra "has the live error even MOVED"
                    // clause. That clause is load-bearing there and would be
                    // actively harmful here, and the two projects that decide it
                    // pull in opposite directions:
                    //  - f_TestSUP98's reduced pre-solve converges through a
                    //    661-iteration excursion in which best-so-far is EXACTLY
                    //    frozen while the live error swings by 2.4e+02. Only the
                    //    range clause saves it, so the pre-solve needs it.
                    //  - f_/j_CASHNK's full solve is the one place in the corpus
                    //    where this watch does real work: best-so-far is exactly
                    //    frozen (6.936 / 0.0007399) for 19 consecutive windows
                    //    while the live error runs a perfect limit cycle
                    //    (lo/hi bit-identical every window, range/lo = 108 / 5e+04).
                    //    A range clause would read that oscillation as movement and
                    //    never fire, costing those two projects their rescue -
                    //    they converge either way, but at 10001 iterations instead
                    //    of ~1000. Measured, not reasoned.
                    // So the discriminator that works HERE is best-so-far alone:
                    // 588 C drops 2.1e-3 per window and survives; CASHNK drops
                    // EXACTLY 0 and fires.
                    //
                    // COST, stated plainly: CASHNK now fires at 999 rather than
                    // 641 / 561, so its hand-off to the phase-extinction retry is
                    // ~360-440 iterations later. That is the price of removing a
                    // false positive on a converging project, and it is paid on
                    // the two projects that were already the corpus's cheapest
                    // beneficiaries of this watch.
                    //
                    // The threshold is pa_OptimaTol's own default and no new
                    // number enters the file; `best` is updated on every genuine
                    // fall, even one too small to count, so the reference the next
                    // window is measured against never drifts upward on noise.
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
        // Applied after each solve() - see point 1 above.
        // Logged, not silent: this watch changes a solve's OUTCOME (it forces a
        // failure so the retry tiers can start), and until 2026-09-03 the only
        // way to tell whether it had fired was to A/B the whole run against a
        // build with the field off. The pre-solve watch has always said so in
        // its own discard message; this one now does too. info level, once per
        // solve that actually stalls - reset() clears the flag before each.
        auto applyStall = [stallWatch]( Optima::Result& r ) {
            // Deliberately does NOT fold the early-trend stop: that is done once,
            // at the primary solve's own call site, and folding it here would
            // mark every RETRY as failed too, since the flag stays set.
            if( stallWatch->stalled || stallWatch->timedOut )
            {
                ipm_logger->info( "CalculateEquilibriumStateOptima: full solve abandoned after {}"
                                  " - pa_OptimaStallWindow={}",
                                  stallWatch->timedOut ? "the wall-clock budget"
                                                       : "a window with no meaningful progress",
                                  stallWatch->window );
                // ... and into the TRACE as a DECIDE, not only into the log. This site abandons a
                // solve and hands it to the retry tiers, i.e. it changes both the outcome and the
                // cost - CLAUDE.md s4's rule for where a DECIDE belongs - and it was the one such
                // site with no record, so the freeze's `# dec` census could not see it at all and
                // no corpus-scale question about this guard was answerable without scraping logs.
                // Opened by work item 42b (plan v5 s136.8): the guard spends pa_OptimaStallWindow
                // iterations PER WINDOW deciding, and a stall costs at least TWO windows - the first
                // survives and rolls the reference forward, the second declares it. This record is
                // what established that: on LimBrine it reads window=500 run=500 iters=1000, so the
                // 1087-iteration cold solve is ONE abandoned solve spanning two windows plus an
                // 87-iteration re-solve, not two abandoned solves as the arithmetic alone suggested.
                // Whether that shape is corpus-wide or one fixture's quirk needs this record to ask.
                // `iters` is Optima's own count for the abandoned solve; `run` is how far into the
                // current window it got; `kind` separates the two ways this fires, because a
                // wall-clock timeout is not a stall and must never be pooled with one (and
                // pa_OptimaMaxSeconds is off by default, so `kind=timeout` should be absent).
                native_trace_decide( "optimastall kind=%s window=%ld run=%ld iters=%ld"
                                     " besterr=%.6e bestcomp=%.6e",
                                     stallWatch->timedOut ? "timeout" : "stall",
                                     (long)stallWatch->window, (long)stallWatch->run,
                                     (long)r.iterations, stallWatch->bestErr, stallWatch->bestComp );
                r.succeeded = false;
            }
        };

        Optima::Solver solver;
        solver.setOptions( options );
        // Needed only for ROP's own fallback below (Reaktoro's own retry
        // re-solves from the ORIGINAL state, not the diverged one) -
        // saved unconditionally since it's cheap (one State copy) and
        // keeps this line's placement obviously correct (before the first
        // solve mutates `state` in place) regardless of which mode runs.
        const Optima::State initialState = state;

        // Species this call's phase-extinction retry (below) deliberately
        // fixed at the floor, so the phase-assemblage stability check further
        // down can exempt them - see that exemption's own comment for why.
        // "Present" threshold: NOT bare pm.DSM (native's own PhaseSelect()
        // threshold, ipm_chemical.cpp) - confirmed a real false-positive
        // source by direct testing on tools/Cu-Pourbaix: native GEMS3K
        // represents an absent phase as exactly 0, but Optima's box-
        // constrained solve can NEVER reach exactly 0 (every species has a
        // dcFloor lower bound, for the log-based Hessian's well-posedness)
        // - a phase pinned exactly at dcFloor is Optima's own numeric
        // representation of "absent", the same role 0 plays natively, not
        // a real (if small) precipitated amount. Confirmed by direct trace
        // on every one of the Cu-Pourbaix control-condition sweep's 5
        // targets: C(cr)/Cu3(CO3)2(OH)2(c) reported "present but unstable"
        // with YF==dcFloor exactly (1e-13) while this project's own
        // pa_p->DS default is far smaller (1e-20) than GEMS3K's own
        // numerical floor - a mismatch between this ONE project's DSM
        // setting and dcFloor, not a real wrong assemblage. Reuses the
        // exact same threshold*10 margin already validated for this
        // purpose in DetectPhaseCollapseAndReseed() above.
        const double presenceThreshold = std::max( pm.DSM, dcFloor * 10. );

        // Phase-selection state, shared by the compaction probe just below and
        // by the phase-selection loop after the solve, because a phase pinned
        // by the probe must be READMISSIBLE by that loop - which needs both its
        // lifecycle flag and its original box.
        //   phSelState[k]: 0 = untouched, 1 = deactivated, 2 = readmitted
        //                  (terminal - a readmitted phase is never dropped
        //                  again, which is what bounds the loop).
        std::vector<char> phSelState( (size_t)std::max( pm.FI, 1L ), 0 );
        std::vector<double> savedLo( (size_t)std::max( L, 1L ), 0. );
        std::vector<double> savedHi( (size_t)std::max( L, 1L ), 0. );
        // Iterations spent in the compaction probe below - folded into
        // optimaIterTotal after the primary solve, since that is declared
        // there (from result.iterations) and this runs before it.
        long int probeIterations = 0;

        // ---- PSSC-equivalent phase compaction (pa_OptimaPhaseCompaction, default 0 = off) ----
        // A short, deliberately unconverged PROBE solve, purely to obtain a
        // dual good enough to say which phases are absent - then pin those at
        // the floor so the real solve below runs on the reduced active set.
        // See BASE_PARAM::OptimaPhaseCompaction (ms_multi.h) for why this
        // cannot be done after the fact and why misclassification is
        // self-correcting rather than silent.
        if( !reaktoroMode && pa_p->OptimaPhaseCompaction > 0 && L > 0 )
        {
            Optima::Options probeOpts = options;
            probeOpts.maxiters = pa_p->OptimaPhaseCompaction;
            Optima::Solver probeSolver;
            probeSolver.setOptions( probeOpts );
            Optima::State probeState = state;
            stallWatch->reset();
            Optima::Result probeResult = probeSolver.solve( problem, probeState );
            probeIterations += probeResult.iterations;

            // Which phase holds the solvent - never compacted away, whatever a
            // short probe happens to say about it (the solvent-collapse trap
            // this file documents at length is exactly a probe-length state
            // reporting water as absent).
            long int aqPhIdxC = -1;
            {
                long int j0a = 0;
                for( long int k = 0; k < pm.FIs; k++ )
                {
                    if( pm.LO >= j0a && pm.LO < j0a + pm.L1[k] ) { aqPhIdxC = k; break; }
                    j0a += pm.L1[k];
                }
            }
            long int nPinnedPh = 0, nPinnedSp = 0, j0c = 0;
            for( long int k = 0; k < pm.FI; k++ )
            {
                const long int j1c = j0c + pm.L1[k];
                if( k != aqPhIdxC && j1c <= L )
                {
                    // Same exemption as the stability check: a phase under any
                    // genuine kinetic restriction is the caller's business, not
                    // this pass's (mirrors native's KinConstrDC/KinConstrPh).
                    bool kinC = false;
                    for( long int j = j0c; j < j1c; j++ )
                        if( pm.DUL[j] < 1e6 || pm.DLL[j] > 0.0 ) { kinC = true; break; }
                    double totC = 0.;
                    for( long int j = j0c; j < j1c; j++ )
                        totC += std::max( probeState.x[j], 0. );
                    if( !kinC && totC < presenceThreshold )
                    {
                        for( long int j = j0c; j < j1c; j++ )
                        {
                            savedLo[(size_t)j] = problem.xlower[j];
                            savedHi[(size_t)j] = problem.xupper[j];
                            problem.xlower[j] = dcFloor;
                            problem.xupper[j] = dcFloor;
                            state.x[j] = dcFloor;
                            nPinnedSp++;
                        }
                        phSelState[(size_t)k] = 1;
                        nPinnedPh++;
                    }
                }
                j0c = j1c;
            }
            ipm_logger->warn( "CalculateEquilibriumStateOptima: phase compaction probe ({} iterations) "
                               "pinned {} of {} phases ({} of {} species) as absent",
                               probeResult.iterations, nPinnedPh, pm.FI, nPinnedSp, L );
            native_trace_decide( "compaction probeiters=%ld pinnedph=%ld ofph=%ld "
                                 "pinnedsp=%ld ofsp=%ld",
                                 (long)probeResult.iterations, (long)nPinnedPh,
                                 (long)pm.FI, (long)nPinnedSp, (long)L );
        }

        // ---- Tiered Hessian: a short attempt with the cheap Hessian first ----
        // pa_OptimaFDHessianDelay = N runs up to N iterations with the FD
        // PartiallyExact loop suppressed. If that converges, its state is adopted
        // and the primary solve below re-runs from it STILL CHEAP, confirming an
        // already-converged point in a handful of iterations (see the *fdSuppress
        // assignment at the end of this block for why it is not re-confirmed with
        // the exact Hessian). If it does NOT converge the state is left untouched
        // and the primary solve starts from the original seed with FD on - the N
        // iterations are spent and that trajectory is discarded, which is required
        // rather than tidy: see fdSuppress's declaration.
        //
        // Worth it because the FD loop costs a full activity evaluation per basic
        // variable per iteration and most of the corpus does not need it: measured
        // at N = 1500, 11x total wall time on the 17-project subset and 2.4x on the
        // 301-point solvus sweep, every G identical and the limb accuracy unchanged.
        // N must exceed the slowest cheap-Hessian convergence in the corpus
        // (f_/j_GEOTHERM, ~1100) or those projects restart and lose the win.
        long int fdCheapIterations = 0;
        bool cheapAttemptWon = false;
        // STRUCTURAL SKIP: some phase models are known to need the exact columns,
        // and for those the attempt is not merely wasted - see the f_CASHNK
        // measurement below. Decide it from the phase models rather than from the
        // iteration count, once, before spending anything.
        //
        // Two families, both established by section 21's independent FD-on/FD-off
        // A/B over all three corpora and reproduced exactly by this field's own
        // 25-project sweep:
        //   * multisite / reciprocal solid solutions (Berman, CEF, Modified
        //     Bragg-Williams). pa_PhaseHessianFloor's exact block only covers
        //     end-members above a presence threshold, but in a multisite model an
        //     end-member at low amount still carries real curvature through its
        //     SITE fractions - so that filter discards exactly the columns that
        //     matter and the FD loop is the only place the curvature exists.
        //   * fluid EoS phases (PR78, PRSV, CG, SRK, CORK, ...).
        // Note it is NOT "non-ideal": Van Laar is non-ideal and is one of the
        // biggest winners here (f_/j_Solvus 145/108 -> 68 iterations).
        //
        // Across gems3k + gems3k-fail + gems3k-psina exactly six projects match
        // (f_/j_CASHNK, j_TiQ_PRSV, CASH+CsSr, CASH+_G_csh_sol, CSHSnplus) and
        // they are precisely the FD-sensitive set section 21 identified - no
        // over-trigger, and no project outside it regressed under the delay.
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
            ipm_logger->warn( "CalculateEquilibriumStateOptima: pa_OptimaFDHessianDelay ignored - "
                              "this system has a multisite or fluid-EoS phase, which needs the "
                              "finite-difference Hessian columns from the first iteration" );
        if( kFDDelay > 0 && kFDHessian && !fdRequiredByModel )
        {
            Optima::Options cheapOpts = options;
            cheapOpts.maxiters = kFDDelay;
            Optima::Solver cheapSolver;
            cheapSolver.setOptions( cheapOpts );
            Optima::State cheapState = state;
            *fdSuppress = true;
            stallWatch->reset();
            // The early-trend trigger must NOT fire inside this attempt, and
            // getting that wrong is silent and expensive. cheapOpts is a COPY of
            // options, so it carries the same convergence.check hook; a trend
            // stop here would set the flag, and the fold at the primary solve
            // below would then mark THAT solve failed however well it went.
            // Measured before this guard: j_10TH_G_seawater went OK 229 -> ERR
            // 1226 for exactly this reason. The cheap attempt is itself already
            // a discarded probe - there is nothing for the repair loop to look
            // at until the real trajectory has run.
            const bool trendWasArmed = *earlyTrendArmed;
            *earlyTrendArmed = false;
            Optima::Result cheapResult = cheapSolver.solve( problem, cheapState );
            *earlyTrendArmed = trendWasArmed;
            *earlyTrend = false;
            *fdSuppress = false;
            fdCheapIterations = cheapResult.iterations;
            cheapAttemptWon = ( cheapResult.succeeded && !stallWatch->stalled && !stallWatch->timedOut );
            if( cheapAttemptWon )
                state = cheapState;
            else
                ipm_logger->warn( "CalculateEquilibriumStateOptima: cheap-Hessian attempt did not "
                                  "converge in {} iterations - discarding it and restarting with the "
                                  "finite-difference Hessian (pa_OptimaFDHessianDelay = {})",
                                  cheapResult.iterations, kFDDelay );
            // If the cheap attempt WON, the primary solve below re-runs from its
            // state with the FD loop STILL suppressed - it only has to confirm a
            // point that is already converged, in a handful of iterations.
            // Deliberately not "confirm it with the exact Hessian": that was tried
            // and it FAILS. On the 301-point solvus sweep, 142 of 301 points came
            // back BAD_GEM_AOP with a KKT stationarity residual of ~0.5 at H+, at
            // phase amounts identical to the passing run's to six figures - the FD
            // columns evaluated AT a converged point are near-singular, and the few
            // steps Optima then takes leave a state its own masked error test still
            // accepts while our unmasked KKT check correctly rejects it. Matches
            // Reaktoro's published rule exactly: the exact Hessian is for when the
            // quasi-Newton method FAILS, not for double-checking when it succeeds.
            *fdSuppress = cheapAttemptWon;
        }

        // ---- WORK ITEM 22: is there anything for the repair loop to find? ----
        // GEMS3K_PRESOLVE_PHSTAB_PROBE=<path>: run the phase-assemblage stability
        // check ONCE on the state this call is about to start from, before the
        // primary solve, and print what it would have said. Pure instrumentation -
        // nothing here feeds the solve, and it is not reached at all when the
        // variable is unset.
        //
        // The question it exists to answer (plan v5 section 104.6, work item 22):
        // pa_OptimaEarlyStabilityAt buys its win by spending N iterations to get a
        // dual good enough for the phase-selection repair loop to act on. Across
        // every row the corpus-wide flip moved, ONE predicate separates the wins
        // from the losses - did that loop then find anything? Where it did, the cap
        // paid for itself; where it did not, the cap cost exactly the probe. So the
        // switch worth having is a cheaper way to ask the same question, and the
        // obvious candidate is to ask it of the state we already hold.
        //
        // WHY THAT IS NOT CIRCULAR ON A WARM LEG, which is where the losses are.
        // A warm call starts from a previous converged answer, so pm.U[], pm.Y_la[]
        // and pm.Gamma[] already describe a real solved state and the scan costs one
        // activity-coefficient refresh rather than N iterations. Contrast the dual
        // gate (section 104.6(a)): on a warm leg dwrel is bit-exactly ZERO, because
        // nothing has moved yet - the least informative run produces the most
        // confident "settled". This scan has the opposite property, since a warm
        // start's dual is informative precisely because it came from a solve.
        //
        // WHAT TO READ OFF IT. `pre` is what the scan says at entry, `warm` says
        // whether the dual it read was a previous answer or a cold seed. The
        // predicate is confirmed if `pre` finds a violation exactly on the rows
        // where the post-solve `phasesel` DECIDE record deactivates something, and
        // finds none on 07PSIna_G_simple_0_0_0_150_0 SOP and j_CASHNK SOP - the two
        // rows where the cap costs the probe and buys nothing.
        //
        // NOT SUFFICIENT ON ITS OWN, and section 104.6 says why: section 101.6 has
        // this same loop deactivating CASH+Sn at a 0.044 gap and being WRONG where
        // 0.00796 here is RIGHT, so a scan that merely finds a candidate does not
        // establish that acting on it is safe. The gap column is printed so that
        // question can be asked of the whole corpus rather than of two points.
        {
            static FILE* preDbg = []() -> FILE* {
                const char* fn = std::getenv( "GEMS3K_PRESOLVE_PHSTAB_PROBE" );
                return fn ? fopen( fn, "a" ) : nullptr; }();
            if( preDbg && L > 0 )
            {
                // The scan reads derived arrays, not `state` - refresh them from
                // the seed this call will actually start from, or it scores
                // whatever the previous call happened to leave behind. Both are
                // recomputed after the solve on every path, so this cannot leak
                // into the reported answer.
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
                            preViol, preAbsent, &preRanked, &preCensus );
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
        stallWatch->reset();
        Optima::Sensitivity sensitivity( dims );
        // Arm the cap here and nowhere earlier - see earlyCapArmed's declaration
        // for why late arming is what scopes it to the first attempt.
        *earlyCapArmed = ( earlyCapN > 0 );
        Optima::Result result = wantSens ? solver.solve( problem, state, sensitivity )
                                         : solver.solve( problem, state );
        if( g_fdDiagFloorCols > 0 )
        {
            native_trace_decide( "fddiagfloor hits=%ld cols=%ld", g_fdDiagFloorHits, g_fdDiagFloorCols );
            g_fdDiagFloorHits = g_fdDiagFloorCols = 0;
        }
        // Returning true from convergence.check means "stop", and Optima reports
        // that as SUCCESS - so an early-trend stop has to be folded back to a
        // failure HERE, before anything reads result.succeeded. Getting this
        // wrong is silent: earlyCapHit below tests !result.succeeded, so folding
        // it later (in applyStall, where the stall watch does its own fold)
        // leaves the flag false and the safety net never fires.
        if( *earlyTrend || *earlyCap ) result.succeeded = false;
        // Did pa_OptimaEarlyStabilityAt - rather than the problem itself - stop
        // the first attempt? Consumed by the early-probe safety net further down.
        // Covers both rules: the positive form's iteration cap and the negative
        // form's phase-decay trend trigger, each of which now sets its own flag in
        // the convergence hook. The cap arm used to be inferred instead, from
        // "the solve failed AND it ran at least N iterations" - which was the only
        // thing available while the cap was a maxiters budget, and which could not
        // distinguish the cap stopping the attempt from the attempt exhausting a
        // budget that happened to be larger.
        const bool earlyCapHit = !result.succeeded && ( *earlyCap || *earlyTrend );
        // ---- ITEM 24: CAN THE PROBE'S COST BE RECOVERED, NOT ONLY SHORTENED? -
        // The early-probe safety net far below re-solves from `initialState`, so a
        // probe that finds nothing costs exactly N - measured at N = 200/100/50/25
        // with no exception (plan v5 section 105.3) - and the run total is
        // `N + full`. Those N iterations are not WRONG, they are THROWN AWAY:
        // section 103.1's finding that the re-solve reproduces the unarmed solve
        // exactly is the same statement said approvingly. Resuming the net from the
        // probe's own end state instead would make the total `~full`, i.e. the
        // probe free rather than merely short. Item 23 took the warm probe from 200
        // to 25; this asks for 0.
        //
        // CAPTURED HERE, from the primary solve, and deliberately NOT read off
        // `state` at the net: every retry tier between the two reassigns `state`,
        // and some of them pin species at the floor - which is the artefact
        // unpinForNet() exists to undo. The state captured here was produced on the
        // ORIGINAL box, so it is consistent with the problem the net restores, and
        // no species sits at a pin the net has just lifted.
        //
        // NOT THE DEFAULT WHILE IT IS UNMEASURED, and section 83 is why this is a
        // measurement and not an obvious fix: switching an Optima trajectory
        // mid-flight has failed before - "the first cheap iterations move the
        // iterate somewhere the exact columns cannot recover from, so the
        // trajectory has to be discarded rather than corrected". A resumed re-solve
        // may not simply subtract N; it may cost more than the restart does, or
        // converge somewhere the restart would not have. GEMS3K_OPTIMA_NET_RESUME=1
        // selects the resume; unset or 0 is the shipped restart. One getenv in a
        // static initialiser, and the State is copied only on a call where the
        // probe actually fired.
        //
        // MODE 2 IS THE RESUME WITH A FALLBACK, and it exists because mode 1 LOSES
        // ANSWERS. Measured 2026-09-10 on f_Solvus_G_test3 armed at -5, AOP, nine
        // nudges at 1e-15: the restart converges on 8 of 9 draws, the bare resume on
        // only 5 of 9 - the resumed iterate lands somewhere the solver cannot always
        // recover from, which is section 83's failure mode arriving in a milder form
        // than a corrupted trajectory. A lost answer outranks any iteration saving
        // (CLAUDE.md s2), so the bare resume is not defaultable at any measured win.
        // Mode 2 tries the resume and, if it does NOT converge, re-solves once more
        // from `initialState` exactly as mode 0 would have - so its worst case is
        // BOUNDED at one extra solve on a run that was already failing, and its best
        // case is every win mode 1 has. That is the same shape as the guard that made
        // pa_OptimaStallWindow defaultable.
        //
        // MODE 3 DEFERS THAT FALLBACK UNTIL THE WHOLE LADDER HAS FAILED, and it exists
        // because mode 2 pays for the guard on rows that never needed it. A re-solve
        // that STALLS is not the same event as a run that LOSES ITS ANSWER: on
        // f_CASHNK the resumed re-solve stalls and the phase-extinction tier converges
        // it anyway, so mode 2's extra solve buys nothing and costs 1000 iterations;
        // on f_Solvus_G_test3 the resumed re-solve also stalls and the tier does NOT
        // rescue it, and there the extra solve is the answer. Both are
        // `netresolve outcome=stalled`, so the outcome cannot tell them apart - but
        // the END OF THE LADDER can, because by then the tier has had its turn. So
        // mode 3 asks the question once, at the only point where it is answerable:
        // did this call, with every rescue it has, still fail? Same bounded worst case
        // as mode 2 - one extra solve on a run that has already failed - paid on
        // strictly fewer rows.
        // Resolved in optima_net_resume_mode() (ms_multi.h) so the trace's EFF
        // line and this call site cannot disagree - the failure mode s95.4 cost
        // a ladder measurement to, applied before it can happen again.
        const int netResumeMode = optima_net_resume_mode();
        const long int earlyProbeIters = earlyCapHit ? (long int)result.iterations : 0;
        std::unique_ptr<Optima::State> earlyProbeState;
        if( earlyCapHit && netResumeMode >= 1 && netResumeMode <= 3 )
            earlyProbeState.reset( new Optima::State( state ) );
        // DID THE PROBLEM FAIL, or did WE stop the attempt? Every failure-only
        // retry tier below is gated on `!result.succeeded`, and the fold-back above
        // makes that true even though the solve was progressing perfectly well -
        // so without this distinction the early stop hands a deliberately truncated
        // state to repair machinery whose entire premise is that the solve failed.
        //
        // MEASURED, and it is what made the 2026-09-09 default flip red where the
        // freeze could not see it (`solvus.aop` 7 of 61 points lost, plus
        // `solvus.stallwatch`'s control). At 588 C on the solvus sweep the cap fires
        // at iteration 200 with the dual settled to 1.476e-15 - a textbook correct
        // early look - and then:
        //
        //   reached iteration 200 with the dual settled ... (max|dw|/max|w| = 1.476e-15)
        //   phase-extinction retry - fixing 3 species of a vanishing interchangeable phase
        //   pNP=0 succeeded=true iterations=3            <- and the point is LOST
        //
        // The phase-extinction tier is built for a genuinely stalled twin-phase
        // system, where a redundant phase is heading for extinction and the solver
        // has given up. Handed a state stopped at 200, it reads a phase that is
        // merely still SHRINKING as vanishing, pins three species at the floor, and
        // "converges" in 3 iterations to an answer the stability check rejects.
        // Worse, it sets result.succeeded = true, so pa_OptimaEarlyStabilityAt's own
        // safety net - the one thing that could have restored the answer - is
        // skipped. A false failure propagating into a tier that then reports a false
        // success is precisely how a lost answer gets past every net.
        //
        // So when WE stopped the attempt, only the phase-selection repair loop may
        // act: it is the tier this field exists to reach, it is not gated on
        // `!result.succeeded` (it runs on every call), and its own verdict is
        // recomputed from the current state rather than assuming a failure. If it
        // finds nothing, the net re-solves at the full budget exactly as designed.
        //
        // A LAMBDA AND NOT A bool, and that distinction is not stylistic: `result` is
        // REASSIGNED by every retry tier as the ladder proceeds, so each
        // `!result.succeeded` test below is deliberately a question about the state
        // reached by the tiers ABOVE it, not about the primary solve. Snapshotting
        // this into a bool here was tried first and silently changed behaviour at the
        // SHIPPED default, where the whole mechanism is supposed to be unreachable:
        // j_CASHNK went AOP OK 1001 -> FAIL and HOP OK 4586 -> FAIL with
        // pa_OptimaEarlyStabilityAt = 0, because a tier that an earlier retry had
        // already fixed still saw "failed" and ran anyway. earlyCapHit IS a snapshot,
        // correctly - it is a statement about how the primary attempt ended and
        // nothing below changes that.
        auto genuineFailure = [&]() -> bool { return !result.succeeded && !earlyCapHit; };
        // Disarm both for every retry: the trigger exists to give the repair loop
        // an early look ONCE. Left armed, a phase that keeps dissolving - or a
        // dual that stays settled past N - would end every attempt in turn and the
        // run could never finish.
        *earlyTrendArmed = false;
        *earlyCapArmed   = false;
        // A DECIDE record for the trace, so which form fired, when, and on what
        // dual movement is recoverable from the data rather than from the log
        // (GEMS3K/CLAUDE.md s4: the trace is mode-attributed by position, the log
        // is not). Emitted for the cap and the trend alike, and only when one of
        // them actually acted - it is a decision that changes both the answer and
        // the cost, which is the bar for a DECIDE site.
        if( earlyCapHit )
            native_trace_decide( "earlystop form=%s fired=1 at=%ld dwrel=%.6e",
                                 ( *earlyCap ? "cap" : "trend" ),
                                 (long int)result.iterations, *earlyStopDwRel );
        else if( earlyCapN > 0 && *earlyStopDwRel >= 0. )
            // The cap reached N and the dual gate refused - the case that keeps
            // f_Solvus and CSHSnplus at their baselines. Traced because "the
            // mechanism did not fire" and "the mechanism was never reached" are
            // otherwise the same observable, and they are not the same thing.
            native_trace_decide( "earlystop form=cap fired=0 at=%ld dwrel=%.6e",
                                 earlyCapN, *earlyStopDwRel );
        // NOTE: no maxiters restore here any more. The cap stopped being a budget
        // clamp on 2026-09-09 (see its arming site above), so `options` has
        // carried the full budget all along and every retry below already has it.
        if( wantSens )
        {
            // Captured from the PRIMARY solve only. A retry re-solves a
            // different problem (reseeded, or with phases pinned), so its
            // derivatives would not describe the state finally returned unless
            // the retry itself succeeded - out of scope for this plumbing step.
            optima_dndb_rows = L + R;
            optima_dndb_cols = N;
            optima_dndb.assign( (size_t)(( L + R ) * N), 0. );
            for( long int j = 0; j < L + R; j++ )
                for( long int i = 0; i < N; i++ )
                    optima_dndb[ (size_t)( j * N + i ) ] = sensitivity.xc( j, i );
        }
        // Every retry below runs with the exact columns, whatever the first
        // attempt used: reaching a retry is itself the escalation trigger.
        *fdSuppress = false;
        applyStall( result );
        // Optima's own iteration count is otherwise invisible to any
        // caller (native's GEM_CalcTime()/GEM_Iterations() convention is
        // the only cross-solver-mode reporting path GEMS3K has) - summed
        // across every solve() call this function makes (primary plus any
        // retry), into pm.ITG below, since Optima's whole Newton loop is
        // architecturally one phase, unlike native's separate FIA(MBR)/IPM
        // split (pm.ITF is left at 0 for this reason, not left unset).
        long int optimaIterTotal = result.iterations + probeIterations + fdCheapIterations + dimReduceIters;

        // The phase-extinction tier, as a re-runnable unit rather than a block.
        // It is ASSEMBLED at its own site below (AOP/SOP only - ROP never
        // reaches it) and is left empty for ROP, so an empty check is a mode
        // test and not a defensive one.
        //
        // WHY IT HAS TO BE CALLABLE TWICE, measured on j_CASHNK HOP at a cap of
        // 200. An unarmed call rescues that project through
        //   primary solve stalls at 4500 -> phase-extinction tier -> OK at 4501.
        // pa_OptimaEarlyStabilityAt's own safety net re-solves at the full budget
        // when the probe finds nothing, and that re-solve reproduces the unarmed
        // primary solve EXACTLY - including its stall at 4500. What it did not
        // reproduce is the ladder that follows a stall, because the net is the
        // last rung: the tiers all sit above it and have already had their turn,
        // on the truncated state. So the row was lost at 4700 iterations (200 +
        // 4500) where the unarmed call returns OK at 4501 - not because the
        // re-solve differed from an unarmed solve, but because nothing was left
        // to rescue it. Handing the tier the full-budget state closes that, and
        // the state it then judges is the solver's own final word, so it is
        // judged on the ordinary (non-early-stop) path.
        std::function<void(bool)> runExtinctionTier;

        if( reaktoroMode )
        {
            // TWO SEPARATE REORDERING ATTEMPTS, both tried and REVERTED
            // the same session - do not try a third without a genuinely
            // new hypothesis, see below for why. The efficiency goal both
            // attempts shared: skip a wasted ~200-iteration toggle attempt
            // on systems where the trap, not a box-bound issue, is the
            // real problem, by trying the water retry first whenever the
            // trap signature is already present after the first solve.
            //
            // Attempt 1 reused the SAME `Optima::Solver` C++ object across
            // all retry attempts and broke `j_CASHNK` (silently accepted a
            // wrong `pH=0` answer as "OK"). Hypothesis: `Optima::Solver`
            // is not a pure, stateless function of `(problem,state)`, so
            // reuse across calls made call ORDER change the outcome.
            //
            // Attempt 2 gave every retry its OWN freshly constructed
            // `Optima::Solver` instance (confirmed first, via the full
            // 25-project sweep, that fresh instances alone are a no-op in
            // the already-validated toggle-first order) - and hit the
            // EXACT SAME `j_CASHNK` regression anyway. This DISPROVES the
            // attempt-1 hypothesis: solver-object statefulness was not
            // the (or not the whole) cause. The real explanation is more
            // fundamental and cannot be engineered around by an instance-
            // freshness trick: this is a non-convex NLP, and the two
            // retries' different trial states are different starting
            // points for Newton's own iteration - different starting
            // points can converge to DIFFERENT local KKT-stationary
            // points, one of which (on `j_CASHNK`) happens to be a real
            // but WRONG stationary point that this function's own KKT/
            // mass-balance/target checks don't catch (a legitimate local
            // optimum can still satisfy first-order conditions while being
            // physically wrong - the same reason native GEMS3K's own PSSC
            // phase-selection pass exists at all, per this file's own
            // "Key internals" notes). Reordering which retry runs FIRST
            // changes which local optimum gets found - toggle-first
            // happens to avoid the bad one on every system tested so far,
            // water-first does not. There is no reason to expect a THIRD
            // ordering choice (or any other instance-hygiene trick) to be
            // safe without the same full-sweep re-validation this file has
            // now applied twice - re-validate the FULL 25-project sweep,
            // not just the two systems this investigation started from,
            // before ever touching this order again.
            //
            // Reverted to the toggle-first order below (validated correct
            // across the full sweep, twice now) - fresh Optima::Solver
            // instances per attempt are KEPT (harmless, arguably better
            // hygiene, confirmed not to change any validated result).
            auto tryToggleRetry = [&]()
            {
                // Reaktoro's OWN fallback, ported verbatim from
                // EquilibriumSolver.cpp (read from source, not guessed):
                // toggle backtracksearch.apply_min_max_fix_and_accept and
                // re-solve from the ORIGINAL (pre-first-solve) state.
                Optima::Options toggled = options;
                toggled.backtracksearch.apply_min_max_fix_and_accept =
                    !toggled.backtracksearch.apply_min_max_fix_and_accept;
                Optima::Solver solverToggle;
                solverToggle.setOptions( toggled );
                state = initialState;
                stallWatch->reset();
                result = solverToggle.solve( problem, state );
                applyStall( result );
                optimaIterTotal += result.iterations;
            };
            // Detects the "aqueous solvent collapses to its floor" trap
            // signature on the CURRENT state - shares its formula with
            // tryWaterRetry() below (kept as a separate, duplicated inline
            // computation rather than routed through the shared
            // DetectSolventCollapseAndReseed() helper, per the same
            // subtle-behavior-change reasoning already on record for the
            // pre-reordering version of this code: that helper computes no
            // candidate value at all once the phase already dominates,
            // which would silently drop the validated "still try water
            // when !succeeded even if not currently trapped" case below).
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
            // GEMS3K-specific retry, NOT part of Reaktoro's own mechanism -
            // see this method's own declaration comment (ms_multi.h) and
            // GEMS3K's CLAUDE.md for the full history. Returns whether it
            // actually ran a solve (false if the current state doesn't
            // warrant one, e.g. water already exceeds the computed bound).
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
                        Optima::Result waterResult = solverWater.solve( problem, waterState );
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
            // Matches the exact condition validated across the full
            // 25-project sweep: trigger on EITHER `!result.succeeded` OR
            // the trap signature (not the trap signature alone).
            // The trap-signature arm keeps running after an early stop: unlike the
            // tiers gated on failure alone it is a POSITIVE observation about the
            // state in hand (the solvent collapsed), not an inference from "the
            // solve failed", so a truncated trajectory does not make it a false
            // positive - and a collapsed solvent at iteration N is real whatever
            // stopped the attempt.
            if( genuineFailure() || solventTrapped() )
                tryWaterRetry();

            // Generalization added 2026-08-24 (GEMS3K's CLAUDE.md, "Why a
            // single, general... initial approximation formula is not
            // simply available"): the "phase reported absent" cliff in
            // PrimalChemicalPotentials() (YF[k]<=pm.DSM) is generic to
            // EVERY multicomponent phase, not aqueous-specific - checked
            // here in ADDITION to (not replacing) the aqueous-specific
            // retry just above, which is left untouched to avoid any
            // behavior change to its own already-validated path.
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
                ipm_logger->warn( "CalculateEquilibriumStateOptima: general phase-collapse retry (ROP) - "
                                   "{} species reseeded", reseeds.size() );
                Optima::State retryState = state; // build on whatever the prior attempt(s) left, not the original seed
                for( const auto& rs : reseeds )
                    retryState.x[rs.first] = rs.second;
                Optima::Solver solverPhase; // fresh instance - see the note above tryToggleRetry's own solver
                solverPhase.setOptions( options );
                stallWatch->reset();
                Optima::Result retryResult = solverPhase.solve( problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                state = retryState;
                result = retryResult;
            }
        }
        else
        {
        // Solvent-seed retry, only if the plain solve above didn't already
        // succeed AND the aqueous solvent (pm.LO) ended up not dominating
        // its own phase by mass - the precise, confirmed signature of the
        // "aqueous solvent collapses to its floor" failure (GEMS3K's
        // CLAUDE.md, 2026-08-23/24). Root cause, found this session by
        // diffing Optima's own per-iteration output table (options.output,
        // a stock feature also used on the Reaktoro side, see
        // debug-optima-vs-reaktoro/reaktoro_flowline.py) against Reaktoro's
        // identical trace on the same chemical system
        // (Resources/gems3k/j_Flowline_G_series1_...): AutoInitialApproximation()'s
        // raw LP-simplex cold start ITSELF - not anything Optima's Newton
        // iteration does - seeded the solvent at its numerical floor
        // (~1e-5) while assigning essentially all of the H/O mass to two
        // trace redox species instead (H2(aq)~214mol, O2(aq)~102mol): a
        // legitimate LP-feasible vertex but a chemically absurd starting
        // point, since a real aqueous solvent should dominate its own
        // phase by mass in virtually every real system. Once seeded this
        // way, the ideal-mixing Hessian diagonal for the solvent,
        // (Xf-Xw)/Xw^2 (both the analytic block above and Reaktoro's own
        // identical formula, EquilibriumHessian.cpp, confirmed by reading
        // its source this session), diverges as Xw->0, fabricating an
        // enormous artificial resistance to the solvent ever growing back
        // out of a near-zero seed - the same class of trap as the
        // pure-phase 1/X[j] bug already fixed elsewhere in this objective,
        // just afflicting the solvent this time. Native GEMS3K's own MBR
        // pass (which gets the exact same AIA seed) never has this problem
        // because MBR's linear system is pure mass-balance/stoichiometry
        // (MassBalanceResiduals(), no PrimalChemicalPotentials()
        // dependency at all) - it can freely move mass between species
        // before IPM's gradient-driven main loop ever runs; Optima's
        // objective has no such gradient-independent recovery phase, so a
        // bad seed here is never corrected by anything else in the
        // pipeline on its own.
        //
        // Gated as a RETRY on the FINAL state's solvent dominance, not on
        // result.succeeded (that gate was tried first and does not work:
        // the trap this fix targets produces a clean, fully-converged
        // Optima result - KKT residual ~1e-9, mass balance satisfied - AT
        // the wrong answer, so result.succeeded is TRUE precisely in the
        // case that needs retrying, and checking it first meant the retry
        // never fired). Also not an unconditional up-front seed change: an
        // earlier version of this fix always reseeded the solvent before
        // ever calling Optima once, whenever it didn't already dominate
        // its own phase - this closed j_Flowline cleanly, but the
        // full-sweep regression check below caught a real cost on two
        // OTHER, unrelated projects (Resources/gems3k/
        // o_Solvus_G_series2_.../t_Solvus_G_series2_...): both already
        // recover from the SAME kind of bad AIA seed on their own, without
        // ever needing this correction, but the large, deliberately-
        // generous jump this fix uses (Xw -> total of every other species
        // in its phase, the magnitude actually needed to avoid the
        // curvature trap - a gentler geometric-mean-scaled jump was tried
        // and confirmed NOT large enough to unstick j_Flowline) introduces
        // a large primal-feasibility violation that cost those two systems
        // enough extra Newton iterations to miss their (unchanged)
        // iteration budget, turning a clean pre-existing success into a
        // real-but-marginal non-convergence (confirmed: Optima's own
        // reported residual was 1.05e-6, only ~2 orders of magnitude short
        // of the 1e-8 tolerance, not a runaway - see the dated CLAUDE.md
        // entry for the full trace). A retry-on-trap structure gives every
        // project that doesn't hit this specific failure mode (the
        // overwhelming majority, including both Solvus projects) the exact
        // same cost and behavior as before this fix existed, and pays the
        // extra-iterations cost of a second solve only on the systems that
        // actually end up trapped.
        {
            double waterSeed = 0.;
            // Guard against a no-op retry: if the reseed value would not
            // actually exceed what the ORIGINAL (pre-primary-solve) seed
            // already had for the solvent, re-running from it can only
            // reproduce the same trajectory at double the cost - confirmed
            // on j_TiQ_PRSV, where the seed's own up-front collapse
            // detection (above) already set pm.Y[pm.LO] to the same bound
            // this post-solve check recomputes, so the "retry" re-solved an
            // identical problem from an identical start and got the
            // identical (already-correct) answer for free. ROP's own
            // tryWaterRetry() already carries this exact guard
            // (initialState.x[pm.LO]) - ported here rather than duplicated
            // independently, see GEMS3K/CLAUDE.md, 2026-08-26.
            if( DetectSolventCollapseAndReseed( state.x.data(), problem.xupper[pm.LO], waterSeed )
                && waterSeed > initialState.x[pm.LO] )
            {
                ipm_logger->warn( "CalculateEquilibriumStateOptima: solvent-collapse retry - "
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
                Optima::Result retryResult = solverWater.solve( problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                state = retryState;
                result = retryResult;
            }

            // Generalization added 2026-08-24 (GEMS3K's CLAUDE.md, "Why a
            // single, general... initial approximation formula is not
            // simply available"): the "phase reported absent" cliff in
            // PrimalChemicalPotentials() (YF[k]<=pm.DSM) is generic to
            // EVERY multicomponent phase, not aqueous-specific - solid
            // solutions (and any other multi-species phase) can hit the
            // exact same trap if a sparse LP-simplex seed leaves their
            // TOTAL amount collapsed, independent of and in addition to
            // the aqueous-specific check just above (which is left
            // untouched, not replaced, to avoid any behavior change to
            // its own already-validated path). Checked on whatever state
            // the water-specific retry above left things at (its own
            // no-op case included), excluding the aqueous phase (already
            // handled its own way above).
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
                ipm_logger->warn( "CalculateEquilibriumStateOptima: general phase-collapse retry - "
                                   "{} species reseeded", reseeds.size() );
                // Base this retry on the CURRENT state (whatever the
                // water-specific retry above already left things at, if
                // it ran), not the original AIA seed - overwriting back
                // to the original seed here would silently discard that
                // retry's own fix if both a solvent collapse and another
                // phase's collapse were present simultaneously.
                Optima::State retryState = state;
                for( const auto& rs : reseeds )
                    retryState.x[rs.first] = rs.second;
                Optima::Solver solverPhase; // fresh instance - see the note above tryToggleRetry's own solver
                solverPhase.setOptions( options );
                stallWatch->reset();
                Optima::Result retryResult = solverPhase.solve( problem, retryState );
                applyStall( retryResult );
                optimaIterTotal += retryResult.iterations;
                state = retryState;
                result = retryResult;
            }

            // Fourth tier, added 2026-08-26: PHASE-EXTINCTION DEACTIVATION -
            // the nearest thing this solver has to native's PSSC pass, and the
            // fix for a failure mode that no amount of seeding, curvature or
            // step control can reach. Fires only when everything above has
            // still not converged, so a project that converges normally pays
            // nothing at all for it.
            //
            // The failure mode, traced end to end on f_CASHNK (GEMS3K's
            // CLAUDE.md, 2026-08-26): that project defines TWO solution phases
            // with identical end-members and identical G0, so exactly one of
            // them is redundant and the solver is free to put the mass in
            // either. Whichever slot loses decays geometrically toward
            // extinction - about one order of magnitude per 2000 iterations,
            // needing ~22000 to reach the floor from 1e-2 - because an
            // interior-point method can only approach a bound asymptotically,
            // never arrive. Native simply removes the phase.
            //
            // But the decay itself is NOT what blocks convergence, which is why
            // this needs deactivation rather than more iterations: once those
            // end-members are sitting at the floor, their reduced gradients
            // s[j] = g[j] + (Wx^T w)[j] are large (O(100-800), since the
            // chemical potential of a species in a near-empty phase carries a
            // log of a vanishing mole fraction) AND change sign between
            // iterations. Optima classifies a variable as "lower unstable" -
            // and therefore MASKS its optimality error to zero
            // (Optima/Stability.cpp, Optima/ResidualErrors.cpp) - only when
            // x[j]==xlower[j] AND s[j]>0. So each sign flip toggles ex[j]
            // between 0 and |s[j]|, and ||ex||max oscillates over three orders
            // of magnitude forever while the primal is EXACTLY frozen and mass
            // balance is perfect (||ew||max ~ 1e-14). Measured on f_CASHNK: the
            // composition is identical to every printed digit across 6000
            // iterations, the basic-variable set does not change even once, and
            // Error still swings 2.14 -> 828.6 against a 1e-8 tolerance.
            //
            // Fixing the box (xlower == xupper) is what breaks the toggle:
            // with a degenerate box, x[j] equals BOTH bounds, so one of
            // is_lower_unstable / is_upper_unstable holds for either sign of
            // s[j] and the variable is masked unconditionally. Measured on
            // f_CASHNK: 10000 iterations and no convergence becomes 57
            // iterations and succeeded=true, at a bit-identical G and pH. Note
            // a merely tightened one-sided bound does NOT work and was tried
            // first - it leaves xlower<xupper, so the toggle survives.
            //
            // Deliberately fixes only problem.xlower/xupper, NOT pm.DUL/pm.DLL.
            // That distinction is the safety net: the KKT check below exempts a
            // degenerate box (correctly - a fixed variable is an equality
            // constraint and carries no sign condition), while the
            // phase-stability check keys off pm.DUL/pm.DLL and so still runs in
            // full on the deactivated phase. If this tier ever removes a phase
            // that genuinely should be present, that check reports "absent but
            // stable - should be present" and the result is rejected rather
            // than silently returned.
            //
            // This tier RUNS after an early stop - unlike the toggle and water tiers,
            // which do not (see genuineFailure's declaration) - but its twin test
            // gains a clause when it does. See kTwinStillFalling below.
            runExtinctionTier = [&]( bool onEarlyStopPath )
            {
                // The aqueous phase index is recomputed here rather than
                // captured: the enclosing block's own aqueousPhaseIdx dies
                // before the safety net that calls this a second time, and a
                // [&] capture of it would dangle.
                long int aqueousPhaseIdx = -1;
                {
                    long int j0aq = 0;
                    for( long int k = 0; k < pm.FIs; k++ )
                    {
                        if( pm.LO >= j0aq && pm.LO < j0aq + pm.L1[k] ) { aqueousPhaseIdx = k; break; }
                        j0aq += pm.L1[k];
                    }
                }
                // Detection is deliberately STRUCTURAL, not a magnitude
                // threshold, because both obvious magnitude-based detectors
                // were tried on f_CASHNK and are wrong for this case:
                //   - "phase total below pm.DSM (or the presence threshold)"
                //     never fires. The redundant phase is at ~1e-5 when the
                //     solver gives up - it is HEADING for extinction, not
                //     there, and that is the whole problem.
                //   - native's own PhaseSelect() criterion (logSI < -DFM,
                //     "present but unstable, remove it") cannot fire either:
                //     two phases built from identical end-members sit on the
                //     SAME tangent plane, so the redundant one has logSI ~ 0,
                //     i.e. it is exactly as stable as its twin. Native does not
                //     remove it - it simply converges anyway, tolerating a
                //     vestigial 3.4e-6, because its Dikin criterion is not
                //     sensitive to the masking discontinuity described above.
                // What IS exact here: two solution phases with identical
                // stoichiometry AND identical standard potentials are
                // thermodynamically interchangeable, so a state with one of
                // them vanishing is the same state as one with it removed.
                // That identity is a property of the project, not a tuned
                // number. The one thing it must not do is break a genuine
                // two-limb split of such a twin pair - which is exactly the
                // Solvus projects' feldspars - so it additionally requires
                // one twin to be smaller than the other by kTwinRatio. A real
                // miscibility gap has limbs of comparable magnitude (the
                // Solvus pair runs ~0.09 vs ~0.16 even 16 C from its critical
                // point); a redundant twin is orders of magnitude down and
                // still falling. Solvus never reaches this code in any case,
                // since the tier is gated on the solve having failed.
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
                // Leal 2014 §2.3.3's two-clause unstable-phase test: below threshold
                // AND decreasing. Opt-in via pa_MbTrendPhaseDecay (default 0 = off);
                // see that field's comment in ms_multi.h for the measurement and for
                // why it is not on by default. Applies IN ADDITION to the twin
                // criterion below - it does not replace it, even though it was
                // measured to subsume it on j_CASHNK, because one case is not enough
                // to retire an exact criterion in favour of a tuned one.
                const long int trendDecayN = pa_p->MbTrendPhaseDecay;
                auto decaying = [&]( long int k ) -> bool
                {
                    if( trendDecayN <= 0 ) return false;
                    if( k < 0 || k >= (long int)phDec->size() ) return false;
                    const double tot = (k < (long int)phTot.size()) ? phTot[k] : 0.;
                    const double pk  = (*phMax)[k];
                    // Leal's two clauses: small, AND monotonically falling for a
                    // sustained run. The 100x drop from its own peak is what
                    // separates "heading for extinction" from "settled small".
                    return (*phDec)[k] >= trendDecayN && pk > 0. && tot < pk * 1.0e-2;
                };
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
                    if( decaying( ka ) )
                    {
                        for( long int j = ph0[(size_t)ka]; j < ph1[(size_t)ka]; j++ )
                            deactivated.push_back( j );
                        ipm_logger->warn( "CalculateEquilibriumStateOptima: phase {} "
                                          "deactivated by trend (consecutive decreases={}, "
                                          "total={:.3e}, peak={:.3e})",
                                          ka, (*phDec)[ka], phTot[(size_t)ka], (*phMax)[ka] );
                        continue;
                    }
                    for( long int kb = 0; kb < pm.FIs; kb++ )
                    {
                        if( kb == ka || kb == aqueousPhaseIdx || pm.L1[kb] <= 1 ) continue;
                        if( phTot[(size_t)kb] <= phTot[(size_t)ka] * kTwinRatio ) continue;
                        if( !interchangeable( ka, kb ) ) continue;
                        // THE "STILL FALLING" CLAUSE - the second half of the test
                        // this tier's own comment above states and the code never
                        // made ("a redundant twin is orders of magnitude down AND
                        // STILL FALLING"). The ratio clause is the first half; until
                        // 2026-09-09 it was the whole test, and that was sound only
                        // because of the invariant the same comment relies on -
                        // "Solvus never reaches this code in any case, since the tier
                        // is gated on the solve having failed".
                        //
                        // pa_OptimaEarlyStabilityAt BREAKS that invariant: its
                        // fold-back makes a perfectly healthy attempt look failed, so
                        // a converging Solvus point does reach here, and the ratio -
                        // a magnitude judgement on a deliberately truncated state -
                        // is then read at a moment when it means nothing.
                        //
                        // MEASURED on the three cases that exercise this, all at a
                        // cap of 200, printing the twin's own counters at the instant
                        // the tier chose to act:
                        //
                        //   case                 total     peak      ratio     falls
                        //   j_CASHNK AOP  RIGHT  4.27e-12  2.832     5.18e+11  70
                        //   j_CASHNK HOP  wrong  3.50e-04  3.50e-04  6.32e+03   0
                        //   Solvus 588 C  WRONG  2.74e-10  4.544e+01 1.59e+11   0
                        //
                        // The RATIO does not separate them - Solvus's 1.59e+11 is
                        // within a factor of 3 of the right case's 5.18e+11, both
                        // enormous. The consecutive-fall count separates them
                        // completely: 70 against 0 and 0. Solvus's twin HAS collapsed
                        // from a peak of 45.4, so "it has not grown yet" is NOT the
                        // story - it is REBOUNDING at iteration 200, the two limbs
                        // still sorting themselves out, which is exactly the "genuine
                        // two-limb split" the comment above says this tier must never
                        // break. j_CASHNK's redundant twin has fallen monotonically
                        // for 70 straight evaluations and is not coming back.
                        //
                        // No threshold is needed at a 70-versus-0 separation, so the
                        // clause is the minimal one - "it went down last evaluation" -
                        // and carries no tuned constant. It is applied ONLY on the
                        // early-stop path: on a genuine failure the state is the
                        // solver's own final word and the ratio means what it always
                        // meant, so that path is left exactly as it was.
                        if( onEarlyStopPath && (*phDec)[(size_t)ka] <= 0 )
                        {
                            ipm_logger->info( "CalculateEquilibriumStateOptima: phase {} looks like a "
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
                    ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-extinction retry - "
                                       "fixing {} species of a vanishing interchangeable phase at the floor",
                                       deactivated.size() );
                    // The kTwinRatio clause actually firing - plan v5 section 95.2
                    // measured its margin at 9.5x and nothing recorded when it acts.
                    native_trace_decide( "extinction fixedsp=%ld",
                                         (long)deactivated.size() );
                    Optima::State retryState = state;
                    for( long int j : deactivated )
                    {
                        problem.xlower[j] = dcFloor;
                        problem.xupper[j] = dcFloor;
                        retryState.x[j]   = dcFloor;
                        extinctFixed[(size_t)j] = 1;
                    }
                    Optima::Solver solverExt; // fresh instance - see the note above tryToggleRetry's own solver
                    solverExt.setOptions( options );
                    stallWatch->reset();
                    Optima::Result retryResult = solverExt.solve( problem, retryState );
                    applyStall( retryResult );
                    optimaIterTotal += retryResult.iterations;
                    state = retryState;
                    result = retryResult;
                }
            };
            if( genuineFailure() || earlyCapHit )
                runExtinctionTier( earlyCapHit );
            // Fifth tier, added 2026-08-27: PSSC-EQUIVALENT PHASE-SELECTION
            // LOOP - an ACTION attached to the phase-assemblage stability
            // detector, which until now could only ever REPORT a wrong
            // assemblage and reject the result. This is the "Phase E" item
            // that has been deferred since the prototype, and the scope note
            // at the top of this file ("no phase-stability (PSSC) pass") is
            // no longer true of the AOP/SOP path because of it.
            //
            // WHY THIS IS THE REMAINING GAP, measured rather than assumed
            // (plan-v5 section 20, GEMS3K's CLAUDE.md 2026-08-27):
            // relevance.cpp counts what native's own converged answer leaves
            // at the numerical floor, and the ABSOLUTE COUNT OF ABSENT PHASES
            // tracks this solver's difficulty far better than problem size -
            // f_Flowline 15 absent phases and 1.15x native, mid_1 117 and
            // 34x, f_TestSUP98 132 and never converges. The mechanism is not
            // a correlation: on f_TestPNTDB the frozen optimality residual
            // is held ENTIRELY by phases that are not there (compreignacite
            // at x=2.4e-13 with s=-3576, then becquerelite, Na-Weeksite,
            // ... all at or within a few 1e-6 of the floor). Native drops
            // them - that is exactly what PhaseSelectionSpeciationCleanup()
            // does on every native call with pa_PC=2 - while this solver has
            // carried every absent phase as a live box-constrained unknown
            // for the whole solve, contributing a row to the KKT system and
            // an ex_j that Optima's bound mask only zeroes when x==xlower
            // AND s>0, which is precisely not the case for these.
            //
            // STRUCTURE - iterative and post-solve, matching native's own
            // PSSC loop rather than a cheap up-front pre-pass. The stability
            // index is a function of the DUAL potentials, so it does not
            // exist before a solve has produced one; there is nothing a
            // pre-pass could read. Bounded at kMaxPhaseSelectLoops, the same
            // 5 native uses before it gives up with W08IPM.
            //
            // TERMINATION is structural, not a tuned iteration count. Each
            // phase carries a lifecycle: untouched -> deactivated ->
            // readmitted, and only those two transitions are ever allowed,
            // so a phase can change at most twice and the loop cannot
            // oscillate a phase in and out indefinitely. Readmission is what
            // makes deactivation safe to attempt at all: if dropping a phase
            // was wrong, the very next pass sees it as "absent but stable"
            // and puts it back, with its original box restored and a seed
            // taken from the bulk composition. A phase that has been
            // readmitted is never dropped again by this loop - if it still
            // disagrees after that, the final check reports the violation
            // and the result is rejected, exactly as before this tier
            // existed.
            //
            // Deactivation fixes ONLY problem.xlower/xupper at dcFloor, never
            // pm.DUL/pm.DLL - the same deliberate distinction the
            // phase-extinction tier above makes, and for the same reason: a
            // degenerate box is what makes Optima mask the variable's
            // optimality error unconditionally (Optima/Stability.cpp -
            // x==xlower AND x==xupper, so one of is_lower_unstable /
            // is_upper_unstable holds for EITHER sign of s), while leaving
            // pm.DUL/pm.DLL untouched means the phase is NOT exempted from
            // the final stability check. A phase this loop removes wrongly is
            // therefore still caught, both by the next pass here and by that
            // check. Note this differs from the phase-extinction tier, whose
            // deactivations ARE exempted - those are proved interchangeable
            // with a present twin and so report logSI ~ 0 tautologically,
            // which is a different situation from a phase removed on its own
            // stability index.
            {
                const int kMaxPhaseSelectLoops = 5;
                // Same threshold, same reason, as the final check's - see
                // the comment at its own call site. It must be identical or
                // this loop chases a violation that check would not report.
                for( int psLoop = 0; psLoop < kMaxPhaseSelectLoops; psLoop++ )
                {
                    // Evaluate the assemblage on the CURRENT state, through
                    // exactly the pipeline the final check uses - the
                    // stability index reads pm.Gamma[] (LINK_UX_MODE),
                    // pm.Y_la[] (CalculateConcentrations, from the committed
                    // dual) and pm.YF[]/pm.YFA[] (refreshed inside the
                    // helper). Anything cheaper here would classify against a
                    // different state than the one finally reported.
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
                                psViol, psAbsent, &psRanked, &psCensus );
                    if( psRanked.empty() )
                        break;                       // assemblage is self-consistent

                    // Walk the ranking rather than stopping at its head. An
                    // unactionable violation - the solvent phase, or one in a
                    // phSelState this tier may not touch - says nothing about
                    // the violations RANKED BEHIND IT, and stopping there made
                    // "nothing to repair" and "nothing repairable at the top of
                    // the list" indistinguishable from outside (plan v5 section
                    // 91, where the loop broke on the first of 289 entries it
                    // could not act on). Costs one scan of an already-computed
                    // vector; the state is not re-evaluated between entries,
                    // because acting on ANY of them invalidates it and the loop
                    // re-solves and re-evaluates from the top.
                    //
                    // What the ranking is worth is bounded by psCensus.clamped:
                    // a logSI saturated at StabilityIndexes()' overflow guard is
                    // not a driving force, and on an unconverged state most
                    // large systems saturate. Reported below, once, when it
                    // dominates - the loop still acts, because a clamped
                    // "absent but stable" is still a phase this tier removed
                    // and can put back.
                    if( psCensus.clamped * 2 > psCensus.scanned && psCensus.scanned > 0 )
                        ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
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
                            // PRESENT but UNSTABLE -> drop it from the assemblage.
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
                            ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                               "deactivating phase {} (present but unstable, logSI gap {})",
                                               psLoop, kBad, psViol );
                            // T5's "164x cliff" was this firing after the primary solve
                            // had spent its whole budget - plan v5 section 95.1.
                            native_trace_decide( "phasesel loop=%ld deactivate=%ld logsigap=%.6e",
                                                 (long)psLoop, (long)kBad, psViol );
                        }
                        else if( psAbsent && phSelState[(size_t)kBad] == 1 )
                        {
                            // ABSENT but STABLE, and it was deactivated by the
                            // compaction probe or by this loop -> that removal was
                            // wrong. Restore the original box and seed the phase off
                            // the floor, bounding each end-member by its OWN limiting
                            // IC (the same per-end-member formula as
                            // DetectPhaseCollapseAndReseed() above, which exists
                            // because a single phase-wide min lets one trace IC cap
                            // the whole phase straight back onto the phase-absence
                            // cliff this is trying to lift it off).
                            //
                            // Readmits EVERY phase that now looks wrongly removed,
                            // not just the worst one. With compaction enabled a single
                            // probe can pin dozens of phases at once, and readmitting
                            // one per loop against a bound of kMaxPhaseSelectLoops
                            // could not undo that. pm.Falp[]/pm.YF[] were just
                            // refreshed by WorstPhaseStabilityViolation(), so this
                            // costs nothing beyond the scan. Deactivation stays
                            // one-at-a-time deliberately - removing a phase is the
                            // side that can be wrong in a way the solver cannot see.
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
                                ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                                   "READMITTING {} wrongly removed phase(s) (worst: phase {}, "
                                                   "absent but stable, logSI gap {})",
                                                   psLoop, nReadmit, kBad, psViol );
                            }
                        }
                        if( !acted )
                            psSkipped++;
                    }   // for psIdx - fall through to the next-worst violation
                    if( acted && psSkipped > 0 )
                        ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
                                           "skipped {} unactionable violation(s) ranked ahead of the one "
                                           "acted on (worst was phase {}, gap {})",
                                           psLoop, psSkipped, psRanked[0].k, psRanked[0].viol );
                    if( !acted )
                    {
                        // Nothing in the WHOLE ranking is actionable - which is
                        // now a measured statement about every violation found,
                        // where before it was a statement about the first one.
                        ipm_logger->warn( "CalculateEquilibriumStateOptima: phase-selection loop {} - "
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
                    Optima::Result psResult = solverPS.solve( problem, psState );
                    applyStall( psResult );
                    optimaIterTotal += psResult.iterations;
                    state = psState;
                    result = psResult;
                }
            }

        }
        } // else (!reaktoroMode) - AOP/SOP's own solvent-collapse retry

        // ---- A NET MUST RESTORE THE PROBLEM, NOT ONLY THE STATE -------------
        // Both safety nets below re-solve from `initialState` so the re-solve
        // reproduces what an unarmed call would have done. That assignment
        // restores the STATE; it does not restore the PROBLEM. Any phase the
        // compaction probe or the phase-selection loop pinned at the floor
        // stayed pinned through the "full budget" re-solve, so the re-solve
        // inherited an assemblage decision taken on a deliberately truncated
        // state instead of starting from the original box - which is the one
        // thing the net exists to undo. Latent until 2026-09-10, when
        // pa_OptimaEarlyStabilityAt = 0 became AUTO and armed the early net on
        // every multisite system (see optima_earlystability_at, ms_multi.h);
        // it is reachable there now.
        //
        // Only phSelState == 1 is un-pinned, and that is the whole scope:
        //   0  never touched - nothing to restore;
        //   1  pinned by the compaction probe or by the phase-selection loop,
        //      originals in savedLo/savedHi - this is the artefact;
        //   2  already readmitted by the loop itself, boxes restored there;
        //   the phase-extinction tier saves nothing BY DESIGN - its pin is the
        //   rescue the net is trying to give the re-solve, not an artefact of
        //   a truncated look, so un-pinning it would undo the fix.
        //
        // RE-PINNED IF THE NET DOES NOT WIN, because the boxes are read again
        // after the nets: the KKT check treats a degenerate box as an equality
        // constraint and exempts the species from the one-sided sign test. Left
        // un-pinned while the pre-net state still holds that phase at the floor,
        // it would manufacture a "at lower bound with a driving force to grow"
        // residual on a result that was reported OK before the net ran - a net
        // that can only ever spend one solve would have started losing answers.
        // Same discipline as the extinction rescue just below: restore unless it
        // converges.
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
                ipm_logger->info( "CalculateEquilibriumStateOptima: safety net - restoring the "
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


        // ---- LAST RESORT: the stall watch's own safety net ------------------
        // pa_OptimaStallWindow exists to hand control to the retry tiers above
        // SOONER, not to decide the verdict. If it fired and none of them
        // recovered, the one thing left worth knowing is whether the watch was
        // simply WRONG - so re-solve once with it disarmed, from the same
        // starting state the primary solve used.
        //
        // WHY THIS MAKES ARMING THE WATCH A ONE-WAY BET, which is the point.
        // Without it, switching the field on can only be justified per project,
        // because a watch that fires on a run that WOULD have converged costs an
        // answer - and that is precisely what happened at 588 C and is why the
        // compiled default was reverted to 0 (see ms_multi.h, and plan v5
        // section 54 for the rule change that removed that particular case).
        // With it, the worst the watch can do is spend one extra solve on a run
        // that was already failing, while the best it can do is large and
        // measured: f_/j_CASHNK converge in ~1001 iterations with the watch and
        // 10001 without, same G/pH/Vs to every digit.
        //
        // Cost is bounded by one ordinary solve and is paid ONLY when the watch
        // fired AND every retry tier failed - so it is free on every project
        // that converges today, and free on the projects the watch actually
        // helps (on CASHNK a retry succeeds, so this never runs).
        //
        // Deliberately gated on `stalled` and NOT on `timedOut`:
        // pa_OptimaMaxSeconds is a wall-clock budget the caller set on purpose,
        // and silently spending another solve after it expires would defeat it.
        //
        // NOT gated to AOP/SOP, and that is deliberate rather than an oversight
        // of the standing "should this be ROP-gated?" question: the watch itself
        // is not mode-gated, so ROP needs the same net. In practice ROP cannot
        // reach this at a window of 500 or more, because it leaves
        // Optima::Options at library defaults (maxiters = 200) and a window
        // longer than the budget can never complete - so this only becomes live
        // for ROP if someone sets a small window explicitly, which is exactly
        // when they would want the net.
        if( !result.succeeded && stallWatch->stalled && stallWatch->window > 0 )
        {
            ipm_logger->info( "CalculateEquilibriumStateOptima: every retry failed after a stall"
                               " - re-solving once with pa_OptimaStallWindow disarmed" );
            const long int savedWindow = stallWatch->window;
            stallWatch->window = 0;      // the check lambda reads this on every call
            stallWatch->reset();
            unpinForNet();               // the PROBLEM too, not only the state - see above
            Optima::Solver solverNS;     // fresh instance, per the retry-ordering note above
            solverNS.setOptions( options );
            Optima::State  nsState  = initialState;
            Optima::Result nsResult = solverNS.solve( problem, nsState );
            optimaIterTotal += nsResult.iterations;
            stallWatch->window = savedWindow;
            if( !nsResult.succeeded )
                repinAfterNet();         // the net did not win - leave the boxes as they were
            if( nsResult.succeeded )
            {
                netUnpinned.clear();     // this state was solved on the restored boxes; keep them
                state  = nsState;
                result = nsResult;
                ipm_logger->info( "CalculateEquilibriumStateOptima: the disarmed re-solve converged"
                                   " in {} iterations - the stall was a false positive",
                                   nsResult.iterations );
            }
        }

        // ---- LAST RESORT: pa_OptimaEarlyStabilityAt's own safety net --------
        // Exactly the same shape, and the same argument, as the stall net above.
        // That field caps the FIRST attempt so the phase-selection repair loop
        // gets its turn before the primary solve has converged on an assemblage
        // it will then have to correct - worth 121x on
        // Resources/gems3k/f_Solvus_G_Test1 (58 iterations if the doomed phase
        // is pinned out up front, against 7007 if the loop only sees it after
        // the fact). But the cap applies to EVERY call, including calls whose
        // assemblage is fine and which simply need their budget: that project's
        // own cold solve needs 823 iterations, so at a cap of 500 it failed
        // outright and the repair loop could not rescue it - there was no
        // violation to repair, only a budget shortfall. That is what kept the
        // field default-off and unusable (plan v5 section 45).
        //
        // Treating the capped attempt as a PROBE removes that: if the cap is
        // what stopped it and nothing downstream recovered, re-solve once from
        // the original state at the full budget, so the only cost of a probe
        // that found nothing is the probe itself. Bounded by one ordinary
        // solve, paid only on a run that is already failing.
        if( !result.succeeded && earlyCapHit )
        {
            ipm_logger->info( "CalculateEquilibriumStateOptima: the pa_OptimaEarlyStabilityAt probe"
                               " found nothing to repair - re-solving once at the full budget" );
            Optima::Solver solverEP;     // fresh instance, per the retry-ordering note above
            solverEP.setOptions( options );   // the full budget; the cap is disarmed by now
            // NOT disarmed, unlike the stall net's own re-solve above - tried on
            // 2026-09-09 and REVERTED the same hour, on a misreading worth recording
            // because the instrument invites it. Instrumenting this re-solve on
            // j_CASHNK HOP printed `it=4500 ok=1 stalled=1`, which reads as "it
            // converged and we threw the answer away". It is not: returning true from
            // convergence.check means STOP, and Optima reports a hook-stop as SUCCESS
            // (see the fold-back note at the primary solve). So ok=1 meant the watch
            // stopped it, the existing fold was right, and the answer was never in
            // hand. Disarming and re-measuring settled it - the re-solve then runs
            // 10200 iterations and still does not converge, twice the cost for the
            // same FAIL. `succeeded` on a solve carrying a convergence hook is not a
            // statement about convergence unless the hook is known not to have fired.
            stallWatch->reset();
            unpinForNet();               // the PROBLEM too, not only the state - see above
            // ITEM 24: `initialState` is the shipped restart, which pays N + full;
            // `*earlyProbeState` resumes the same problem from where the probe
            // stopped and should pay ~full. See the capture site for why the choice
            // is armed there rather than here, and for section 83's counter-
            // precedent. Which one ran is a decision that changes the cost and can
            // change the answer, so it goes in the TRACE, not only in the log.
            const bool epResumed = (bool)earlyProbeState;
            Optima::State  epState  = epResumed ? *earlyProbeState : initialState;
            // dx is the distance the probe moved the iterate, max|x_probe - x_0|
            // over the primal. It is in the record because "the resume is a no-op"
            // and "the resume never happened" are otherwise the same observable -
            // a dx of 0 would mean the start point never changed, which is a
            // plumbing failure and not a measurement (GEMS3K/CLAUDE.md s5).
            double epStartDx = 0.;
            for( long int j = 0; j < L; j++ )
                epStartDx = std::max( epStartDx, std::fabs( epState.x[j] - initialState.x[j] ) );
            native_trace_decide( "netresume from=%s probeat=%ld dx=%.6e",
                                 epResumed ? "probe" : "initial", earlyProbeIters, epStartDx );
            Optima::Result epResult = solverEP.solve( problem, epState );
            optimaIterTotal += epResult.iterations;
            const bool epStalled = ( stallWatch->stalled || stallWatch->timedOut );
            if( epStalled ) epResult.succeeded = false;
            // WHAT STOPPED THE RE-SOLVE decides whether the probe's cost is
            // RECOVERABLE at all (item 24). A re-solve that CONVERGED spent a
            // convergence time, and starting it further along the same trajectory
            // can subtract from that; a re-solve that STALLED spent a BUDGET - the
            // stall window - and where it started cannot change what a budget
            // costs. Traced rather than logged so the census is recoverable from
            // the freeze's own `# dec` lines over the whole corpus instead of from
            // a scrape of 43 log files (GEMS3K/CLAUDE.md s4).
            native_trace_decide( "netresolve outcome=%s it=%ld",
                                 epResult.succeeded ? "converged"
                                                    : ( epStalled ? "stalled" : "failed" ),
                                 (long int)epResult.iterations );
            // THE FALLBACK (mode 2). A resume that did not converge has cost the run
            // its own iterations and bought nothing, so hand the net the start point
            // it would have used unarmed and let it have its ordinary attempt. Fresh
            // solver instance, per the retry-ordering note above; the stall watch is
            // reset so the second attempt is judged on its own progress and not on
            // the first one's. Bounded: one extra solve, only on a run that has
            // already failed twice, and only where the resume was tried at all.
            if( netResumeMode == 2 && epResumed && !epResult.succeeded )
            {
                ipm_logger->info( "CalculateEquilibriumStateOptima: the resumed re-solve did not"
                                   " converge in {} iterations - falling back to the restart",
                                   epResult.iterations );
                Optima::Solver solverFB;
                solverFB.setOptions( options );
                stallWatch->reset();
                Optima::State  fbState  = initialState;
                Optima::Result fbResult = solverFB.solve( problem, fbState );
                optimaIterTotal += fbResult.iterations;
                const bool fbStalled = ( stallWatch->stalled || stallWatch->timedOut );
                if( fbStalled ) fbResult.succeeded = false;
                native_trace_decide( "netfallback outcome=%s it=%ld",
                                     fbResult.succeeded ? "converged"
                                                        : ( fbStalled ? "stalled" : "failed" ),
                                     (long int)fbResult.iterations );
                // Taken UNCONDITIONALLY, not only when it converged: mode 0's own
                // state is what everything below this point is written against - the
                // extinction tier's "this is the solver's final word at the full
                // budget" argument included - and a failed resume is not that word.
                epState  = fbState;
                epResult = fbResult;
            }
            if( epResult.succeeded )
            {
                netUnpinned.clear();     // this state was solved on the restored boxes; keep them
                state  = epState;
                result = epResult;
                ipm_logger->info( "CalculateEquilibriumStateOptima: the full-budget re-solve"
                                   " converged in {} iterations ({} the probe's {} iterations)",
                                   epResult.iterations,
                                   epResumed ? "resumed from" : "restarted, discarding",
                                   earlyProbeIters );
            }
            else if( runExtinctionTier )
            {
                // THE RE-SOLVE REPRODUCED THE UNARMED PRIMARY SOLVE, INCLUDING
                // ITS FAILURE - so it is owed the rescue an unarmed call gets.
                // See runExtinctionTier's own declaration for the measurement:
                // on j_CASHNK HOP the unarmed call stalls at 4500 and the
                // phase-extinction tier converges it in one more iteration,
                // while the armed call spent 200 iterations on the probe, ran
                // every tier on the truncated state (where the tier correctly
                // DECLINES - the twin has not started falling yet), and then had
                // nothing left below the net. The re-solve was never the defect.
                //
                // Judged on the ordinary path (onEarlyStopPath = false), not the
                // early-stop one, and that is the whole point: this state is the
                // solver's own final word at the full budget, so the "still
                // falling" clause the truncated state needs does not apply and
                // the ratio test means what it has always meant.
                //
                // Bounded: one tier attempt on a run that has already failed
                // twice, and it cannot make the outcome worse - the pre-net
                // state and result are restored unless the tier converges.
                const Optima::State  savedState  = state;
                const Optima::Result savedResult = result;
                state  = epState;
                result = epResult;
                runExtinctionTier( false );
                if( result.succeeded )
                    ipm_logger->info( "CalculateEquilibriumStateOptima: the full-budget re-solve"
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

            // ---- MODE 3: THE DEFERRED FALLBACK ---------------------------
            // Asked here and nowhere earlier, because here is the only place the
            // question is answerable: the resumed re-solve AND every rescue below
            // it have now had their turn, so `!result.succeeded` means this call
            // has genuinely lost its answer rather than merely stalled on the way
            // to being rescued. See the capture site for the pair of projects that
            // makes the distinction necessary. Structure deliberately mirrors the
            // block above rather than factoring it: this ladder's ordering carries
            // several separately-argued invariants and a refactor of it would put
            // them all at risk to save twenty lines.
            if( netResumeMode == 3 && epResumed && !result.succeeded )
            {
                ipm_logger->info( "CalculateEquilibriumStateOptima: the resumed re-solve and its"
                                   " rescue both failed - falling back to the restart the unarmed"
                                   " call would have used" );
                unpinForNet();           // repinAfterNet() may have put the boxes back
                stallWatch->reset();
                Optima::Solver solverFB;
                solverFB.setOptions( options );
                Optima::State  fbState  = initialState;
                Optima::Result fbResult = solverFB.solve( problem, fbState );
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
                    ipm_logger->info( "CalculateEquilibriumStateOptima: the restart converged in {}"
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

        ipm_logger->info( "CalculateEquilibriumStateOptima: pNP={} reaktoroMode={} succeeded={} iterations={} nConditions={}",
                           pm.pNP, reaktoroMode, result.succeeded, result.iterations, R );

        for( long int j = 0; j < L; j++ )
        {
            pm.Y[j] = state.x[j];
            pm.X[j] = state.x[j];
        }
        // Optima's Lagrange multipliers ye use the opposite sign convention
        // from GEMS3K's own dual U[] - negate here rather than resign every
        // place U[] is used below.
        for( long int i = 0; i < N; i++ )
            pm.U[i] = -state.ye[i];

        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );

        // Real post-solve concentration/output pipeline - the same one
        // native's own CalculateEquilibriumState() calls internally
        // (ipm_main.cpp) before returning. THIS WAS MISSING ENTIRELY in
        // every earlier version of this function: packDataBr() (node.cpp,
        // always called after either solver, native or Optima) reads
        // pmm->VXc/FVOL[]/FWGT[]/IC directly - fields that are populated
        // ONLY by CalculateConcentrations() (ipm_chemical2.cpp; also sets
        // pm.pH/pe/Eh, superseding the hand-rolled closed-form computation
        // an earlier version of this code had here). Without this call,
        // CNode->Vs/mPS/vPS/pH/Eh after an Optima solve were silently
        // whatever CalculateConcentrations() last left them as - typically
        // stale values from GEM_init()'s own initial load, having no
        // relationship at all to Optima's actual converged Y. Confirmed
        // directly on Resources/gems3k/j_Flowline_G_series1_...: the
        // "Vs 4.25x too small" symptom this file's own investigation notes
        // (GEMS3K's CLAUDE.md, 2026-08-23) used as its primary evidence of
        // a wrong SOLUTION was, at least in significant part, actually
        // evidence of this missing call - state.x[Albite] (traced directly)
        // reached a real, non-trivial converged value while the reported
        // Vs stayed byte-identical to a completely unrelated run, which is
        // only possible if Vs was never actually recomputed from that Y at
        // all. Placed here (after the post-solve X/activity-coefficient
        // refresh, using the final committed pm.U[] dual) so every value it
        // computes reflects Optima's actual solution, not an intermediate
        // one.
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );

        // General KKT/complementary-slackness stationarity check - applies
        // uniformly to every ordinary species, whether R==0 or R>0, unlike
        // the control-condition-specific achievedValueFn check below (which
        // only verifies the R virtual titrant slots) and unlike
        // CheckMassBalanceResiduals() (which only verifies the LINEAR
        // constraint A*Y=B - many different Y satisfy the same mass
        // balance, so that check alone cannot tell a true Gibbs-energy
        // minimum from an arbitrary feasible point). Not theoretical:
        // confirmed necessary by a direct native-vs-Reaktoro cross-check on
        // Resources/gems3k/j_Flowline_G_series1_..., where Optima reported
        // succeeded=true with mass balance satisfied, yet the resulting
        // system volume was 4.25x off from both native GEMS3K and an
        // independent Reaktoro solve of the identical G0/V0 data (GEMS3K's
        // CLAUDE.md, 2026-08-23, "plain AOP as a general drop-in solver").
        // For an interior species (away from both box bounds), true
        // stationarity requires its reduced gradient
        // g[j] = F[j] - sum_i U[i]*A[i,j] to vanish; for a species sitting
        // at a bound, complementary slackness only requires the correctly
        // signed g[j] (non-negative at the lower bound, non-positive at the
        // upper bound) - a species incorrectly left at its lower bound when
        // it should have entered (e.g. a mineral that should have
        // precipitated) shows up here as a strongly negative g[j], exactly
        // the failure mode traced on j_Flowline.
        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
        double maxKKTResidual = 0.;
        long int worstKKTSpecies = -1;
        std::vector<double> gradSaved( (size_t)L, 0. );   // for pa_OptimaZeroAbsent below
        const double kktTol = std::max( pa_p->GAS, 1e-300 );
        for( long int j = 0; j < L; j++ )
        {
            // The gradient Optima actually minimised is pm.F[j] plus, for
            // single-component ("pure") phase species, the log-barrier
            // contribution -tau/X[j] (see the objective lambda above). This
            // check must verify stationarity of THAT objective, not of the
            // barrier-free one, or it reports a spurious violation of size
            // tau/X[j] for every pure-phase species sitting near its floor.
            // Harmless while tau was 1e-16 against a 1e-13 floor (1e-3, just
            // at kktTol); an O(1) false failure once tau is coupled to the
            // floor as Reaktoro couples it. Found by A/B-ing that coupling
            // (CLAUDE.md 2026-08-25, plan-v5 Phase A / A.0).
            double gradJ = pm.F[j];
            if( j >= pm.Ls )
                gradJ -= kLogBarrierTau / std::max( pm.X[j], dcFloor );
            for( long int i = 0; i < N; i++ )
                gradJ -= pm.U[i] * pm.A[ i + j*N ];
            gradSaved[(size_t)j] = gradJ;
            double resid;
            // A KINETICALLY FIXED species (xlower==xupper by construction,
            // from a DUL<1e6 restriction pinning it to a single value - see
            // problem.xupper[j]'s own construction above) is an EQUALITY
            // constraint on this variable, not a one-sided inequality bound -
            // its multiplier is unrestricted in sign, standard KKT theory for
            // a fixed variable. Treating it as "at lower bound" or "at upper
            // bound" and applying the corresponding one-sided sign test is
            // wrong: a real, deliberate exclusion (e.g. a kinetically
            // suppressed mineral, DUL=DLL=0) that still carries a genuine
            // thermodynamic driving force to grow - exactly what it means to
            // be correctly, forcibly excluded - was being reported as a KKT
            // violation. Confirmed on o_/t_Kaolinite: Quartz (DUL=DLL=0,
            // xlower=xupper=dcFloor after the box-construction fix above)
            // flagged a spurious residual of 2.48 even once the composition
            // matched native's own answer exactly. See GEMS3K/CLAUDE.md,
            // 2026-08-26, item 5.
            if( problem.xupper[j] <= problem.xlower[j] + std::max( dcFloor, problem.xlower[j]*1e-6 ) )
                resid = 0.;                        // fixed (degenerate box): no sign test applies
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
        bool kktOk = maxKKTResidual <= kktTol;

        // DUAL DETERMINACY (2026-09-15, plan v5 s123.6). The sign test above reads each
        // reduced gradient at the dual Optima returned. When the INTERIOR species (the ones
        // whose g = 0 fixes the dual) span fewer than N directions, the dual is free along
        // the rest, and every bound-active species' gradient depends on where along them the
        // solver happened to stop - an arbitrary choice, not chemistry. The right question is
        // then whether SOME dual in that free set satisfies every bound condition. Measured on
        // 07PSIna_G_simple_1_0_1_25_0 SHP (C-Ca-H-Mg-O, no redox buffer, Eh undetermined): 15
        // interior species span rank 5 of N = 6, the free direction is the redox one, and
        // H2(aq) at the floor reads g = +0.240 / -0.512 / -2.895 depending only on the warm
        // start, while a feasible interval exists in all three ([-0.66,160], [1.40,162],
        // [7.93,168]). Optima's own ex is zeroed for any bound-pinned variable whatever the
        // sign of s (ResidualErrors.cpp), so it calls all three converged.
        // Runs only when the plain test failed; handles ONE free direction exactly (the redox
        // case). More than one is recorded and left failing until its reach is measured.
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
            // The same interval with NO tolerance slack. The point is chosen from this one when it is
            // not empty: an end of the slack interval puts the bounding species at a residual of
            // exactly kktTol, where rounding decides the verdict (first build: 07PSIna_G_simple_1
            // SHP read after=1.000e-03 and resolved=0 on one warm start, resolved=1 on another).
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
            // Phase 3 WP1: the same verdict onto the CERT line. CertDualFreeDirs() counts the
            // free directions on every path, but only this block SEARCHES them, and only once
            // the plain sign test has failed - so dual_resolved stays -1 ("not searched") on
            // every row that passed, which is what it has to mean.
            native_cert_dualfree_set( resolved );
            native_trace_decide( "dualfree rank=%ld of=%ld dirs=%ld tlo=%.4g thi=%.4g tstar=%.4g resolved=%d "
                                 "before=%.3e after=%.3e lo=%s hi=%s",
                                 (long)rank, (long)N, (long)nDirs, tLo, tHi, tStar, resolved,
                                 before, residAfter, dcName( jLo ).c_str(), dcName( jHi ).c_str() );
        }

        double worstStabilityViol = 0.;
        bool worstStabilityWasAbsent = false;
        long int worstStabilityPhase = WorstPhaseStabilityViolation(
                    presenceThreshold, dcFloor,
                    extinctFixed.empty() ? nullptr : extinctFixed.data(),
                    worstStabilityViol, worstStabilityWasAbsent );
        const bool stabilityOk = ( worstStabilityPhase < 0 );

        // Commit each active condition's titrant into pm.B[] - without
        // this, CheckMassBalanceResiduals() below (and every downstream
        // consumer of pm.B[], including TNode::packDataBr()'s CNode->bIC[]
        // output) would compare Y against the ORIGINAL, un-titrated bulk
        // composition, not the effective one Optima actually solved
        // against. The caller sees the titrated bulk composition on
        // output.
        // Sign: Optima's linear equality is Aex*x_ext=be, i.e.
        // A*Y + sum_k stoich_k*xi_k = B (the titrant column sits on the
        // SAME side as A*Y, not folded into the RHS) - so the effective
        // bulk composition Y actually balances against is
        // B - sum_k stoich_k*xi_k, hence subtract here, not add.
        // Units: this runs while the system is still in pa_DG's INTERNAL scale
        // (ScaleSystemToInternal() above: every amount multiplied by
        // ScFact = pa_DG / sum_i B[i]), and RescaleSystemFromInternal() below
        // brings pm.B/X/Y back but knows nothing of the titrant. So pm.B is
        // committed with the internal value, and titrantAmount is stored in
        // REAL moles - xi / ScFact. Measured 2026-09-14 on the xgems-jupyter
        // Cu-DH project (pa_DG = 1000, sum B = 166.737 mol, ScFact = 5.99747):
        // the value previously returned was the internal one, and at
        // pH 2 / Eh 0, pH 8 / 0.2 and pH 12 / -0.4 the returned bIC moved by
        // exactly -stoich*xi/5.99747 on both the H and the Zz row.
        for( long int k = 0; k < R; k++ )
        {
            for( const auto& rc : conditions[k].stoich )
                pm.B[ rc.first ] -= rc.second * state.x[L+k];
            conditions[k].titrantAmount = state.x[L+k] / ScFact;
        }

        // Verify each active condition actually reached its target - not
        // just that Optima's own KKT residual reports solved. At a
        // bound-active titrant slot, the box constraint's own dual
        // absorbs the KKT residual, so the fixed-gradient condition this
        // mechanism relies on can be unsatisfied even while Optima
        // reports succeeded=true.
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
            if( !conditions[k].targetMet )
            {
                allTargetsMet = false;
                targetMissBuf += conditions[k].name + ": target=" + std::to_string(conditions[k].target)
                               + " achieved=" + std::to_string(conditions[k].achievedValue) + "; ";
            }
        }

        // Mass-balance trustworthiness check, reusing GEMS3K's own
        // per-IC tolerance (pm.DHBM), exactly as the native solver's own
        // testMulti()/CheckMassBalanceResiduals() path does - after the
        // pm.B[] commit above, so this compares against the effective
        // (post-titration) bulk composition, not the pre-titration input.
        long int massBalanceBadIC = CheckMassBalanceResiduals( pm.Y );

        // Report Optima's own total iteration count via the same
        // NumIterFIA/NumIterIPM channel native AIA/SIA already uses
        // (TNode::GEM_CalcTime()/GEM_Iterations()) - previously left at 0
        // unconditionally, which made every AOP/ROP timing comparison
        // report an iteration count of 0 regardless of how much work
        // Optima actually did. pm.ITF (FIA/MBR-equivalent) stays 0 - this
        // solver has no such separate phase - pm.ITG carries the total.
        pm.ITG = optimaIterTotal;

        // ---- pa_OptimaZeroAbsent: report a correctly absent species as EXACTLY
        // ZERO rather than at the box's numerical floor (0 = off; default 1 since
        // 2026-09-03).
        //
        // Placed HERE, after every trustworthiness check above has already run
        // and produced its verdict on the un-zeroed state, and gated on that
        // verdict being clean - so this can never turn an accepted answer into
        // a rejected one, and never fires on a solve already headed for BAD.
        //
        // See BASE_PARAM::OptimaZeroAbsent (ms_multi.h) for why the twenty
        // orders of magnitude between native's pa_DcMin truncation inside GX()
        // and this path's dcFloor are worth closing at all, and for what this
        // does NOT do (native's species leave the PROBLEM; these leave only the
        // ANSWER).
        if( pa_p->OptimaZeroAbsent == 2 && result.succeeded && allTargetsMet
            && massBalanceBadIC < 0 && kktOk && stabilityOk )
        {
            // VALUE 2 - keep the amounts Optima returned, repair only a failing balance (owner 2026-09-15,
            // plan v5 s123.8; BASE_PARAM::OptimaZeroAbsent). Same acceptance gate as the zeroing below, so
            // it never touches a solve already headed for BAD.
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
                        // a DEFAULT SEED (ICIsDefaultSeed, ms_multi.h) neither triggers nor vetoes the repair - its
                        // amount is a placeholder; with nothing marked there are none (the same exemption value 1's
                        // rebalance test applies)
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

            // ZERO ONLY IN ABSENT PHASES, THEN REBALANCE (owner decisions 2026-09-15,
            // plan v5 s123.7). Until then this zeroed every species sitting on the floor,
            // and most of those were DISSOLVED species of the PRESENT aqueous phase
            // (Optima 16/7/5/21 in present phases vs 3/3/4/8 in absent ones on
            // FeNaCl_FyGt_Precip, 07PSIna_G_iron, Al-species, Cu-Pourbaix) - whose small
            // amounts are real information, not absence. Now a phase is zeroed as a whole,
            // and only if EVERY member passes (a)-(d); species of a present phase are never
            // touched. Reaktoro, for comparison, zeroes nothing (species keep its 1e-16
            // floor and x is written back as Optima returned it).
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
                    // (a) not kinetically REQUIRED to be present. DLL > 0 is a
                    //     caller's deliberate retention floor - the gibbsite case
                    //     in optima_regression's metastability suite - and zeroing
                    //     it would silently discard the constraint.
                    if( pm.DLL[j] > 0. ) { absent = false; break; }
                    // (b) the bound it sits on must be the NUMERICAL floor, not a
                    //     real constraint.
                    const double lo = problem.xlower[j];
                    if( lo > dcFloor * ( 1. + 1e-9 ) ) { absent = false; break; }
                    // (c) actually sitting on it - same expression the KKT check
                    //     above uses to classify a variable as bound-active.
                    if( pm.Y[j] > lo + std::max( dcFloor, lo * 1e-6 ) ) { absent = false; break; }
                    // (d) correctly there: a reduced gradient not below -kktTol (the
                    //     verdict's own tolerance; gradSaved is already shifted along a
                    //     free dual direction when the determinacy check above resolved
                    //     one), or a DEGENERATE box (DUL of exactly 0 - a hard kinetic
                    //     exclusion, the o_/t_Kaolinite quartz case).
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

            // REBALANCE: zeroing removes ~dcFloor per species, negligible against a major
            // IC and not against a trace one (plan v5 s123.2: FeNaCl_FyGt_Precip AOP Fe
            // |C|/(B*DHBM) 4.4e-3 -> 6.2e3). If any IC or the charge row ends past both its
            // own residual before zeroing and its tolerance, put the removed amount back onto
            // the present carriers with MassBalanceReproject() - the native repair
            // (pa_MbReproject), rank-revealing pivots, all N rows. If that still leaves an
            // IC past its limit, the zeroing is undone: the worst case is the floor values.
            // A DEFAULT SEED (ICIsDefaultSeed: numerically trace and not marked of interest while
            // others are) neither triggers the rebalance nor the undo - its amount is a placeholder
            // (owner 2026-09-15). With nothing marked there are none and every IC is conserved.
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
                // SELF-GATE. Dropping ~dcFloor from each of several hundred
                // species is utterly negligible against a major IC and NOT
                // negligible against a trace one - the trace-IC sensitivity
                // this branch has measured repeatedly. So the zeroed state is
                // put back through the SAME per-IC test the un-zeroed state
                // just passed, and is kept only if it also passes.
                //
                // MEASURED 2026-09-14e (plan v5 s123.2): that test CANNOT see a
                // trace-IC loss. CheckMassBalanceResiduals()' cutoff is the ABSOLUTE
                // min(DHBM*1e10, 1e-2) = 1e-3 mol at pa_DHB = 1e-13, while the zeroing
                // removes ~1e-13 mol. On FeNaCl_FyGt_Precip AOP it takes Fe (B = 4.8e-4)
                // from |C|/(B*DHBM) = 4.4e-3 to 6.2e3 and the charge residual from
                // 5.1e-19 to 3.2e-13, and is kept; on 07PSIna_G_iron Fe (B = 3e-9) goes
                // to 4.3e8. G, Vs, Ms unchanged. The CERT trace record shows it.
                if( CheckMassBalanceResiduals( pm.Y ) >= 0 )
                {
                    for( long int j = 0; j < L; j++ ) pm.Y[j] = Ysave[(size_t)j];
                    CheckMassBalanceResiduals( pm.Y );   // restore pm.C[] to the accepted state
                    ipm_logger->info( "CalculateEquilibriumStateOptima: pa_OptimaZeroAbsent - "
                                       "reverted, zeroing {} absent species would break the mass "
                                       "balance", nZeroed );
                }
                else
                {
                    // Commit, and recompute everything derived from the primal
                    // so the reported phase amounts, volume and concentrations
                    // describe the state actually being returned.
                    for( long int j = 0; j < L; j++ ) pm.X[j] = pm.Y[j];
                    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                    CalculateActivityCoefficients( LINK_UX_MODE );
                    CalculateConcentrations( pm.X, pm.XF, pm.XFA );
                    ipm_logger->info( "CalculateEquilibriumStateOptima: pa_OptimaZeroAbsent - "
                                       "{} of {} species reported as exactly zero", nZeroed, L );
                }
            }
        }

        if( !result.succeeded && pa_p->DW )
        {
            // DW already gates exactly this decision for the native
            // solver's own "MBR iterations exceeded" case (ipm_main.cpp);
            // reused as-is here rather than adding a parallel flag.
            Error( "E90IPM: Optima solver: ", "Optima::Solver::solve() did not converge (pa_DW forces this to a hard error)" );
        }
        else if( !result.succeeded || !allTargetsMet || massBalanceBadIC >= 0 || !kktOk || !stabilityOk )
        {
            // Soft failure: a state was produced, but it is not fully
            // trustworthy - mirrors AIA/SIA's own BAD_GEM_* semantics
            // (testMulti() already reads pm.MK below, via TNode::GEM_run()).
            pm.MK = 2;
            std::string buf = "Optima solve ";
            buf += result.succeeded ? "reported success" : "did not converge";
            if( !allTargetsMet )
                buf += std::string("; control condition(s) not met: ") + targetMissBuf;
            if( massBalanceBadIC >= 0 )
                buf += "; mass balance residual exceeds tolerance for IC " + char_array_to_string(pm.SB[massBalanceBadIC],3);
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
        // Same read-only system-definition check as native's converged answer (ipm_main.cpp):
        // a fragile definition is a property of the system, not of the solver that met it.
        // Runs while still in pa_DG's internal scale; the check reports real moles itself.
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
        Error( xcpt.title, xcpt.mess );
    }
    catch( std::exception& e )
    {
        // Optima reports its own failures by THROWING std::runtime_error (the errorif/assert macros,
        // optima/Optima/Exception.hpp), not TError. Without this block such an exception skipped the
        // rescale above - leaving MULTI in pa_DG's internal units - and, more seriously, escaped
        // CalculateEquilibriumStateHOP()'s catch(TError&), so HOP/SHP lost native's converged answer and
        // ended T_ERROR_GEM instead of restoring it as BAD_GEM_HOP. Same cleanup, rethrown as TError so
        // every caller sees one failure type. Precedent: GEMS4R's TEqulibrate (equlibrate.cpp:94-118)
        // did exactly this for Reaktoro. Docs/REVIEW-2026-09-29-gems4r.md s3.3.
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

    // pm.FX is the system's total Gibbs energy as a STORED field, and until here the Optima path
    // never wrote it. That it is native-only was already known - ipm_main.cpp's CERT block says so
    // in as many words ("pm.FX and pm.Falp are refreshed by the native path only, so a trace writer
    // reading them would print stale values on Optima rows, the DC_G0() shape") - and the response
    // was defensive: keep FX out of the trace and the certificate. What that left unaudited is that
    // packDataBr() (node.cpp:560, `CNode->Gs = pmm->FX`) ALREADY publishes it, on every call, to
    // every consumer of the documented interface: DATABR.Gs, GEM_to_MT()'s p_Gs, and the Gs field of
    // every exported -dbr file. The instrument was protected and the product was not.
    //
    // WHAT WAS PUBLISHED. MultiConstInit() (ipm_simplex.cpp:909, under a comment reading "???????")
    // seeds pm.FX = 7777777. as an unset marker; the native descent overwrites it (ipm_main.cpp:3700,
    // :3791) and Optima did not. So DATABR.Gs was exactly 7777777. on AOP, SOP, HOP AND SHP -
    // measured 2026-09-22 on seven projects spanning sorption, solid solutions, seawater and large
    // multiphase, all four modes, every project, while TNode::Get_GibbsEnergy() (a RECOMPUTATION via
    // TotalGibbsEnergy()) returned the right value and agreed with native's to eleven digits. HOP and
    // SHP publish the sentinel too although native runs FIRST, because the Optima leg's own
    // MultiConstInit() re-seeds FX after native's descent has set it - so on this field the hybrid
    // modes are strictly worse than the native leg they are built on, which is worth remembering
    // whenever "HOP >= native everywhere" is quoted: like the wall-time result of work item 37, that
    // is a statement about the ANSWER.
    //
    // WHY NO GATE CAUGHT IT. mode_compare reads Get_GibbsEnergy() and freeze.sh parses that, so the
    // freeze's G column is the recomputation and no scored column reads DATABR.Gs at all. The freeze
    // is therefore blind to the defect AND to this fix - the same structural blindness as work item
    // 33's one-ULP DUL/DLL regression, and the second finding in three days that only an instrument
    // standing where a CALLER stands can reach. It was found by disbelieving a warm-standard answer
    // digest that read an identical positive G on two unrelated chemistries.
    //
    // WHY HERE. After RescaleSystemFromInternal(), which divides pm.FX by ScFact itself
    // (ipm_simplex.cpp) - that is how native's stored FX, written in internal units inside the
    // descent, comes out in real units and matches the recomputation. Assigning BEFORE the rescale
    // would therefore have published G/ScFact, and ScFact is 5.66-6.02 on the projects checked, not
    // 1: the placement is load-bearing, not cosmetic. Placed last for the second reason too - every
    // amount the post-solve tiers still move (pa_OptimaZeroAbsent, the extinction tiers above) is
    // committed by now.
    //
    // WHY THIS VALUE. TotalGibbsEnergy() rather than GX(0) because it is exactly what
    // TNode::Get_GibbsEnergy() calls, so the stored field and the recomputation agree BY
    // CONSTRUCTION rather than by coincidence - two accessors for one quantity disagreeing silently
    // is the whole defect being fixed here.
    //
    // WHY SAVED AND RESTORED. TotalGibbsEnergy() IS NOT A PURE READ - it opens with
    // `for(i) pm.X[i] = pm.Y[i]` and then TotalPhasesAmounts(pm.X, pm.XF, pm.XFA), so it writes
    // three arrays packDataBr() goes on to publish. Today that is value-neutral, because X == Y is
    // an invariant here: the post-solve tiers resync X from Y at every branch that touches Y, and
    // RescaleSystemFromInternal() divides X and Y by the same ScFact. But a fix whose correctness
    // rests on an invariant maintained in four places elsewhere in this function will break silently
    // the first time one of them changes, and the symptom would be a published composition quietly
    // reverting to a pre-zeroing vector. Save and restore by copy instead - the same treatment, and
    // for the same reason, as CertCurvMin() in ipm_main.cpp, which also writes solver state and
    // restores it by copy rather than by recomputation. Three vectors once per solve is nothing
    // against the solve itself.
    {
        std::vector<double> Xsave( pm.X,   pm.X   + pm.L );
        std::vector<double> XFsave( pm.XF, pm.XF  + pm.FI );
        std::vector<double> XFAsave( pm.XFA, pm.XFA + pm.FIs );
        pm.FX = TotalGibbsEnergy();
        for( long int j = 0; j < pm.L; j++ )    pm.X[j]   = Xsave[(size_t)j];
        for( long int k = 0; k < pm.FI; k++ )   pm.XF[k]  = XFsave[(size_t)k];
        for( long int k = 0; k < pm.FIs; k++ )  pm.XFA[k] = XFAsave[(size_t)k];
    }

    NumIterFIA = pm.ITF;
    NumIterIPM = pm.ITG;
    pm.t_end = clock();
    pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
    return pm.t_elap_sec;
}

// HYBRID: native selects the species (its own cold IPM/MBR/PSSC pipeline),
// then Optima finishes, warm-started from native's converged primal AND
// dual. Dispatched from TNode::GEM_run() (node.cpp) for NEED_GEM_HOP - see
// NODECODECH's own comment in databr.h for why this is a separate
// caller-selected mode rather than something AOP does internally (the
// user's explicit "switch between gems native and optima solvers, not use
// them together" direction, 2026-08-23 - a caller asking for HOP asks for
// both, by name, and gets both).
//
// Formerly this two-leg orchestration lived inline in node.cpp, with only
// the NATIVE call wrapped in a try/catch - a thrown TError from the OPTIMA
// leg propagated straight past that (there was no catch around it) to
// GEM_run()'s own OUTER catch(TError&), which never calls packDataBr(), so
// a converged native answer was silently discarded whenever the warm
// Optima leg on top of it could not finish (confirmed: native solves
// 07PSIna_G_vcomplex @ 80 C in 693 iterations; the full-dimension warm
// Optima leg on top of it does not finish in any reasonable time, and the
// mode reported nothing at all). Moved here, with the Optima leg's own
// try/catch, so that case degrades to a soft BAD_GEM_HOP carrying native's
// own answer instead - the guarantee this method's own header comment
// states: HOP is never worse than a plain native solve. See
// GEMS3K/CLAUDE.md and Docs/gems3k-optima-plan-v5.md, section 64.6.
double TMultiBase::CalculateEquilibriumStateHOP( long int& NumIterFIA, long int& NumIterIPM,
                                                 bool warmNative )
{
    long int fiaN = 0, ipmN = 0;
    bool nativeOk = true;
    double calcTime = 0.;
    std::vector<double> Ysave, Usave;
    double FXsave = kTotalGibbsEnergyUnset;   // native's total G, external units (restored below)

    // SHP (warmNative) starts the NATIVE leg warm instead of cold. HOP as
    // built runs it cold at EVERY call, which is right for a single
    // equilibrium and wasteful in a sweep, where the previous point's
    // converged state is already on this node: measured over a 301-point
    // temperature sweep of one TNode (j_Solvus_G_series1, 400-700 C,
    // debug-optima-vs-reaktoro/hop_sweep.cpp), stepping it through native's
    // cold path costs 81743 iterations / 1171 ms and through its warm path
    // 14301 / 241 ms - 5.7x fewer iterations, 4.9x less wall time, same G
    // to every digit. HOP pays the first of those at every step.
    //
    // Same detector as the Optima leg's own warm-start guard above: pm.U[]
    // identically zero (or non-finite) means no solver has ever run on this
    // instance, so there is nothing to warm-start FROM and the .dbr file's
    // stored speciation is NOT a warm start - it is a cold start from data
    // frozen at that file's own state point. Reused verbatim rather than
    // reinvented, and deliberately not "has HOP run before": native AIA/SIA
    // and the Optima path both write pm.U[], and handing any of their
    // converged states to a warm native leg is legitimate.
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
            ipm_logger->warn( "CalculateEquilibriumStateHOP: a warm native leg (SHP) was requested "
                              "but this node carries no previous solution (pm.U[] is all zero) - "
                              "running the native leg cold for this call. Use HOP for a first "
                              "solve, or SHP only on a node that has already solved." );
        }
    }

    try
    {
        // pm.pNP = 1 is native's ordinary warm (SIA) start. NOT -1, which is
        // native's own documented "warm, but first raise every zeroed-off
        // species back to a trace amount" convention (DC_RaiseZeroedOff, the
        // pm.pNP <= -1 branch of InitalizeGEM_IPM_Data). That was tried,
        // because pa_OptimaZeroAbsent commits EXACT zeros into pm.Y[] and
        // hence into xDC, so the state a previous SHP call leaves behind can
        // carry more zeros than native's own answer would - and native's SIA
        // then has to re-insert every marginal species by hand. It is a real
        // effect and it is project-specific, not general: over a 301-point
        // temperature sweep of j_Solvus_G_series1, raising cut SHP from 77068
        // iterations (and one failure) to 38239 (and none), while over an
        // 81-point sweep of j_10TH_G_seawater it did the opposite, 3903 ->
        // 8722. Same non-monotone-across-projects signature as every other
        // knob on this branch, so it is not shipped and no field was added
        // for it; the measurement is here so it is not re-derived.
        pm.pNP = wantWarmNative ? 1 : 0;
        calcTime = CalculateEquilibriumState( fiaN, ipmN );
        // Snapshot the converged state BEFORE the Optima leg can touch it.
        // pm.Y[]/pm.U[] are exactly what the Optima leg's warm start reads
        // and what its objective callback then mutates in place every
        // iteration, so this is the only point at which they still
        // describe native's own answer.
        Ysave.assign( pm.Y, pm.Y + pm.L );
        Usave.assign( pm.U, pm.U + pm.N );
        FXsave = pm.FX;
    }
    catch( TError& werr )
    {
        // COLD FALLBACK for the warm native leg, and the reason SHP is safe
        // to offer at all. Native's own SIA is documented to REFUSE states
        // its cold path returns - on ten projects the cold path returns an
        // answer failing native's own mass-balance test, and SIA is the only
        // path that checks it (plan-v5 sections 29.2 and 60.5). On such a
        // project a warm native leg fails at essentially every step, so
        // without this SHP would be strictly worse than HOP there. With it,
        // the leg is simply retried cold and the mode degrades to exactly
        // HOP; the price is one wasted native attempt, which the same
        // measurements put at tens to low hundreds of iterations.
        if( wantWarmNative )
        {
            ipm_logger->warn( "CalculateEquilibriumStateHOP: the warm native leg failed ({}: {}); "
                              "retrying it cold for this node", werr.title, werr.mess );
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
                ipm_logger->warn( "CalculateEquilibriumStateHOP: the native leg failed cold too "
                                   "({}: {}); falling back to a cold Optima solve for this node",
                                   nerr2.title, nerr2.mess );
                fiaN = ipmN = 0;
                calcTime = 0.;
            }
        }
        else
        {
            // No assemblage to hand over. Degrade to a plain cold Optima solve
            // rather than failing outright - the same "a caller shouldn't have
            // to know" reasoning as the USE_OPTIMA_SOLVER fallback in
            // node.cpp - and say so, because the result is then an AOP result
            // wearing a HOP status.
            nativeOk = false;
            ipm_logger->warn( "CalculateEquilibriumStateHOP: the native leg failed ({}: {}); falling "
                               "back to a cold Optima solve for this node", werr.title, werr.mess );
            fiaN = ipmN = 0;
            calcTime = 0.;
        }
    }

    // Warm-start Optima from native's converged primal AND dual. Warm, not
    // cold, is the whole point: measured native -> warm Optima 7/7 where
    // native -> native SIA fails (plan-v5 section 29.2), and the seeded
    // dual is what makes a correct (x,y) pair cost O(1) iterations rather
    // than several hundred (section 27).
    pm.pNP = nativeOk ? 1 : 0;
    long int fiaO = 0, ipmO = 0;
    // Everything accumulated into calcTime up to here belongs to the native leg
    // (including a cold retry after a failed warm one, which re-assigns rather
    // than adds). Work item 38 - see HopLegSplit in ms_multi.h for why the sum
    // reported below is not enough on its own.
    const double timeNativeLeg = calcTime;

    // Tell the Optima leg it is running on top of native's own assemblage,
    // which is what lets the dimension reduction - otherwise cold-start-only
    // - run in front of it WHEN THE PROJECT ASKS FOR IT (pa_OptimaDimReduce
    // > 0; AUTO does not reach it, because on a project whose warm
    // verification already costs 2 iterations a reduction is a 14x loss -
    // measured, see the gate's own comment and ms_multi.h). Set only when
    // native actually produced an assemblage: with nativeOk == false this
    // degenerates to a plain cold AOP call, which reaches the reduction
    // through pm.pNP == 0 on its own. RAII because
    // CalculateEquilibriumStateOptima() has many exits including thrown
    // Error(...), and this flag must not leak into a later call.
    struct HopLegGuard {
        TMultiBase* m;
        bool prev;
        HopLegGuard( TMultiBase* mm, bool on ) : m(mm), prev(mm->optima_hop_leg)
            { m->optima_hop_leg = on; }
        ~HopLegGuard() { m->optima_hop_leg = prev; }
    } hopLegGuard( this, nativeOk );

    // Work item 38: the per-leg record is written by a DESTRUCTOR, not at the end of the
    // function, so that a call which THROWS is recorded too.
    //
    // WHY IT HAD TO MOVE. The first version wrote it just before `return`, and the
    // !nativeOk branch below deliberately lets a failure propagate - so every failed
    // two-leg call emitted nothing. Measured on the CalcColumn transport loop: 22 of 4200
    // HOP solves errored and produced 4178 records, and a consumer summing them was
    // summing only the calls that returned, with nothing in the record saying so. The
    // count was recoverable (it equals the caller's own err count) but unmarked, which is
    // this branch's recurring defect shape: an absent row and a measured one must not look
    // alike. A destructor also means a future early return cannot silently skip it.
    //
    // WHAT A FAILED CALL CAN AND CANNOT REPORT. fiaO/ipmO are real: the Optima path's own
    // catch(TError&) sets them before re-throwing, and those iterations are work actually
    // spent - the same convention NumIterFIA/NumIterIPM already use for a discarded attempt
    // above. timeOptima is NOT recoverable, because `calcTime +=` never completed; it is
    // left at 0 and `failed=1` is what says the 0 is not a measurement. Do not later
    // "improve" this by timing the throw - see the TIMES note below.
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
            // A DECIDE record so the ITERATION split reaches a freeze the same way every
            // other solver choice does, rather than needing its own harness. Costs
            // nothing when GEMS3K_NATIVE_TRACE_FILE is unset.
            //
            // The two TIMES are deliberately NOT in this record and must not be added to
            // it. Iteration counts are deterministic for a fixed input; wall times are
            // not, so a `# dec` line carrying them would differ between any two freezes
            // of identical code and turn a column whose whole job is to flag a mechanism
            // behaving differently into noise. Same reason `dimreducepass` is excluded by
            // freeze.sh (2026-09-12). The times are on hop_split for a caller that wants
            // them, which is where a timing measurement belongs.
            native_trace_decide( "hop-legsplit warm=%d nativeok=%d failed=%d "
                                 "fian=%ld ipmn=%ld fiao=%ld ipmo=%ld",
                                 warm ? 1 : 0, nok ? 1 : 0, completed ? 0 : 1,
                                 (long)fN, (long)iN, (long)fO, (long)iO );
        }
    } legRecorder{ this, warmNative, nativeOk, fiaN, ipmN, fiaO, ipmO, calcTime, timeNativeLeg };

    if( !nativeOk )
    {
        // Nothing to fall back to - let a failure here propagate exactly
        // as a plain cold AOP call's would (there is no "worse than
        // native" floor to defend when native itself never produced one).
        calcTime += CalculateEquilibriumStateOptima( fiaO, ipmO, false,
                              /* runKinetics */ false );  // the native leg already advanced it
    }
    else
    {
        try
        {
            calcTime += CalculateEquilibriumStateOptima( fiaO, ipmO, false,
                              /* runKinetics */ false );  // the native leg already advanced it
        }
        catch( TError& oerr )
        {
            // Restore native's converged primal and dual, and re-establish
            // every quantity derived from the primal - the Optima leg's own
            // objective callback may have mutated pm.X[]/pm.XF[]/pm.XFA[]/
            // activity coefficients/concentrations in place before it
            // failed, so restoring Y alone would leave those stale. Same
            // recompute sequence pa_OptimaZeroAbsent's own commit path
            // already uses, above.
            for( long int j = 0; j < pm.L; j++ ) pm.Y[j] = pm.X[j] = Ysave[(size_t)j];
            for( long int i = 0; i < pm.N; i++ ) pm.U[i] = Usave[(size_t)i];
            // pm.FX too: the Optima leg re-initialises it to kTotalGibbsEnergyUnset and fails before
            // assigning it, so without this every BAD_GEM_HOP/SHP row published the sentinel through
            // packDataBr() -> DATABR.Gs (seen on GEMS4R's Opa-CI; REVIEW-2026-09-29-gems4r.md s4).
            // Native's own value, saved after its leg rescaled, so already in external units.
            pm.FX = FXsave;
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            CalculateConcentrations( pm.X, pm.XF, pm.XFA );
            // Soft failure, mirroring CalculateEquilibriumStateOptima()'s
            // own BAD_GEM_* convention (testMulti() reads pm.MK below, via
            // TNode::GEM_run()) - the Optima leg could not finish, but
            // native's own answer is real and is what gets reported here,
            // not lost.
            pm.MK = 2;
            std::string buf = std::string("Optima leg failed after a successful native solve (")
                             + oerr.title + oerr.mess
                             + "); reporting native's own converged state instead";
            setErrorMessage( 21, "W21IPM: HOP: ", buf.c_str() );
            // fiaO/ipmO were already set by CalculateEquilibriumStateOptima()'s
            // OWN catch(TError&) block before it re-threw - real iterations
            // spent on the discarded attempt, not zero. Left as they are:
            // only the STATE that attempt produced is discarded, not the
            // reported cost of having made it.
            ipm_logger->warn( "CalculateEquilibriumStateHOP: Optima leg failed after a successful "
                               "native solve ({}: {}) - restoring and reporting native's own answer "
                               "as BAD_GEM_HOP", oerr.title, oerr.mess );
        }
    }

    // Report the TRUE total cost of both legs (or, on the restore path, the
    // cost actually incurred - the discarded Optima attempt's iterations
    // are real work spent, so they are still counted; only its STATE is
    // discarded).
    NumIterFIA = fiaN + fiaO;
    NumIterIPM = ipmN + ipmO;

    legRecorder.completed = true;       // see HopLegRecorder above - the record is written
    return calcTime;                    // by its DESTRUCTOR, on this path and on a throw alike
}

#endif // USE_OPTIMA_SOLVER
