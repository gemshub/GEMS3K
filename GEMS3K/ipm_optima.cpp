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
#include <cstdlib>
#include <functional>
#include <sstream>

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
bool TMultiBase::LPGibbsDual( std::vector<double>& yOut )
{
    const long int N = pm.N;
    const long int L = pm.L;
    if( N <= 0 || L <= 0 || pm.G0 == nullptr )
        return false;
    auto aFn = [this,N]( long int i, long int j ) { return pm.A[ i + j*N ]; };
    std::vector<double> nDummy;
    std::vector<double> y( (size_t)N, 0. );
    if( !TwoPhaseSimplexMinSum( N, L, aFn, pm.B, nDummy, pm.G0, y.data() ) )
        return false;
    for( long int i = 0; i < N; i++ )
        if( !std::isfinite( y[(size_t)i] ) )
            return false;
    yOut.swap( y );
    return true;
}

long int TMultiBase::WorstPhaseStabilityViolation( double presenceThreshold,
                                                   const char* exemptSpecies,
                                                   double& violOut, bool& wasAbsentOut )
{
    const BASE_PARAM *pa_p = base_param();
    const long int L = pm.L;

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
        if( kinConstrPh ) continue;

        const double logSI = pm.Falp[k];
        const bool present = ( pm.YF[k] >= presenceThreshold );
        // Same thresholds/logic as native's own PhaseSelect()
        // (ipm_chemical.cpp): a stable phase (logSI > DF) not
        // currently in the assemblage, or an unstable one
        // (logSI < -DFM) that IS, both indicate the wrong assemblage
        // was reached.
        if( !present && logSI > pa_p->DF )
        {
            const double viol = logSI - pa_p->DF;
            if( viol > worstStabilityViol )
            { worstStabilityViol = viol; worstStabilityPhase = k; worstStabilityWasAbsent = true; }
        }
        else if( present && logSI < -pa_p->DFM )
        {
            const double viol = -pa_p->DFM - logSI;
            if( viol > worstStabilityViol )
            { worstStabilityViol = viol; worstStabilityPhase = k; worstStabilityWasAbsent = false; }
        }
    }
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
bool TMultiBase::OptimaReducedPreSolve( long int maxPasses, double dcFloor,
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
    const double kLogBarrierTau = pa_p->LogBarrierTau;
    const double kPhaseHessianFloor = pa_p->PhaseHessianFloor;
    const bool hasAq = HasAqueousPhase();

    // Box bounds, built exactly as the full path builds them - including the
    // "pm.DUL[j] < 1e6 without a `> 0.` guard" convention, so a DUL of exactly
    // zero stays the hard exclusion GEMS3K means it to be.
    std::vector<double> xlo( (size_t)L ), xhi( (size_t)L );
    for( long int j = 0; j < L; j++ )
    {
        xlo[(size_t)j] = std::max( pm.DLL[j], dcFloor );
        xhi[(size_t)j] = ( pm.DUL[j] < 1e6 )
                          ? std::max( pm.DUL[j], dcFloor ) : std::max( pm.SMols, 1.0 ) * 10.;
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
        const double dimTol = pa_p->OptimaDimReduceTol;
        std::vector<double> lpDual;
        if( dimTol > 0. && pm.G0 != nullptr && LPGibbsDual( lpDual ) )
        {
            long int added = 0;
            for( long int j = 0; j < L; j++ )
            {
                if( act[(size_t)j] ) continue;
                if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;
                double z = 0.;
                for( long int i = 0; i < N; i++ ) z += lpDual[(size_t)i] * pm.A[ i + j*N ];
                if( pm.G0[j] - z < dimTol ) { act[(size_t)j] = 1; added++; }
            }
            ipm_logger->info( "OptimaReducedPreSolve: LP-Gibbs pricing admitted {} extra species "
                               "(tol={} RT)", added, dimTol );
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
    options.maxiters = (unsigned)std::max( 2000L, (long int)pa_p->IIM );
    options.convergence.tolerance = pa_p->OptimaTol;

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
                     kFDHessian, kMoleFracHessian, hasAq, &nxToJ, &jToNx, &xlo, &act, &Fbase]
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
                        for( long int j = j0; j < j1; j++ )
                        {
                            const long int sj = jToNx[(size_t)j];
                            if( sj < 0 ) continue;
                            const double Xj = std::max( pm.X[j], dcFloor );
                            if( kMoleFracHessian )
                            {
                                for( long int i = j0; i < j1; i++ )
                                {
                                    const long int si = jToNx[(size_t)i];
                                    if( si >= 0 ) res.fxx(sj,si) = -1.0 / Xf;
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


        Optima::Solver solver;
        solver.setOptions( options );
        Optima::Result result = solver.solve( problem, state );
        iterationsOut += (long int)result.iterations;

        if( !result.succeeded )
        {
            ipm_logger->info( "OptimaReducedPreSolve: pass {} did not converge on {} of {} "
                               "species - discarding the reduced pre-solve", pass, nS, L );
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
        long int readmitted = 0;
        for( long int j = 0; j < L; j++ )
        {
            if( act[(size_t)j] ) continue;
            if( xhi[(size_t)j] <= xlo[(size_t)j] * ( 1. + 1e-9 ) ) continue;  // fixed by its box
            double dual = 0.;
            for( long int i = 0; i < N; i++ )
                dual += pm.U[i] * pm.A[ i + j*N ];
            if( pm.F[j] - dual < 0. ) { act[(size_t)j] = 1; readmitted++; }
        }

        ipm_logger->info( "OptimaReducedPreSolve: pass {} - {} of {} species active, "
                           "{} Optima iterations, {} readmitted",
                           pass, nS, L, result.iterations, readmitted );

        if( readmitted == 0 )
            return true;   // fixed point: the reduced answer satisfies the full KKT conditions
    }

    ipm_logger->warn( "OptimaReducedPreSolve: readmission did not settle within {} passes - "
                       "discarding (an unsettled active set is not a fixed point, and its dual is "
                       "wrong about everything still omitted)", maxPasses );
    return discard();
}

double TMultiBase::CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool reaktoroMode )
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

    pm.t_start = clock();
    pm.t_end = pm.t_start;
    pm.t_elap_sec = 0.0;
    pm.ITF = pm.ITG = 0;
    pm.Ec = pm.MK = pm.PZ = 0;
    setErrorMessage( 0, "", "" );

    if( pa_p->DG > 1e-5 )
    {
        ScFact = SystemTotalMolesIC();
        ScaleSystemToInternal( ScFact );
    }

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
                                            ? pm.DUL[pm.LO] : std::max( pm.SMols, 1.0 ) * 10.;
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
        long int dimReduceIters = 0;
        bool dimReduceDone = false;
        if( !reaktoroMode && R == 0 && pa_p->OptimaDimReduce > 0 )
        {
            long int nActive = 0;
            dimReduceDone = OptimaReducedPreSolve( pa_p->OptimaDimReduce, dcFloor,
                                                   dimReduceIters, nActive );
            if( dimReduceDone )
                ipm_logger->info( "CalculateEquilibriumStateOptima: dimension-reduction pre-solve "
                                   "produced a warm start over {} of {} species in {} iterations",
                                   nActive, L, dimReduceIters );
            // Re-establish the same consistent (Y, X, XF/XFA, activity
            // coefficients) state the seed block above leaves behind. The
            // pre-solve's own objective callback mutated all of those while it
            // ran, and on the discard path pm.Y[] is still the seed - so this
            // is correct whether it succeeded or not, and makes "discarded"
            // mean genuinely discarded.
            TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
            for( long int j = 0; j < L; j++ )
                pm.X[j] = pm.Y[j];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
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
                                 ? std::max( pm.DUL[j], dcFloor ) : std::max( pm.SMols, 1.0 ) * 10.;
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
        const double titrantBound = std::max( pm.SMols, 1.0 ) * 2.0;
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
        const long int kFDDelay = ( !reaktoroMode && pa_p->OptimaFDHessianDelay > 0 )
                                  ? pa_p->OptimaFDHessianDelay : 0;
        const bool kMoleFracHessian = ( pa_p->OptimaMoleFracHessian != 0 );
        problem.f = [this, L, R, dcFloor, &fixedGrad, kLogBarrierTau, kPhaseHessianFloor, kFDHessian, kMoleFracHessian, fdSuppress, hasAq, phLast, phDec, phMax]
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
                        for( long int j = j0; j < j1; j++ )
                        {
                            const double Xj = std::max( pm.X[j], dcFloor );
                            if( kMoleFracHessian )
                            {
                                for( long int i = j0; i < j1; i++ )
                                    res.fxx(j,i) = -1.0 / Xf;
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
                        pm.X[i] = Xi + h;
                        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                        CalculateActivityCoefficients( LINK_UX_MODE );
                        PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );
                        for( long int j = 0; j < L; j++ )
                            res.fxx(j,i) = ( pm.F[j] - Fbase[j] ) / h;
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
            // pa_OptimaEarlyStabilityAt: cap the FIRST attempt only, so the
            // phase-selection repair loop below gets its turn before the primary
            // solve has converged on an assemblage it will then have to correct.
            // Restored to the full budget immediately after that first solve, so
            // every retry keeps the budget it had. See the field's comment in
            // ms_multi.h for the measurement that motivates it.
            if( pa_p->OptimaEarlyStabilityAt > 0 )
                options.maxiters = (unsigned)std::min( (long int)options.maxiters,
                                                       pa_p->OptimaEarlyStabilityAt );
            options.convergence.tolerance = pa_p->OptimaTol;
            // Trust region on per-variable Newton-step growth - a field this
            // branch added to its local Optima checkout. Default 0.0 = off;
            // measured HARMFUL at every nonzero value tried (GEMS3K/CLAUDE.md
            // 2026-08-24), kept only as re-runnable infrastructure.
            options.backtracksearch.max_step_ratio = pa_p->OptimaMaxStepRatio;
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

        // ---- Stall/freeze limit (pa_OptimaStallWindow, default 500; 0 = off) ----
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
            long int run = 0;
            bool stalled = false;
            // Wall-clock budget (pa_OptimaMaxSeconds). Deliberately NOT cleared by
            // reset(): the budget covers the whole call - primary solve plus every
            // retry - not each attempt separately, so the deadline is set once.
            double maxSeconds = 0.;
            std::chrono::steady_clock::time_point started;
            bool timedOut = false;
            void reset() {
                bestErr = bestComp = std::numeric_limits<double>::infinity();
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
        if( stallWatch->window > 0 || stallWatch->maxSeconds > 0. )
            options.convergence.check =
                [stallWatch]( Optima::ConvergenceCheckArgs const& args ) -> bool
                {
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
                    bool improved = false;
                    if( e    < stallWatch->bestErr  ) { stallWatch->bestErr  = e;    improved = true; }
                    if( comp < stallWatch->bestComp ) { stallWatch->bestComp = comp; improved = true; }
                    if( improved ) { stallWatch->run = 0; return false; }
                    if( ++stallWatch->run < stallWatch->window ) return false;
                    stallWatch->stalled = true;
                    return true;   // stop now; folded back to a failure below
                };
        // Applied after each solve() - see point 1 above.
        auto applyStall = [stallWatch]( Optima::Result& r ) {
            if( stallWatch->stalled || stallWatch->timedOut ) r.succeeded = false;
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
            Optima::Result cheapResult = cheapSolver.solve( problem, cheapState );
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

        std::vector<char> extinctFixed( (size_t)std::max(L,1L), 0 );
        stallWatch->reset();
        Optima::Sensitivity sensitivity( dims );
        Optima::Result result = wantSens ? solver.solve( problem, state, sensitivity )
                                         : solver.solve( problem, state );
        if( pa_p->OptimaEarlyStabilityAt > 0 )
        {
            // First attempt is over - give everything downstream the real budget.
            options.maxiters = (unsigned)std::max( 2000L, (long int)pa_p->IIM );
            solver.setOptions( options );
        }
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

            if( !result.succeeded )
                tryToggleRetry();
            // Matches the exact condition validated across the full
            // 25-project sweep: trigger on EITHER `!result.succeeded` OR
            // the trap signature (not the trap signature alone).
            if( !result.succeeded || solventTrapped() )
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
            if( !result.succeeded )
            {
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
            }
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
                    const long int kBad = WorstPhaseStabilityViolation(
                                presenceThreshold,
                                extinctFixed.empty() ? nullptr : extinctFixed.data(),
                                psViol, psAbsent );
                    if( kBad < 0 )
                        break;                       // assemblage is self-consistent
                    if( kBad == aqueousPhaseIdx )
                        break;                       // never remove or reseed the solvent phase here

                    long int jb = 0;
                    for( long int k = 0; k < kBad; k++ )
                        jb += pm.L1[k];
                    const long int je = jb + pm.L1[kBad];
                    if( je > L )
                        break;

                    bool acted = false;
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
                        break;   // nothing this loop is allowed to do - let the final check report it

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
        const bool kktOk = maxKKTResidual <= kktTol;

        double worstStabilityViol = 0.;
        bool worstStabilityWasAbsent = false;
        long int worstStabilityPhase = WorstPhaseStabilityViolation(
                    presenceThreshold,
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
        for( long int k = 0; k < R; k++ )
        {
            conditions[k].titrantAmount = state.x[L+k];
            for( const auto& rc : conditions[k].stoich )
                pm.B[ rc.first ] -= rc.second * conditions[k].titrantAmount;
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

    if( pa_p->DG > 1e-5 )
        RescaleSystemFromInternal( ScFact );

    NumIterFIA = pm.ITF;
    NumIterIPM = pm.ITG;
    pm.t_end = clock();
    pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
    return pm.t_elap_sec;
}

#endif // USE_OPTIMA_SOLVER
