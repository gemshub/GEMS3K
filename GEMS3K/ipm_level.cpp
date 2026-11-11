//-------------------------------------------------------------------
/// \file ipm_level.cpp
/// The leveled solution: a small dense simplex for the linear Gibbs problem, and the column-generation
/// seed that adds mixed phases at their best compositions. Used by the native start and by Optima.
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
//
// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include "gems3k_timing.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include "ms_multi.h"

extern thread_local bool g_tpdStrongSearch;   // defined in ipm_main.cpp

// pa_OptimaLevelStart: the column-generation LP's last round, and the certified start it gives.
thread_local bool g_cgConverged = false;
thread_local long int g_cgPseudo = 0;
// Two used compositions of one phase lie far apart: the phase wants to split, which one copy cannot hold.
thread_local bool g_cgGap = false;
thread_local std::vector<double> g_cgDual;

// Small self-contained two-phase primal simplex (dense tableau, Bland's rule for both
// entering and leaving variables, the standard anti-cycling guarantee): solves
// minimize sum_j c_j x_j  s.t.  A*x = b (A: nRows x nCols, via the `a` accessor), x >= 0.
// Runs once per seed, not per iteration, so a dense O(rows*cols) tableau is adequate.
// Returns false without modifying xOut when the system is infeasible or the iteration
// cap/numerics fail (not distinguished; the caller falls back to a simpler seed), and true
// with xOut sized nCols on success.
// `cost` (optional, null = all ones) is the objective; `yOut` (optional) receives the optimal
// dual of the equality rows, i.e. y with c_j - sum_i y_i a_ij >= 0 for every column. With
// cost = pm.G0 that dual is a first-order estimate of pm.U[] (used by OptimaReducedPreSolve()).
bool TwoPhaseSimplexMinSum( long int nRows, long int nCols,
                             const std::function<double(long int,long int)>& a,
                             const double* b, std::vector<double>& xOut,
                             const double* cost, double* yOut )
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

// pa_OptimaCgSeed (value = TPD tolerance in RT): a cold seed for Optima by column generation.
// The species Gibbs LP (min sum_c cost_c n_c, A n = b, n >= 0; species columns priced at
// G0 + fDQF) gives a vertex and its dual u; every non-ideal condensed solution phase is then
// searched for its minimum tangent-plane distance against u (NativeTpdPhase()); each
// composition y with TPD < -tol becomes a pseudo-compound column a_c = sum_j y_j A[:,j],
// cost_c = TPD + a_c.u = sum_j y_j (G0_j + fDQF_j + ln y_j + lnGam_j(y)); the LP is re-solved
// until no phase prices negative (at most 20 rounds). The species amounts returned are the
// species columns plus each pseudo-compound's amount spread over its composition, so a solid
// solution's mass is not put into a single end-member.
// In plain words: a first guess that already considers mixed phases at sensible compositions.
bool TMultiBase::ColumnGenerationSeed( std::vector<double>& nOut, double tol )
{
    gems3k_timing::Scope timing( gems3k_timing::Level );
    const long int N = pm.N, L = pm.L;
    if( N <= 0 || L <= 0 || !pm.G0 || !pm.U || !pm.A || !pm.B )
        return false;
    struct Col { long int j; long int k; long int jb; std::vector<double> y; std::vector<double> a; double cost; };
    std::vector<Col> cols;
    for( long int j = 0; j < L; j++ )
    {
        Col c; c.j = j; c.k = -1; c.jb = j;
        c.a.resize( (size_t)N );
        for( long int i = 0; i < N; i++ ) c.a[(size_t)i] = pm.A[ i + j*N ];
        c.cost = pm.G0[j] + ( pm.fDQF ? pm.fDQF[j] : 0. );
        cols.push_back( c );
    }
    std::vector<double> Usave( pm.U, pm.U + N );
    std::vector<double> x, y( (size_t)N, 0. ), cost;
    int rounds = 0; long int added = 0, lastAdded = 0;
    bool solved = false;
    for( ; rounds < 20; rounds++ )
    {
        cost.resize( cols.size() );
        for( size_t c = 0; c < cols.size(); c++ ) cost[c] = cols[c].cost;
        auto aFn = [&cols]( long int i, long int c ) { return cols[(size_t)c].a[(size_t)i]; };
        std::fill( y.begin(), y.end(), 0. );
        if( !TwoPhaseSimplexMinSum( N, (long int)cols.size(), aFn, pm.B, x, cost.data(), y.data() ) )
            break;
        solved = true;
        for( long int i = 0; i < N; i++ ) pm.U[i] = y[(size_t)i];
        lastAdded = 0;
        long int jb = 0;
        for( long int k = 0; k < pm.FIs; k++ )
        {
            const long int n = pm.L1[k];
            const char ph = pm.PHC[k];
            if( n > 1 && ph != PH_AQUEL && ph != PH_GASMIX && ph != PH_PLASMA && ph != PH_FLUID && ph != PH_SORPTION
                && ph != PH_POLYEL && ph != PH_ADSORPT && ph != PH_IONEX )
            {
                std::vector<double> yb;
                const double tpd = NativeTpdPhase( k, jb, yb );
                if( tpd < -tol && tpd > -1e299 )
                {
                    bool dup = false;
                    for( const Col& c : cols )
                        if( c.k == k )
                        {
                            double d = 0.;
                            for( long int a = 0; a < n; a++ ) d = std::max( d, std::fabs( c.y[(size_t)a] - yb[(size_t)a] ) );
                            if( d < 1e-4 ) { dup = true; break; }
                        }
                    if( !dup )
                    {
                        Col c; c.j = -1; c.k = k; c.jb = jb; c.y = yb; c.a.assign( (size_t)N, 0. );
                        for( long int a = 0; a < n; a++ )
                            for( long int i = 0; i < N; i++ ) c.a[(size_t)i] += yb[(size_t)a] * pm.A[ i + (jb+a)*N ];
                        double au = 0.;
                        for( long int i = 0; i < N; i++ ) au += c.a[(size_t)i] * y[(size_t)i];
                        c.cost = tpd + au;
                        cols.push_back( c );
                        lastAdded++;
                    }
                }
            }
            jb += n;
        }
        added += lastAdded;
        if( !lastAdded ) break;
    }
    for( long int i = 0; i < N; i++ ) pm.U[i] = Usave[(size_t)i];
    if( !solved || x.size() != cols.size() )
        return false;
    nOut.assign( (size_t)L, 0. );
    for( size_t c = 0; c < cols.size(); c++ )
    {
        if( x[c] <= 0. ) continue;
        if( cols[c].k < 0 ) nOut[(size_t)cols[c].j] += x[c];
        else for( size_t a = 0; a < cols[c].y.size(); a++ ) nOut[(size_t)(cols[c].jb + (long int)a)] += x[c] * cols[c].y[a];
    }
    double bScale = 1., maxResid = 0.;
    for( long int i = 0; i < N; i++ ) bScale = std::max( bScale, std::fabs( pm.B[i] ) );
    for( long int i = 0; i < N; i++ )
    {
        double sum = 0.;
        for( long int j = 0; j < L; j++ ) sum += pm.A[ i + j*N ] * nOut[(size_t)j];
        maxResid = std::max( maxResid, std::fabs( sum - pm.B[i] ) );
    }
    long int nPseudoUsed = 0;
    for( size_t c = 0; c < cols.size(); c++ ) if( cols[c].k >= 0 && x[c] > 0. ) nPseudoUsed++;
    g_cgGap = false;
    for( size_t c = 0; c < cols.size() && !g_cgGap; c++ )
        for( size_t d = c+1; d < cols.size() && cols[c].k >= 0 && x[c] > 0.; d++ )
            if( cols[d].k == cols[c].k && x[d] > 0. )
            {
                double dy = 0.;
                for( size_t a = 0; a < cols[c].y.size(); a++ ) dy = std::max( dy, std::fabs( cols[c].y[a] - cols[d].y[a] ) );
                if( dy > 0.1 ) { g_cgGap = true; break; }
            }
    g_cgConverged = ( lastAdded == 0 );
    g_cgPseudo = nPseudoUsed;
    g_cgDual = y;
    if( g_tpdStrongSearch && base_param()->TpdStrongBasins > 0 )
        native_trace_decide( "tpdscreen basins=%ld", (long)base_param()->TpdStrongBasins );
    native_trace_decide( "cgseed rounds=%d columns_added=%ld pseudo_in_vertex=%ld mb_resid=%.3e tol=%.1e",
                         rounds + 1, (long)added, (long)nPseudoUsed, maxResid, tol );
    if( maxResid > std::max( 1e-6, 1e-8 * bScale ) )
        return false;
    return true;
}

// pa_OptimaFastAcceptCut: checks a finished answer for a missing phase that would form.
double TMultiBase::AbsentPhaseWorstDrive( long int& worstPhase )
{
    worstPhase = -1;
    double worst = 0., tot = 0.;
    for( long int k = 0; k < pm.FI; k++ ) tot += pm.XF[k];
    long int jb = 0;
    for( long int k = 0; k < pm.FI; k++ )
    {
        const long int n = pm.L1[k];
        const char ph = pm.PHC[k];
        bool heldOut = pm.DUL != nullptr;   // every member held out by the input (upper limit 0): it cannot form
        for( long int j = jb; j < jb + n && heldOut; j++ ) heldOut = pm.DUL[j] <= 0.;
        const bool absent = pm.XF[k] <= 1e-9 * tot && !heldOut;
        if( absent && n == 1 )
        {
            double au = 0.;
            for( long int i = 0; i < pm.N; i++ ) au += pm.A[ i + jb*pm.N ] * pm.U[i];
            const double g = pm.G0[jb] + ( pm.fDQF ? pm.fDQF[jb] : 0. ) - au;
            if( g < worst ) { worst = g; worstPhase = k; }
        }
        else if( absent && k < pm.FIs && ph != PH_AQUEL && ph != PH_GASMIX && ph != PH_PLASMA && ph != PH_FLUID
                 && ph != PH_SORPTION && ph != PH_POLYEL && ph != PH_ADSORPT && ph != PH_IONEX )
        {
            std::vector<double> yb;
            g_tpdStrongSearch = true;
            const double t = NativeTpdPhase( k, jb, yb );
            g_tpdStrongSearch = false;
            if( t < worst ) { worst = t; worstPhase = k; }
        }
        jb += n;
    }
    return worst;
}
