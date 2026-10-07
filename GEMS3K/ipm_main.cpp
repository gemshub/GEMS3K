//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_main.cpp
/// Implementation of parts of the Interior Points Method (IPM) module
/// for convex programming Gibbs energy minimization
/// Uses: JAMA/C++ Linear Algebra Package based on the Template
/// Numerical Toolkit (TNT) - an interface for scientific computing in C++,
/// (c) Roldan Pozo, NIST (USA), http://math.nist.gov/tnt/download.html
//
// Copyright (c) 1992-2012  D.Kulik, S.Dmitrieva, K.Chudnenko
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

#include <cstdlib>
#include <cstdarg>
#include "ms_multi.h"
#include "jama_lu.h"
#include "jama_cholesky.h"
#include "kinetics.h"
#include "v_service.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <chrono>
#include <mutex>
#include <set>
#include <map>
#include <string>

// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> TMultiBase::ipm_logger = spdlog::stdout_color_mt("ipm");

#define uDDtrace false

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Solver event trace - see the declaration in ms_multi.h. The file is opened in a
/// function-local static (thread-safe initialisation); concurrent fprintf() calls are
/// individually locked, so lines from different threads may interleave but are not corrupted.
FILE* native_trace_file()
{
    static FILE* fp = []() -> FILE*
    {
        const char* fn = std::getenv( "GEMS3K_NATIVE_TRACE_FILE" );
        return fn ? fopen( fn, "a" ) : nullptr;
    }();
    return fp;
}

/// Nesting depth of native_trace_quiet(). While > 0, native_trace_run_header() and
/// native_trace_run_result() write nothing (DECIDE and event records are unaffected). Used for
/// inner solves that belong to one GEM_run() call (e.g. pa_ColdRetryNudges), since the RUN/KEY
/// records are attributed to modes by position.
static thread_local int native_trace_quiet_depth = 0;
void native_trace_quiet( bool on )
{
    native_trace_quiet_depth += on ? 1 : -1;
}

/// Per-iteration IPM descent record (GEMS3K_IPM_PROBE); see the call site in
/// InteriorPointsMethod().
FILE* ipm_probe_file()
{
    static FILE* fp = []() -> FILE*
    {
        const char* fn = std::getenv( "GEMS3K_IPM_PROBE" );
        return fn ? fopen( fn, "a" ) : nullptr;
    }();
    return fp;
}

/// One DECIDE record - see the declaration in ms_multi.h.
/// pa_IpmAugmentedKKT: fallbacks to the normal equations in the current InteriorPointsMethod()
/// call (reset at its entry); only the first emits DECIDE "ipmkkt-fallback". thread_local
/// because nodes may be solved concurrently on separate threads.
static thread_local long int s_ipmKktFallbacks = 0;

void native_trace_decide( const char* fmt, ... )
{
    FILE* ntf = native_trace_file();
    if( !ntf )
        return;                       // zero cost when the trace is not enabled
    fputs( "DECIDE ", ntf );
    va_list ap;
    va_start( ap, fmt );
    vfprintf( ntf, fmt, ap );
    va_end( ap );
    fputc( '\n', ntf );
    fflush( ntf );
}

/// Complete run configuration, once per GEM_run() call, for every solver mode. Three lines:
/// RUN (mode, T, P, system shape), BULK (the bulk composition with IC names) and SET (every
/// BASE_PARAM field in force), plus EFF (below). Written into the GEMS3K_NATIVE_TRACE_FILE trace
/// from TNode::GEM_run(), so it covers AOP/SOP/ROP/HOP/SHP as well as native. pm.B[] is read
/// in caller units (after unpackDataBr(), before any internal rescaling).
void native_trace_run_header( const MULTI& pm, const BASE_PARAM* pa, long int mode )
{
    // One call, one certificate: reset the free-dual search result here, since only the
    // header is guaranteed to run before the dispatch that may set it.
    native_cert_dualfree_reset();
    FILE* ntf = native_trace_file();
    if( !ntf || !pa || native_trace_quiet_depth > 0 )
        return;

    const char* mname;
    switch( mode )
    {
      case 1:  mname = "AIA";  break;
      case 5:  mname = "SIA";  break;
      case 10: mname = "AOP";  break;
      case 14: mname = "SOP";  break;
      case 18: mname = "ROP";  break;
      case 22: mname = "HOP";  break;
      case 26: mname = "SHP";  break;
      default: mname = "?";    break;
    }

    fprintf( ntf, "RUN   mode=%s(%ld) pNP=%ld TK=%.4f Pbar=%.6e"
                  " N=%ld L=%ld Ls=%ld FI=%ld FIs=%ld LO=%ld\n",
             mname, (long)mode, (long)pm.pNP, pm.T, pm.P,
             (long)pm.N, (long)pm.L, (long)pm.Ls,
             (long)pm.FI, (long)pm.FIs, (long)pm.LO );

    fprintf( ntf, "BULK " );
    for( long int i = 0; i < pm.N; i++ )
    {
        // pm.SB[] is a fixed-width packed char array, so the name comes back
        // blank-padded; trim it or the line stops parsing on whitespace.
        std::string icn = char_array_to_string( pm.SB[i], MAXICNAME );
        while( !icn.empty() && icn.back() == ' ' ) icn.pop_back();
        fprintf( ntf, " %s=%.10e", icn.c_str(), pm.B[i] );
    }
    fprintf( ntf, "\n" );

    // One line, key=value, every field of BASE_PARAM in declaration order. Keep this list in
    // step with BASE_PARAM (ms_multi.h): a field missing here is invisible in the trace.
    fprintf( ntf, "SET  "
             " pa_PC=%d pa_PD=%d pa_PRD=%d pa_PSM=%d pa_DP=%d pa_DW=%d pa_DT=%d"
             " pa_PLLG=%d pa_PE=%d pa_IIM=%d"
             " pa_DG=%.6e pa_DHB=%.6e pa_DS=%.6e pa_DK=%.6e pa_DF=%.6e pa_DFM=%.6e"
             " pa_DFYw=%.6e pa_DFYaq=%.6e pa_DFYid=%.6e pa_DFYr=%.6e pa_DFYh=%.6e"
             " pa_DFYc=%.6e pa_DFYs=%.6e pa_DB=%.6e pa_AG=%.6e pa_DGC=%.6e"
             " pa_GAR=%.6e pa_GAH=%.6e pa_GAS=%.6e pa_DNS=%.6e pa_XwMin=%.6e"
             " pa_ScMin=%.6e pa_DcMin=%.6e pa_PhMin=%.6e pa_ICmin=%.6e"
             " pa_EPS=%.6e pa_IEPS=%.6e pa_DKIN=%.6e"
             " pa_PSTALL=%d pa_OptimaTol=%.6e pa_LogBarrierTau=%.6e"
             " pa_OptimaMaxStepRatio=%.6e pa_PhaseHessianFloor=%.6e"
             " pa_OptimaStallWindow=%ld pa_OptimaMaxSeconds=%.6e"
             " pa_OptimaFDHessian=%ld pa_OptimaMoleFracHessian=%ld"
             " pa_OptimaPhaseCompaction=%ld pa_OptimaFDHessianDelay=%ld"
             " pa_OptimaDcFloor=%.6e pa_MbClassRule=%.6e pa_MbTrendPhaseDecay=%ld"
             " pa_OptimaEarlyStabilityAt=%ld pa_OptimaDimReduce=%ld"
             " pa_OptimaDimReduceTol=%.6e pa_MbPivotSplit=%ld pa_OptimaZeroAbsent=%ld"
             " pa_OptimaReadmitSeed=%.6e pa_IpmStallWindow=%d pa_MbReproject=%d"
             " pa_DeterminacyWarn=%.6e pa_ColdRetryNudges=%ld"
             " pa_OptimaPreSolveFirstIters=%ld pa_LpDualFillout=%ld"
             " pa_FilloutBudget=%.6e pa_StabTPD=%ld"
             " pa_IpmAugmentedKKT=%ld pa_IpmLoopTweaks=%ld"
             " pa_OptimaLineSearch=%.6e pa_OptimaFDDiagFloor=%ld"
             " pa_OptimaLSStallEscape=%ld pa_OptimaLSWindow=%ld pa_OptimaLSRejectWorse=%ld"
             " pa_OptimaTpdAccept=%.6e pa_OptimaCgSeed=%.6e pa_OptimaColdRetry=%ld pa_OptimaFinish=%ld pa_OptimaAcceptRepair=%ld\n",
             (int)pa->PC, (int)pa->PD, (int)pa->PRD, (int)pa->PSM, (int)pa->DP,
             (int)pa->DW, (int)pa->DT, (int)pa->PLLG, (int)pa->PE, (int)pa->IIM,
             pa->DG, pa->DHB, pa->DS, pa->DK, pa->DF, pa->DFM,
             pa->DFYw, pa->DFYaq, pa->DFYid, pa->DFYr, pa->DFYh,
             pa->DFYc, pa->DFYs, pa->DB, pa->AG, pa->DGC,
             pa->GAR, pa->GAH, pa->GAS, pa->DNS, pa->XwMin,
             pa->ScMin, pa->DcMin, pa->PhMin, pa->ICmin,
             pa->EPS, pa->IEPS, pa->DKIN,
             (int)pa->PSTALL, pa->OptimaTol, pa->LogBarrierTau,
             pa->OptimaMaxStepRatio, pa->PhaseHessianFloor,
             (long)pa->OptimaStallWindow, pa->OptimaMaxSeconds,
             (long)pa->OptimaFDHessian, (long)pa->OptimaMoleFracHessian,
             (long)pa->OptimaPhaseCompaction, (long)pa->OptimaFDHessianDelay,
             pa->OptimaDcFloor, pa->MbClassRule, (long)pa->MbTrendPhaseDecay,
             (long)pa->OptimaEarlyStabilityAt, (long)pa->OptimaDimReduce,
             pa->OptimaDimReduceTol, (long)pa->MbPivotSplit,
             (long)pa->OptimaZeroAbsent, pa->OptimaReadmitSeed,
             (int)pa->IpmStallWindow, (int)pa->MbReproject, pa->DeterminacyWarn,
             (long)pa->ColdRetryNudges, (long)pa->OptimaPreSolveFirstIters,
             (long)pa->LpDualFillout, pa->FilloutBudget, (long)pa->StabTPD,
             (long)pa->IpmAugmentedKKT, (long)pa->IpmLoopTweaks,
             pa->OptimaLineSearch, (long)pa->OptimaFDDiagFloor,
             (long)pa->OptimaLSStallEscape, (long)pa->OptimaLSWindow,
             (long)pa->OptimaLSRejectWorse,
             pa->OptimaTpdAccept, pa->OptimaCgSeed, (long)pa->OptimaColdRetry, (long)pa->OptimaFinish,
             (long)pa->OptimaAcceptRepair );

    // ---- EFF: the effective value of every auto-gated setting -------------
    // The SET line is literal (what the project file carries). For a three-valued field a
    // configured 0 means "decide from the problem", so each such field is resolved here
    // through the same function the solver calls, and printed with its gate inputs: `_cfg` is
    // the configured value, `_eff` what the solver will use. Any new auto-gated field must be
    // added here. The header runs before the solve, so for the pre-solve it reports the pass
    // count and `reached=` what can be decided now (it also needs a cold leg, or a HOP leg
    // under an explicit positive setting; pNP is on the RUN line).
    {
        const long int drCfg = (long)pa->OptimaDimReduce;
        const long int drEff = optima_dimreduce_passes( drCfg, (long)pm.L );
        const char* drGate = ( drCfg > 0 ) ? "EXPLICIT"
                           : ( drCfg < 0 ) ? "OFF" : "AUTO";
        // A warm leg never takes the cold-start path, and ROP (referenceMode) skips the
        // pre-solve entirely.
        const char* drReached = ( drEff <= 0 )      ? "no(off)"
                              : ( mode == 18 )      ? "no(ROP)"
                              : ( pm.pNP != 0 && drCfg <= 0 ) ? "no(warm,AUTO)"
                              : "maybe";
        // pa_OptimaEarlyStabilityAt, AUTO-gated on the presence of a multisite solid-solution
        // model (optima_earlystability_at()); the multisite phase count is printed as the gate input.
        const long int esCfg  = (long)pa->OptimaEarlyStabilityAt;
        const long int esMulti = optima_multisite_phase_count( pm.sMod, pm.FIs );
        // AUTO is leg-dependent, so the leg is printed as a gate input. It is derived from the
        // mode, not from pm.pNP: on HOP/SHP this header runs before the native leg, when pm.pNP
        // is still 0, while the Optima leg is warm. SOP, HOP and SHP have a warm Optima leg;
        // AOP and ROP are cold.
        const bool esWarm = ( mode == 14 || mode == 22 || mode == 26 );
        const long int esEff  = optima_earlystability_at( esCfg, esMulti, esWarm );
        const char* esGate = ( esCfg > 0 ) ? "EXPLICIT-CAP"
                           : ( esCfg < 0 ) ? "EXPLICIT-TREND"
                           : esWarm        ? "AUTO-WARM" : "AUTO-COLD";
        // A native AIA/SIA call never reaches either form.
        const char* esReached = ( esEff == 0 ) ? "no(off)"
                              : ( mode == 1 || mode == 5 ) ? "no(native)"
                              : "maybe";
        fprintf( ntf, "EFF  "
                 " pa_OptimaDimReduce_cfg=%ld pa_OptimaDimReduce_eff=%ld"
                 " dimreduce_gate=%s dimreduce_L=%ld dimreduce_minDC=%ld"
                 " dimreduce_reached=%s"
                 " pa_OptimaEarlyStabilityAt_cfg=%ld pa_OptimaEarlyStabilityAt_eff=%ld"
                 " earlystability_gate=%s earlystability_multisitePh=%ld"
                 " earlystability_autoCap=%ld earlystability_autoWarmCap=%ld"
                 " earlystability_warmLeg=%d earlystability_reached=%s"
                 " optima_netresume_eff=%d optima_netresume_src=%s"
                 " presolve_firstbudget_eff=%ld presolve_passbudget_full=%ld\n",
                 drCfg, drEff, drGate, (long)pm.L,
                 (long)kOptimaDimReduceAutoMinDC, drReached,
                 esCfg, esEff, esGate, esMulti,
                 (long)kOptimaEarlyStabilityAutoCap,
                 (long)kOptimaEarlyStabilityAutoWarmCap, (int)esWarm, esReached,
                 optima_net_resume_mode(),
                 std::getenv( "GEMS3K_OPTIMA_NET_RESUME" ) ? "env" : "default",
                 optima_presolve_pass_budget( (long)pa->OptimaPreSolveFirstIters, (long)pa->IIM, true ),
                 optima_presolve_pass_budget( (long)pa->OptimaPreSolveFirstIters, (long)pa->IIM, false ) );
    }
    fflush( ntf );
}

// ---------------------------------------------------------------------------
// The outcome KEY - the regime the solve landed in
// ---------------------------------------------------------------------------
// Emitted once per GEM_run() call, after the dispatch: the present-phase name list with an
// order-independent 64-bit hash of it, pH, pe, ionic strength, and each present phase's
// amount and molar volume. Lets results be grouped by the assemblage they reached (e.g. to
// look up settings by regime). The fluid-root type of a cubic-EoS phase is not classified;
// the molar volume it would be derived from is emitted instead. No new MULTI or BASE_PARAM
// member.

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Report-only certificate instruments; see the declarations in ms_multi.h.

/// Separator-safe form of a species or phase name for a trace payload: the DECIDE records join
/// names with ',' inside fields separated by ' ', ':' and '=', and real names contain all four.
/// Same substitution as the KEY/CERT writers, kept separate so their output is unchanged.
static std::string trace_safe_name( std::string s )
{
    while( !s.empty() && ( s.back() == ' ' || s.back() == '\t' || s.back() == '\0' ) ) s.pop_back();
    std::string t; bool sp = false;
    for( char c : s )
    {
        if( c == ' ' || c == '\t' ) { sp = true; continue; }
        if( sp && !t.empty() ) t += '_';
        sp = false;
        t += ( c == ',' || c == ':' || c == '=' ) ? '/' : c;
    }
    return t.empty() ? std::string( "-" ) : t;
}

/// Optima's free-dual search result for the CERT record - see the declaration.
static int native_cert_dualfree_resolved = -1;
void native_cert_dualfree_reset() { native_cert_dualfree_resolved = -1; }
void native_cert_dualfree_set( int resolved ) { native_cert_dualfree_resolved = resolved; }
int  native_cert_dualfree_get() { return native_cert_dualfree_resolved; }

/// Smallest eigenvalue of a small dense symmetric block, by cyclic Jacobi sweeps (the same
/// sweep, convergence test and rotation as ipm_optima.cpp's SymEigFloorInPlace(), which returns
/// the floored matrix instead). Returns false and leaves `lmin` (and `lmaxOut`, if given)
/// untouched if the rotations do not settle. `lmaxOut` (optional) returns the largest
/// eigenvalue, for a condition number.
static bool CertSymEigMin( std::vector<double> a, int n, double& lmin, double* lmaxOut = nullptr )
{
    if( n < 1 ) return false;
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
                if( fabs(apq) < 1e-300 ) continue;
                const double app = a[(size_t)p*n+p], aqq = a[(size_t)q*n+q];
                const double theta = ( aqq - app ) / ( 2.*apq );
                const double t = ( theta >= 0. ? 1. : -1. ) /
                                 ( fabs(theta) + sqrt( theta*theta + 1. ) );
                const double c = 1./sqrt( t*t + 1. ), sn = t*c;
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
                }
            }
    }
    if( !converged ) return false;
    lmin = a[0];
    double lmax = a[0];
    for( int i = 1; i < n; i++ )
    {
        lmin = std::min( lmin, a[(size_t)i*n+i] );
        lmax = std::max( lmax, a[(size_t)i*n+i] );
    }
    if( lmaxOut ) *lmaxOut = lmax;
    return true;
}

void TMultiBase::CertPrimalPotentials( std::vector<double>& F ) const
{
    F.assign( (size_t)pm.L, 0. );
    if( !pm.X || !pm.XF || !pm.G0 || !pm.fDQF || !pm.F0 || !pm.DCCW || !pm.L1 )
        return;

    // The phase loop of PrimalChemicalPotentials(), with its skip tests and its
    // NonLogTerm / logXw / logYFk bookkeeping kept as locals so no pm.* scalar moves.
    double NonLogTerm = 0., logXw = 0., logYFk = 0.;
    long int j = 0;
    for( long int k = 0; k < pm.FI; k++ )
    {
        const long int i = j + pm.L1[k];
        const double Yf = pm.XF[k];
        double YFk = 0.;
        if( pm.FIs && k < pm.FIs && pm.XFA )
            YFk = pm.XFA[k];

        if( pm.L1[k] == 1L && Yf < pm.PhMinM ) { j = i; continue; }
        if( Yf <= pm.DSM || ( pm.PHC[k] == PH_AQUEL &&
            ( Yf <= pm.DSM || pm.X[pm.LO] <= pm.XwMinM ) ) ) { j = i; continue; }
        if( Yf >= 1e6 ) { j = i; continue; }   // PrimalChemicalPotentials() throws here; a
                                               // report-only field skips the phase instead

        NonLogTerm = 0.;
        if( ( pm.PHC[k] == PH_AQUEL && YFk >= pm.XwMinM )
         || ( pm.PHC[k] == PH_SORPTION && YFk >= pm.ScMinM )
         || ( pm.PHC[k] == PH_POLYEL && YFk >= pm.ScMinM ) )
        {
            logXw = log( YFk );
            if( k >= pm.FIs || pm.sMod[k][SPHAS_TYP] != SM_AQPITZ )
                NonLogTerm = 1. - YFk / Yf;
        }
        if( pm.L1[k] > 1 )
            logYFk = log( Yf );

        for( ; j < i; j++ )
        {
            if( pm.X[j] < std::min( pm.DcMinM, pm.lowPosNum ) )
                continue;
            const double Gj = pm.G0[j] + pm.fDQF[j] + pm.F0[j];   // rebuilt, not pm.G[j]
            const double lx = log( pm.X[j] );
            switch( pm.DCCW[j] )
            {
              case DC_SINGLE:       F[(size_t)j] = Gj;                                       break;
              case DC_ASYM_SPECIES: F[(size_t)j] = Gj + lx - logXw + NonLogTerm;             break;
              case DC_ASYM_CARRIER: F[(size_t)j] = Gj + lx - logYFk + NonLogTerm + 1.0
                                                   - 1.0/( 1.0 - NonLogTerm );               break;
              case DC_SYMMETRIC:    F[(size_t)j] = Gj + lx - logYFk;                         break;
              default:              F[(size_t)j] = 0.;                                       break;
            }
        }
        j = i;
    }
}

/// Is species j scorable for the reduced-gradient tests - present at the answer, and with a
/// potential CertPrimalPotentials() actually filled in. A species left at 0 there is either
/// below pm.DcMinM or in a skipped phase, and its F carries no information.
static inline bool cert_species_scorable( const MULTI& pm, const std::vector<double>& F, long int j )
{
    return pm.X[j] > pm.DcMinM && F[(size_t)j] != 0.;
}

/// The numerical floor on a species amount, as CalculateEquilibriumStateOptima() builds its
/// box lower bounds: pa_OptimaDcFloor when set, else pa_DHB.
static inline double cert_dc_floor( const MULTI& pm, const BASE_PARAM* pa )
{
    const double f = pa ? ( pa->OptimaDcFloor > 0. ? pa->OptimaDcFloor : std::max( pa->DHB, 1e-300 ) )
                        : 1e-300;
    return std::max( f, pm.DcMinM );
}

/// The box a species sits in, as the Optima KKT check classifies it: 0 interior, -1 at the
/// lower bound, +1 at the upper bound, 2 degenerate (DUL <= DLL, a kinetically fixed species -
/// an equality constraint, no sign test). The lower-bound tolerance is the numerical floor
/// max(DLL, dcFloor), not pm.DcMinM, so a floor-pinned species is not read as interior.
static inline int cert_box_state( const MULTI& pm, long int j, double dcFloor )
{
    const double lo = pm.DLL ? pm.DLL[j] : 0.;
    const double hi = pm.DUL ? pm.DUL[j] : 1e6;
    const double tol = std::max( dcFloor, fabs(lo)*1e-6 );
    if( hi <= lo + tol )              return 2;
    if( pm.X[j] <= lo + tol )         return -1;
    if( pm.X[j] >= hi * ( 1.-1e-6 ) ) return 1;
    return 0;
}

double TMultiBase::CertKktMax( const std::vector<double>& F, long int& worstJ ) const
{
    worstJ = -1;
    if( !pm.X || !pm.A || !pm.U || pm.N <= 0 || (long int)F.size() != pm.L )
        return -1.;

    // pm.N, not pm.NR: NR drops the last IC row while the aqueous phase is absent, and the
    // Optima path's KKT check sums over all N rows, so one index set serves every mode.
    const double dcFloor = cert_dc_floor( pm, base_param() );
    double worst = -1.;
    for( long int j = 0; j < pm.L; j++ )
    {
        if( !cert_species_scorable( pm, F, j ) )
            continue;
        double s = F[(size_t)j];
        for( long int i = 0; i < pm.N; i++ )
            s -= pm.U[i] * pm.A[ i + j*pm.N ];
        double resid;
        switch( cert_box_state( pm, j, dcFloor ) )
        {
          case 2:  resid = 0.;                  break;   // degenerate box: unrestricted sign
          case -1: resid = std::max( -s, 0. );  break;   // at lower bound: s should be >= 0
          case 1:  resid = std::max(  s, 0. );  break;   // at upper bound: s should be <= 0
          default: resid = fabs( s );           break;   // interior: s should vanish
        }
        if( resid > worst ) { worst = resid; worstJ = j; }
    }
    return worst;
}

long int TMultiBase::CertDualFreeDirs( const std::vector<double>& F, long int& rank, long int& nIC ) const
{
    rank = 0;
    nIC = pm.N;
    if( !pm.X || !pm.A || pm.N <= 0 || (long int)F.size() != pm.L )
        return -1;

    const long int N = pm.N;
    const double dcFloor = cert_dc_floor( pm, base_param() );
    std::vector<double> Q, r( (size_t)N );
    auto orthogonalize = [&]( const std::vector<double>& basis, long int nb ) {
        for( long int c = 0; c < nb; c++ )
        {
            double d = 0.;
            for( long int i = 0; i < N; i++ ) d += basis[(size_t)(c*N + i)] * r[(size_t)i];
            for( long int i = 0; i < N; i++ ) r[(size_t)i] -= d * basis[(size_t)(c*N + i)];
        }
        double nrm = 0.;
        for( long int i = 0; i < N; i++ ) nrm += r[(size_t)i] * r[(size_t)i];
        return sqrt( nrm );
    };
    for( long int j = 0; j < pm.L && rank < N; j++ )
    {
        if( !cert_species_scorable( pm, F, j ) )
            continue;
        if( cert_box_state( pm, j, dcFloor ) != 0 )   // only the INTERIOR species fix the dual
            continue;
        double nrm0 = 0.;
        for( long int i = 0; i < N; i++ ) { r[(size_t)i] = pm.A[ i + j*N ]; nrm0 += r[(size_t)i]*r[(size_t)i]; }
        nrm0 = sqrt( nrm0 );
        if( !( nrm0 > 0. ) ) continue;
        const double nrm = orthogonalize( Q, rank );
        if( nrm < 1e-8 * nrm0 ) continue;
        Q.resize( (size_t)((rank+1)*N) );
        for( long int i = 0; i < N; i++ ) Q[(size_t)(rank*N + i)] = r[(size_t)i] / nrm;
        rank++;
    }
    return N - rank;
}

/// cond(M) = lambda_max/lambda_min of a small dense symmetric M, via CertSymEigMin(), set to
/// 1e300 when the sweep does not converge or either eigenvalue is non-positive (on a PSD Gram
/// that is rounding on a near-singular matrix: "no bound").
static double CertCondFromGram( const std::vector<double>& M, long int n )
{
    if( n <= 0 ) return 1e300;
    double lmin = 0., lmax = 0.;
    if( !CertSymEigMin( M, (int)n, lmin, &lmax ) ) return 1e300;
    if( !( lmin > 0. ) || !( lmax > 0. ) ) return 1e300;
    return std::min( lmax / lmin, 1e300 );
}

/// M rescaled by its own symmetric Jacobi diagonal, d_i = sqrt(M_ii) where M_ii > 0, else 1
/// (a structurally empty row is left undivided).
static std::vector<double> CertJacobiScale( const std::vector<double>& M, long int n )
{
    std::vector<double> d( (size_t)n, 1. );
    for( long int i = 0; i < n; i++ )
        if( M[(size_t)(i*n+i)] > 0. ) d[(size_t)i] = sqrt( M[(size_t)(i*n+i)] );
    std::vector<double> Ms( M.size() );
    for( long int i = 0; i < n; i++ )
        for( long int k = 0; k < n; k++ )
            Ms[(size_t)(i*n+k)] = M[(size_t)(i*n+k)] / ( d[(size_t)i] * d[(size_t)k] );
    return Ms;
}

void TMultiBase::CertRank( CertRankReport& r ) const
{
    r = CertRankReport();
    if( !pm.X || !pm.A || pm.N <= 0 )
        return;

    const long int N = pm.N;

    // Present species: pm.X[j] > pm.DcMinM (as in cert_species_scorable(), without its F test).
    std::vector<long int> presentIdx;
    presentIdx.reserve( (size_t)pm.L );
    for( long int j = 0; j < pm.L; j++ )
        if( pm.X[j] > pm.DcMinM ) presentIdx.push_back( j );
    const long int pres = (long int)presentIdx.size();

    r.of = N;
    r.pres = pres;
    if( pres <= 0 ) return;

    // Row scale for `rank` and `sv_ratio`: each IC row divided by its own max |entry| over the
    // present columns, so a trace IC's row does not read as dependent merely because it is
    // small. A row with no present column keeps scale 1.
    std::vector<double> rowScale( (size_t)N, 1. );
    for( long int i = 0; i < N; i++ )
    {
        double m = 0.;
        for( long int jc = 0; jc < pres; jc++ )
            m = std::max( m, fabs( pm.A[ i + presentIdx[(size_t)jc]*N ] ) );
        if( m > 0. ) rowScale[(size_t)i] = m;
    }

    // rank: modified Gram-Schmidt over the row-scaled present columns, with the same acceptance
    // test as CertDualFreeDirs(). That routine ranks the interior species, unscaled, to ask
    // whether the dual is fixed; this ranks every present species, row-scaled.
    {
        std::vector<double> Q, col( (size_t)N ), resid( (size_t)N );
        long int rank = 0;
        for( long int jc = 0; jc < pres && rank < N; jc++ )
        {
            const long int j = presentIdx[(size_t)jc];
            double nrm0 = 0.;
            for( long int i = 0; i < N; i++ )
            {
                col[(size_t)i] = pm.A[ i + j*N ] / rowScale[(size_t)i];
                nrm0 += col[(size_t)i]*col[(size_t)i];
            }
            nrm0 = sqrt( nrm0 );
            if( !( nrm0 > 0. ) ) continue;
            resid = col;
            for( long int c = 0; c < rank; c++ )
            {
                double d = 0.;
                for( long int i = 0; i < N; i++ ) d += Q[(size_t)(c*N+i)] * resid[(size_t)i];
                for( long int i = 0; i < N; i++ ) resid[(size_t)i] -= d * Q[(size_t)(c*N+i)];
            }
            double nrm = 0.;
            for( long int i = 0; i < N; i++ ) nrm += resid[(size_t)i]*resid[(size_t)i];
            nrm = sqrt( nrm );
            if( nrm < 1e-8 * nrm0 ) continue;
            Q.resize( (size_t)((rank+1)*N) );
            for( long int i = 0; i < N; i++ ) Q[(size_t)(rank*N+i)] = resid[(size_t)i] / nrm;
            rank++;
        }
        r.rank = rank;
    }

    // sv_ratio / sv_ratio_raw: sigma_min/sigma_max of A_present, scaled and raw, from the
    // eigenvalues of the N x N Gram A_present A_present^T (squared singular values).
    auto buildGram = [&]( bool scaled )
    {
        std::vector<double> M( (size_t)(N*N), 0. );
        for( long int jc = 0; jc < pres; jc++ )
        {
            const long int j = presentIdx[(size_t)jc];
            for( long int i = 0; i < N; i++ )
            {
                const double ai = pm.A[ i + j*N ] / ( scaled ? rowScale[(size_t)i] : 1. );
                if( ai == 0. ) continue;
                for( long int k = i; k < N; k++ )
                    M[(size_t)(i*N+k)] += ai * ( pm.A[ k + j*N ] / ( scaled ? rowScale[(size_t)k] : 1. ) );
            }
        }
        for( long int i = 0; i < N; i++ )
            for( long int k = 0; k < i; k++ )
                M[(size_t)(i*N+k)] = M[(size_t)(k*N+i)];
        return M;
    };
    {
        double lmin, lmax;
        if( CertSymEigMin( buildGram( true ), (int)N, lmin, &lmax ) && lmin >= 0. && lmax > 0. )
            r.sv_ratio = sqrt( lmin / lmax );
        if( CertSymEigMin( buildGram( false ), (int)N, lmin, &lmax ) && lmin >= 0. && lmax > 0. )
            r.sv_ratio_raw = sqrt( lmin / lmax );
    }

    // chg_res / chg_span: is the charge row, restricted to the present columns, already a linear
    // combination of the element rows over the same columns? Modified Gram-Schmidt over vectors
    // indexed by the present columns, same 1e-8 relative test, first charge row Z = N - E. Each
    // row is divided by the same rowScale[i] as above, so a trace IC's row is not swamped.
    const long int Z = N - pm.E;
    if( pm.E > 0 && Z > 0 )
    {
        std::vector<double> Q, col( (size_t)pres ), resid( (size_t)pres );
        long int erank = 0;
        for( long int i = 0; i < Z; i++ )
        {
            double nrm0 = 0.;
            for( long int jc = 0; jc < pres; jc++ )
            {
                col[(size_t)jc] = pm.A[ i + presentIdx[(size_t)jc]*N ] / rowScale[(size_t)i];
                nrm0 += col[(size_t)jc]*col[(size_t)jc];
            }
            nrm0 = sqrt( nrm0 );
            if( !( nrm0 > 0. ) ) continue;
            resid = col;
            for( long int c = 0; c < erank; c++ )
            {
                double d = 0.;
                for( long int jc = 0; jc < pres; jc++ ) d += Q[(size_t)(c*pres+jc)] * resid[(size_t)jc];
                for( long int jc = 0; jc < pres; jc++ ) resid[(size_t)jc] -= d * Q[(size_t)(c*pres+jc)];
            }
            double nrm = 0.;
            for( long int jc = 0; jc < pres; jc++ ) nrm += resid[(size_t)jc]*resid[(size_t)jc];
            nrm = sqrt( nrm );
            if( nrm < 1e-8 * nrm0 ) continue;
            Q.resize( (size_t)((erank+1)*pres) );
            for( long int jc = 0; jc < pres; jc++ ) Q[(size_t)(erank*pres+jc)] = resid[(size_t)jc] / nrm;
            erank++;
        }
        double cnrm0 = 0.;
        for( long int jc = 0; jc < pres; jc++ )
        {
            col[(size_t)jc] = pm.A[ Z + presentIdx[(size_t)jc]*N ] / rowScale[(size_t)Z];
            cnrm0 += col[(size_t)jc]*col[(size_t)jc];
        }
        cnrm0 = sqrt( cnrm0 );
        if( cnrm0 > 0. )
        {
            resid = col;
            for( long int c = 0; c < erank; c++ )
            {
                double d = 0.;
                for( long int jc = 0; jc < pres; jc++ ) d += Q[(size_t)(c*pres+jc)] * resid[(size_t)jc];
                for( long int jc = 0; jc < pres; jc++ ) resid[(size_t)jc] -= d * Q[(size_t)(c*pres+jc)];
            }
            double rnrm = 0.;
            for( long int jc = 0; jc < pres; jc++ ) rnrm += resid[(size_t)jc]*resid[(size_t)jc];
            r.chg_res = sqrt( rnrm ) / cnrm0;
        }
        else
            r.chg_res = 0.;   // charge row is identically zero over the present columns: trivially in span
        r.chg_span = ( r.chg_res < 1e-10 ) ? 1 : 0;
    }

    // cond_ipm / cond_mbr: condition number of A_p diag(w) A_p^T for the two weights the solver
    // stages apply - IPM's (WeightMultipliers(false)) and MBR's (WeightMultipliers(true)). A local
    // reconstruction of that function's arithmetic at pm.X[j], not a call to it and not a read of
    // pm.W[] (live solver scratch, which a report-only path must not write). The 1.34e120 clamp
    // and the "BOTH_LIM takes the min after squaring" rule for the MBR shape are the same.
    auto weightAt = [&]( long int j, bool square ) -> double
    {
        const char rlc = pm.RLC ? pm.RLC[j] : (char)NO_LIM;
        const double lo = pm.DLL ? pm.DLL[j] : 0.;
        const double hi = pm.DUL ? pm.DUL[j] : 0.;
        auto clampSq = []( double w )
        {
            if( fabs(w) > 1.34e120 ) w = signbit(w) ? -1.34e120 : 1.34e120;
            return w*w;
        };
        switch( rlc )
        {
          case UPPER_LIM:
          {
              const double w1 = hi - pm.X[j];
              return square ? clampSq( w1 ) : std::max( w1, 0. );
          }
          case BOTH_LIM:
          {
              const double w1 = pm.X[j] - lo;
              const double w2 = hi - pm.X[j];
              if( square ) return std::min( clampSq( w1 ), clampSq( w2 ) );
              const double w = std::min( w1, w2 );
              return w < 0. ? 0. : w;
          }
          case NO_LIM:
          case LOWER_LIM:
          default:
          {
              const double w1 = pm.X[j] - lo;
              return square ? clampSq( w1 ) : std::max( w1, 0. );
          }
        }
    };
    auto buildWeightedGram = [&]( bool square )
    {
        std::vector<double> M( (size_t)(N*N), 0. );
        for( long int jc = 0; jc < pres; jc++ )
        {
            const long int j = presentIdx[(size_t)jc];
            const double w = weightAt( j, square );
            if( w == 0. ) continue;
            for( long int i = 0; i < N; i++ )
            {
                const double ai = pm.A[ i + j*N ];
                if( ai == 0. ) continue;
                const double wai = w*ai;
                for( long int k = i; k < N; k++ )
                    M[(size_t)(i*N+k)] += wai * pm.A[ k + j*N ];
            }
        }
        for( long int i = 0; i < N; i++ )
            for( long int k = 0; k < i; k++ )
                M[(size_t)(i*N+k)] = M[(size_t)(k*N+i)];
        return M;
    };
    {
        const std::vector<double> Mi = buildWeightedGram( false );
        r.cond_ipm     = CertCondFromGram( Mi, N );
        r.cond_ipm_jac = CertCondFromGram( CertJacobiScale( Mi, N ), N );
        const std::vector<double> Mm = buildWeightedGram( true );
        r.cond_mbr     = CertCondFromGram( Mm, N );
        r.cond_mbr_jac = CertCondFromGram( CertJacobiScale( Mm, N ), N );
    }
}

double TMultiBase::CertCurvMin( long int& worstPhase )
{
    worstPhase = -1;
    const double kNone = 1e300;
    if( !pm.X || !pm.XF || !pm.L1 || pm.FIs <= 0 )
        return kNone;

    // Restore by copy, not by recomputation: CalculateActivityCoefficients(LINK_UX_MODE) is not
    // idempotent (it accumulates pm.lnGmo and blends pm.F0 through the smoothing factor
    // pm.FitVar[3]), so re-running it at the original composition does not reproduce the state.
    // The refresh at the original X is still run to re-seed the TSolMod objects' own
    // composition; the arrays are then overwritten by the saved copies, so a later warm call
    // reads a bit-identical state.
    struct CertSave { double* p; std::vector<double> v; };
    std::vector<CertSave> saved;
    auto keep = [&saved]( double* p, long int n ) {
        if( p && n > 0 ) saved.push_back( { p, std::vector<double>( p, p + n ) } ); };
    keep( pm.X, pm.L );      keep( pm.XF, pm.FI );     keep( pm.XFA, pm.FIs );
    keep( pm.G, pm.L );      keep( pm.lnGam, pm.L );   keep( pm.lnGmo, pm.L );
    keep( pm.Gamma, pm.L );  keep( pm.F0, pm.L );      keep( pm.fDQF, pm.L );
    keep( pm.Wx, pm.L );     keep( pm.FitVar, 5 );
    const std::vector<double> Xsave( pm.X, pm.X + pm.L );
    auto restoreBase = [&saved]() {
        for( const CertSave& c : saved ) std::copy( c.v.begin(), c.v.end(), c.p ); };

    // Every FD column starts from the same base state, and F at the base is computed once:
    // CalculateActivityCoefficients(LINK_UX_MODE) is history-dependent, so the base is restored
    // between columns.
    std::vector<double> Fbase;
    CertPrimalPotentials( Fbase );

    double lmin = kNone;
    bool perturbed = false;
    try
    {
        long int p0 = 0;
        for( long int k = 0; k < pm.FIs; k++ )
        {
            const long int p1 = p0 + pm.L1[k];
            const long int nEnd = p1 - p0;
            const bool isAq = ( pm.LO >= p0 && pm.LO < p1 );
            // pm.XF[k] > pm.DSM is the presence test the KEY record uses, so curv_min is scored
            // over the phases KEY calls present.
            if( isAq || nEnd <= 1 || p1 > pm.L || pm.XF[k] <= pm.DSM ) { p0 = p1; continue; }

            // Present end-members only, on the same test as the pa_PhaseHessianFloor site.
            const double dcFloor = cert_dc_floor( pm, base_param() );
            double phTot = 0.;
            for( long int j = p0; j < p1; j++ ) phTot += std::max( pm.X[j], 0. );
            std::vector<long int> pres;
            for( long int j = p0; j < p1; j++ )
                if( pm.X[j] > std::max( dcFloor * 1e3, phTot * 1e-6 ) )
                    pres.push_back( j );
            const int nP = (int)pres.size();
            if( nP <= 1 ) { p0 = p1; continue; }

            std::vector<double> Fpert, fxx( (size_t)nP*nP, 0. );
            for( int c = 0; c < nP; c++ )
            {
                const long int jc = pres[(size_t)c];
                const double Xi = Xsave[(size_t)jc];
                const double h = fabs(Xi) * 1e-7;
                if( !( h > 0. ) ) continue;
                pm.X[jc] = Xi + h;
                perturbed = true;
                TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
                CalculateActivityCoefficients( LINK_UX_MODE );
                CertPrimalPotentials( Fpert );
                for( int rr = 0; rr < nP; rr++ )
                    fxx[(size_t)rr*nP + c] = ( Fpert[(size_t)pres[(size_t)rr]]
                                             - Fbase[(size_t)pres[(size_t)rr]] ) / h;
                restoreBase();          // back to the exact base state before the next column
            }
            std::vector<double> blk( (size_t)nP*nP );
            for( int a = 0; a < nP; a++ )
                for( int b = 0; b < nP; b++ )
                    blk[(size_t)a*nP+b] = 0.5 * ( fxx[(size_t)a*nP+b] + fxx[(size_t)b*nP+a] );
            double lk = 0.;
            if( CertSymEigMin( blk, nP, lk ) && lk < lmin )
            {
                lmin = lk;
                worstPhase = k;
            }
            p0 = p1;
        }
    }
    catch( ... )
    {
        // A solution model that throws on a perturbed composition must not turn a completed
        // solve into a failed call (this runs after packDataBr()): report nothing and restore.
        lmin = kNone;
        worstPhase = -1;
    }

    if( perturbed )
    {
        // pm.X is at the base state here; the refresh re-seeds the TSolMod objects, whose own
        // internal composition no array copy can restore.
        restoreBase();
        try
        {
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );   // re-seed the TSolMod objects
        }
        catch( ... ) { /* the copies below are what make the state exact */ }
        restoreBase();
    }
    return lmin;
}

double TMultiBase::CertStabTPD( long int& worstPhase, long int& nScanned, long int& nDisagree )
{
    worstPhase = -1; nScanned = 0; nDisagree = 0;
    const double kNone = 1e300, kTol = 1e-6;
    if( !pm.X || !pm.XF || !pm.L1 || !pm.U || !pm.A || !pm.lnGam || pm.FIs <= 0 )
        return kNone;
    const BASE_PARAM* pa = base_param();
    if( !pa || pa->StabTPD < 1 )
        return kNone;

    // The same save set and restore-by-copy as CertCurvMin().
    struct CertSave { double* p; std::vector<double> v; };
    std::vector<CertSave> saved;
    auto keep = [&saved]( double* p, long int n ) {
        if( p && n > 0 ) saved.push_back( { p, std::vector<double>( p, p + n ) } ); };
    keep( pm.X, pm.L );      keep( pm.XF, pm.FI );     keep( pm.XFA, pm.FIs );
    keep( pm.G, pm.L );      keep( pm.lnGam, pm.L );   keep( pm.lnGmo, pm.L );
    keep( pm.Gamma, pm.L );  keep( pm.F0, pm.L );      keep( pm.fDQF, pm.L );
    keep( pm.Wx, pm.L );     keep( pm.FitVar, 5 );
    auto restoreBase = [&saved]() {
        for( const CertSave& c : saved ) std::copy( c.v.begin(), c.v.end(), c.p ); };

    const double dcFloor = cert_dc_floor( pm, pa );
    const long int N = pm.N;
    double sumXF = 0.;
    for( long int k = 0; k < pm.FI; k++ ) sumXF += std::max( pm.XF[k], 0. );
    double worst = kNone;
    bool perturbed = false;
    long int p0 = 0;
    for( long int k = 0; k < pm.FIs; k++ )
    {
        const long int p1 = p0 + pm.L1[k];
        const long int n = pm.L1[k];
        const char ph = pm.PHC[k];
        if( n <= 1 || p1 > pm.L || ph == PH_AQUEL || ph == PH_SORPTION || ph == PH_POLYEL
            || ph == PH_ADSORPT || ph == PH_IONEX ) { p0 = p1; continue; }
        double xmax = 0.;
        for( long int j = p0; j < p1; j++ ) xmax = std::max( xmax, pm.X[j] );
        // A trace phase counts as absent too: a phase holding under 1e-6 of the system's total
        // phase amount carries no material share.
        const bool absent = pm.XF[k] <= pm.DSM || xmax <= dcFloor * 1e3 || pm.XF[k] < 1e-6 * sumXF;
        if( !absent ) { p0 = p1; continue; }

        std::vector<double> c0( (size_t)n );
        for( long int j = p0; j < p1; j++ )
        {
            double au = 0.;
            for( long int i = 0; i < N; i++ ) au += pm.A[ i + j*N ] * pm.U[i];
            c0[(size_t)(j-p0)] = au - ( pm.G0[j] + pm.fDQF[j] );
        }
        // TPD(y) = sum_j y_j (ln y_j + lnGam_j(y) - c_j), recorded at every composition evaluated
        // (starts and substitution iterates alike). A negative value anywhere certifies that the
        // phase lowers G; no stationarity is needed.
        double tpdMin = 1e300;
        // lnTM(y) and the substitution step W(y)/sum W(y); the phase is put at 1 mol so every
        // presence gate inside CalculateActivityCoefficients() passes.
        auto evalY = [&]( const std::vector<double>& y, std::vector<double>& ynew ) -> double {
            restoreBase();
            perturbed = true;
            for( long int j = p0; j < p1; j++ ) pm.X[j] = std::max( y[(size_t)(j-p0)], 1e-300 );
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            double mx = -1e300;
            std::vector<double> lw( (size_t)n );
            for( long int a = 0; a < n; a++ )
            {
                lw[(size_t)a] = c0[(size_t)a] - pm.lnGam[p0+a];
                mx = std::max( mx, lw[(size_t)a] );
            }
            double s = 0.;
            for( double v : lw ) s += exp( v - mx );
            ynew.resize( (size_t)n );
            for( long int a = 0; a < n; a++ ) ynew[(size_t)a] = exp( lw[(size_t)a] - mx ) / s;
            double tpd = 0.;
            for( long int a = 0; a < n; a++ )
                if( y[(size_t)a] > 0. ) tpd += y[(size_t)a] * ( log( y[(size_t)a] ) - lw[(size_t)a] );
            tpdMin = std::min( tpdMin, tpd );
            return mx + log( s );
        };
        std::vector<double> y0( (size_t)n, 1.0 / n ), tmp;
        double s0 = 0.;
        for( long int j = p0; j < p1; j++ ) s0 += std::max( pm.X[j], 0. );
        if( s0 > 0. )
            for( long int j = p0; j < p1; j++ ) y0[(size_t)(j-p0)] = std::max( pm.X[j], 0. ) / s0;
        try
        {
            (void)evalY( y0, tmp );
            std::vector<std::vector<double>> starts;
            for( long int a = 0; a < n; a++ )
            {
                std::vector<double> v( (size_t)n, 1e-6 / std::max( 1L, n - 1 ) );
                v[(size_t)a] = 1. - 1e-6;
                starts.push_back( v );
            }
            {
                const double mx = *std::max_element( c0.begin(), c0.end() );
                double s = 0.;
                std::vector<double> v( (size_t)n );
                for( long int a = 0; a < n; a++ ) { v[(size_t)a] = exp( c0[(size_t)a] - mx ); s += v[(size_t)a]; }
                for( double& x : v ) x /= s;
                starts.push_back( v );
            }
            starts.push_back( y0 );
            // Denser starts: the centroid, a 19-point grid for a binary, edge midpoints (capped at
            // 64 starts) otherwise.
            starts.push_back( std::vector<double>( (size_t)n, 1.0 / n ) );
            if( n == 2 )
                for( int g = 1; g < 20; g++ ) starts.push_back( { g / 20., 1. - g / 20. } );
            else
                for( long int a = 0; a < n && starts.size() < 64; a++ )
                    for( long int b2 = a + 1; b2 < n && starts.size() < 64; b2++ )
                    {
                        std::vector<double> v( (size_t)n, 1e-6 ); v[(size_t)a] = 0.5; v[(size_t)b2] = 0.5;
                        double t = 0.; for( double x : v ) t += x; for( double& x : v ) x /= t;
                        starts.push_back( v );
                    }
            for( auto y : starts )
            {
                double lt = 0.;
                bool conv = false;
                for( int it = 0; it < 300; it++ )
                {
                    lt = evalY( y, tmp );
                    double d = 0.;
                    for( long int a = 0; a < n; a++ ) d = std::max( d, fabs( tmp[(size_t)a] - y[(size_t)a] ) );
                    y = tmp;
                    if( d < 1e-11 ) { conv = true; break; }
                }
                (void)conv; (void)lt;
            }
            nScanned++;
            if( tpdMin < 1e300 )
            {
                if( tpdMin < worst ) { worst = tpdMin; worstPhase = k; }
                // Disagreement against the solver's own single-point index (pm.Falp, log10; <= 0
                // reads "stable"). StabilityIndexes() writes the sentinel -1 for a zero-amount
                // non-ideal phase, which also reads "stable".
                const double falp = pm.Falp ? pm.Falp[k] : 0.;
                if( falp <= kTol && tpdMin < -kTol ) nDisagree++;
            }
        }
        catch( ... ) { /* a model that throws at a trial composition is skipped, never fatal */ }
        restoreBase();
        p0 = p1;
    }
    if( perturbed )
    {
        try
        {
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );   // re-seed the TSolMod objects
        }
        catch( ... ) {}
        restoreBase();
    }
    return worst;
}

// The composition search of CertStabTPD() for one absent non-ideal phase k (species p0 ..
// p0+L1[k]-1): returns min TPD (RT per mole of phase; < 0: the phase lowers G against the
// current dual pm.U) and the composition where it was found. Used by the TPD acceptance of the
// Optima path, the column-generation seed, and by PhaseSelectionSpeciationCleanup() when
// GEMS3K_NATIVE_TPD_INSERT is set (in place of the single-point index Falp). Same
// save/restore-by-copy as CertStabTPD(); returns 1e300 if the model throws or nothing was
// evaluated.
// In plain words: searches for the composition at which a missing mixed phase would be most
// stable, and says whether it should form.
double TMultiBase::NativeTpdPhase( long int k, long int p0, std::vector<double>& ybest )
{
    const long int n = pm.L1[k], p1 = p0 + n, N = pm.N;
    ybest.assign( (size_t)n, 1.0 / std::max( 1L, n ) );
    if( n <= 1 || p1 > pm.L || !pm.X || !pm.U || !pm.A || !pm.lnGam )
        return 1e300;
    struct Save { double* p; std::vector<double> v; };
    std::vector<Save> saved;
    auto keep = [&saved]( double* p, long int m ) { if( p && m > 0 ) saved.push_back( { p, std::vector<double>( p, p + m ) } ); };
    keep( pm.X, pm.L );      keep( pm.XF, pm.FI );     keep( pm.XFA, pm.FIs );
    keep( pm.G, pm.L );      keep( pm.lnGam, pm.L );   keep( pm.lnGmo, pm.L );
    keep( pm.Gamma, pm.L );  keep( pm.F0, pm.L );      keep( pm.fDQF, pm.L );
    keep( pm.Wx, pm.L );     keep( pm.FitVar, 5 );
    auto restoreBase = [&saved]() { for( const Save& c : saved ) std::copy( c.v.begin(), c.v.end(), c.p ); };

    std::vector<double> c0( (size_t)n );
    for( long int j = p0; j < p1; j++ )
    {
        double au = 0.;
        for( long int i = 0; i < N; i++ ) au += pm.A[ i + j*N ] * pm.U[i];
        c0[(size_t)(j-p0)] = au - ( pm.G0[j] + pm.fDQF[j] );
    }
    double tpdMin = 1e300;
    auto evalY = [&]( const std::vector<double>& y, std::vector<double>& ynew ) {
        restoreBase();
        for( long int j = p0; j < p1; j++ ) pm.X[j] = std::max( y[(size_t)(j-p0)], 1e-300 );
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );
        double mx = -1e300;
        std::vector<double> lw( (size_t)n );
        for( long int a = 0; a < n; a++ ) { lw[(size_t)a] = c0[(size_t)a] - pm.lnGam[p0+a]; mx = std::max( mx, lw[(size_t)a] ); }
        double s = 0.;
        for( double v : lw ) s += exp( v - mx );
        ynew.resize( (size_t)n );
        for( long int a = 0; a < n; a++ ) ynew[(size_t)a] = exp( lw[(size_t)a] - mx ) / s;
        double tpd = 0.;
        for( long int a = 0; a < n; a++ )
            if( y[(size_t)a] > 0. ) tpd += y[(size_t)a] * ( log( y[(size_t)a] ) - lw[(size_t)a] );
        if( tpd < tpdMin ) { tpdMin = tpd; ybest = y; }
    };
    try
    {
        std::vector<std::vector<double>> starts;
        for( long int a = 0; a < n; a++ )
        {
            std::vector<double> v( (size_t)n, 1e-6 / std::max( 1L, n - 1 ) );
            v[(size_t)a] = 1. - 1e-6;
            starts.push_back( v );
        }
        {
            const double mx = *std::max_element( c0.begin(), c0.end() );
            double s = 0.;
            std::vector<double> v( (size_t)n );
            for( long int a = 0; a < n; a++ ) { v[(size_t)a] = exp( c0[(size_t)a] - mx ); s += v[(size_t)a]; }
            for( double& x : v ) x /= s;
            starts.push_back( v );
        }
        starts.push_back( std::vector<double>( (size_t)n, 1.0 / n ) );
        if( n == 2 )
            for( int g = 1; g < 20; g++ ) starts.push_back( { g / 20., 1. - g / 20. } );
        else
            for( long int a = 0; a < n && starts.size() < 64; a++ )
                for( long int b2 = a + 1; b2 < n && starts.size() < 64; b2++ )
                {
                    std::vector<double> v( (size_t)n, 1e-6 ); v[(size_t)a] = 0.5; v[(size_t)b2] = 0.5;
                    double t = 0.; for( double x : v ) t += x; for( double& x : v ) x /= t;
                    starts.push_back( v );
                }
        std::vector<double> tmp;
        for( auto y : starts )
            for( int it = 0; it < 300; it++ )
            {
                evalY( y, tmp );
                double d = 0.;
                for( long int a = 0; a < n; a++ ) d = std::max( d, fabs( tmp[(size_t)a] - y[(size_t)a] ) );
                y = tmp;
                if( d < 1e-11 ) break;
            }
    }
    catch( ... ) { tpdMin = 1e300; }
    try
    {
        restoreBase();
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );   // re-seed the TSolMod objects at the base state
    }
    catch( ... ) {}
    restoreBase();
    return tpdMin;
}

void native_trace_run_result( const MULTI& pm, long int mode, long int status, TMultiBase* mb )
{
    FILE* ntf = native_trace_file();
    if( !ntf || native_trace_quiet_depth > 0 )
        return;

    // The requested mode, captured by the caller before the dispatch.
    const char* mname;
    switch( mode )
    {
      case 1:  mname = "AIA";  break;
      case 5:  mname = "SIA";  break;
      case 10: mname = "AOP";  break;
      case 14: mname = "SOP";  break;
      case 18: mname = "ROP";  break;
      case 22: mname = "HOP";  break;
      case 26: mname = "SHP";  break;
      default: mname = "?";    break;
    }

    // FNV-1a over the present-phase names, order-independent by XOR-folding each phase's hash:
    // the assemblage is a set.
    unsigned long long akey = 0ull;
    long int nPresent = 0;

    // Phase names - full, and unique within the system (names=2). pm.SF[k] is MAXSYMB characters
    // of phase-class code and padding, then the MAXPHNAME-character name: the internal run of
    // blanks becomes one underscore and the trailing padding is trimmed, so each phase is one
    // token. A repeated name carries its occurrence number, "~2", "~3", counted over all FI
    // phases in system order, so a phase keeps the same label whichever phases are present.
    std::vector<std::string> keyName( (size_t)pm.FI );
    {
        std::map<std::string, long int> seen;
        for( long int k = 0; k < pm.FI; k++ )
        {
            std::string pn = char_array_to_string( pm.SF[k], MAXSYMB + MAXPHNAME );
            while( !pn.empty() && pn.back() == ' ' ) pn.pop_back();
            std::string t; bool sp = false;
            for( char c : pn )
            {
                if( c == ' ' || c == '\t' ) { sp = true; continue; }
                if( sp && !t.empty() ) t += '_';
                sp = false;
                // ',' ':' '=' are this line's own separators; some phase names contain a comma.
                t += ( c == ',' || c == ':' || c == '=' ) ? '/' : c;
            }
            const long int occ = ++seen[t];
            keyName[(size_t)k] = occ > 1 ? t + "~" + std::to_string( occ ) : t;
        }
    }

    fprintf( ntf, "KEY   mode=%s status=%ld pH=%.6f pe=%.6f IS=%.6e phases=",
             mname, (long)status, pm.pH, pm.pe, pm.IC );
    for( long int k = 0; k < pm.FI; k++ )
    {
        if( pm.XF[k] <= pm.DSM )
            continue;
        const std::string& pn = keyName[(size_t)k];
        unsigned long long h = 1469598103934665603ull;
        for( char c : pn ) { h ^= (unsigned char)c; h *= 1099511628211ull; }
        akey ^= h;
        // amount and molar volume: the second is the raw material for the fluid-root
        // axis this deliberately does not classify (see above). FVOL is cm3.
        const double vmol = ( pm.XF[k] > 0. && pm.FVOL != nullptr ? pm.FVOL[k] / pm.XF[k] : 0. );
        fprintf( ntf, "%s%s:%.6e:%.6e", ( nPresent ? "," : "" ), pn.c_str(), pm.XF[k], vmol );
        nPresent++;
    }
    fprintf( ntf, " nph=%ld akey=%016llx names=2\n", (long)nPresent, akey );

    // CERT record - the answer certificate, computed on the returned amounts pm.X for every mode:
    //   mb_rel  worst |C_i| / (B_i * DHBM) over the ordinary ICs [0, N-E); > 1 fails the relative
    //           test MBR applies (an IC with B_i = 0 and a non-zero residual reads 1e300)
    //   mb_abs  worst |C_i| over the same range
    //   chg_abs worst |C_i| over the charge row(s) [N-E, N)
    //   mb_pass mb_rel <= 1 (pa_DT's absolute floor is not applied)
    // Report-only fields, not part of mb_pass:
    //   kkt_max        worst sign-aware reduced-gradient residual over the present species, in
    //                  RT, at the dual pm.U[] (compare with pa_GAS; read dual_dirs first).
    //   curv_min       smallest eigenvalue of any present multicomponent non-aqueous phase's
    //                  symmetrised FD curvature block; < 0 = converged inside a spinodal;
    //                  1e300 = no phase qualified.
    //   dual_dirs      free dual directions = N - rank of the interior species' stoichiometry;
    //                  > 0 means the dual is not determined by the answer, so kkt_max on the
    //                  bound-active species depends on where the solver stopped.
    //                  dual_rank/dual_of carry the rank and N.
    //   dual_resolved  -1 no free-dual search ran, 0 searched and failed, 1 searched and resolved
    //                  (only the Optima path searches).
    // Total G and pm.Falp are not recorded: they are refreshed by the native path only. Runs
    // only when the trace file is open.
    double certKkt = -1., certCurv = 1e300;
    long int certKktJ = -1, certCurvK = -1, certRank = 0, certNic = pm.N, certDirs = -1;
    // stab_ss / stab_ph / stab_n / stab_dis - pa_StabTPD: min TPD in RT over the absent
    // multicomponent phases scanned (1e300 = none scanned or pa_StabTPD = 0), the phase carrying
    // it, how many were scanned, and how many the single-point index misclassifies.
    double certStab = 1e300;
    long int certStabK = -1, certStabN = 0, certStabDis = 0;
    if( mb )
    {
        std::vector<double> certF;
        mb->CertPrimalPotentials( certF );
        certKkt  = mb->CertKktMax( certF, certKktJ );
        certDirs = mb->CertDualFreeDirs( certF, certRank, certNic );
        certCurv = mb->CertCurvMin( certCurvK );
        certStab = mb->CertStabTPD( certStabK, certStabN, certStabDis );
    }
    // Separator-safe, as the KEY record's names: ' ', ',', ':' and '=' are this line's separators.
    auto certSafeName = []( const char* raw, size_t len ) {
        std::string s0 = char_array_to_string( raw, (int)len );
        while( !s0.empty() && ( s0.back() == ' ' || s0.back() == '\0' ) ) s0.pop_back();
        std::string t; bool sp = false;
        for( char c : s0 )
        {
            if( c == ' ' || c == '\t' ) { sp = true; continue; }
            if( sp && !t.empty() ) t += '_';
            sp = false;
            t += ( c == ',' || c == ':' || c == '=' ) ? '/' : c;
        }
        return t.empty() ? std::string( "-" ) : t;
    };
    const std::string certKktName  = certKktJ  >= 0 ? certSafeName( pm.SM[certKktJ], MAXDCNAME )
                                                    : std::string( "-" );
    const std::string certCurvName = certCurvK >= 0 ? certSafeName( pm.SF[certCurvK], MAXSYMB + MAXPHNAME )
                                                    : std::string( "-" );
    const std::string certStabName = certStabK >= 0 ? certSafeName( pm.SF[certStabK], MAXSYMB + MAXPHNAME )
                                                    : std::string( "-" );

    if( pm.X && pm.B && pm.A && pm.N > 0 )
    {
        const long int Z = pm.N - pm.E;
        long int iRel = -1, iAbs = -1, iChg = -1, iSeed = -1;
        double rel = 0., absr = 0., chg = 0., relSeed = 0.;
        for( long int i = 0; i < pm.N; i++ )
        {
            double c = pm.B[i];
            for( long int j = 0; j < pm.L; j++ )
                c -= pm.A[i + j*pm.N] * pm.X[j];
            const double a = fabs( c );
            if( i >= Z )
            {
                if( a > chg ) { chg = a; iChg = i; }
                continue;
            }
            if( a > absr ) { absr = a; iAbs = i; }
            const double bar = pm.B[i] * pm.DHBM;
            const double r = bar > 0. ? a / bar : ( a > 0. ? 1e300 : 0. );
            // A default seed (ICIsDefaultSeed()) is reported apart as mb_seed_rel and does not
            // enter mb_pass.
            if( mb && mb->ICIsDefaultSeed( i ) )
            {
                if( r > relSeed ) { relSeed = r; iSeed = i; }
                continue;
            }
            if( r > rel ) { rel = r; iRel = i; }
        }
        auto icName = [&pm]( long int i ) {
            if( i < 0 ) return std::string( "-" );
            std::string s = char_array_to_string( pm.SB[i], MAXICNAME );
            while( !s.empty() && ( s.back() == ' ' || s.back() == '\0' ) ) s.pop_back();
            return s.empty() ? std::string( "-" ) : s;
        };
        // mb_rel_b: the bulk amount of the IC carrying mb_rel. mb_seed_rel / mb_seed_ic: worst
        // relative residual over the default seeds, which mb_pass leaves out (0 and "-" when
        // nothing is marked of interest). New fields are appended at the end of the line, so
        // existing readers are unaffected.
        fprintf( ntf, "CERT  mode=%s status=%ld mb_rel=%.3e mb_rel_ic=%s mb_rel_b=%.3e mb_abs=%.3e"
                      " mb_abs_ic=%s chg_abs=%.3e chg_ic=%s mb_pass=%d mb_seed_rel=%.3e mb_seed_ic=%s"
                      " kkt_max=%.3e kkt_species=%s curv_min=%.6e curv_phase=%s"
                      " dual_dirs=%ld dual_rank=%ld dual_of=%ld dual_resolved=%d"
                      " stab_ss=%.6e stab_ph=%s stab_n=%ld stab_dis=%ld\n",
                 mname, (long)status, rel, icName( iRel ).c_str(), iRel >= 0 ? pm.B[iRel] : 0.,
                 absr, icName( iAbs ).c_str(), chg, icName( iChg ).c_str(), rel <= 1. ? 1 : 0,
                 relSeed, icName( iSeed ).c_str(),
                 certKkt, certKktName.c_str(), certCurv, certCurvName.c_str(),
                 (long)certDirs, (long)certRank, (long)certNic, native_cert_dualfree_get(),
                 certStab, certStabName.c_str(), (long)certStabN, (long)certStabDis );
    }

    // RANK record: the geometry and conditioning of the present species' stoichiometry at the
    // returned answer, for every mode.
    //   rank/of/pres   numerical rank of the present species' stoichiometry (row-scaled) against
    //                  pm.N, and how many species were present.
    //   sv_ratio(_raw) sigma_min/sigma_max of the present columns, row-scaled and unscaled.
    //   chg_res/span   whether the charge row, over present columns, is implied by the element
    //                  rows at this answer (chg_span=1: it carries no independent information).
    //   cond_ipm(_jac) condition number of A_p diag(w) A_p^T for IPM's weight, plain and after
    //   cond_mbr(_jac) symmetric Jacobi scaling; likewise for MBR's weight.
    // Report-only; computed only when the trace file is open.
    if( mb && pm.X && pm.A && pm.N > 0 )
    {
        TMultiBase::CertRankReport rk;
        mb->CertRank( rk );
        fprintf( ntf, "RANK  mode=%s status=%ld rank=%ld of=%ld pres=%ld sv_ratio=%.3e"
                      " sv_ratio_raw=%.3e chg_res=%.3e chg_span=%d cond_ipm=%.3e cond_ipm_jac=%.3e"
                      " cond_mbr=%.3e cond_mbr_jac=%.3e\n",
                 mname, (long)status, (long)rk.rank, (long)rk.of, (long)rk.pres,
                 rk.sv_ratio, rk.sv_ratio_raw, rk.chg_res, rk.chg_span,
                 rk.cond_ipm, rk.cond_ipm_jac, rk.cond_mbr, rk.cond_mbr_jac );
    }
    fflush( ntf );
}

/// Worst relative and worst absolute mass-balance residual of the amount vector `amt`, recomputed
/// locally (pm.C[] is MBR's scratch and one update stale), over the ordinary IC range [0, Z) -
/// the charge row is excluded, as in MBR's convergence loops - each with the IC that carries it.
/// `rel` is normalised so that rel > 1 means "fails the relative test MBR applies".
static void native_trace_mb_of( const MULTI& pm, const double* amt,
                                long int& iRelOut, double& relOut,
                                long int& iAbsOut, double& absOut )
{
    const long int Z = pm.N - pm.E;
    iRelOut = iAbsOut = -1; relOut = absOut = 0.;
    for( long int i = 0; i < Z; i++ )
    {
        double c = pm.B[i];
        for( long int j = 0; j < pm.L; j++ )
            c -= pm.A[i + j*pm.N] * amt[j];
        const double a = fabs( c );
        if( a > absOut ) { absOut = a; iAbsOut = i; }
        const double bar = pm.B[i] * pm.DHBM;
        const double r = bar > 0. ? a / bar : ( a > 0. ? 1e300 : 0. );
        if( r > relOut ) { relOut = r; iRelOut = i; }
    }
}
static void native_trace_mb( const MULTI& pm, long int& iRelOut, double& relOut,
                                              long int& iAbsOut, double& absOut )
{
    native_trace_mb_of( pm, pm.Y, iRelOut, relOut, iAbsOut, absOut );
}

/// One STAGE line for MBR/IPM entry and exit. `what` is "enter" or "exit";
/// eRet < 0 omits the return-code field. The residual pair is only meaningful
/// once a primal exists, hence `withResidual`.
///
/// `binding` names which of MBR's two tests would reject this state, using the
/// project's own pa_DT exactly as the convergence branches do: with DT == 0
/// only the relative test exists; with DT != 0 an IC must exceed BOTH (the
/// absolute cutoff is a floor under the relative bar - see the long comment at
/// that branch). "ok" means the state passes, i.e. MBR would call it converged.
static void native_trace_stage( const MULTI& pm, const BASE_PARAM* pa_p,
                                const char* stage, const char* what,
                                long int eRet, bool withResidual )
{
    FILE* fp = native_trace_file();
    if( !fp ) return;
    long int nz = 0;
    for( long int j = 0; j < pm.L; j++ )
        if( pm.Y[j] == 0. ) nz++;
    fprintf( fp, "STAGE %-3s %-5s K2=%ld ITF=%ld ITG=%ld IT=%ld zeroDC=%ld",
             stage, what, (long)pm.K2, (long)pm.ITF, (long)pm.ITG, (long)pm.IT, (long)nz );
    if( eRet >= 0 )
        fprintf( fp, " eRet=%ld", (long)eRet );
    if( withResidual )
    {
        long int ir, ia; double r, a;
        native_trace_mb( pm, ir, r, ia, a );
        const double e = fabs( (double)pa_p->DT );
        const double absCut = ( e < 2. ) ? pm.DHBM : pow( 10., -e );
        const bool relFail = ( r > 1. );
        const bool absFail = ( a > absCut );
        const char* binding = !pa_p->DT ? ( relFail ? "rel" : "ok" )
                                        : ( ( relFail && absFail ) ? "rel+abs" : "ok" );
        fprintf( fp, " worstRelIC=%s rel=%.6e worstAbsIC=%s abs=%.6e absCut=%.3e DT=%d binding=%s",
                 ir >= 0 ? char_array_to_string( pm.SB[ir], MAXICNAME ).c_str() : "-", r,
                 ia >= 0 ? char_array_to_string( pm.SB[ia], MAXICNAME ).c_str() : "-", a,
                 absCut, (int)pa_p->DT, binding );
    }
    fprintf( fp, "\n" );
    fflush( fp );
}

/// One MBRX line naming the exit path MassBalanceRefinement() took. The return code does not
/// distinguish "every IC passed" from "gave up, and the strict check did not apply": the
/// post-loop guard
///
///     if( pa_p->DW && ( WhereCalledFrom == 0L || pm.pNP ) )
///
/// is false for a cold call's second MBR (WhereCalledFrom = pm.K2 >= 1, pm.pNP = 0), which then
/// returns iRet = 0 on a state its own per-IC test rejects.
static void native_trace_mbr_exit( const MULTI& pm, const BASE_PARAM* pa_p,
                                   const char* reason, long int whereFrom,
                                   long int it1, long int iRet, bool restored )
{
    FILE* fp = native_trace_file();
    if( !fp ) return;
    const bool strictApplies = ( pa_p->DW && ( whereFrom == 0L || pm.pNP ) );
    fprintf( fp, "MBRX  reason=%s from=%ld IT1=%ld DP=%d iRet=%ld restoredBest=%d"
                 " DW=%d pNP=%ld strictCheckApplies=%d\n",
             reason, (long)whereFrom, (long)it1, (int)pa_p->DP, (long)iRet,
             (int)restored, (int)pa_p->DW, (long)pm.pNP, (int)strictApplies );
    fflush( fp );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Call to GEM IPM calculation of equilibrium state in MULTI
/// (with already scaled GEM problem)
// pa_MbReproject: repairs an unsatisfied mass balance on the answer by projecting Y back onto
// A.Y = b with one solve over a set of N species. The set is rank-revealing: species in
// decreasing amount, each kept only if its column raises the rank (the most abundant species
// alone can be exactly linearly dependent, e.g. H2O = H+ + OH-). All N rows take part, the
// charge row included, so the correction cannot introduce a charge imbalance.
bool TMultiBase::MassBalanceReproject( double* amt, bool keepPartial )
{
    const long int N = pm.N, L = pm.L;
    if( N < 1 || L < N || !pm.A || !pm.B || !amt ) return false;

    long int i1 = -1, i2 = -1; double relOld = 0., absOld = 0.;
    native_trace_mb_of( pm, amt, i1, relOld, i2, absOld );
    if( !( relOld > 0. ) ) return false;

    // The state on entry, for a full revert if the whole attempt fails to improve.
    std::vector<double> Xorig( amt, amt + L );

    // Iterated projection: when a component has to be clamped, the next pass sorts it last, so
    // the rank-revealing selection picks a different carrier for that direction. Each pass must
    // strictly improve the worst relative residual or it is undone and the loop stops.
    // Bounded: the extra passes are kept only if they bring the answer under its own tolerance;
    // otherwise the single-pass result is kept, since a partial repair still moves the answer.
    // With keepPartial the multi-pass result is kept anyway (used where amt only seeds the next
    // solver iterations and is not the returned answer).
    const long int maxPass = 8;
    long int nClamped = 0, nPass = 0;
    double relCur = relOld, absCur = absOld;
    bool improved = false;
    std::vector<double> Xsingle; double relSingle = 0., absSingle = 0.;
    std::vector<long int> pivKept, pivSingle;       // carrier species of the repair that is kept

    for( long int pass = 0; pass < maxPass; pass++ )
    {
    std::vector<double> Xpass( amt, amt + L );      // undo target for THIS pass alone

    std::vector<double> C( (size_t)N );
    for( long int i = 0; i < N; i++ )
    {
        double c = pm.B[i];
        for( long int j = 0; j < L; j++ ) c -= pm.A[i + j*N] * amt[j];
        C[(size_t)i] = c;
    }

    // Rank-revealing pivot by modified Gram-Schmidt against an orthonormal basis of
    // the columns accepted so far. O(L*N^2), and only reached on a run that has
    // already failed its own mass-balance test.
    std::vector<long int> ord( (size_t)L );
    for( long int j = 0; j < L; j++ ) ord[(size_t)j] = j;
    std::sort( ord.begin(), ord.end(),
               [amt]( long int a, long int b ) { return amt[a] > amt[b]; } );

    std::vector<long int> piv;
    std::vector<double> Q, r( (size_t)N );
    for( size_t t = 0; t < ord.size() && (long int)piv.size() < N; t++ )
    {
        const long int j = ord[t];
        double nrm0 = 0.;
        for( long int i = 0; i < N; i++ )
        { r[(size_t)i] = pm.A[i + j*N]; nrm0 += r[(size_t)i]*r[(size_t)i]; }
        nrm0 = sqrt( nrm0 );
        if( !( nrm0 > 0. ) ) continue;
        const long int k = (long int)piv.size();
        for( long int c = 0; c < k; c++ )
        {
            double d = 0.;
            for( long int i = 0; i < N; i++ ) d += Q[(size_t)(c*N + i)] * r[(size_t)i];
            for( long int i = 0; i < N; i++ ) r[(size_t)i] -= d * Q[(size_t)(c*N + i)];
        }
        double nrm = 0.;
        for( long int i = 0; i < N; i++ ) nrm += r[(size_t)i]*r[(size_t)i];
        nrm = sqrt( nrm );
        if( nrm < 1e-8 * nrm0 ) continue;            // dependent on the set so far
        Q.resize( (size_t)((k+1)*N) );
        for( long int i = 0; i < N; i++ ) Q[(size_t)(k*N + i)] = r[(size_t)i] / nrm;
        piv.push_back( j );
    }
    if( (long int)piv.size() < N ) break;            // rank(A) < N - nothing to do

    // Ap * dy = C, Gaussian elimination with partial pivoting (row-major).
    bool singular = false;
    std::vector<double> M( (size_t)N*N ), rhs( C ), dy( (size_t)N, 0. );
    for( long int i = 0; i < N; i++ )
        for( long int c = 0; c < N; c++ )
            M[(size_t)(i*N + c)] = pm.A[i + piv[(size_t)c]*N];
    for( long int k = 0; k < N; k++ )
    {
        long int pk = k; double best = fabs( M[(size_t)(k*N + k)] );
        for( long int i = k+1; i < N; i++ )
        { const double v = fabs( M[(size_t)(i*N + k)] ); if( v > best ) { best = v; pk = i; } }
        if( !( best > 1e-300 ) ) { singular = true; break; }
        if( pk != k )
        {
            for( long int c = k; c < N; c++ )
                std::swap( M[(size_t)(k*N + c)], M[(size_t)(pk*N + c)] );
            std::swap( rhs[(size_t)k], rhs[(size_t)pk] );
        }
        for( long int i = k+1; i < N; i++ )
        {
            const double f = M[(size_t)(i*N + k)] / M[(size_t)(k*N + k)];
            if( f == 0. ) continue;
            for( long int c = k; c < N; c++ )
                M[(size_t)(i*N + c)] -= f * M[(size_t)(k*N + c)];
            rhs[(size_t)i] -= f * rhs[(size_t)k];
        }
    }
    if( singular ) break;
    for( long int k = N-1; k >= 0; k-- )
    {
        double sum = rhs[(size_t)k];
        for( long int c = k+1; c < N; c++ ) sum -= M[(size_t)(k*N + c)] * dy[(size_t)c];
        dy[(size_t)k] = sum / M[(size_t)(k*N + k)];
    }

    // Feasibility: a component that would go negative is clamped, and the acceptance test below
    // decides (a clamped step no longer satisfies Ap.dy = C exactly, so it is kept only if it
    // still strictly improves the worst relative residual). The number of clamped components
    // separates an ill-conditioned but solved projection from a carrier that ran out.
    long int nClampedPass = 0;
    for( long int c = 0; c < N; c++ )
    {
        const double x = amt[piv[(size_t)c]];
        if( x + dy[(size_t)c] < 0. ) { dy[(size_t)c] = -x; nClampedPass++; }
    }

    for( long int c = 0; c < N; c++ ) amt[piv[(size_t)c]] += dy[(size_t)c];

    double relNew = 0., absNew = 0.;
    native_trace_mb_of( pm, amt, i1, relNew, i2, absNew );
    if( !( relNew < relCur ) )               // this PASS did not help - undo this pass
    {
        for( long int j = 0; j < L; j++ ) amt[j] = Xpass[(size_t)j];
        break;
    }
    relCur = relNew; absCur = absNew;
    nClamped += nClampedPass; nPass++; improved = true; pivKept = piv;
    if( pass == 0 )                  // remember exactly what the old single pass gave
    { Xsingle.assign( amt, amt + L ); relSingle = relNew; absSingle = absNew; pivSingle = piv; }

    if( nClampedPass == 0 ) break;   // a clean solve leaves nothing for a further pass
    if( relCur < 1. ) break;         // the answer now passes its own mass-balance test
    }   // pass loop

    // Short of a full repair, fall back to the single-pass result (see "Bounded" above),
    // unless keepPartial.
    if( improved && !( relCur < 1. ) && !Xsingle.empty() && keepPartial )
    {
        if( nPass > 1 )
            native_trace_decide( "mbreproject-partialkept species=%ld passes=%ld "
                                 "relsingle=%.3e relkept=%.3e",
                                 (long)N, (long)nPass, relSingle, relCur );
    }
    else if( improved && !( relCur < 1. ) && !Xsingle.empty() )
    {
        // Diagnostic only - this record gates nothing; the discard below is unconditional.
        // Counts, for a discarded multi-pass repair, the ICs left with a single carrier and the
        // rank, compared with the single-pass result. O(N*L), only when the trace is open.
        if( nPass > 1 && native_trace_file() )
        {
            const double thr = std::min( pm.lowPosNum, pm.DcMinM );
            auto singleCarrierICs = [&]( const double* x ) {
                long int ns = 0;
                for( long int i = 0; i < N; i++ )
                {
                    long int c = 0;
                    for( long int j = 0; j < L && c < 2; j++ )
                        if( pm.A[i + j*N] != 0. && x[j] > thr ) c++;
                    if( c == 1 ) ns++;
                }
                return ns;
            };
            const long int singleS = singleCarrierICs( Xsingle.data() );
            const long int singleK = singleCarrierICs( amt );
            native_trace_decide( "mbreproject-partialdiscarded species=%ld passes=%ld "
                                 "relsingle=%.3e relkept=%.3e singlecarrier=%ld>%ld",
                                 (long)N, (long)nPass, relSingle, relCur, singleS, singleK );
        }
        for( long int j = 0; j < L; j++ ) amt[j] = Xsingle[(size_t)j];
        relCur = relSingle; absCur = absSingle; nPass = 1; pivKept = pivSingle;
    }

    if( !improved )
    {
        for( long int j = 0; j < L; j++ ) amt[j] = Xorig[(size_t)j];
        // DECIDE record of a repair that was computed and then reverted.
        native_trace_decide( "mbreproject-reverted species=%ld clamped=%ld relbefore=%.3e "
                             "absbefore=%.3e", (long)N, (long)nClamped, relOld, absOld );
        return false;
    }
    const double relNew = relCur, absNew = absCur;

    // Keep both amount vectors and everything derived from them consistent with
    // the repaired state - pm.pH, FVOL, IC and the rest come from
    // CalculateConcentrations. Whichever array was repaired, the other follows it.
    for( long int j = 0; j < L; j++ ) { pm.X[j] = amt[j]; pm.Y[j] = amt[j]; }
    TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
    CalculateConcentrations( pm.X, pm.XF, pm.XFA );

    // Names for the log: the worst element before and after, and the species that carried the repair.
    long int j1 = -1, j2 = -1; double relChk = 0., absChk = 0.;
    native_trace_mb_of( pm, amt, j1, relChk, j2, absChk );
    auto icName = []( const MULTI& m, long int i ) {
        return i >= 0 ? char_array_to_string( m.SB[i], MAXICNAME ) : std::string( "-" ); };
    std::string carriers;
    for( size_t t = 0; t < pivKept.size() && t < 10; t++ )
        carriers += ( t ? ", " : "" ) + char_array_to_string( pm.SM[pivKept[t]], MAXDCNAME );
    if( pivKept.size() > 10 ) carriers += ", ...";
    if( !( relNew < 1. ) )
        gems_logger->warn( "Mass balance is still off after repair (worst element {0}, {1:.1e}x its tolerance). "
                           "Try: a warm restart from this result (SIA, SOP or SHP), or raise the amount of the smallest "
                           "element {0} (xGEMS: Material.min_amount or the min_amount argument). The test allows a residual of "
                           "pa_DHB = {2:.0e} times the element's own amount, so a very small element is checked very strictly: "
                           "an element that is only a placeholder can be raised until the check passes, as long as the amount "
                           "stays negligible next to the real components. The 'less material than the floor' warning "
                           "(Optima modes) gives the smallest amount that clears the solver's floor.",
                           icName( pm, j1 ), relNew, pm.DHBM );
    else
        gems_logger->debug( "pa_MbReproject: mass balance repaired in {} pass(es) over {} species ({}), {} "
                            "component(s) clamped at zero. Worst relative residual {:.3e}x ({}) -> {:.3e}x "
                            "({}) its own tolerance; worst absolute {:.3e} ({}) -> {:.3e} ({}) mol.",
                            nPass, N, carriers, nClamped, relOld, icName( pm, i1 ), relNew, icName( pm, j1 ),
                            absOld, icName( pm, i2 ), absNew, icName( pm, j2 ) );
    // DECIDE record of a repair that fired.
    native_trace_decide( "mbreproject species=%ld passes=%ld clamped=%ld relbefore=%.3e "
                         "relafter=%.3e absbefore=%.3e absafter=%.3e",
                         (long)N, (long)nPass, (long)nClamped, relOld, relNew, absOld, absNew );
    return true;
}

// Element classes - see the declarations in ms_multi.h. One place computes them, so the
// zeroing's rebalance test and the CERT record agree on which IC is a default seed.
bool TMultiBase::ICIsNumericalTrace( long int i ) const
{
    if( !pm.B || i < 0 || i >= pm.N - pm.E ) return false;      // charge rows are never trace
    const BASE_PARAM* pa = base_param();
    if( !pa || !( pa->DHB > 0. ) ) return false;
    double sumB = 0.;                                            // SystemTotalMolesIC()'s sum: ordinary ICs
    for( long int k = 0; k < pm.N - pm.E; k++ ) sumB += pm.B[k];
    if( !( sumB > 0. ) ) return false;
    const double floorInt = pa->OptimaDcFloor > 0. ? pa->OptimaDcFloor : pa->DHB;
    // the IC's amount in the units the floor is expressed in; B/sum(B) is the same on internal and real amounts
    const double bInt = pa->DG > 1e-5 ? pm.B[i] * ( pa->DG / sumB ) : pm.B[i];
    return bInt * pa->DHB < floorInt;
}

bool TMultiBase::ICIsDefaultSeed( long int i ) const
{
    if( elementsOfInterest.empty() || !ICIsNumericalTrace( i ) ) return false;
    std::string s = char_array_to_string( pm.SB[i], MAXICNAME );
    s.erase( s.find_last_not_of( std::string( " \t\0", 3 ) ) + 1 );
    return std::find( elementsOfInterest.begin(), elementsOfInterest.end(), s ) == elementsOfInterest.end();
}

// SubFloorElementCheck - see the declaration in ms_multi.h. One warning per call, naming every sub-floor element.
void TMultiBase::SubFloorElementCheck( double dcFloor ) const
{
    const long int N = pm.N, L = pm.L;
    if( N < 1 || L < 1 || !pm.A || !pm.B || !( dcFloor > 0. ) ) return;
    const double toReal = pm.SizeFactor > 0. ? 1. / pm.SizeFactor : 1.;
    std::string report;
    long int nWarn = 0;
    double floorHint = std::numeric_limits<double>::infinity(), bulkHint = 0.;
    for( long int i = 0; i < N; i++ )
    {
        if( !( pm.B[i] > 0. ) ) continue;
        if( pm.ICC && ( pm.ICC[i] == IC_CHARGE || pm.ICC[i] == IC_VOLUME ) ) continue;
        double need = 0., stoich = 0.;
        long int carriers = 0;
        for( long int j = 0; j < L; j++ )
        {
            const double a = pm.A[i + j*N];
            if( !( a > 0. ) ) continue;
            need += a * std::max( pm.DLL ? pm.DLL[j] : 0., dcFloor );
            stoich += a;
            carriers++;
        }
        if( carriers == 0 || need <= pm.B[i] ) continue;
        std::string name = char_array_to_string( pm.SB[i], MAXICNAME );
        name.erase( name.find_last_not_of( std::string( " \t\0", 3 ) ) + 1 );
        std::string tag;
        if( !elementsOfInterest.empty() )
            tag = ICIsDefaultSeed( i ) ? " [default seed]" : " [of interest]";
        report += fmt::format( "{}{}{}: {:.2g} mol, but its {} species must each hold at least {:.2g} mol, "
                               "{:.2g} mol of it in total", report.empty() ? "" : "; ", name, tag,
                               pm.B[i] * toReal, carriers, dcFloor * toReal, need * toReal );
        floorHint = std::min( floorHint, 0.1 * pm.B[i] / stoich );
        bulkHint = std::max( bulkHint, 10. * need * toReal );
        nWarn++;
    }
    if( !nWarn ) return;
    // The suggested amount as a fraction of everything in the system, so the user can judge that it is insignificant.
    double totalReal = 0.;
    for( long int i = 0; i < N; i++ )
        if( pm.B[i] > 0. && !( pm.ICC && ( pm.ICC[i] == IC_CHARGE || pm.ICC[i] == IC_VOLUME ) ) )
            totalReal += pm.B[i] * toReal;
    const double fraction = totalReal > 0. ? bulkHint / totalReal : 0.;
    gems_logger->warn(
        "Optima: {0} element(s) have less material than the solver's floor amount allows. Their mass "
        "balance is then fixed afterwards, so their amounts and pH/Eh may be wrong (native modes are not "
        "affected). To avoid it, raise the amount of each of these elements to about {1:.1g} mol or more: that "
        "is the smallest amount that clears the floor, and it is only {4:.0e} of the {5:.3g} mol in the whole system. "
        "If the amount is a minimum that the program added for an element you did not specify, raise that minimum "
        "(xGEMS: Material.min_amount or the min_amount argument of equilibrate/setB/clear, default 1e-11 mol); "
        "or remove the element if it is a placeholder; or set pa_OptimaDcFloor to {2:.1g} or lower (in the "
        "project's -ipm file; not in GEMS). Elements: {3}",
        nWarn, bulkHint, floorHint, report, fraction, totalReal );
}

// EnergyDeterminacyCheck: is each present phase's amount actually fixed by the energy?
// Warns when a present phase's amount is not fixed to pa_DeterminacyWarn (0 = skipped).
//
// Model. Moving present phase k by t mol while keeping A.x = b costs at least
//     E(t) = (1/2) t^2 / c_k,
// c_k the phase's compliance - the cheapest way to move k with every other species free to
// compensate under A dn = 0. Answers whose energies differ by less than the energy resolution
// eps_G are indistinguishable, so k's amount is fixed only to dn_k = sqrt( 2 eps_G c_k ).
//
// Curvature. A species of a multi-component phase has ideal curvature 1/x_j; a pure phase has
// none of its own, so pure phases are free variables.
//
// With S the active multi-component species (weights X_S), P the active pure phases,
//     K0 = [ A_S X_S A_S'   A_P ]        v_k = [ A_S X_S g_S ]      w_k = g_S' X_S g_S
//          [ A_P'           0   ]              [ g_P         ]
// for phase k's indicator g = (g_S, g_P), eliminating the border of the least-energy KKT system
// gives c_k = w_k - v_k' K0^-1 v_k (>= 0). One factorisation of K0, size N + |P|, then one solve
// per present phase. c_k == 0 means k is pinned by mass balance alone.
//
// When K0 is singular there are two causes:
//  (a) A_P z = 0 - pure phases with dependent stoichiometries: moving along z costs nothing, so
//      every pure phase with z_q != 0 has an amount the energy does not fix at all. A phase q is
//      degenerate iff dropping it does not lower rank(A_P); these are named, a maximal
//      independent subset stays in K0, and every other phase is still checked.
//  (b) A_S' y = A_P' y = 0 - a redundant IC row over the active species (e.g. the charge row as
//      the valence sum of the element rows); such rows are dropped before assembly.
// Both selections are greedy Gram-Schmidt on the unweighted stoichiometry at 1e-9 relative.
// A pivot below m*DBL_EPSILON of its original column scale is then reported as
// `determinacy-singular` and the check makes no statement.
//
// Energy resolution: eps_G = DBL_EPSILON * sum_j |x_j mu_j| (RT units), the rounding floor of G,
// mu_j = sum_i a_ij u_i. Species pinned at a kinetic bound (DLL/DUL) are excluded.
// Read-only: nothing here writes solver state.
void TMultiBase::EnergyDeterminacyCheck()
{
    const long int N = pm.N, L = pm.L, FI = pm.FI, FIs = pm.FIs;
    if( N < 1 || L < 1 || FI < 1 || !pm.A || !pm.X || !pm.U || !pm.L1 ) return;
    const bool probe = getenv( "GEMS3K_DETERMINACY_PROBE" ) != nullptr;
    const double warnRel = (double)base_param()->DeterminacyWarn;
    if( !( warnRel > 0. ) && !probe ) return;       // off: costs nothing

    std::vector<long int> phaseOf( (size_t)L, -1 );
    for( long int k = 0, jb = 0; k < FI; jb += pm.L1[k], k++ )
        for( long int j = jb; j < jb + pm.L1[k] && j < L; j++ ) phaseOf[(size_t)j] = k;

    std::vector<char> act( (size_t)L, 0 );
    double gabs = 0.;
    for( long int j = 0; j < L; j++ )
    {
        const double x = pm.X[j];
        if( !( x > 0. ) || phaseOf[(size_t)j] < 0 ) continue;
        if( pm.DUL && pm.DUL[j] < 1e6 && x >= pm.DUL[j] * ( 1. - 1e-9 ) ) continue;
        if( pm.DLL && pm.DLL[j] > 0. && x <= pm.DLL[j] * ( 1. + 1e-9 ) ) continue;
        act[(size_t)j] = 1;
        gabs += fabs( x * DC_DualChemicalPotential( pm.U, pm.A + j*N, pm.NR, j ) );
    }
    const double epsG = std::numeric_limits<double>::epsilon() * gabs;
    if( !( epsG > 0. ) ) return;

    // Greedy rank-revealing selection (see the header comment): keep vecs[i] iff its Gram-Schmidt
    // residual against the vectors already kept exceeds 1e-9 of its own norm.
    auto independentSubset = []( const std::vector<std::vector<double>>& vecs ) {
        std::vector<std::vector<double>> Q;
        std::vector<char> keep( vecs.size(), 0 );
        for( size_t i = 0; i < vecs.size(); i++ )
        {
            std::vector<double> r = vecs[i];
            double n0 = 0.;
            for( double e : r ) n0 += e*e;
            n0 = sqrt( n0 );
            if( !( n0 > 0. ) ) continue;
            for( int pass = 0; pass < 2; pass++ )          // twice is enough (Kahan-Parlett)
                for( const auto& q : Q )
                {
                    double d = 0.;
                    for( size_t t = 0; t < r.size(); t++ ) d += q[t] * r[t];
                    for( size_t t = 0; t < r.size(); t++ ) r[t] -= d * q[t];
                }
            double n1 = 0.;
            for( double e : r ) n1 += e*e;
            n1 = sqrt( n1 );
            if( n1 > 1e-9 * n0 )
            {
                for( double& e : r ) e /= n1;
                Q.push_back( std::move( r ) );
                keep[i] = 1;
            }
        }
        return keep;
    };

    std::vector<long int> actCols, Pall, rows;
    for( long int j = 0; j < L; j++ )
        if( act[(size_t)j] )
        {
            actCols.push_back( j );
            if( phaseOf[(size_t)j] >= FIs ) Pall.push_back( j );
        }

    // (b) IC rows: those touching an active species, minus any redundant over the active species.
    {
        std::vector<long int> cand;
        std::vector<std::vector<double>> rv;
        for( long int r = 0; r < N; r++ )
        {
            std::vector<double> row( actCols.size() );
            bool nz = false;
            for( size_t c = 0; c < actCols.size(); c++ )
                if( ( row[c] = pm.A[r + actCols[c]*N] ) != 0. ) nz = true;
            if( nz ) { cand.push_back( r ); rv.push_back( std::move( row ) ); }
        }
        const std::vector<char> keep = independentSubset( rv );
        for( size_t i = 0; i < cand.size(); i++ ) if( keep[i] ) rows.push_back( cand[i] );
        if( probe && rows.size() < cand.size() )
            fprintf( stderr, "DETPROBE redundant-rows=%ld of %ld\n", (long)( cand.size() - rows.size() ), (long)cand.size() );
    }
    const long int n = (long int)rows.size();
    if( n < 1 ) return;
    auto A = [&]( long int r, long int j ) { return pm.A[rows[(size_t)r] + j*N]; };
    auto trimmedPhaseName = [&]( long int k ) {
        // feeds comma-joined DECIDE payloads - see trace_safe_name()
        return trace_safe_name( char_array_to_string( pm.SF[k] + MAXSYMB, MAXPHNAME ) );
    };

    // (a) Pure phases: a maximal independent subset P goes into K0; a phase is DEGENERATE iff
    // removing it does not lower rank(A_P), i.e. it lies on some null vector of A_P.
    std::vector<long int> P;
    std::vector<char> degenerate( (size_t)FI, 0 );
    long int nDegenerate = 0;
    {
        auto colsOf = [&]( long int skip ) {
            std::vector<std::vector<double>> cv;
            for( size_t q = 0; q < Pall.size(); q++ )
            {
                if( (long int)q == skip ) continue;
                std::vector<double> col( (size_t)n );
                for( long int r = 0; r < n; r++ ) col[(size_t)r] = A( r, Pall[q] );
                cv.push_back( std::move( col ) );
            }
            return cv;
        };
        const std::vector<char> keep = independentSubset( colsOf( -1 ) );
        long int rank = 0;
        for( size_t q = 0; q < Pall.size(); q++ ) if( keep[q] ) { P.push_back( Pall[q] ); rank++; }
        if( rank < (long int)Pall.size() )
        {
            std::string names;
            for( size_t q = 0; q < Pall.size(); q++ )
            {
                long int rq = 0;
                for( char f : independentSubset( colsOf( (long int)q ) ) ) rq += f;
                const long int k = phaseOf[(size_t)Pall[q]];
                if( rq == rank && !degenerate[(size_t)k] )
                {
                    degenerate[(size_t)k] = 1;
                    nDegenerate++;
                    names += ( names.empty() ? "" : "," ) + trimmedPhaseName( k );
                }
            }
            native_trace_decide( "determinacy-degenerate purephases=%ld rank=%ld phases=%s",
                                 (long)Pall.size(), rank, names.c_str() );
            if( probe ) fprintf( stderr, "DETPROBE degenerate purephases=%ld rank=%ld phases=%s\n",
                                 (long)Pall.size(), rank, names.c_str() );
        }
    }
    const long int p = (long int)P.size(), m = n + p;

    std::vector<double> K( (size_t)(m*m), 0. );
    for( long int j = 0; j < L; j++ )
    {
        if( !act[(size_t)j] || phaseOf[(size_t)j] >= FIs ) continue;
        for( long int r = 0; r < n; r++ )
        {
            const double ar = A( r, j );
            if( ar != 0. ) for( long int c = 0; c < n; c++ ) K[(size_t)(r*m + c)] += pm.X[j] * ar * A( c, j );
        }
    }
    for( long int q = 0; q < p; q++ )
        for( long int r = 0; r < n; r++ )
            K[(size_t)(r*m + n+q)] = K[(size_t)((n+q)*m + r)] = A( r, P[(size_t)q] );

    // LU with partial pivoting; a pivot below m*DBL_EPSILON of its ORIGINAL column's scale is
    // singular (see the header comment - both structural causes were removed above).
    std::vector<double> colScale( (size_t)m, 0. );
    for( long int c = 0; c < m; c++ )
        for( long int i = 0; i < m; i++ ) colScale[(size_t)c] = std::max( colScale[(size_t)c], fabs( K[(size_t)(i*m + c)] ) );
    const double pivTol = (double)m * std::numeric_limits<double>::epsilon();
    std::vector<long int> perm( (size_t)m );
    for( long int i = 0; i < m; i++ ) perm[(size_t)i] = i;
    for( long int c = 0; c < m; c++ )
    {
        long int pc = c; double best = 0.;
        for( long int i = c; i < m; i++ )
        { const double v = fabs( K[(size_t)(i*m + c)] ); if( v > best ) { best = v; pc = i; } }
        if( !( best > 1e-300 ) || !( best > pivTol * colScale[(size_t)c] ) )
        {
            native_trace_decide( "determinacy-singular species=%ld purephases=%ld", (long)n, (long)p );
            if( probe ) fprintf( stderr, "DETPROBE singular=1 n=%ld p=%ld col=%ld relpivot=%.3e\n", (long)n, (long)p, (long)c,
                                 colScale[(size_t)c] > 0. ? best / colScale[(size_t)c] : 0. );
            return;
        }
        if( pc != c )
        {
            for( long int cc = 0; cc < m; cc++ ) std::swap( K[(size_t)(c*m + cc)], K[(size_t)(pc*m + cc)] );
            std::swap( perm[(size_t)c], perm[(size_t)pc] );
        }
        for( long int i = c+1; i < m; i++ )
        {
            const double f = ( K[(size_t)(i*m + c)] /= K[(size_t)(c*m + c)] );
            if( f != 0. ) for( long int cc = c+1; cc < m; cc++ ) K[(size_t)(i*m + cc)] -= f * K[(size_t)(c*m + cc)];
        }
    }

    long int nWarn = 0; double worstRel = 0.; long int kWorst = -1;
    std::string listed, listedDegenerate;
    long int nListedDegenerate = 0, nListedOther = 0;
    std::vector<double> v( (size_t)m ), z( (size_t)m );
    for( long int k = 0, jb = 0; k < FI; jb += pm.L1[k], k++ )
    {
        const long int je = std::min( jb + pm.L1[k], L );
        double w = 0., nk = 0.; bool present = false;
        std::fill( v.begin(), v.end(), 0. );
        for( long int j = jb; j < je; j++ )
        {
            if( !act[(size_t)j] ) continue;
            present = true; nk += pm.X[j];
            if( k < FIs )
            { w += pm.X[j]; for( long int r = 0; r < n; r++ ) v[(size_t)r] += pm.X[j] * A( r, j ); }
        }
        if( !present ) continue;
        if( degenerate[(size_t)k] )
        {
            // Not fixed at all (cause (a) above): no c_k is computed. Listed first, so the
            // 8-name cap never hides it.
            if( warnRel > 0. )
            {
                if( nListedDegenerate++ < 8 )
                    listedDegenerate += ( listedDegenerate.empty() ? "" : ", " ) + trimmedPhaseName( k )
                                     + " (interchangeable)";
                nWarn++;
            }
            continue;
        }
        if( k >= FIs )
            for( long int q = 0; q < p; q++ ) if( phaseOf[(size_t)P[(size_t)q]] == k ) v[(size_t)(n+q)] = 1.;

        for( long int i = 0; i < m; i++ ) z[(size_t)i] = v[(size_t)perm[(size_t)i]];
        for( long int i = 0; i < m; i++ )
            for( long int c = 0; c < i; c++ ) z[(size_t)i] -= K[(size_t)(i*m + c)] * z[(size_t)c];
        for( long int i = m-1; i >= 0; i-- )
        {
            for( long int c = i+1; c < m; c++ ) z[(size_t)i] -= K[(size_t)(i*m + c)] * z[(size_t)c];
            z[(size_t)i] /= K[(size_t)(i*m + i)];
        }
        double vKv = 0.;
        for( long int i = 0; i < m; i++ ) vKv += v[(size_t)i] * z[(size_t)i];
        const double ck  = std::max( w - vKv, 0. );
        const double rel = sqrt( 2. * epsG * ck ) / nk;
        const std::string name = trimmedPhaseName( k );

        if( probe )
            fprintf( stderr, "DETPROBE phase=%s class=%c pure=%d n=%.6e rel=%.3e c=%.3e epsG=%.3e\n",
                     name.c_str(), pm.PHC ? pm.PHC[k] : '?', (int)( k >= FIs ), nk, rel, ck, epsG );
        if( warnRel > 0. && rel >= warnRel )
        {
            // rel >= 1: the energy cannot even say whether the phase is present; named as such.
            if( nListedOther++ < 8 )
                listed += ( listed.empty() ? "" : ", " ) + name
                        + ( rel >= 1. ? std::string( " (presence)" ) : fmt::format( " ({:.2g} %)", 100. * rel ) );
            nWarn++;
            if( rel > worstRel ) { worstRel = rel; kWorst = k; }
        }
    }

    if( nWarn > 0 )
    {
        if( !listedDegenerate.empty() )
        {
            listed = listedDegenerate + ( listed.empty() ? "" : ", " ) + listed;
            worstRel = std::numeric_limits<double>::infinity();
        }
        const std::string worstTxt = nDegenerate > 0
            ? fmt::format( "whose amount is not fixed AT ALL: {} present pure phase(s) have stoichiometries that "
                           "compensate one another exactly (e.g. the same substance entered twice), so only their "
                           "combination is determined", nDegenerate )
            : worstRel >= 1.
            ? std::string( "whose amount lies inside its own energy resolution, so even its PRESENCE is undetermined" )
            : fmt::format( "fixed only to +-{:.2g} %", 100. * worstRel );
        if( nDegenerate > 0 ) kWorst = -1;
        gems_logger->warn(
            "Some phase amounts are not fixed by the energy: {} phase(s) are determined only to worse than "
            "{:.0f} % (worst: {}, {}). Other values are equally valid, so treat these amounts as "
            "approximate. Phases: {}",
            nWarn, 100. * warnRel,
            kWorst >= 0 ? trimmedPhaseName( kWorst ) : std::string( "an interchangeable phase" ), worstTxt, listed );
        native_trace_decide( "undetermined phases=%ld threshold=%.0e worst=%.2e list=%s",
                             (long)nWarn, warnRel, worstRel, listed.c_str() );
    }
}

// ExcludeRedundantDCs: a species entered twice is removed from the solve.
//
// Redundant means thermodynamically indistinguishable: identical stoichiometry (all N rows,
// charge included), identical DC class code, and identical standard properties at the current
// T,P (G0 including DQF terms, H0, S0, Cp0, molar volume), in one of two placements:
//  (a) twice in the same multi-component phase - the pair behaves as one species with G0
//      lowered by RT ln 2, so the duplicate changes the answer;
//  (b) as two single-species phases - their amounts are interchangeable and only the sum is
//      determined.
// Not redundant: same formula with different properties (polymorphs); two multi-component
// phases with identical member lists (that is how a miscibility gap is modelled); species of a
// multi-site phase; species of sorption / polyelectrolyte phases; a copy whose end-member (DMc)
// coefficients differ or that appears in the interaction-parameter index (IPx); the solvent.
// A pair with user metastability limits on either copy (DLL > 0 or DUL < 1e6) is reported but
// not changed.
// Removal uses the kinetic-exclusion path: for this call only, every copy after the first gets
// DLL = DUL = 0 with RLC = BOTH_LIM, and any starting amount moves to the kept copy. The caller's
// DATABR dll/dul are never written; RestoreRedundantDCs() restores pm.DLL/DUL/RLC after the
// solve. Warns once per distinct finding per process; a DECIDE record on every call.
std::vector<TMultiBase::RedundantDCHold> TMultiBase::ExcludeRedundantDCs()
{
    std::vector<RedundantDCHold> held;
    const long int N = pm.N, L = pm.L, FI = pm.FI, FIs = pm.FIs;
    if( N < 1 || L < 2 || FI < 1 || !pm.A || !pm.L1 || !pm.G0 || !pm.DCC || !pm.PHC
        || !pm.DUL || !pm.DLL || !pm.RLC )
        return held;

    auto same = []( double a, double b ) {
        return fabs( a - b ) <= 1e-12 * std::max( { fabs( a ), fabs( b ), 1. } );
    };
    auto identical = [&]( long int j1, long int j2 ) {
        if( pm.DCC[j1] != pm.DCC[j2] ) return false;
        for( long int i = 0; i < N; i++ )
            if( pm.A[i + j1*N] != pm.A[i + j2*N] ) return false;
        if( !same( pm.G0[j1], pm.G0[j2] ) ) return false;
        if( pm.H0 && !same( pm.H0[j1], pm.H0[j2] ) ) return false;
        if( pm.S0 && !same( pm.S0[j1], pm.S0[j2] ) ) return false;
        if( pm.Cp0 && !same( pm.Cp0[j1], pm.Cp0[j2] ) ) return false;
        if( pm.Vol && !same( pm.Vol[j1], pm.Vol[j2] ) ) return false;
        return true;
    };
    auto freeBounds = [&]( long int j ) { return !( pm.DLL[j] > 0. ) && !( pm.DUL[j] < 1e6 ); };
    // both feed comma-joined DECIDE payloads - see trace_safe_name()
    auto dcName = [&]( long int j ) {
        return trace_safe_name( char_array_to_string( pm.SM[j], MAXDCNAME ) );
    };
    auto phName = [&]( long int k ) {
        return trace_safe_name( char_array_to_string( pm.SF[k] + MAXSYMB, MAXPHNAME ) );
    };

    std::vector<long int> jb( (size_t)FI + 1, 0 );
    for( long int k = 0; k < FI; k++ ) jb[(size_t)k+1] = jb[(size_t)k] + pm.L1[k];

    std::vector<char> removed( (size_t)L, 0 );
    std::string report, reportKept;
    // kKeep/kDrop are the phases of keep/drop (equal for an in-phase pair); why = the keep rule that decided.
    auto hold = [&]( long int keep, long int drop, const char* kind, long int kKeep, long int kDrop, const char* why ) {
        if( !( freeBounds( keep ) && freeBounds( drop ) ) )
        {
            reportKept += fmt::format( "{}{}:{}={}({})", reportKept.empty() ? "" : ",", kind,
                                       dcName( keep ), dcName( drop ), kKeep == kDrop ? phName( kKeep ) : phName( kKeep ) + "/" + phName( kDrop ) );
            return;
        }
        removed[(size_t)drop] = 1;
        held.push_back( { drop, pm.RLC[drop], pm.DLL[drop], pm.DUL[drop] } );
        if( pm.X ) { pm.X[keep] += pm.X[drop]; pm.X[drop] = 0.; }
        if( pm.Y ) { pm.Y[keep] += pm.Y[drop]; pm.Y[drop] = 0.; }
        pm.DLL[drop] = 0.; pm.DUL[drop] = 0.; pm.RLC[drop] = BOTH_LIM;
        report += fmt::format( "{}{}:{}>{}({}{})", report.empty() ? "" : ",", kind, dcName( drop ), dcName( keep ),
                               kKeep == kDrop ? phName( kKeep ) : phName( kDrop ) + ">" + phName( kKeep ),
                               why && *why ? std::string( ",kept:" ) + why : std::string() );
    };
    // Suspicious: same class, stoichiometry and G0 at this T,P, but another standard property
    // differs. Not interchangeable, so nothing is removed; reported only.
    std::string reportSuspicious;
    auto suspicious = [&]( long int j1, long int j2 ) {
        if( pm.DCC[j1] != pm.DCC[j2] || !same( pm.G0[j1], pm.G0[j2] ) ) return false;
        for( long int i = 0; i < N; i++ )
            if( pm.A[i + j1*N] != pm.A[i + j2*N] ) return false;
        return true;   // caller has already found identical() false
    };

    // (a) within one multi-component phase. Offsets into IPx/DMc exactly as
    // CalculateActivityCoefficients() walks them (single-species non-gas phases carry none).
    long int ipe = 0, jde = 0;
    for( long int k = 0; k < FIs && k < FI; k++ )
    {
        const long int b = jb[(size_t)k], e = jb[(size_t)k+1], n1 = pm.L1[k];
        if( n1 == 1 && !( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_PLASMA || pm.PHC[k] == PH_FLUID ) )
            continue;
        const long int nPar = pm.LsMod ? pm.LsMod[k*3] : 0, maxOrd = pm.LsMod ? pm.LsMod[k*3+1] : 0;
        const long int perDC = pm.LsMdc ? pm.LsMdc[k*3] : 0, nSub = pm.LsMdc ? pm.LsMdc[k*3+1] : 0;
        const long int ipb = ipe, jdb = jde;
        ipe += nPar * maxOrd;
        jde += perDC * n1;
        if( n1 < 2 || nSub > 0 || pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL )
            continue;
        std::vector<char> inIPx( (size_t)n1, 0 );
        if( pm.IPx )
            for( long int t = 0; t < nPar * maxOrd; t++ )
            {
                const long int m = pm.IPx[ipb + t];
                if( m >= 0 && m < n1 ) inIPx[(size_t)m] = 1;
            }
        for( long int j2 = b + 1; j2 < e; j2++ )
        {
            if( removed[(size_t)j2] || j2 == pm.LO || inIPx[(size_t)(j2 - b)] ) continue;
            for( long int j1 = b; j1 < j2; j1++ )
            {
                if( removed[(size_t)j1] || j1 == pm.LO || inIPx[(size_t)(j1 - b)] || !identical( j1, j2 ) ) continue;
                bool sameDMc = true;
                if( pm.DMc )
                    for( long int c = 0; c < perDC && sameDMc; c++ )
                        sameDMc = pm.DMc[jdb + (j1 - b)*perDC + c] == pm.DMc[jdb + (j2 - b)*perDC + c];
                if( !sameDMc ) continue;
                hold( j1, j2, "inphase", k, k, "" );
                break;
            }
        }
    }

    // (b) single-species phases of the same class. Which copy to keep:
    //   1. "pure":   a pure phase (k >= FIs) over a single-species solution phase (k < FIs);
    //   2. "name":   otherwise the phase whose name the species symbol abbreviates (letters in
    //                order, first letter matching, case-insensitive: "Brc" -> Brucite);
    //   3. "first":  otherwise the first listed.
    auto abbreviates = []( const std::string& sym, const std::string& name ) {
        if( sym.empty() || name.empty() || tolower( (unsigned char)sym[0] ) != tolower( (unsigned char)name[0] ) )
            return false;
        size_t q = 0;
        for( char c : name )
            if( q < sym.size() && tolower( (unsigned char)c ) == tolower( (unsigned char)sym[q] ) ) q++;
        return q == sym.size();
    };
    for( long int k2 = 1; k2 < FI; k2++ )
    {
        if( pm.L1[k2] != 1 || pm.PHC[k2] == PH_SORPTION || pm.PHC[k2] == PH_POLYEL ) continue;
        const long int j2 = jb[(size_t)k2];
        if( removed[(size_t)j2] || ( k2 < FIs && pm.LsMdc && pm.LsMdc[k2*3] > 0 ) ) continue;
        for( long int k1 = 0; k1 < k2; k1++ )
        {
            if( pm.L1[k1] != 1 || pm.PHC[k1] != pm.PHC[k2] ) continue;
            const long int j1 = jb[(size_t)k1];
            if( removed[(size_t)j1] || ( k1 < FIs && pm.LsMdc && pm.LsMdc[k1*3] > 0 ) ) continue;
            if( !identical( j1, j2 ) )
            {
                if( suspicious( j1, j2 ) )
                    reportSuspicious += fmt::format( "{}{}/{}({})", reportSuspicious.empty() ? "" : ",",
                                                     phName( k1 ), phName( k2 ), dcName( j1 ) );
                continue;
            }
            long int keepK = k1, dropK = k2; const char* why = "first";
            const bool pure1 = k1 >= FIs, pure2 = k2 >= FIs;
            const bool name1 = abbreviates( dcName( j1 ), phName( k1 ) ), name2 = abbreviates( dcName( j2 ), phName( k2 ) );
            if( pure1 != pure2 )      { why = "pure"; if( pure2 ) { keepK = k2; dropK = k1; } }
            else if( name1 != name2 ) { why = "name"; if( name2 ) { keepK = k2; dropK = k1; } }
            hold( jb[(size_t)keepK], jb[(size_t)dropK], "purephase", keepK, dropK, why );
            break;
        }
    }

    if( !reportSuspicious.empty() )
        native_trace_decide( "redundant-dc-suspicious list=%s", reportSuspicious.c_str() );
    if( report.empty() && reportKept.empty() && reportSuspicious.empty() )
        return held;
    if( !held.empty() && pm.pNP && pm.X && pm.Y )
    {   // warm start: the moved amounts must be reflected in the phase totals already derived
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
    }
    if( !report.empty() )
        native_trace_decide( "redundant-dc removed=%ld list=%s", (long)held.size(), report.c_str() );
    if( !reportKept.empty() )
        native_trace_decide( "redundant-dc-constrained list=%s", reportKept.c_str() );

    static std::mutex warnedMutex;
    static std::set<std::string> warned;
    const std::string key = report + "|" + reportKept + "|" + reportSuspicious;
    {
        std::lock_guard<std::mutex> lock( warnedMutex );
        if( !warned.insert( key ).second ) return held;
    }
    if( !report.empty() )
        gems_logger->warn(
            "Redundant species: {} duplicate(s) were removed from the solve and are reported as 0 (the "
            "copy kept carries the amount). Removed>kept: {}. Try: delete the duplicates from the project.",
            held.size(), report );
    if( !reportKept.empty() )
        gems_logger->warn(
            "Redundant species kept because a copy has metastability limits (DLL/DUL): {}. Their amounts "
            "are not fixed separately unless the limits differ. Try: remove one copy.", reportKept );
    if( !reportSuspicious.empty() )
        gems_logger->warn(
            "Suspicious species pairs (nothing removed): same stoichiometry and G0 but different H0, S0, "
            "Cp0 or V0. Try: check the data of these entries: {}", reportSuspicious );
    return held;
}

void TMultiBase::RestoreRedundantDCs( const std::vector<RedundantDCHold>& held )
{
    for( const auto& h : held )
    {
        pm.RLC[h.j] = h.rlc;
        pm.DLL[h.j] = h.dll;
        pm.DUL[h.j] = h.dul;
    }
}

// StrandedElementCheck: an element that can only live in one multi-component phase and holds
// that phase open (e.g. trace elements whose only carriers are aqueous, in a system with very
// little water). Read-only, on the converged answer: IC i (not charge/volume, B[i] > 0) whose
// carriers - species with a(i,j) != 0 not excluded by DUL = 0 - all belong to one
// multi-component phase k. Then
//   share_i = sum_{carriers in k} X[j] / XF[k]    (fraction of phase k's moles that carry i)
//   trace_k = XF[k] / sum_k XF[k]                 (phase k's size relative to the system)
// Warn when share_i >= kStrandedShareWarn and trace_k <= kStrandedTraceWarn.
// GEMS3K_STRANDED_PROBE prints every confined element with both numbers. An element held only
// by pure phases is not flagged. The remedy is a host phase or removing the element.
void TMultiBase::StrandedElementCheck()
{
    const long int N = pm.N, L = pm.L, FI = pm.FI;
    if( N < 1 || L < 1 || FI < 1 || !pm.A || !pm.X || !pm.XF || !pm.L1 || !pm.B ) return;
    const bool probe = getenv( "GEMS3K_STRANDED_PROBE" ) != nullptr;
    static const double kStrandedShareWarn = 1e-2, kStrandedTraceWarn = 1e-6;

    std::vector<long int> phaseOf( (size_t)L, -1 );
    for( long int k = 0, jb = 0; k < FI; jb += pm.L1[k], k++ )
        for( long int j = jb; j < jb + pm.L1[k] && j < L; j++ ) phaseOf[(size_t)j] = k;
    double total = 0.;
    for( long int k = 0; k < FI; k++ ) if( pm.XF[k] > 0. ) total += pm.XF[k];
    if( !( total > 0. ) ) return;

    auto trimmed = []( std::string s ) { s.erase( s.find_last_not_of( " \t" ) + 1 ); return s; };
    std::string report;
    long int nWarn = 0;
    for( long int i = 0; i < N; i++ )
    {
        if( !( pm.B[i] > 0. ) ) continue;
        if( pm.ICC && ( pm.ICC[i] == IC_CHARGE || pm.ICC[i] == IC_VOLUME ) ) continue;
        long int k = -1; bool confined = true; double carried = 0.;
        for( long int j = 0; j < L && confined; j++ )
        {
            if( pm.A[i + j*N] == 0. ) continue;
            if( pm.DUL && pm.DUL[j] < 1e6 && !( pm.DUL[j] > 0. ) ) continue;    // excluded species
            const long int kj = phaseOf[(size_t)j];
            if( k < 0 ) k = kj; else if( kj != k ) confined = false;
            if( pm.X[j] > 0. ) carried += pm.X[j];
        }
        if( !confined || k < 0 || pm.L1[k] < 2 ) continue;
        const double nk = pm.XF[k] > 0. ? pm.XF[k] : 0.;
        const double share = nk > 0. ? carried / nk : 1.;
        const double trace = nk / total;
        // icName is compared against elementsOfInterest, so it stays the trimmed name; only the
        // copy that goes into the DECIDE payload is made separator-safe.
        const std::string icName = trimmed( char_array_to_string( pm.SB[i], MAXICNAME ) );
        const std::string icNameSafe = trace_safe_name( icName );
        const std::string phName = trace_safe_name( char_array_to_string( pm.SF[k] + MAXSYMB, MAXPHNAME ) );
        if( probe )
            fprintf( stderr, "STRANDPROBE ic=%s phase=%s class=%c b_internal=%.3e phase_mol=%.3e share=%.3e trace=%.3e\n",
                     icName.c_str(), phName.c_str(), pm.PHC ? pm.PHC[k] : '?', pm.B[i],
                     pm.SizeFactor > 0. ? nk / pm.SizeFactor : nk, share, trace );
        if( share >= kStrandedShareWarn && trace <= kStrandedTraceWarn )
        {
            // nk is in pa_DG's internal scale here (GibbsEnergyMinimization runs rescaled); report real moles
            const double nkReal = pm.SizeFactor > 0. ? nk / pm.SizeFactor : nk;
            // with elements marked of interest, name which kind this is: the remedy differs
            std::string tag;
            if( !elementsOfInterest.empty() )
                tag = std::find( elementsOfInterest.begin(), elementsOfInterest.end(), icName ) != elementsOfInterest.end()
                      ? " [of interest]" : " [default seed]";
            report += fmt::format( "{}{} in {} ({:.2g} mol, {:.0f} % of it){}", report.empty() ? "" : "; ",
                                   icNameSafe, phName, nkReal, 100. * share, tag );
            nWarn++;
        }
    }
    if( !nWarn ) return;
    native_trace_decide( "stranded-element n=%ld list=%s", nWarn, report.c_str() );

    static std::mutex warnedMutex;
    static std::set<std::string> warned;
    std::string key;   // element + phase names only: the amounts change along a sweep
    for( size_t a = 0, b; a < report.size(); a = b + 2 )
    {
        b = report.find( "; ", a ); if( b == std::string::npos ) b = report.size();
        const std::string item = report.substr( a, b - a );
        key += item.substr( 0, item.find( " (" ) ) + ";";
    }
    {
        std::lock_guard<std::mutex> lock( warnedMutex );
        if( !warned.insert( key ).second ) return;
    }
    gems_logger->warn(
        "Fragile system: {} element(s) exist only in one solution phase, which is present here only as a "
        "trace made mostly of that element. Results can be lost or change under tiny input changes, and "
        "pH/Eh/ionic strength are not meaningful. Try: add a phase that can host the element, or remove "
        "the element.{} Element in phase: {}", nWarn,
        elementsOfInterest.empty() ? "" : " An element tagged [default seed] is only a placeholder, so "
        "removing it is the direct fix.", report );
}

void TMultiBase::GibbsEnergyMinimization()
{
  bool IAstatus;
  Reset_uDD( 0L, uDDtrace); // Experimental - added 06.05.2011 KD
  pm.CondNum = 0.;      // reset condition-number diagnostics for this GEM call;
  pm.CondNumDiag = 0.;  // updated to the worst case seen across internal linear solves below
  pm.SolveTimeMs = 0.;      // per-phase timing: reset, accumulated across internal linear
  pm.CondNumTimeMs = 0.;    // solves below (SolveTimeMs) and its diagnostics-only subset
  pm.SolveCallCount = 0;    // (CondNumTimeMs); lets per-system overhead be measured directly,
                            // from a single run, instead of differencing two separate builds

  TNode::ipmlog_file->debug(" GEMIPM TC={}", pm.TCc);

  // One CALL header per equilibrium calculation, so a trace file holding several
  // calls can be split, and so the thresholds every later PSSC line is compared
  // against are recorded once rather than repeated per phase.
  if( FILE* ntf = native_trace_file() )
  {
      const BASE_PARAM* pa0 = base_param();
      fprintf( ntf, "CALL  pNP=%ld N=%ld L=%ld Ls=%ld FI=%ld TK=%.4f Pbar=%.6e"
                    " PC=%d DF=%.3e DFM=%.3e PRD=%d DSM=%.3e DcMinM=%.3e DHBM=%.3e"
                    " DT=%d IIM=%d DP=%d\n",
               (long)pm.pNP, (long)pm.N, (long)pm.L, (long)pm.Ls, (long)pm.FI,
               pm.T, pm.P, (int)pa0->PC, pa0->DF, pa0->DFM, (int)pa0->PRD,
               pm.DSM, pm.DcMinM, pm.DHBM, (int)pa0->DT, (int)pa0->IIM, (int)pa0->DP );
      fflush( ntf );
  }

#ifndef NDEBUG
  // DATABR values as received, before any internal processing (debug level), to tell a
  // caller-side problem from a solver-side one.
  if(gems_logger->should_log(spdlog::level::debug)) {
      gems_logger->debug("GibbsEnergyMinimization() entry - bulk composition B[]:");
      for(int i = 0; i < pm.N; i++)
          gems_logger->debug("  B[{}] = {:.6e}", i, pm.B[i]);
      gems_logger->debug("GibbsEnergyMinimization() entry - DC limits DUL[]/DLL[]:");
      for(int j = 0; j < pm.L; j++)
          gems_logger->debug("  j={} DUL={:.6e} DLL={:.6e}", j, pm.DUL[j], pm.DLL[j]);
  }
#endif

FORCED_AIA:
    GEM_IPM_Init();
   if( pm.pNP )
   {
      if( pm.ITaia <=30 )       // Foolproof
           pm.IT = 30;
       else
           pm.IT = pm.ITaia;  // Setting number of iterations for the smoothing parameter
      Set_DC_limits( true );  // Experimental location for setting AMRs 29.07.15
   }

   IAstatus = GEM_IPM_InitialApproximation( );
   if( FILE* ntf = native_trace_file() )
   {
       long int nz = 0;
       for( long int j = 0; j < pm.L; j++ )
           if( pm.Y[j] == 0. ) nz++;
       fprintf( ntf, "INIT  pNP=%ld IAstatus=%d zeroDC=%ld of %ld\n",
                (long)pm.pNP, (int)IAstatus, (long)nz, (long)pm.L );
       fflush( ntf );
   }
   if( IAstatus == false )
   {
      //Wrapper call for the IPM iteration sequence
      GEM_IPM( -1 );
      if( !pm.pNP )
          pm.ITaia = pm.IT;

       // calculation of demo data for gases
       for( long int ii=0; ii<pm.N; ii++ )
           pm.U_r[ii] = pm.U[ii]*pm.RT;
       GasParcP();  // do we really need it?
   }
   pm.IT = pm.ITG;

   // testing results
   if( pm.MK == 2 )
   {	if( pm.pNP )
        {                     // SIA mode failed
            pm.pNP = 0;
            pm.MK = 0;
            Reset_uDD( 0L, uDDtrace );  // resetting u divergence detector
            goto FORCED_AIA;  // Trying again with AIA set after bad SIA
         }
        else {                 // AIA mode failed
           if( nCNud <= 0L )   // Generic AIA IPM failure, except the case of u divergence
               Error( pm.errorCode ,pm.errorBuf );
           // Now trying again with AIA down to cnr-2 IPM iteration, no PhaseSelection()
           // and possibly cleanup only for species with much lower activity than concentration
           pm.ITG = 0;
           goto FORCED_AIA;
       }
   }
   pm.FitVar[0] = bfc_mass();  // getting total mass of solid phases in the system

   // Exact-zero census of the answer, into the trace: an interior-point method cannot produce
   // an exact zero, so each zero here was written by an explicit removal step (see the PSSC
   // lines).
   if( FILE* ntf = native_trace_file() )
   {
       long int nzX = 0, nzY = 0, nzPh = 0;
       for( long int j = 0; j < pm.L; j++ )
       {
           if( pm.X[j] == 0. ) nzX++;
           if( pm.Y[j] == 0. ) nzY++;
       }
       for( long int k = 0; k < pm.FI; k++ )
           if( pm.XF[k] == 0. ) nzPh++;
       fprintf( ntf, "ANSWER MK=%ld PZ=%ld K2=%ld ITF=%ld ITG=%ld zeroDC_X=%ld zeroDC_Y=%ld"
                     " of %ld zeroPH=%ld of %ld\n",
                (long)pm.MK, (long)pm.PZ, (long)pm.K2, (long)pm.ITF, (long)pm.ITG,
                (long)nzX, (long)nzY, (long)pm.L, (long)nzPh, (long)pm.FI );
       fflush( ntf );
   }


   // ---- Mass-balance check of the answer this call returns: warn only. Native's cold path
   // does not check the state it returns (a cold call's second MBR is exempt from the strict
   // guard), so the fact is reported instead of changing the verdict. Not NDEBUG-gated.
   // Suppressed when pm.MK/pm.PZ already mark the solution as bad. One O(N*L) pass.
   if( !pm.MK && !pm.PZ )
   {
       long int iRel = -1, iAbs = -1; double rel = 0., absr = 0.;
       native_trace_mb_of( pm, pm.X, iRel, rel, iAbs, absr );
       const double dtExp  = fabs( (double)base_param()->DT );
       const double absCut = ( dtExp < 2. ) ? pm.DHBM : pow( 10., -dtExp );
       // Same per-IC rule the convergence branches apply: with DT == 0 only the
       // relative test exists; with DT != 0 an IC must exceed BOTH.
       bool fails = !base_param()->DT ? ( rel > 1. )
                                      : ( rel > 1. && absr > absCut );
       // pa_MbReproject: try to repair the state before reporting it. Restores pm.X unless the
       // worst relative residual strictly improves.
       if( fails && iRel >= 0 && base_param()->MbReproject )
       {
           if( MassBalanceReproject( pm.X ) )
           {
               native_trace_mb_of( pm, pm.X, iRel, rel, iAbs, absr );
               fails = !base_param()->DT ? ( rel > 1. )
                                         : ( rel > 1. && absr > absCut );
           }
       }

       if( fails && iRel >= 0 )
           gems_logger->warn(
               "Mass balance not met: element {} is {:.1e}x its tolerance (worst absolute {:.1e} mol, "
               "element {}). The answer is returned unchanged, but a warm (SIA) re-solve would reject it. "
               "Try: set pa_DT (GEMS: Pa_DPV[2]) to -6 or lower so trace elements are judged by an absolute 10^DT mol, or reduce pa_DHB (GEMS: Pa_DHB).",
               char_array_to_string( pm.SB[iRel], MAXICNAME ), rel, absr,
               iAbs >= 0 ? char_array_to_string( pm.SB[iAbs], MAXICNAME ) : std::string( "-" ) );
   }

   if( !pm.MK && !pm.PZ )
   {
       EnergyDeterminacyCheck();
       StrandedElementCheck();
   }

   if( pm.MK || pm.PZ ) // no good solution
       /*TProfil::pm->*/testMulti();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Main sequence of GEM IPM algorithm implementation.
///  Main place for implementation of diagnostics and setup
///  of IPM precision and convergence
///  rLoop is the index of the primal solution refinement loop (for tracing)
///   or -1 if this is main GEM_IPM call
//
void TMultiBase::GEM_IPM( long int /*rLoop*/ )
{
    long int i, j, eRet, status=0; long int csRet=0;
// bool CleanAfterIPM = true;
    const BASE_PARAM *pa_p = base_param();

#ifdef GEMITERTRACE
  to_text_file( "MultiDumpB.txt" );   // Debugging
#endif

    pm.W1=0; pm.K2=0;         // internal counters and indicators
    pm.Ec = pm.MK = pm.PZ = 0;    // Return codes
    // Per solve, across this call's phase-selection passes (see insBudgetTried).
    insBudgetTried.assign( static_cast<size_t>( pm.FI ), 0 );
    if(!nCNud && !cnr )
        setErrorMessage( 0, "" , "");  // empty error info
 //   if( TProfil::pm->pa.p.PLLG == 0 )  // Disabled by DK 11.05.2011
 //       TProfil::pm->pa.p.PLLG = 20;  // Changed 28.04.2010 KD

//    if( pm.pULR && pm.PLIM )
//        Set_DC_limits( DC_LIM_INIT );
//          Set_DC_limits( false );

mEFD:  // Mass balance refinement (formerly EnterFeasibleDomain())
     native_trace_stage( pm, pa_p, "MBR", "enter", -1, true );
     eRet = MassBalanceRefinement( pm.K2 ); // Here the MBR() algorithm is called
     native_trace_stage( pm, pa_p, "MBR", "exit", eRet, true );

#ifndef NDEBUG
   if(gems_logger->should_log(spdlog::level::debug)) {
       gems_logger->debug("After MassBalanceRefinement (eRet={}) - species amounts:", eRet);
       for(int j = 0; j < pm.L; j++)
           gems_logger->debug("  j={} Y={:.6e} X={:.6e}", j, pm.Y[j], pm.X[j]);
   }
#endif

#ifdef GEMITERTRACE
to_text_file( "MultiDumpC.txt" );   // Debugging
#endif

// STEPWISE (2)  - stop point to examine output from EFD()
STEP_POINT("After FIA");

    switch( eRet )
    {
     case 0:  // OK
         break;
     case 5:  // Initial Lagrange multiplier for metastability broken for DC
     case 4:  // Initial mass balance broken for IC
     case 3:  // too small step length in descent algorithm
     case 2:  // max number of iterations has been exceeded in MassBalanceRefinement()
     case 1: // degeneration in R matrix  for MassBalanceRefinement()
                 if( pm.pNP )
                 {   // bad SIA mode - trying the AIA mode
                pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
                TNode::ipmlog_file->trace("ITF={} ITG={}  IT={}  ! PIA->AIA on E04IPM", pm.ITF, pm.ITG, pm.IT);
                goto FORCED_AIA;
   	         }
   	         else
                         Error( pm.errorCode ,pm.errorBuf );
              break;
    }

   // calling the MainIPMDescent() minimization algorithm
   native_trace_stage( pm, pa_p, "IPM", "enter", -1, true );
   eRet = InteriorPointsMethod( status/*, pm.K2*/ );
   native_trace_stage( pm, pa_p, "IPM", "exit", eRet, true );

#ifndef NDEBUG
   if(gems_logger->should_log(spdlog::level::debug)) {
       gems_logger->debug("After InteriorPointsMethod (eRet={}, status={}) - species amounts:", eRet, status);
       for(int j = 0; j < pm.L; j++)
           gems_logger->debug("  j={} Y={:.6e} X={:.6e}", j, pm.Y[j], pm.X[j]);
   }
#endif

#ifdef GEMITERTRACE
to_text_file( "MultiDumpD.txt" );   // Debugging
#endif

  // STEPWISE (3)  - stop point to examine output from IPM()
   STEP_POINT("After IPM");

// Diagnostics of IPM results
   switch( eRet )
   {
   case 0:  // OK
 #ifndef NDEBUG
       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("Before CalculateActivityCoefficients:");
           for(int j = 0; j < pm.L; j++)
               gems_logger->debug("  j={} Y={:.6e} X={:.6e} W={:.6e} F={:.6e} lnGam={:.6e}",
                                  j, pm.Y[j], pm.X[j], pm.W[j], pm.F[j], pm.lnGam[j]);
       }
#endif

       CalculateActivityCoefficients(LINK_PP_MODE);
#ifndef NDEBUG
       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("After CalculateActivityCoefficients:");
           for(int j = 0; j < pm.L; j++)
               gems_logger->debug("  j={} Y={:.6e} X={:.6e} W={:.6e} F={:.6e} lnGam={:.6e}",
                                  j, pm.Y[j], pm.X[j], pm.W[j], pm.F[j], pm.lnGam[j]);
       }
#endif
       break;
     case 2:  // max number of iterations has been exceeded in InteriorPointsMethod()
     case 1: // degeneration in R matrix  for InteriorPointsMethod()
         if( pm.pNP )
         {   // bad PIA mode - trying the AIA mode
                pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
                TNode::ipmlog_file->trace("ITF={} ITG={}  IT={}   ! PIA->AIA on E06IPM", pm.ITF, pm.ITG, pm.IT);
                goto FORCED_AIA;
         }

         TNode::ipmlog_file->trace("ITF={} ITG={}  IT={}   AIA: DX->1e-4, DHBM->1e-6 on E06IPM", pm.ITF, pm.ITG, pm.IT);
         Error( pm.errorCode, pm.errorBuf );
         break;
     case 3:  // bad CalculateActivityCoefficients() status in SIA mode
     case 4: // Mass balance broken after DualTh recover of DC amounts
         if( pm.pNP )
         {   // bad SIA mode - trying the AIA mode
                pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
	        goto FORCED_AIA;
         }
         Error( pm.errorCode, pm.errorBuf );
         break;
     case 5: // Divergence in dual solution approximation
                // no or only partial cleanup can be done
         pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
         if( pm.pNP )
         {   // bad SIA mode - trying the AIA mode 
              goto FORCED_AIA;
         }
         goto FORCED_AIA; // even if in AIA, start over and go until r-1 then finish and do MBR
         break;
   }

   // Here the entry to new PSSC() module controlled by PC >= 2
   if( pa_p->PC >= 2 && nCNud <= 0 ) // only if there is no divergence in the dual solution
   {
       long int ps_rcode, k_miss, k_unst, cleanupStatus = 1;
       if( pa_p->PC > 2 )
           cleanupStatus = 0; // in this case separate SpeciationCleanup() is called

#ifndef NDEBUG
       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("Before PhaseSelectionSpeciationCleanup - species amounts:");
           for(int i = 0; i < pm.N; i++)
               gems_logger->debug("  x[{}] = {:.6e}", i, pm.X[i]); // or whatever the species array is
       }
#endif
       ps_rcode = PhaseSelectionSpeciationCleanup( k_miss, k_unst, cleanupStatus );

       // pa_MbReproject, second call site: PSSC's speciation cleanup can zero or raise single
       // species without looking at the mass balance and leaves the repair to the next MBR, which
       // on a warm call carries the strict guard. Applied only when no phase was inserted
       // (k_miss < 0); with an insertion MBR must re-equilibrate. A correction that would undo a
       // phase elimination fails the method's own acceptance test and is reverted. PSSC works on
       // pm.Y, so pm.Y is repaired (pm.X and derived values re-synchronised on success).
       // A partial repair is kept: pm.Y here is the start of the next iterations, not the answer.
       if( k_miss < 0 && base_param()->MbReproject )
           MassBalanceReproject( pm.Y, true );

#ifndef NDEBUG
       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("After PhaseSelectionSpeciationCleanup - species amounts:");
           for(int i = 0; i < pm.N; i++)
               gems_logger->debug("  x[{}] = {:.6e}", i, pm.X[i]);
           gems_logger->debug("cleanupStatus={}, k_miss={}, k_unst={}", cleanupStatus, k_miss, k_unst);
       }

       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("After PSSC - Y vs W:");
           for(int j = 0; j < pm.L; j++)
               if(pm.Y[j] > min(pm.lowPosNum, pm.DcMinM))
                   gems_logger->debug("  j={} Y={:.6e} W={:.6e} ratio={:.3f}",
                                      j, pm.Y[j], pm.W[j],
                                      pm.W[j] > 0 ? pm.Y[j]/pm.W[j] : 0.0);
       }
#endif

  // STEPWISE (3)  - stop point to examine output from SpeciationCleanup()
  STEP_POINT("After PSSC()");

       switch( ps_rcode )  // analyzing return code of PSSC()
       {
             case 1:   // IPM solution is final and consistent, no phases were inserted
                       pm.PZ = 0;
                       break;
             case 0:   // some phases were inserted and a new IPM loop is needed
                      TNode::ipmlog_file->trace("ITF={} ITG={} K2={} k_miss={}  k_unst={}  ! (new Selekt loop)",
                                    pm.ITF, pm.ITG, pm.K2, k_miss, k_unst);
                      pm.PZ = 1;
                      goto mEFD;
             default:
             case -1:  // the IPM solution is inconsistent after 5 phase insertion loops
             {
                 std::string pmbuf = std::to_string(k_miss)+ ": ";
                 if(k_miss >=0 )
                     pmbuf += name_for_message(pm.SF[k_miss],20);
                 std::string pubuf = std::to_string(k_unst)+ ": ";
                 if(k_unst >=0 )
                    pubuf += name_for_message(pm.SF[k_unst],20);

                 std::string buf = " Computed phase assemblage remains inconsistent after 5 phase selection loops.\n"
                                    " Problematic phase(s): ";
                 buf +=  pmbuf + "  " + pubuf +"\n";
                 setErrorMessage( 8, "W08IPM: PSSC():", buf.c_str() );
                 if( pm.pNP )
                 {   // bad SIA mode - there are inconsistent phases after 5 attempts. Attempting AIA mode
                         pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
                         TNode::ipmlog_file->trace("ITF={} ITG={}  IT={} k_miss={}  k_unst={} ! PIA->AIA on E08IPM (Selekt)",
                                                  pm.ITF, pm.ITG, pm.IT, k_miss, k_unst);
                         goto FORCED_AIA;
                 }
                 else
                 { pm.PZ = 2; // IPM solution could not be improved in PhaseSelect() -
                                //   therefore, some inconsistent phases remain
                   // return;
                 }
             }
       } // end switch

   }
   else if( pa_p->PC == 1 && nCNud <= 0 ) // Old PhaseSelect() mode PC = 1
   {
       //=================== calling old Phase Selection algorithm =====================
        long int ps_rcode, k_miss, k_unst, RaiseZeroDCs = 0;

        ps_rcode = PhaseSelect( k_miss, k_unst, RaiseZeroDCs );

        if( (ps_rcode == 0 || ps_rcode == 1) && pa_p->PRD != 0 )
            // This block is calling the separate cleanup speciation function
        {
           double AmThExp, AmountThreshold, ChemPotDiffCutoff = 1e-2;
     //      long int eRet;

           AmThExp = abs(pa_p->PRD);
           if( AmThExp < 4.)
               AmThExp = 4.;
           AmountThreshold = pow(10,-AmThExp);
           if( pa_p->GAS > 1e-6 )
                ChemPotDiffCutoff = pa_p->GAS;
           for( j=0; j<pm.L; j++ )
               pm.XY[j]=pm.Y[j];    // Storing a copy of speciation vector
           csRet = SpeciationCleanup( AmountThreshold, ChemPotDiffCutoff );
           if( csRet == 1 || csRet == 2 || csRet == -1 || csRet == -2 )
           {  //  Significant cleanup has been done - mass balance refinement is necessary
              TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
              for( j=0; j<pm.L; j++ )
                  pm.X[j]=pm.Y[j];
              TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
              CalculateConcentrations( pm.X, pm.XF, pm.XFA );  // also ln activities (DualTh)

 // STEPWISE (3)  - stop point to examine output from SpeciationCleanup()
  STEP_POINT("After Cleanup");
           }
           if( csRet != 0 )  {   // Cleanup removed something
      //         for( j=0; j<pm.L; j++ )   // restoring the Y vector
      //             pm.Y[j]=pm.XY[j];
               pm.W1 = 1;
     //          goto mEFD;   // Forced after cleanup ( check what to do with pm.K2 )
           }
        } // end cleanup

        CalculateActivityCoefficients( LINK_PP_MODE);

        switch( ps_rcode )
           {
              case 1:   // IPM solution is final and consistent, no phases were inserted
                        pm.PZ = 0;
                        break;
              case 0:   // some phases were inserted and a new IPM loop is needed
                        TNode::ipmlog_file->trace("ITF={} ITG={}  IT={}  K2={} k_miss={}  k_unst={}  ! (new Selekt loop)",
                                     pm.ITF, pm.ITG, pm.IT, pm.K2, k_miss, k_unst);
                        pm.PZ = 1;
                        goto mEFD;
              default:
              case -1:  // the IPM solution is inconsistent after 3 phase selection loops
              {
                  std::string pmbuf = std::to_string(k_miss)+ ": ";
                  if(k_miss >=0 )
                       pmbuf += name_for_message(pm.SF[k_miss],20);
                  std::string pubuf = std::to_string(k_unst)+ ": ";
                  if(k_unst >=0 )
                      pubuf += name_for_message(pm.SF[k_unst],20);

                  std::string buf = " Computed phase assemblage remains inconsistent after 3 phase selection loops.\n"
                                    " Problematic phase(s): ";
                  buf +=  pmbuf + "  " + pubuf +"\n";
                  setErrorMessage( 8, "W09IPM: Phase Selection:", buf.c_str() );
                  if( pm.pNP )
                  {   // bad SIA mode - there are inconsistent phases after 3 attempts. Attempting AIA mode
                      pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
                      TNode::ipmlog_file->trace("ITF={} ITG={}  IT={} k_miss={}  k_unst={}  ! PIA->AIA on E08IPM (Selekt)",
                                               pm.ITF, pm.ITG, pm.IT, k_miss, k_unst);
                      goto FORCED_AIA;
                  }
                  else
                  { pm.PZ = 2; // IPM solution could not be improved in PhaseSelect() -
                                 //   therefore, some inconsistent phases remain
                    // return;
                  }
              }
        } // end switch
   }  // end old PhaseSelect
/*
   if( pa->p.PRD != 0 && !( pa->p.PC == 2 || pa->p.PC == 1 ) && nCNud <= 0 )
   {    // This block is calling the selarate cleanup speciation function
      double AmThExp, AmountThreshold, ChemPotDiffCutoff = 1e-2;
//      long int eRet;

      AmThExp = (double)abs(pa->p.PRD);
      if( AmThExp < 4.)
          AmThExp = 4.;
      AmountThreshold = pow(10,-AmThExp);
      if( pa->p.GAS > 1e-6 )
           ChemPotDiffCutoff = pa->p.GAS;
      for( j=0; j<pm.L; j++ )
          pm.XY[j]=pm.Y[j];    // Storing a copy of speciation vector
      csRet = SpeciationCleanup( AmountThreshold, ChemPotDiffCutoff );
      if( csRet == 1 || csRet == 2 || csRet == -1 || csRet == -2 )
      {  //  Significant cleanup has been done - mass balance refinement is necessary
         TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
         for( j=0; j<pm.L; j++ )
             pm.X[j]=pm.Y[j];
         TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
         CalculateConcentrations( pm.X, pm.XF, pm.XFA );  // also ln activities (DualTh)

// STEPWISE (3)  - stop point to examine output from SpeciationCleanup()
   STEP_POINT("After Cleanup");

      }
      if( csRet != 0 )  {   // Cleanup removed something
 //         for( j=0; j<pm.L; j++ )   // restoring the Y vector
 //             pm.Y[j]=pm.XY[j];
          pm.W1 = 1;
//          goto mEFD;   // Forced after cleanup ( check what to do with pm.K2 )
      }
   } // end cleanup
*/
   if( pa_p->PC == 3 )
        XmaxSAT_IPM2();  // Install upper limits to xj of surface species (questionable)!

  if( nCNud <= 0 )
  {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
       if(gems_logger->should_log(spdlog::level::debug)) {
           gems_logger->debug("W and F arrays before second MBR:");
           for(int j = 0; j < pm.L; j++)
               if(pm.Y[j] > min(pm.lowPosNum, pm.DcMinM))
                   gems_logger->debug("  j={} Y={:.6e} W={:.6e} F={:.6e}", j, pm.Y[j], pm.W[j], pm.F[j]);
       }
#endif
    native_trace_stage( pm, pa_p, "MBR", "enter", -1, true );
    eRet = MassBalanceRefinement( pm.K2 ); // Mass balance improvement in all normal cases
    native_trace_stage( pm, pa_p, "MBR", "exit", eRet, true );
    switch( eRet )
    {
      case 0:  // OK - refinement of concentrations and activity coefficients
        for( j=0; j<pm.L; j++ )
           pm.X[j]=pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );  // also ln activities (DualTh)
        if( pm.PD >= 2 )
        {
           CalculateActivityCoefficients( LINK_UX_MODE);
        }
        break;
     case 5:  // Cleaned-up Lagrange multiplier for metastability broken for DC
     case 4:  // Cleaned-up mass balance broken for IC
     case 3:  // too small step length in MB refinement algorithm after cleanup
     case 2:  // max number of iterations has been exceeded in MassBalanceRefinement()
     case 1: // degeneration of R matrix in MassBalanceRefinement() after cleanup
             if( pm.pNP )
             {   // bad SIA mode - trying the AIA mode
                pm.MK = 2;   // Set to check in CalculateEquilibriumState() later on
                TNode::ipmlog_file->trace("ITF={} ITG={}  IT={}  ! PIA->AIA on E04IPM", pm.ITF, pm.ITG, pm.IT);
                goto FORCED_AIA;
             }
             else
               Error( pm.errorCode ,pm.errorBuf );
          break;
   }
 }

   pm.t_end = clock();
   pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
// STEPWISE (4) Stop point after PhaseSelect()
   STEP_POINT("Before Refine()");
//   if( pm.MK == 2 )
//       goto FORCED_AIA;

   TNode::ipmlog_file->trace("ITF={} ITG={}  IT={} MBPRL={}  {}",
       pm.ITF, pm.ITG, pm.IT, pm.W1, ( pm.pNP ? " Ok after PIA": " Ok after AIA"));

FORCED_AIA:   // Finish
   pm.FI1 = 0;  // Recomputing the number of non-zeroed-off phases
   pm.FI1s = 0;
   for( i=0; i<pm.FI; i++ )
   {
       if( pm.YF[i] >= min( pm.PhMinM, 1e-22 ) )  // Check 1e-22 !!!!!
       {
            pm.FI1++;
            if( i < pm.FIs )
                pm.FI1s++;
       }
   }
   for( i=0; i<pm.L; i++)
      pm.G[i] = pm.G0[i];
   // At pm.MK == 1, normal return after successful improvement of mass balance precision
   pm.t_end = clock();   // Fix pure runtime
   pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);

#ifdef GEMITERTRACE
to_text_file( "MultiDumpE.txt" );   // Debugging
#endif
}

// ------------------------------------------------------------------------------------------------------
/// Finding out whether the automatic initial approximation is necessary for
/// launching the IPM algorithm.
/// Uses a modified simplex method with two-side constraints (Karpov ea 1997)
/// \return
/// false - OK for IPM
/// true  - OK solved (pure phases only in the system)
//
bool TMultiBase::GEM_IPM_InitialApproximation(  )
{
    long int i, j, k, NN, eCode=-1L;
    double minB;//, sfactor;
    const BASE_PARAM *pa_p = base_param();

#ifdef GEMITERTRACE
to_text_file( "MultiDumpA.txt" );   // Debugging
#endif

// Scaling the IPM numerical controls for the system total amount and minimum b(IC)
    NN = pm.N - pm.E;    // Charge is not checked!
    minB = pm.B[0]; // pa->p.DB;
    for(i=0;i<NN;i++)
    {
        if( pm.B[i] < pa_p->DB )
        {
           if( eCode < 0  )
    	   {
              eCode = i;  // Error state is activated
              pm.PZ = 3;
              std::string buf = "Too small input amount of independent component ";
                          buf += name_for_message(pm.SB[i],3)+" = "+std::to_string(pm.B[i])
                          +" mol. Try: raise its bulk amount, or remove it from the system.";
              setErrorMessage( 20, "W20IPM: IPM Main Descent:", buf.c_str());
    	   }
           else
           {
              addErrorMessage((std::string(", ")+name_for_message(pm.SB[i],3)+" = "+std::to_string(pm.B[i])).c_str());
           }
           pm.B[i] = pa_p->DB;
        }
        if( pm.B[i] < minB )
           minB = pm.B[i];      // Looking for the smallest IC amount
    }

    if( eCode >= 0 )
    {
       /*TProfil::pm->*/testMulti();
       pm.PZ = 0;
       setErrorMessage( -1, "" , "");
    }

    /*sfactor =*/ RescaleToSize( false ); //  replacing calcSfactor();

   bool AllPhasesPure = true;   // Added by DK on 09.03.2010
   // checking if all phases are pure
   for( k=0; k < pm.FI; k++ )
       if( pm.L1[k] > 1 )
           AllPhasesPure = false;
   if( AllPhasesPure == true )  // Provisional
       pm.pNP = 0;  // call SolveSimplex() also in SIA mode!

   // Analyzing if the Simplex LPP approximation is necessary
    if( !pm.pNP  )
    {
        // Preparing to call SolveSimplex() - "cold start"
//    	pm.FitVar[4] = 1.0; // Debugging: no smoothing
//        pm.FitVar[4] = pa->p.AG;  //  initializing the smoothing parameter
        pm.ITaia = 0;             // resetting the previous number of AIA iterations
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        //      pm.IC = 0.0;  For reproducibility of simplex LPP-based IA?
        pm.PCI = 1.0;
        pm.logCDvalues[0] = pm.logCDvalues[1] = pm.logCDvalues[2] = pm.logCDvalues[3] =
              pm.logCDvalues[4] = log( pm.PCI );  // reset CD sampler array

        // Cleaning vectors of activity coefficients
        for( j=0; j<pm.L; j++ )
        {
            if( noZero(pm.lnGmf[j]) )
                pm.lnGam[j] = pm.lnGmf[j]; // setting up fixed act.coeff. for SolveSimplex()
            else pm.lnGam[j] = 0.;
            pm.Gamma[j] = 1.;
        }
        for( j=0; j<pm.L; j++)
            pm.G[j] = pm.G0[j] + pm.lnGam[j];  // Provisory cleanup 4.12.2009 DK
        if( pm.LO )
        {
           CalculateConcentrations( pm.X, pm.XF, pm.XFA );  // cleanup for aq phase?
           pm.IC = 0.0;  // Important for the LPP-based AIA reproducibility
           if( pm.E && pm.FIat > 0 )
           {
              for( k=0; k<pm.FIs; k++ )
              {
                 long int ist;
                 if( pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION )
                     for( ist=0; ist<pm.FIat; ist++ ) // loop over surface types
                     {
                        pm.XetaA[k][ist] = 0.0;
                        pm.XetaB[k][ist] = 0.0;
                        pm.XpsiA[k][ist] = 0.0;
                        pm.XpsiB[k][ist] = 0.0;
                        pm.XpsiD[k][ist] = 0.0;
                        pm.XcapD[k][ist] = 0.0;
                     }  // ist
              }  // k
           } // FIat
        }   // LO
//        if( pa->p.PC == 2 )
//           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
        pm.IT = 0; pm.ITF += 1; // Assuming SolveSimplex() time equal to one iteration of MBR()
//        pm.PCI = 0.0;
     // Calling the simplex LP approximation here
        AutoInitialApproximation( );
// experimental 15.03.10 (probably correct)
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        for( j=0; j< pm.L; j++ )
            pm.X[j] = pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );

        // Calculation of mass-balance residuals and DC concentrations in phases
        MassBalanceResiduals( pm.N, pm.L, pm.A, pm.X, pm.B, pm.C);
        CalculateConcentrations( pm.X, pm.XF, pm.XFA );  // also ln activities (DualTh)
//  STEPWISE (0) - stop point for examining results from LPP-based IA
STEP_POINT( "End Simplex" );
        if( AllPhasesPure )     // bugfix DK 09.03.2010   was if(!pm.FIs)
        {                       // no multi-component phases!
            pm.W1=0; pm.K2=0;               // set internal counters
            pm.Ec = pm.MK = pm.PZ = 0;
#ifdef GEMITERTRACE
to_text_file( "MultiDumpLP.txt" );   // Debugging
#endif

   pm.t_end = clock();
   pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
           pm.FI1 = 0;
           pm.FI1s = 0;
           for( i=0; i<pm.FI; i++ )
           if( pm.YF[i] > 1e-18 )
           {
             pm.FI1++;
             if( i < pm.FIs )
                pm.FI1s++;
           }
           return true; // If so, the GEM problem is already solved !
        }
        // Setting default trace amounts to DCs that were zeroed off
        // pa_FilloutBudget needs the simplex solution (which species the LP zeroed), and
        // DC_RaiseZeroedOff() overwrites it; snapshot only when the budget is on.
        std::vector<double> yLpFill;
        if( FilloutBudgetValue() > 0. )
            yLpFill.assign( pm.Y, pm.Y + pm.L );
        DC_RaiseZeroedOff( 0, pm.L );
        // pa_FilloutBudget: the class constants can ask for more of an element than the system
        // holds. Applied after any fill-out, so it bounds whatever was written.
        ApplyFilloutBudget( yLpFill );
        // this operation greatly affects the accuracy of mass balance!
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        for( j=0; j< pm.L; j++ )
            pm.X[j] = pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        //        if( pa->p.PC == 2 )
        //           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
        if( pm.PD == 2 /* && pm.Lads==0 */ )
        {
            pm.FitVar[4] = -1.0;   // To avoid smoothing when F0[j] is calculated first time
            CalculateActivityCoefficients( LINK_UX_MODE);
            pm.FitVar[4] = 1.0;
        }
#ifdef GEMITERTRACE
to_text_file( "MultiDumpAA.txt" );   // Debugging
#endif
    }
    else  // Taking previous GEMIPM result as an initial approximation
    {
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        for( j=0; j< pm.L; j++ )
            pm.X[j] = pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );

pm.PCI = 1.; // SD 05/05/2010 for smaller number of iterations for systems with adsorbtion

         pm.logCDvalues[0] = pm.logCDvalues[1] = pm.logCDvalues[2] = pm.logCDvalues[3] =
         pm.logCDvalues[4] = log( pm.PCI );  // reset CD sampler array
     if( pm.PD >= 2 /* && pm.Lads==0 */ )
        {
                pm.FitVar[4] = -1.0;   // To avoid smoothing when F0[j] is calculated first time
                CalculateActivityCoefficients( LINK_UX_MODE );
                pm.FitVar[4] = 1.0;
                // if( pm.PD >= 3 )
                   // CalculateActivityCoefficients( LINK_PP_MODE );  // Temporarily disabled (DK 06.07.2009)
        }

        if( pm.pNP <= -1 )
        {  // With raising species and phases zeroed off by SolveSimplex()
           // Setting default trace amounts of DCs that were zeroed off
           DC_RaiseZeroedOff( 0, pm.L );
        }
     }

// STEPWISE (1) - stop point to see IA from old solution or raised LPP IA
STEP_POINT("Before FIA");

    return false;
}

// ------------------- ------------------ ----------------
/// Calculation of a feasible IPM approximation, refinement of the mass balance.
//
/// Algorithm: see Karpov, Chudnenko, Kulik 1997 Amer.J.Sci. vol 297 p. 798-799
/// (Appendix B)
//
/// Control: MaxResidualRatio, 0 (deactivated), > DHBM and < 1 - accuracy for
///     "trace" independent components (max residual for i should not exceed
///     B[i]*MaxResidualRatio)
//
/// \param WhereCalledFrom, 0 - at entry after automatic LPP-based IA;
///                         1 - at entry in SIA (start without SolveSimplex()
///                         2 - after post-IPM cleanup
///                         3 - additional (after PhaseSelection)
/// \return  0 -  OK,
///          1 -  no SLE solution at the specified precision pa.p.DHB
///          2  - used up more than pa.p.DP iterations
///          3  - too small step length (< 1e-6), no descent possible
///          4  - error in Initial mass balance residuals (debugging)
///          5  - error in MetastabilityLagrangeMultiplier() (debugging)
//
long int TMultiBase::MassBalanceRefinement( long int WhereCalledFrom )
{
    long int IT1;
    long int I, J, Z,  N, sRet, iRet=0, j, jK;
    double LM;//, pmp_PCI;
    const BASE_PARAM *pa_p = base_param();

    ErrorIf( !pm.MU || !pm.W, "MassBalanceRefinement()",
                              "Error of memory allocation for pm.MU or pm.W." );

    // calculation of total mole amounts of phases
    TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );

//    if( pm.PLIM )
//        Set_DC_limits(  DC_LIM_INIT );

    // Adjustment of primal approximation according to kinetic constraints
    // Now returns <0 (OK) or index of DC that caused a problem
    jK = MetastabilityLagrangeMultiplier();
    if( jK >= 0 )
    {  // Experimental
        std::string buf = "(EFD("+std::to_string(WhereCalledFrom);
                    buf += ")) Invalid initial Lagrange multiplier for metastability-constrained DC ";
                    buf += name_for_message( pm.SM[jK], MAXDCNAME);
        setErrorMessage( 17, "E17IPM: Mass Balance Refinement: ", buf.c_str());
        native_trace_mbr_exit( pm, pa_p, "metastability", WhereCalledFrom, -1, 5, false );
        return 5;
    }


//----------------------------------------------------------------------------
// BEGIN:  main loop
    // AbsMbCutoff_stall computed once before the loop —
    // mirrors the convergence check logic for the combined mode
    double AbsMbCutoff_stall = pm.DHBM;
    // Track previous max residual to detect oscillation
    double prev_maxResidual = std::numeric_limits<double>::max();
    long int stalledIter = 0;

    // Best-residual tracking. Every early exit below (LM too small, stall detector, pa_p->DP
    // exhausted) returns whatever pm.Y holds at that moment, which need not be the best state
    // MBR visited. bestY snapshots pm.Y at the iteration with the smallest ResidualMetric(),
    // before that iteration's update. It can only replace the returned state by one with a
    // strictly smaller residual, and only on exits that are already non-ideal.
    // In plain words: if the mass-balance step gives up, return the best point it found, not
    // the last one.
    std::vector<double> bestY( pm.L );
    double bestResidual = std::numeric_limits<double>::max();
    bool haveBest = false, trRestored = false;

    // Worst normalised mass-balance residual over all N ICs (> 1: some IC is outside its
    // tolerance), with the same tolerance logic as the per-IC checks below but always over
    // the full range, so states from different iterations are comparable.
    auto ResidualMetric = [&]() -> double
    {
        double worst = 0.;
        for( long int ii=0; ii<pm.N; ii++ )
        {
            double denom = pm.B[ii] * pm.DHBM;
            if( pa_p->DT )
                denom = std::max( denom, AbsMbCutoff_stall );
            if( denom <= 0. )
                denom = ( pm.DHBM > 0. ? pm.DHBM : 1e-16 );
            worst = std::max( worst, fabs(pm.C[ii]) / denom );
        }
        return worst;
    };

    for( IT1=0; IT1 < pa_p->DP; IT1++, pm.ITF++ )
    {
        // get size of task
        pm.NR=pm.N;
        if( pm.LO )
        {   if( pm.YF[0] < pm.DSM && pm.YFA[0] < pm.XwMinM ) // fixed 30.08.2009 DK
                 pm.NR= pm.N-1;
        }
        N=pm.NR;
       // Calculation of mass-balance residuals in IPM
       MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C);
       // Testing mass balance residuals
       Z = pm.N - pm.E;
       if( pa_p->MbClassRule > 0. )
       {
           // pa_MbClassRule (default 0 takes the branches below unchanged): a relative
           // threshold for trace ICs and an absolute one for major ICs, classified by ratio to
           // the largest IC (scale-invariant under pa_DG rescaling).
           double maxB = 0.;
           for( I=0; I<Z; I++ )
               if( pm.B[I] > maxB ) maxB = pm.B[I];
           double AbsMbCutoff_cls;
           {
               const double e = fabs( (double)pa_p->DT );
               // |DT| >= 2 names the major absolute cutoff; otherwise fall back so the switch
               // works alone.
               AbsMbCutoff_cls = ( e >= 2. ) ? pow( 10., -e ) : pm.DHBM * 1e5;
               if( e >= 2. ) AbsMbCutoff_stall = AbsMbCutoff_cls;
           }
           const double traceB = pa_p->MbClassRule * maxB;
           for( I=0; I<Z; I++ )
           {
               const bool isTrace = ( pm.B[I] < traceB );
               if( isTrace ? ( fabs(pm.C[I]) > pm.B[I] * pm.DHBM )
                           : ( fabs(pm.C[I]) > AbsMbCutoff_cls ) )
                   break;
           }
       }
       else if( !pa_p->DT )
       {   // relative balance accuracy for all ICs
           for( I=0;I<Z;I++ )
             if( fabs(pm.C[I]) > pm.B[I] * pm.DHBM )
               break;
       }
       else { // combined balance accuracy - an absolute floor under the relative test
           // An IC counts as not converged only if it exceeds both thresholds, i.e. the bar is
           // max(absolute, relative) per IC: a trace IC is judged on the absolute cutoff, a
           // major one keeps its relative bar. Only projects that set pa_DT take this branch
           // (|DT| >= 2 names the absolute cutoff as 10^-|DT|). This makes the test physical;
           // it does not remove a degeneracy the test is detecting.
           double AbsMbAccExp, AbsMbCutoff;
           AbsMbAccExp = abs( pa_p->DT );
           if( AbsMbAccExp < 2. )  // If DT is set to 1 or -1 then DHBM is used also as the absolute cutoff
               AbsMbCutoff = pm.DHBM;
           else
           {
               AbsMbCutoff = pow( 10, -AbsMbAccExp );
               AbsMbCutoff_stall = pow( 10., -AbsMbAccExp );
           }
           for( I=0;I<Z;I++ )
              if( fabs( pm.C[I]) > AbsMbCutoff && fabs(pm.C[I]) > pm.B[I] * pm.DHBM )
                  break;
       }
       // Best-residual snapshot: pm.C[] is this iteration's residual for the pre-update pm.Y.
       {
           double curResidual = ResidualMetric();
           if( curResidual < bestResidual )
           {
               bestResidual = curResidual;
               for( j=0; j<pm.L; j++ )
                   bestY[j] = pm.Y[j];
               haveBest = true;
           }
       }

       if( I == Z ) // balance residuals OK
       { // very experimental - updating activity coefficients
           for( j=0; j< pm.L; j++ )
               pm.X[j] = pm.Y[j];
           TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
           //        if( pa->p.PC == 2 )
           //           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
           if( pm.PD == 2 /* && pm.Lads==0 */ )
           {
//               CalculateConcentrations( pm.X, pm.XF, pm.XFA );
               CalculateActivityCoefficients( LINK_UX_MODE);
           }
           if(iRet==1) {
               iRet=0;  // no SLE solution on internal iterations SD 02/2026
           }
           native_trace_mbr_exit( pm, pa_p, "converged", WhereCalledFrom, IT1, iRet, trRestored );
           return iRet;       // mass balance refinement finished OK
       }

       WeightMultipliers( true );  // creating R matrix

       // Assembling and solving the system of linearized equations
       sRet = MakeAndSolveSystemOfLinearEquations( N, true );

       if( sRet == 1 )  // error: no SLE solution!
       {
           iRet = 1;
           std::string buf = "(EFD("+std::to_string(WhereCalledFrom)+ ")";
           buf += "Degeneration in R matrix (fault in SLE solver).\n"
                  "Mass balance cannot be improved, not possible to proceed.";
           setErrorMessage( 5, "E05IPM: Mass Balance Refinement: " , buf.c_str() );
           gems_logger->warn(buf);
       }
       else
           iRet = 0;  // reset: matrix was singular in a previous iteration but recovered —
               // stale Uefd was used for MU descent direction in the failed iteration,
               // subsequent iterations may succeed as Y evolves away from singularity

      // SOLVED: solution of linear matrix has been obtained
         //          pm.PCI = calcDikin( N, true);
      /*pmp_PCI =*/ DikinsCriterion( N, true);  // calc of MU values and Dikin criterion

      LM = StepSizeEstimate( true ); // Estimation of the MBR() iteration step size LM

      if( LM < 1e-6 ) //SD 1e-6 was
      {  // Experimental
          iRet = 3;
          std::string buf = "(MBR("+std::to_string(WhereCalledFrom);
                      buf += ")): Too small LM step size - cannot converge (check pa_DG, GEMS: Pa_DG)";
          setErrorMessage( 3, "E03IPM: Mass Balance Refinement", buf.c_str() );
          break;
       }
      if( LM > 1.)
         LM = 1.;
      ipm_logger->trace("LM {}", LM);

      // calculation of new primal solution approximation
      // from step size LM and the gradient vector MU
      for(J=0;J<pm.L;J++)
            pm.Y[J] += LM * pm.MU[J];

      // Stall detection: exit early if residuals are no longer improving
      // for 3 consecutive iterations. Gated behind pa_p->PSTALL (default 1,
      // preserving prior behavior) -- with it disabled, a stalled run no
      // longer exits early and instead keeps iterating until pa_p->DP is
      // exhausted, surfacing as a real "Maximum allowed number of MBR
      // iterations exceeded" failure below instead of a silent iRet=0.
      if( pa_p->PSTALL )
      {
          double maxDeltaY = 0.0;
          for( J=0; J<pm.L; J++ )
              maxDeltaY = std::max( maxDeltaY, fabs(pm.MU[J]) * LM );

          // Stall detection: compute max residual across all failing active ICs
          // Start from I — ICs before I already passed the convergence check
          double cur_maxResidual = 0.0;
          bool any_failing = false;
          bool degenerate_cause = false;
          for( long int II=I; II<pm.N; II++ )
          {
              bool ic_failing;
              double ic_tol = pm.B[II] * pm.DHBM;
              if( !pa_p->DT )
                  ic_failing = fabs(pm.C[II]) > ic_tol;
              else
                  ic_failing = fabs(pm.C[II]) > AbsMbCutoff_stall ||
                               fabs(pm.C[II]) > ic_tol;
              if( ic_failing )
              {
                  if( pm.B[II] >= pm.DcMinM )
                  {
                      // Non-negligible IC is failing — not purely degenerate
                      any_failing = true;
                      cur_maxResidual = std::max( cur_maxResidual, fabs(pm.C[II]) );
                      degenerate_cause = false;  // permanently cleared
                  }
                  else if( !any_failing )
                  {
                      // Negligible IC failing, no non-negligible failures yet
                      degenerate_cause = true;
                  }

                  //if( pm.B[II] < 1e-15 ) // need to remove hardcoded value
                  //{
                  //  gems_logger->warn("MBR({}): B[{}]={:.3e}", WhereCalledFrom, II, pm.B[II]);
                  //  iRet = 0; // treat as degenerate but not physically negligible, to avoid triggering AIA fallback
                  //}
              }
          }

          bool stalled = false;
          if( any_failing )
          {
              // Stall if residuals are not improving AND steps are tiny
              // relative to the residual magnitude
              if( maxDeltaY < cur_maxResidual * pm.DHBM )
                  stalled = true;
              // Also stall if residuals are oscillating (not decreasing)
              if( cur_maxResidual >= prev_maxResidual * 0.999 )
                  stalledIter++;
              else
              {
                  stalledIter = 0;  // residuals improving, reset
                  prev_maxResidual = cur_maxResidual;
              }
          }

          if( stalled || stalledIter >= 10 )
          {
              gems_logger->debug("MBR({}): stall at IT1={} stalledIter={} "
                                     "maxDeltaY={:.3e} curRes={:.3e} prevRes={:.3e}",
                                     WhereCalledFrom, IT1, stalledIter,
                                     maxDeltaY, cur_maxResidual, prev_maxResidual);
              // Reset iRet if stall is caused only by physically degenerate ICs:
              // - degeneracy is due to negligible bulk amount, not numerical failure
              // - active ICs have converged sufficiently
              // - returning iRet=1 would incorrectly trigger AIA fallback
              // if( iRet == 1 && degenerate_cause )
                  iRet = 0;
              break;
          }
      }

// STEPWISE (5) Stop point at end of iteration of FIA()
STEP_POINT("FIA Iteration");

}  /* End loop on IT1 */
//----------------------------------------------------------------------------

    // Best-residual restore. Every path reaching this point is a non-ideal exit (LM too small,
    // stall detector, or pa_p->DP exhausted). If the in-loop snapshot has a strictly smaller
    // residual than the state about to be returned, put it back.
    if( haveBest )
    {
        MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C );
        double curResidualFinal = ResidualMetric();
        if( bestResidual < curResidualFinal )
        {
            gems_logger->debug("MBR({}): restoring best-residual state (best={:.3e} vs current={:.3e})",
                               WhereCalledFrom, bestResidual, curResidualFinal);
            for( j=0; j<pm.L; j++ )
                pm.Y[j] = bestY[j];
            MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C );
            trRestored = true;
        }
    }

    //  Prescribed mass balance precision cannot be reached
                    // Temporary workaround for pathological systems 06.05.2010 DK
   const bool dpExhaustedWithoutStallDetection = !pa_p->PSTALL && IT1 == pa_p->DP;
   if( (pa_p->DW && ( WhereCalledFrom == 0L || pm.pNP )) ||
       dpExhaustedWithoutStallDetection )  // DW behavior plus mandatory failure on PSTALL=0 exhaustion
   {
       iRet = 2;
       std::string buf = "(MBR("+std::to_string(WhereCalledFrom);
                   buf += ")) Maximum allowed number of MBR iterations (";
                   buf += std::to_string(pa_p->DP) +") exceeded. Try: pa_DT (GEMS: Pa_DPV[2]) = -6 or lower for trace elements, or check the species stoichiometry.";
       setErrorMessage( 4, "E04IPM: Mass Balance Refinement: ", buf.c_str());
       native_trace_mbr_exit( pm, pa_p, "budget_strict", WhereCalledFrom, IT1, iRet, trRestored );
       return iRet; // no MBR() solution
   }
   // very experimental - updating activity coefficients after MBR()
   for( j=0; j< pm.L; j++ )
         pm.X[j] = pm.Y[j];
   TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
              //        if( pa->p.PC == 2 )
              //           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
   if( pm.PD == 2 /* && pm.Lads==0 */ )
   {
   //        CalculateConcentrations( pm.X, pm.XF, pm.XFA );
         CalculateActivityCoefficients( LINK_UX_MODE);
   }
   native_trace_mbr_exit( pm, pa_p, "nonideal_lenient", WhereCalledFrom, pa_p->DP, iRet, trRestored );
   return iRet;   // inaccurate MBR() solution
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of chemical equilibrium using the Interior Points
///  Method algorithm (see Karpov et al., 1997, p. 785-786)
///  GEM IPM
/// \return  0, if converged;
///          1, in the case of R matrix degeneration
///          2, (more than max iteration) - no convergence
///              or user's interruption
///          3, CalculateActivityCoefficients() returns bad (non-zero) status
///          4, Mass balance broken  in DualTH (Mol_u)
///          5, Divergence in dual solution u vector has been detected
//
long int TMultiBase::InteriorPointsMethod( long int &status/*, long int rLoop*/ )
{
    bool StatusDivg;
    long int N, IT1,J,Z,iRet,i,  nDivIC;
    double LM=0., LM1=1., FX1,    DivTol;
    // Noise-stall accept, pa_IpmStallWindow (see the field in ms_multi.h). All signals are
    // required, and none references pa_DK; each alone was found unsafe.
    const long int kIpmStallMaxW = 200;
    const double kIpmStallFXTol   = 1.e-9;  // energy flat over the window
    const double kIpmStallCompTol = 1.e-6;  // sumX and maxX flat over the window
    const double kIpmStallSpRel   = 1.e-3;  // every species flat RELATIVE TO ITSELF,
    const double kIpmStallSpNegl  = 1.e-9;  //   unless it is this small a share of the total
    const double kIpmStallIncLo   = 0.35;   // PCI increases on 35-65 % of steps,
    const double kIpmStallIncHi   = 0.65;   //   i.e. it is bouncing, not moving
    std::vector<double> stW_pci, stW_fx, stW_sum, stW_max;
    // One row of pm.X[] per iteration, held in a ring buffer.
    std::vector<double> stW_xf;
    long int stW_n = 0, stW_head = 0;
    const BASE_PARAM *pa_p = base_param();
    ipmKktRescues = 0;
    s_ipmKktFallbacks = 0;

    status = 0;
    if( pm.FIs )
      for( J=0; J<pm.Ls; J++ )
            pm.lnGmo[J] = pm.lnGam[J];

    pm.FX=GX( LM  );  // calculation of G(x)

    if( pm.FIs ) // multicomponent phases are present
      for(Z=0; Z<pm.FIs; Z++)
        pm.YFA[Z]=pm.XFA[Z];

//----------------------------------------------------------------------------
//  Main loop of IPM iterations
    for( IT1 = 0; IT1 < pa_p->IIM; IT1++, pm.IT++, pm.ITG++ )
    {
        StatusDivg = false;
        pm.NR=pm.N;
        if( pm.LO ) // water-solvent is present
        {
            if( pm.YF[0]<pm.DSM && pm.Y[pm.LO]< pm.XwMinM )  // fixed 30.08.2009 DK
                pm.NR=pm.N-1;
        }
        N = pm.NR;

#ifdef GEMITERTRACE
to_text_file( "MultiDumpDC1.txt" );   // Debugging
#endif

        PrimalChemicalPotentials( pm.F, pm.Y, pm.YF, pm.YFA );

        // Saving previous content of the U vector to Uc vector
        if( pm.PCI <= pm.DXM * 10. ) // only at low enough Dikin criterion values
        {
           for(J=0;J<pm.N;J++)
              pm.Uc[J][0] = pm.U[J];
        }
        // Setting weight multipliers for DC
        WeightMultipliers( false );

        // Making and solving the R matrix of IPM linearized equations
        iRet = MakeAndSolveSystemOfLinearEquations( N, false );
        if( iRet == 1 )
        {
            setErrorMessage( 7, "E07IPM: IPM Main Descent: ",
   " Degeneration in R matrix (fault in the linearized system solver).\n"
   " It is not possible to obtain a valid GEM IPM solution.\n"
   " Try: check the system for duplicate species or elements that no species contains.\n" );
          return 1;
        }

   if( !nCNud && base_param()->PLLG )   // disabled if PLLG = 0
   { // Experimental - added 06.05.2011 by DK
      Increment_uDD( pm.ITG, uDDtrace );
//   DivTol = pow( 10., -fabs( (double)TProfil::pm->pa.p.PLLG ) );
      DivTol = (double)base_param()->PLLG;
      if( fabs(DivTol) >= 30000. )
          DivTol = 1e6;  // this is to allow complete tracing in the case of divergence
//      if( pm.ITG )
//          DivTol /= pm.ITG;
//       DivTol -= log(pm.ITG);
//      if( DivTol < 0.3 )
//          DivTol = 0.3;
      // Checking the dual solution for divergence
      if(DivTol < 0. )
         nDivIC = Check_uDD( 0, -DivTol, uDDtrace );
      else
         nDivIC = Check_uDD( 1, DivTol, uDDtrace );

      if( nDivIC )
      { // Printing error message
        char buf[512];
          memset(buf, '\0', 512);
        StatusDivg = true;
        sprintf( buf, "Divergence in dual solution approximation (u) \n at IPM iteration %ld with gen.tolerance %g "
             "for %ld ICs:   %-6.5s", pm.ITG, DivTol, nDivIC, pm.SB[ICNud[0]] );
        setErrorMessage( 14, "W14IPM: IPM Main Descent:", buf);
        for( Z =1; Z<nCNud; Z++ ) {
            addErrorMessage((std::string(" ")+name_for_message(pm.SB[ICNud[Z]],6)).c_str());
        }
      }
   }

// Got the dual solution u vector - calculating the Dikin's Criterion of IPM convergence
   pm.PCI = DikinsCriterion( N, false );

#ifdef GEMITERTRACE
to_text_file( "MultiDumpDC.txt" );   // Debugging
#endif

       if( StatusDivg )
           return 5L;

       // Initial estimate of IPM descent step size LM
       LM = StepSizeEstimate( false );
       LM1 = OptimizeStepSize( LM ); // Finding an optimal value of the descent step size
       FX1 = GX( LM1 ); // Calculation of the total Gibbs energy of the system G(X)
                          // and copying of Y, YF vectors into X,XF, respectively.
       pm.FX=FX1;
       // temporary
       for(i=4; i>0; i-- )
            pm.logCDvalues[i] = pm.logCDvalues[i-1];
       pm.logCDvalues[0] = log( pm.PCI );  // updating CD sampler array

       if( pm.PHC[0] == PH_AQUEL && ( pm.XF[0] < pm.DSM ||
            pm.X[pm.LO] <= pm.XwMinM ))    // fixed 28.04.2010 DK
       {
           pm.XF[0] = 0.;  // elimination of aqueous phase if too little amount
           pm.XFA[0] = 0.;
       }

       // Main IPM iteration done
       // Main calculation of activity coefficients
        if( pm.PD >= 2 )
            status = CalculateActivityCoefficients( LINK_UX_MODE );

if( pm.pNP && status ) // && rLoop < 0  )
{
        setErrorMessage( 18, "E18IPM: IPM Main Descent", "Bad CalculateActivityCoefficients() status in SIA mode");
	return 3L;
}

// STEPWISE (6)  Stop point at IPM() main iteration
STEP_POINT( "IPM Iteration" );

        // Per-iteration descent record, GEMS3K_IPM_PROBE=<path> (PCI, FX, mass-balance columns).
        // Shows when the Dikin criterion has become noise around DXM while the composition has
        // settled. Cost when unset: one null test. O(N*L) per iteration when set.
        if( FILE* ipf = ipm_probe_file() )
        {
            double sumX = 0., minX = 1e300, maxX = 0.;
            for( long int jj = 0; jj < pm.L; jj++ )
            {
                sumX += pm.X[jj];
                if( pm.X[jj] > maxX ) maxX = pm.X[jj];
                if( pm.X[jj] > 0. && pm.X[jj] < minX ) minX = pm.X[jj];
            }
            // Mass-balance residual of the current primal pm.X, normalised so that mbRel > 1
            // means "fails the per-IC test MBR applies" (a physical tolerance, pa_DHB).
            long int iRel, iAbs; double mbRel, mbAbs;
            native_trace_mb_of( pm, pm.X, iRel, mbRel, iAbs, mbAbs );
            fprintf( ipf, "IPMIT %ld PCI=%.10e DXM=%.6e LM=%.10e LM1=%.10e FX=%.14e"
                          " NR=%ld sumX=%.10e minX=%.6e maxX=%.10e mbRel=%.6e\n",
                     (long)IT1, pm.PCI, pm.DXM, LM, LM1, pm.FX,
                     (long)N, sumX, minX, maxX, mbRel );
            fflush( ipf );
        }

        if( pm.PCI <= pm.DXM )  // Dikin criterion satisfied - converged!
            goto CONVERGED;
        if( pa_p->IpmStallWindow > 0 )
        {
            const long int W = std::min( (long int)pa_p->IpmStallWindow, kIpmStallMaxW );
            double sX = 0., xX = 0.;
            for( long int jj = 0; jj < pm.L; jj++ )
            {   sX += pm.X[jj];  if( pm.X[jj] > xX ) xX = pm.X[jj];  }
            stW_pci.push_back( pm.PCI ); stW_fx.push_back( pm.FX );
            stW_sum.push_back( sX );     stW_max.push_back( xX );
            if( (long int)stW_pci.size() > W + 1 )
            {   stW_pci.erase( stW_pci.begin() ); stW_fx.erase( stW_fx.begin() );
                stW_sum.erase( stW_sum.begin() ); stW_max.erase( stW_max.begin() ); }
            // Per-species history, from the same pm.X the scalars above are taken from.
            const size_t nR = (size_t)W + 1, nSp = (size_t)pm.L;
            if( stW_xf.size() != nR * nSp )
            {   stW_xf.assign( nR * nSp, 0. );  stW_n = 0;  stW_head = 0; }
            for( size_t k = 0; k < nSp; k++ )
                stW_xf[ (size_t)stW_head * nSp + k ] = pm.X[k];
            stW_head = (long int)( ( (size_t)stW_head + 1 ) % nR );
            if( (size_t)stW_n < nR ) stW_n++;
            if( (long int)stW_pci.size() == W + 1 )
            {
                auto spreadOK = []( const std::vector<double>& v, double tol ) -> bool
                {   double lo = v[0], hi = v[0];
                    for( size_t q = 1; q < v.size(); q++ )
                    {   if( v[q] < lo ) lo = v[q];  if( v[q] > hi ) hi = v[q]; }
                    const double ref = fabs( v.back() );
                    return ref > 0. ? ( hi - lo ) <= tol * ref : ( hi - lo ) == 0.;
                };
                long int up = 0;
                for( size_t q = 1; q < stW_pci.size(); q++ )
                    if( stW_pci[q] > stW_pci[q-1] ) up++;
                const double inc = (double)up / (double)( stW_pci.size() - 1 );
                // Every species amount flat, under both normalisers: the aggregate sumX/maxX
                // pair is blind to a redistribution at nearly constant total (between phases,
                // or between end-members of one solid solution). One clause catches large
                // absolute motion in a big species, the other large relative motion in a small
                // one. O(L) amortised: the scan bails on the first species still moving.
                bool xFlat = ( (size_t)stW_n == nR );
                if( xFlat )
                {
                    const size_t nw = ( (size_t)stW_head + nR - 1 ) % nR;  // newest row
                    double tot = 0.;
                    for( size_t k = 0; k < nSp; k++ ) tot += stW_xf[ nw * nSp + k ];
                    xFlat = ( tot > 0. );
                    const double negl = kIpmStallSpNegl * tot;
                    for( size_t k = 0; k < nSp && xFlat; k++ )
                    {
                        double lo = stW_xf[k], hi = stW_xf[k];
                        for( size_t q = 1; q < nR; q++ )
                        {   const double v = stW_xf[ q * nSp + k ];
                            if( v < lo ) lo = v;   if( v > hi ) hi = v; }
                        // Second clause: relative to the species' own largest amount over the
                        // window (not the newest value, or a vanishing species would exempt
                        // itself); species below kIpmStallSpNegl of the total are exempt.
                        if( ( hi - lo ) > kIpmStallCompTol * tot ) xFlat = false;
                        else if( hi > negl && ( hi - lo ) > kIpmStallSpRel * hi ) xFlat = false;
                    }
                }
                if( xFlat
                 && spreadOK( stW_fx,  kIpmStallFXTol )
                 && spreadOK( stW_sum, kIpmStallCompTol )
                 && spreadOK( stW_max, kIpmStallCompTol )
                 && inc >= kIpmStallIncLo && inc <= kIpmStallIncHi )
                    goto CONVERGED;   // the criterion is noise; the state is not moving
            }
        }
        if( nCNud > 0L && (IT1 >= cnr-2 && IT1 >= 2 ) )  // finish here because u vector diverges at further IPM iterations
            goto CONDITIONALLY_CONVERGED;
        // Restoring vectors Y and YF from X and XF for the next IPM iteration
        Restore_Y_YF_Vectors();
    } // end of the main IPM cycle
    // DXM was not reached in IPM iterations
    setErrorMessage( 6, "E06IPM: IPM Main Descent: " ,
            "IPM convergence criterion tolerance (pa_DK, GEMS: Pa_DK) could not be reached"
    		" (more than pa_IIM iterations done, GEMS: Pa_IIM).\n"
            " Try: raise pa_IIM (default 7000), or increase pa_DK (default 1e-6) a little.\n" );
    return 2L;  // bad convergence - too many IPM iterations or deterioration of dual solution!
//----------------------------------------------------------------------------
CONVERGED:
// if( !StatusDivg )
   pm.PCI = pm.DXM * 0.999999; // temporary - for smoothing
  return 0L;
CONDITIONALLY_CONVERGED:
   pm.PZ = 5; // Evtl. do something to reconfigure or circumvent PSSC()
  return 0L;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of mass-balance residuals in GEM IPM CSD.
/// \param   N - number of IC in IPM problem
/// \param   L -   number of DC in IPM problem
/// \param   A - DC stoichiometry matrix (LxN)
/// \param   Y - moles  DC quantities in IPM solution (L)
/// \param   B - Input bulk chem. compos. (N)
/// \param   C - mass balance residuals (N)
void TMultiBase::MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                                   double *B, double *C )
{
    long int ii, jj, i;
    for(ii=0;ii<N;ii++)
        C[ii]=B[ii];
    for(jj=0;jj<L;jj++)
     for( i=arrL[jj]; i<arrL[jj+1]; i++ )
     {  ii = arrAN[i];
         C[ii]-=(*(A+jj*N+ii))*Y[jj];
     }

    if(ipm_logger->should_log(spdlog::level::trace)) {
        ipm_logger->trace("MassBalanceResiduals {}", to_string(C, N));
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Diagnostics for a severe break of mass balance (abs.moles)
/// after GEM IPM PhaseSelect(). When pm.X is passed as parameter
/// \return -1 (Ok) or index of the first IC for which the balance is broken
long int
TMultiBase::CheckMassBalanceResiduals(double *Y )
{
    double cutoff;
    long int iRet = -1L;
    std::string buf;

    cutoff = min( pm.DHBM * 1e10, 1e-2 );  // 11.05.2010 DK
    MassBalanceResiduals( pm.N, pm.L, pm.A, Y, pm.B, pm.C);

    for(long int i=0; i<(pm.N - pm.E); i++)
    {
        if( fabs( pm.C[i] ) < cutoff )
            continue;
        if( iRet < 0  )
        {
            iRet = i;  // Error state is activated
            buf = "Mass balance is broken on iteration %ld  for ICs %-3.3s";
            buf +=  std::to_string(pm.ITG)+"  for ICs ";
            buf +=  name_for_message(pm.SB[i],3);
            setErrorMessage( 2, "E02IPM: PSSC(): ", buf.c_str());
        }
        else
        {
            addErrorMessage((std::string(", ")+name_for_message(pm.SB[i],3)).c_str());
        }
    } // i
    return iRet;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Interior Points Method
/// subroutine for unconditional optimization of the descent step length
/// on the interval 0 to LM.
/// uses the "Golden Section" algorithm
/// Formerly called LMD()
/// \return optimal value of LM which provides the largest possible monotonous
/// decrease in G(X)
//
double TMultiBase::OptimizeStepSize( double LM )
{
    double A,B,C,LM1,LM2;
    double FX1,FX2;
    A=0.0;
    B=LM;
    if( LM<2. )
        C=.05*LM;
    else C=.1;
    if( B-A<C)
        goto OCT;
    LM1=A+.382*(B-A);
    LM2=A+.618*(B-A);

    FX1= GX( LM1 );
    FX2= GX( LM2 );

SH1:
    if( fabs(FX1 - FX2) <= ( (fabs(FX1) < fabs(FX2) ? fabs(FX2) : fabs(FX1)) * 1e-13 ) )
          goto OCT; // Important fix of loss of IPM3 convergence (fix by Svitlana SD on 5.Mar.2020)
    if( FX1>FX2 )  // && fabs(FX1 - FX2) > ( (fabs(FX1) < fabs(FX2) ? fabs(FX2) : fabs(FX1)) * 1e-16 ))
        goto SH2;
    else goto SH3;
SH2:
    A=LM1;
    if( B-A<C)
        goto OCT;
    LM1=LM2;
    FX1=FX2;
    LM2=A+.618*(B-A);
    FX2=GX( LM2 );

    goto SH1;
SH3:
    B=LM2;
    if( B-A<C)
        goto OCT;
    LM2=LM1;
    FX2=FX1;
    LM1=A+.382*(B-A);
    FX1=GX( LM1 );
    goto SH1;
OCT:
    LM1=A+(B-A)/2;
    return(LM1);
}

//===================================================================

/// Cleaning the unstable phase with index k >= 0 (if k < 0 only DC will be cleaned)
void TMultiBase::DC_ZeroOff( long int jStart, long int jEnd, long int k )
{
  if( k >=0 )
     pm.YF[k] = 0.;

  for(long int j=jStart; j<jEnd; j++ )
     pm.Y[j] =  0.0;
}

/// Inserting minor quantities of DC which were zeroed off by SolveSimplex().
/// Important for the automatic initial approximation with solution phases
///  (k = -1)  or inserting a solution phase after PhaseSelect() (k >= 0)
//
void TMultiBase::DC_RaiseZeroedOff( long int jStart, long int jEnd, long int k )
{
//  double sfactor = scalingFactor;
//  if( fabs( sfactor ) > 1. )   // can reach 30 at total moles in system above 300000 (DK 11.03.2008)
//	  sfactor = 1.;       // Workaround for very large systems (insertion breaks the EFD convergence)
  if( k >= 0 )
       pm.YF[k] = 0.;

  for(long int j=jStart; j<jEnd; j++ )
  {
     switch( pm.DCC[j] )
     {
       case DC_AQ_PROTON:
       case DC_AQ_ELECTRON:
       case DC_AQ_SPECIES:
       case DC_AQ_SURCOMP:
            if( k >= 0 || pm.Y[j] < pm.DFYaqM )
               pm.Y[j] =  pm.DFYaqM;
           break;
       case DC_AQ_SOLVCOM:
       case DC_AQ_SOLVENT:
            if( k >= 0 || pm.Y[j] < pm.DFYwM )
                pm.Y[j] =  pm.DFYwM;
            break;
       case DC_GAS_H2O:
       case DC_GAS_CO2:
       case DC_GAS_H2:
       case DC_GAS_N2:
       case DC_GAS_COMP:
       case DC_SOL_IDEAL:
case DC_SCM_SPECIES:
            if( k >= 0 || pm.Y[j] < pm.DFYidM )
                  pm.Y[j] = pm.DFYidM;
             break;
       case DC_SOL_MINOR: case DC_SOL_MINDEP:
            if( k >= 0 || pm.Y[j] < pm.DFYhM )
                   pm.Y[j] = pm.DFYhM;
             break;
       case DC_SOL_MAJOR: case DC_SOL_MAJDEP:
            if( k >= 0 || pm.Y[j] < pm.DFYrM )
                  pm.Y[j] =  pm.DFYrM;
             break;
       case DC_SCP_CONDEN:
             if( k >= 0 )
             {                // Added 05.11.2007 DK
                 pm.Y[j] =  pm.DFYsM;
                 break;
             }
             if( pm.Y[j] < pm.DFYcM )
                  pm.Y[j] =  pm.DFYcM;
              break;
                    // implementation for adsorption?
       default:
             if( k >= 0 || pm.Y[j] < pm.DFYaqM )
                   pm.Y[j] =  pm.DFYaqM;
             break;
     }
     if( k >=0 )
     pm.YF[k] += pm.Y[j];
   } // i
}

/// The effective budget: GEMS3K_FILLOUT_BUDGET if set (negative = unset), else the field. The
/// environment override is kept for the CTest case fillout.budget. Not cached, so it can be
/// changed within one process. Read by both the call site (which must snapshot pm.Y before
/// DC_RaiseZeroedOff()) and ApplyFilloutBudget(), so they agree.
double TMultiBase::FilloutBudgetValue() const
{
    const char* v = std::getenv( "GEMS3K_FILLOUT_BUDGET" );
    const double envF = ( v && *v ) ? atof( v ) : -1.;
    return ( envF >= 0. ) ? envF : base_param()->FilloutBudget;
}

/// pa_FilloutBudget: caps how much DC_RaiseZeroedOff()'s class constants may change the mass
/// balance. yLp is pm.Y as the simplex left it, before the raise.
void TMultiBase::ApplyFilloutBudget( const std::vector<double>& yLp )
{
    // GEMS3K_FILLOUT_BUDGET overrides the field (negative = unset, so 0 stays a usable value).
    const double f = FilloutBudgetValue();
    if( !( f > 0. ) || (long int)yLp.size() != (size_t)pm.L || !pm.A || !pm.B )
        return;
    const long int N = pm.N, L = pm.L;
    const long int Zlim = N - pm.E;          // ordinary ICs only, as MBR's own loops scan
    if( Zlim <= 0 )
        return;

    // All of the excess comes from the raise: the LP solution satisfies A n = b exactly
    // (except where the LP supplies none of an IC, when there is nothing to scale).
    std::vector<double> raised( (size_t)Zlim, 0. );
    for( long int j = 0; j < L; j++ )
    {
        const double add = pm.Y[j] - yLp[(size_t)j];
        if( add <= 0. )
            continue;
        for( long int i = 0; i < Zlim; i++ )
        {
            const double a = pm.A[ i + j*N ];
            if( a > 0. )
                raised[(size_t)i] += add * a;
        }
    }

    long int nScaled = 0;
    double worstS = 1.;
    for( long int j = 0; j < L; j++ )
    {
        const double add = pm.Y[j] - yLp[(size_t)j];
        if( add <= 0. )
            continue;
        double s = 1.;
        for( long int i = 0; i < Zlim; i++ )
        {
            const double a = pm.A[ i + j*N ];
            if( a <= 0. || !( pm.B[i] > 0. ) || raised[(size_t)i] <= 0. )
                continue;
            s = std::min( s, f * pm.B[i] / raised[(size_t)i] );
        }
        if( s < 1. )
        {
            pm.Y[j] = yLp[(size_t)j] + add * s;
            nScaled++;
            worstS = std::min( worstS, s );
        }
    }
    if( nScaled > 0 )
        native_trace_decide( "filloutbudget f=%.6e scaled=%ld of %ld worst_s=%.6e",
                             f, nScaled, (long)L, worstS );
}

/// Adjustment of primal approximation according to kinetic constraints
long int TMultiBase::MetastabilityLagrangeMultiplier()
{
    double E = base_param()->DKIN; //1E-8;  Default min value of Lagrange multiplier p
//    E = 1E-30;

    for(long int J=0;J<pm.L;J++)
    {
        if( pm.Y[J] < 0. )   // negative number of moles!
            return J;
        if( pm.Y[J] < min( pm.lowPosNum, pm.DcMinM ))
            continue;
// kg44 why use a switch? Much to complicated! Simply correct all the values that are to big or to small. 
// values that are in the intervall given by DLL and DUL need no change.	
        if(pm.Y[J]<pm.DLL[J])
	{
	  if ((pm.DUL[J] - pm.DLL[J])/2.0 > E) 
	    pm.Y[J]=pm.DLL[J]+E;
	  else 
	    pm.Y[J]=pm.DLL[J]+(pm.DUL[J]-pm.DLL[J])/2.0; // this seems to work better than setting it directly to constraints
	}
        if (pm.Y[J]>pm.DUL[J])
	{
	  if ((pm.DUL[J] - pm.DLL[J])/2.0 > E) 
	    pm.Y[J]=pm.DUL[J]-E;
	  else 
	    pm.Y[J]=pm.DLL[J]+(pm.DUL[J]-pm.DLL[J])/2.0; // this seems to work better than setting it directly to constraints
	}
	  
    }   // J
    return -1L;
}

/// Calculation of weight multipliers for DCs
void TMultiBase::WeightMultipliers( bool square )
{
  long int J;
  double  W1, W2;

  for( J=0; J<pm.L; J++)
  {
    switch( pm.RLC[J] )
    {
      case NO_LIM:
      case LOWER_LIM:
           W1=(pm.Y[J]-pm.DLL[J]);
           if( square )
           {
             if(fabs(W1) > 1.34e120)    //  1.34e154
             {
                 if(signbit(W1))
                      W1 = -1.34e120;
                 else
                     W1 = 1.34e120;
             }
             pm.W[J]= W1 * W1;
           }
           else
             pm.W[J] = max( W1, 0. );
           break;
      case UPPER_LIM:
           W1=(pm.DUL[J]-pm.Y[J]);
           if( square )
           {
               if(fabs(W1) > 1.34e120)
               {
                  if(signbit(W1))
                     W1 = -1.34e120;
                  else
                     W1 = 1.34e120;
               }
               pm.W[J]= W1 * W1;
           }
           else
             pm.W[J] = max( W1, 0.);
           break;
      case BOTH_LIM:
           W1=(pm.Y[J]-pm.DLL[J]);
           W2=(pm.DUL[J]-pm.Y[J]);
           if( square )
           {
              if(fabs(W1) > 1.34e120)
              {
                 if(signbit(W1))
                   W1 = -1.34e120;
                 else
                   W1 = 1.34e120;
               }
             W1 = W1*W1;
             if(fabs(W2) > 1.34e120)
             {
                if(signbit(W2))
                   W2 = -1.34e120;
                else
                   W2 = 1.34e120;
             }
             W2 = W2*W2;
           }
           pm.W[J]=( W1 < W2 ) ? W1 : W2 ;
           if( !square && pm.W[J] < 0. ) pm.W[J]=0.;
           break;
      default: // error
          setErrorMessage( 16, "E16IPM: IPM Main Descent:", "Error in codes of some DC metastability constraints" );
          Error( pm.errorCode, pm.errorBuf );
    }
  } // J
}

#define  a(j,i) ((*(pm.A+(i)+(j)*Na)))

#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
// The following diagnostics (condition-number estimation, per-phase timing) add real overhead
// per linear solve, so they are compiled only with the CMake option
// ENABLE_BENCHMARK_DIAGNOSTICS. Other benchmark-only instrumentation should use the same macro.

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Cheap proxy for conditioning: ratio of largest to smallest |diagonal|
/// entry of the (symmetric) IPM matrix AA (N x N, column-major, AA[i+k*N]).
/// O(N), computed before any decomposition. Returns +inf if the smallest
/// diagonal magnitude is (numerically) zero.
static double DiagRatioConditionProxy( const double* AA, long int N )
{
    if( N <= 0 )
        return 0.;
    double diag_max = 0., diag_min = std::numeric_limits<double>::max();
    for( long int i = 0; i < N; i++ )
    {
        double d = fabs( AA[i + i*N] );
        if( d > diag_max )
            diag_max = d;
        if( d < diag_min )
            diag_min = d;
    }
    if( diag_min < 1e-300 )
        return std::numeric_limits<double>::infinity();
    return diag_max / diag_min;
}

/// A few steps of power iteration on the symmetric matrix AA (N x N,
/// column-major, AA[i+k*N]) to estimate its largest-magnitude eigenvalue.
/// Reuses no factorization (plain matrix-vector products) so it works
/// regardless of whether Cholesky or LU ends up solving the system.
static double PowerIterationMaxEig( const double* AA, long int N, int iters )
{
    if( N <= 0 )
        return 0.;
    std::vector<double> x(N), Ax(N);
    for( long int i = 0; i < N; i++ )
        x[i] = 1. + 0.001 * (double)(i % 7);

    double lambda_max = 0.;
    for( int it = 0; it < iters; it++ )
    {
        for( long int i = 0; i < N; i++ )
        {
            double s = 0.;
            for( long int j = 0; j < N; j++ )
                s += AA[i + j*N] * x[j];
            Ax[i] = s;
        }
        double norm = 0.;
        for( long int i = 0; i < N; i++ )
            norm += Ax[i]*Ax[i];
        norm = sqrt( norm );
        if( norm < 1e-300 )
            return 0.;
        for( long int i = 0; i < N; i++ )
            x[i] = Ax[i] / norm;
        lambda_max = norm;
    }
    return lambda_max;
}

/// A few steps of inverse iteration using an already-computed factorization
/// (JAMA::Cholesky or JAMA::LU, both expose Array1D<double> solve(const
/// Array1D<double>&)) to estimate the smallest-magnitude eigenvalue of the
/// factorized matrix. Cheap: each step is just a triangular solve reusing
/// the factors already paid for by the main solve() call.
template <class Decomp>
static double InverseIterationMinEig( Decomp& decomp, long int N, int iters )
{
    if( N <= 0 )
        return 0.;
    Array1D<double> x(N);
    for( long int i = 0; i < N; i++ )
        x[(int)i] = 1. + 0.001 * (double)(i % 7);

    double inv_norm = 0.;
    for( int it = 0; it < iters; it++ )
    {
        Array1D<double> y = decomp.solve( x );
        double norm = 0.;
        for( long int i = 0; i < N; i++ )
            norm += y[(int)i]*y[(int)i];
        norm = sqrt( norm );
        if( norm < 1e-300 )
            return std::numeric_limits<double>::infinity();
        for( long int i = 0; i < N; i++ )
            x[(int)i] = y[(int)i] / norm;
        inv_norm = norm;
    }
    if( inv_norm < 1e-300 )
        return std::numeric_limits<double>::infinity();
    return 1. / inv_norm;
}
#endif // GEMS3K_BENCHMARK_DIAGNOSTICS

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// pa_IpmAugmentedKKT = 1 or 2: the main-loop solve for pm.U that never forms A^T W A.
/// Both modes solve (A_act^T W A_act + D) u = A_act^T W F over the same species as the
/// normal-equations assembly below (Y > min(lowPosNum, DcMinM)), differing from it only by D
/// and rounding. pm.MU is not written here (DikinsCriterion() recomputes it from pm.U).
/// \return 0 solved, 2 singular - the caller then takes this step by the normal equations
///         (DECIDE ipmkkt-fallback). Never 1: the augmented solve must not fail where the
///         normal equations would have gone on.
long int TMultiBase::SolveIpmAugmented( long int N )
{
    const long int Na = pm.N;   // stride of the a(j,i) macro
    const long int mode = base_param()->IpmAugmentedKKT;
    const double kRescueRel = 1e-12;

    std::vector<long int> act;
    for( long int jj = 0; jj < pm.L; jj++ )
        if( pm.Y[jj] > min( pm.lowPosNum, pm.DcMinM ) )
            act.push_back( jj );
    const long int La = (long int)act.size();

    // D is zero on every row some active species carries, and 1e-12 * max_i d_i on a row none
    // does (d_i = sum_j W_j a_ji^2, the normal matrix's own diagonal). A uniform D would move
    // the weak directions of u in this badly conditioned system.
    std::vector<double> dreg( (size_t)N, 0. );
    for( long int r = 0; r < La; r++ )
    {
        const long int j = act[(size_t)r];
        for( long int i = 0; i < N; i++ )
        {   const double v = a(j,i);
            dreg[(size_t)i] += pm.W[j] * v * v; }
    }
    double dmax = 0.;
    for( long int i = 0; i < N; i++ )
        dmax = std::max( dmax, dreg[(size_t)i] );
    if( !( dmax > 0. ) )
        return 2;                  // nothing active carries any IC: no scale to regularise with
    long int nZero = 0, iZero = -1;
    for( long int i = 0; i < N; i++ )
    {
        if( dreg[(size_t)i] > 0. )
            dreg[(size_t)i] = 0.;
        else
        {   // Zero-row rescue: u_i is undetermined; the floor sets it to 0 instead of failing.
            dreg[(size_t)i] = kRescueRel * dmax;
            if( iZero < 0 ) iZero = i;
            nZero++;
        }
    }
    if( nZero > 0 && ipmKktRescues++ == 0 )
        native_trace_decide( "ipmkkt-zerorow mode=%ld ics=%ld first=%s itg=%ld",
                             (long)mode, (long)nZero,
                             char_array_to_string( pm.SB[iZero], MAXICNAME ).c_str(), (long)pm.ITG );

    if( mode == 1 )
    {
        // Saddle-point form, dense LU with partial pivoting:
        //   [ I          -W A_act ] [ x ]   [ -W F ]
        //   [ -A_act^T   -D       ] [ u ] = [  0   ]
        const long int K = La + N;
        Array2D<double> KKT( K, K, 0.0 );
        Array1D<double> V( K, 0.0 );
        for( long int r = 0; r < La; r++ )
        {
            const long int j = act[(size_t)r];
            KKT[(int)r][(int)r] = 1.0;
            for( long int c = 0; c < N; c++ )
                KKT[(int)r][(int)(La + c)] = -pm.W[j] * a(j,c);
            V[(int)r] = -pm.W[j] * pm.F[j];
        }
        for( long int r = 0; r < N; r++ )
        {
            for( long int c = 0; c < La; c++ )
                KKT[(int)(La + r)][(int)c] = -a(act[(size_t)c], r);
            KKT[(int)(La + r)][(int)(La + r)] = -dreg[(size_t)r];
        }
        JAMA::LU<double> lu( KKT );
        if( !lu.isNonsingular() )
        {
            ipm_logger->warn("SolveIpmAugmented (pa_IpmAugmentedKKT=1): augmented matrix singular, "
                             "L_act={} N={}", La, N);
            return 2;
        }
        Array1D<double> S = lu.solve( V );
        for( long int i = 0; i < N; i++ )
            pm.U[i] = S[(int)(La + i)];
        return 0;
    }

    // mode 2: least squares  min |W^1/2 (A u - F)|^2 + u^T D u  by Householder QR of the
    // (La+N) x N stack [W^1/2 A_act ; D^1/2], right-hand side [W^1/2 F ; 0]. Column-major.
    const long int M = La + N;
    std::vector<double> Q( (size_t)M * (size_t)N, 0. ), b( (size_t)M, 0. );
    for( long int r = 0; r < La; r++ )
    {
        const long int j = act[(size_t)r];
        const double sw = sqrt( std::max( pm.W[j], 0. ) );
        for( long int c = 0; c < N; c++ )
            Q[(size_t)c * M + r] = sw * a(j,c);
        b[(size_t)r] = sw * pm.F[j];
    }
    for( long int i = 0; i < N; i++ )
        Q[(size_t)i * M + La + i] = sqrt( dreg[(size_t)i] );

    std::vector<double> Rd( (size_t)N, 0. ), v( (size_t)M, 0. );
    for( long int k = 0; k < N; k++ )
    {
        double* qk = &Q[(size_t)k * M];
        double nrm = 0.;
        for( long int r = k; r < M; r++ ) nrm += qk[r] * qk[r];
        nrm = sqrt( nrm );
        if( !( nrm > 0. ) )
        {
            ipm_logger->warn("SolveIpmAugmented (pa_IpmAugmentedKKT=2): zero column {} of {}", k, N);
            return 2;
        }
        const double alpha = ( qk[k] > 0. ) ? -nrm : nrm;
        double vn2 = 0.;
        for( long int r = k; r < M; r++ ) { v[(size_t)r] = qk[r]; }
        v[(size_t)k] -= alpha;
        for( long int r = k; r < M; r++ ) vn2 += v[(size_t)r] * v[(size_t)r];
        if( vn2 > 0. )
        {
            for( long int c = k; c < N; c++ )
            {
                double* qc = &Q[(size_t)c * M];
                double s = 0.;
                for( long int r = k; r < M; r++ ) s += v[(size_t)r] * qc[r];
                s *= 2. / vn2;
                for( long int r = k; r < M; r++ ) qc[r] -= s * v[(size_t)r];
            }
            double s = 0.;
            for( long int r = k; r < M; r++ ) s += v[(size_t)r] * b[(size_t)r];
            s *= 2. / vn2;
            for( long int r = k; r < M; r++ ) b[(size_t)r] -= s * v[(size_t)r];
        }
        Rd[(size_t)k] = qk[k];
    }
    double rmax = 0.;
    for( long int k = 0; k < N; k++ ) rmax = std::max( rmax, fabs( Rd[(size_t)k] ) );
    for( long int k = 0; k < N; k++ )
        if( !( fabs( Rd[(size_t)k] ) > 1e-14 * rmax ) )
        {
            ipm_logger->warn("SolveIpmAugmented (pa_IpmAugmentedKKT=2): rank-deficient R at {} of {}", k, N);
            return 2;
        }
    for( long int k = N - 1; k >= 0; k-- )
    {
        double s = b[(size_t)k];
        for( long int c = k + 1; c < N; c++ )
            s -= Q[(size_t)c * M + k] * pm.U[c];
        pm.U[k] = s / Rd[(size_t)k];
    }
    return 0;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Make and Solve a system of linear equations to find the dual vector
/// approximation using a method of Cholesky Decomposition. Good if a
/// square matrix R happens to be symmetric and positive defined.
/// If Cholesky Decomposition does not solve the problem, an attempt is done
/// to solve the SLE using method of LU Decomposition
/// (A = L*U , L is lower triangular ( has elements only on the diagonal and below )
///   U is is upper triangular ( has elements only on the diagonal and above))
/// \param
///    initAppr - Inital approximation point(true) or iteration of IPM (false)
///    N - dimension of the matrix R (number of equations)
/// \return 0  - solved OK;
///         1  - no solution, degenerated or inconsistent system
long int TMultiBase::MakeAndSolveSystemOfLinearEquations( long int N, bool initAppr )
{
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
    auto solve_t0 = std::chrono::steady_clock::now();
    double diag_ms = 0.;   // condition-number-diagnostics-only time, this call
    pm.SolveCallCount++;
#endif

    long int ii, i, jj, kk, k, Na = pm.N;


    if( !initAppr && base_param()->IpmAugmentedKKT > 0 )
    {
        const long int ret = SolveIpmAugmented( N );
        if( ret == 2 )
        {   // singular for the augmented solve: take this step by the normal equations below
            if( s_ipmKktFallbacks++ == 0 )
                native_trace_decide( "ipmkkt-fallback mode=%ld itg=%ld N=%ld",
                                     (long)base_param()->IpmAugmentedKKT, (long)pm.ITG, (long)N );
        }
        else
        {
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
            pm.SolveTimeMs += std::chrono::duration<double, std::milli>(
                                  std::chrono::high_resolution_clock::now() - solve_t0 ).count();
            pm.CondNumTimeMs += diag_ms;
#endif
            return ret;
        }
    }

    Alloc_A_B( N );
    // Diagonal (Jacobi) preconditioning scale factors for the initAppr (MBR)
    // branch only - see the scaling block below the assembly for the rationale.
    // Left empty in the main-IPM (initAppr=false) branch, which is not
    // preconditioned.
    std::vector<double> Dscale;

    // Making the matrix of IPM linear equations
    for( kk = 0; kk < N; kk++)
        for( ii=0; ii < N; ii++ )
            (*(AA+(ii)+(kk)*N)) = 0.;

    for( jj=0; jj < pm.L; jj++ )
    {
        if( pm.Y[jj] > min( pm.lowPosNum, pm.DcMinM ) )
        {
            for( k = arrL[jj]; k < arrL[jj+1]; k++)
                for( i = arrL[jj]; i < arrL[jj+1]; i++ )
                {   ii = arrAN[i];
                    kk = arrAN[k];
                    if( ii >= N || kk >= N )
                        continue;
                    (*(AA+(ii)+(kk)*N)) += a(jj,ii) * a(jj,kk) * pm.W[jj];
                }
        }
    }

    if( initAppr )
        for( ii = 0; ii < N; ii++ )
            BB[ii] = pm.C[ii];
    else {
        for( ii = 0; ii < N; ii++ )
            BB[ii] = 0.;
        for( jj=0; jj < pm.L; jj++ )
            if( pm.Y[jj] > min( pm.lowPosNum, pm.DcMinM ) )
                for( i = arrL[jj]; i < arrL[jj+1]; i++ )
                {   ii = arrAN[i];
                    if( ii >= N )
                        continue;
                    BB[ii] += pm.F[jj] * a(jj,ii) * pm.W[jj];
                }
    }

    // Diagonal (Jacobi) preconditioning of the initAppr (MBR) Schur-complement matrix:
    // A' = D*A*D, B' = D*B with D = diag(1/sqrt(|A_ii|)), which keeps symmetry and positive
    // definiteness and normalises every nonzero diagonal entry to 1. The dual is unscaled
    // (U = D*U') at the unpack site below. A row with a structurally zero diagonal gets
    // Dscale = 1 (it is singular anyway and goes to the existing singular-matrix path).
    // Rescaling cannot repair a genuine near-singularity of A D^-1 A^T (e.g. the H and O rows
    // tied by water's fixed 2:1 stoichiometry).
    if( initAppr )
    {
        Dscale.assign( N, 1. );
        for( ii = 0; ii < N; ii++ )
        {
            double d = fabs( *(AA+(ii)+(ii)*N) );
            if( d > 1e-300 )
                Dscale[ii] = 1. / sqrt( d );
        }
        for( kk = 0; kk < N; kk++ )
            for( ii = 0; ii < N; ii++ )
                (*(AA+(ii)+(kk)*N)) *= Dscale[ii] * Dscale[kk];
        for( ii = 0; ii < N; ii++ )
            BB[ii] *= Dscale[ii];
    }

#ifndef PGf90
    Array2D<double> A( N, N, AA );
    Array1D<double> B( N, BB );
#else
    Array2D<double> A( N, N);
    Array1D<double> B( N );
    for( kk = 0; kk < N; kk++)
        for( ii = 0; ii < N; ii++ )
            A[kk][ii] = (*(AA+(ii)+(kk)*N));
    for( ii = 0; ii < N; ii++ )
        B[ii] = BB[ii];
#endif

#ifndef NDEBUG
        ipm_logger->debug("MakeAndSolveSystemOfLinearEquations\n {} \n {} \n",
                          A.to_string(), B.to_string());
#endif

#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
    // Condition-number diagnostics: cheap diagonal-ratio proxy (always available)
    // and a power-iteration estimate of the largest eigenvalue (reused below for
    // both the Cholesky and LU branches, whichever ends up solving the system).
    // Timed separately from the rest of the solve (diag_ms) so its cost can be
    // attributed per system instead of inferred from noisy whole-call timing.
    auto diag_t0 = std::chrono::steady_clock::now();
    pm.CondNumDiag = std::max( pm.CondNumDiag, DiagRatioConditionProxy( AA, N ) );
    double lambda_max = PowerIterationMaxEig( AA, N, 6 );
    diag_ms += std::chrono::duration<double, std::milli>(
                   std::chrono::steady_clock::now() - diag_t0 ).count();
#endif

    // From here on, the NIST TNT Jama/C++ linear algebra package is used
    //    (credit: http://math.nist.gov/tnt/download.html)
    // this routine constructs the Cholesky decomposition, A = L x LT .
    JAMA::Cholesky<double> chol(A);
 #ifndef NDEBUG
    ipm_logger->debug("Cholesky Decomposition\n{}", chol.to_string());
 #endif

    if( chol.is_spd() )
    {
        B = chol.solve( B );
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
        auto inv_t0 = std::chrono::steady_clock::now();
        double lambda_min = InverseIterationMinEig( chol, N, 6 );
        diag_ms += std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - inv_t0 ).count();
        double cond = ( lambda_min > 1e-300 ) ? lambda_max / lambda_min
                                               : std::numeric_limits<double>::infinity();
        pm.CondNum = std::max( pm.CondNum, cond );
#endif
    }
    else
    {
        // no solution by Cholesky decomposition; Trying the LU Decompositon
        // The LU decompostion with pivoting always exists, even if the matrix is
        // singular, so the constructor will never fail.
        JAMA::LU<double> lu(A);
                // The primary use of the LU decomposition is in the solution
        // of square systems of simultaneous linear equations.
        // This will fail if isNonsingular() returns false.
        if( !lu.isNonsingular() )
        {
            // Singular matrix — log diagnostics to identify the cause
            ipm_logger->warn("MakeAndSolveSystemOfLinearEquations FAILED:");
            ipm_logger->warn("LU Decomposition\n{}", lu.to_string());
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
            pm.CondNum = std::numeric_limits<double>::infinity();
#endif
            for( long int r=0; r<N; r++ )
            {
                double row_norm = 0.0;
                for( long int c=0; c<N; c++ )
                    row_norm += std::abs( *(AA+(r)+(c)*N) );
                if( row_norm < 1e-20 )
                    ipm_logger->trace("  IC[{}] B={:.3e} — zero row, "
                                      "no active species contribute to this IC",
                                      r, pm.B[r]);
            }
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
            pm.SolveTimeMs += std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - solve_t0 ).count();
            pm.CondNumTimeMs += diag_ms;
#endif
            return 1; // Singular matrix - too bad! No solution ...
        }

        B = lu.solve( B );
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
        auto inv_t0 = std::chrono::steady_clock::now();
        double lambda_min = InverseIterationMinEig( lu, N, 6 );
        diag_ms += std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - inv_t0 ).count();
        double cond = ( lambda_min > 1e-300 ) ? lambda_max / lambda_min
                                               : std::numeric_limits<double>::infinity();
        pm.CondNum = std::max( pm.CondNum, cond );
#endif
    }

    if( initAppr )
    {
        // Unscale: the system just solved was A'U' = B' with A' = D*A*D (the
        // preconditioned matrix assembled above), so the true dual is U = D*U'.
        for( ii = 0; ii < N; ii++ )
            pm.Uefd[ii] = B[(int)ii] * Dscale[ii];
    }
    else {
        for( ii = 0; ii < N; ii++ )
            pm.U[ii] = B[(int)ii];
    }
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
    pm.SolveTimeMs += std::chrono::duration<double, std::milli>(
                          std::chrono::steady_clock::now() - solve_t0 ).count();
    pm.CondNumTimeMs += diag_ms;
#endif
    return 0;
}

#undef a

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of MU values (in the vector of direction of descent) and Dikin criterion
/// \param initAppr - use in MassBalanceRefinement() (true) or main iteration of IPM (false)
/// \param N - dimension of the matrix R (number of equations)
double TMultiBase::DikinsCriterion(  long int N, bool initAppr )
{
  long int  J;
  double Mu, PCI=0., qMu;

  for(J=0;J<pm.L;J++)
  {
    if( pm.Y[J] > min( pm.lowPosNum, pm.DcMinM ) )
    {
      if( initAppr )
      {
          Mu = DC_DualChemicalPotential( pm.Uefd, pm.A+J*pm.N, N, J );
          qMu = Mu*pm.W[J];
          if( fabs(qMu) > 1.34e120 )  // workaround of NAN in PCI += qMu*qMu
          {
            if(signbit(qMu))
              qMu = -1.34e120;
            else
              qMu = 1.34e120;
          }
          pm.MU[J] = qMu;
          PCI += qMu*qMu;
//          PCI += sqrt(fabs(qMu));  // Experimental - absolute differences?
      }
      else {
          Mu = DC_DualChemicalPotential( pm.U, pm.A+J*pm.N, N, J );
          Mu -= pm.F[J];
          qMu =  Mu*pm.W[J];
          pm.MU[J] = qMu;
          PCI += fabs(qMu);
//          PCI += qMu*qMu;    // sum of squares (see Chudnenko ea 2001 report)
//          PCI += fabs(qMu*Mu);   // As it was before 2009
      }
    }
    else
      pm.MU[J]=0.; // initializing dual potentials
  }
  if( initAppr )
  {
     if( PCI > pm.lowPosNum  )
     {
         PCI=1./sqrt(PCI);
//         PCI = 1./PCI;
//         PCI = 1./PCI/PCI;
     }
         else PCI=1.; // zero Psi value ?
  }
  else {  // if PCI += qMu * qMu
          ;
//      PCI = sqrt( PCI );
//      PCI *= PCI;
  }

  ipm_logger->trace("DikinsCriterion {}", PCI);
  return PCI;
}

/// Estimation of the descent step length LM
/// \param initAppr - MBR() (true) or iteration of IPM (false)
double TMultiBase::StepSizeEstimate(  bool initAppr )
{
    long int J, Z = -1;
    double LM=1., LM1=1., Mu;

    for(J=0;J<pm.L;J++)
    {
        Mu = pm.MU[J];
        if( pm.RLC[J]==NO_LIM || pm.RLC[J]==LOWER_LIM || pm.RLC[J]==BOTH_LIM )
        {
            if( Mu < 0 && fabs(Mu) > pm.lowPosNum )
            {
                if( Z == -1 )
                {
                    Z = J;
                    LM = (-1)*(pm.Y[Z]-pm.DLL[Z])/Mu;
                }
                else
                {
                    LM1 = (-1)*(pm.Y[J]-pm.DLL[J])/Mu;
                    if( LM > LM1)
                        LM = LM1;
                }
            }
        }
        if( pm.RLC[J]==UPPER_LIM || pm.RLC[J]==BOTH_LIM )
        {
            if( Mu > pm.lowPosNum ) // *100.)
            {
                if( Z == -1 )
                {
                    Z = J;
                    LM = (pm.DUL[Z]-pm.Y[Z])/Mu;
                }
                else
                {
                    LM1=(pm.DUL[J]-pm.Y[J])/Mu;
                    if( LM>LM1)
                        LM=LM1;
                }
            }
        }
    }

    if( initAppr )
    {
        if( Z == -1 )
            LM = pm.PCI;
        else
            LM *= .95;     // Smoothing of final lambda value
    }
    else
    {
        if( Z == -1 ) {
            LM = 1./sqrt(pm.PCI);  // Might cause infinite loop in OptimizeStepSize() if PCI is too low?
        }
        //     LM = min( LM, 10./pm.DX );
        LM = min( LM, 1.0e10 );  // Set an empirical upper limit for LM to prevent freezing
    }
    return LM;
}

/// Restoring primal vectors Y and YF
void TMultiBase::Restore_Y_YF_Vectors()
{
    long int Z, I, JJ = 0;

    for( Z=0; Z<pm.FI ; Z++ )
    {
        if( pm.XF[Z] <= pm.DSM ||
                ( pm.PHC[Z] == PH_SORPTION &&
                  ( pm.XFA[Z] < base_param()->ScMin) ) )
        {
            pm.YF[Z]= 0.;
            if( pm.FIs && Z<pm.FIs )
                pm.YFA[Z] = 0.;
            for(I=JJ; I<JJ+pm.L1[Z]; I++)
            {
                pm.Y[I]=0.;
                pm.lnGam[I] = 0.;
            }
        }
        else
        {
            pm.YF[Z] = pm.XF[Z];
            if( pm.FIs && Z < pm.FIs )
                pm.YFA[Z] = pm.XFA[Z];
            for(I = JJ; I < JJ+pm.L1[Z]; I++)
                pm.Y[I]=pm.X[I];
        }
        JJ += pm.L1[Z];
    } // Z

}

/// Calculation of the system size scaling factor and modified thresholds/cutoffs/insertion values
/// Replaces calcSfactor()
double TMultiBase::RescaleToSize( bool /*standard_size*/ )
{
    double SizeFactor=1.;
    const BASE_PARAM *pa_p = base_param();

    pm.SizeFactor = 1.;
//  re-scaling numeric settings
    pm.DHBM = SizeFactor * pa_p->DHB; // Mass balance accuracy threshold
    pm.DXM =  SizeFactor * pa_p->DK;   // Dikin' convergence threshold
//    pm.DX = pa_p->DK;
    pm.DSM =  SizeFactor * pa_p->DS;   // Cutoff for solution phase amount
// Cutoff amounts for DCs
    pm.XwMinM = SizeFactor * pa_p->XwMin;  // cutoff for the amount of water-solvent
    pm.ScMinM = SizeFactor * pa_p->ScMin;  // cutoff for amount of the sorbent
    pm.DcMinM = SizeFactor * pa_p->DcMin;  // cutoff for Ls set (amount of solution phase component)
    pm.PhMinM = SizeFactor * pa_p->PhMin;  // cutoff for single-comp.phase amount and its DC
  // insertion values before SolveSimplex() (re-scaled to system size)
    pm.DFYwM = SizeFactor * pa_p->DFYw;
    pm.DFYaqM = SizeFactor * pa_p->DFYaq;
    pm.DFYidM = SizeFactor * pa_p->DFYid;
    pm.DFYrM = SizeFactor * pa_p->DFYr;
    pm.DFYhM = SizeFactor * pa_p->DFYh;
    pm.DFYcM = SizeFactor * pa_p->DFYc;
    // Insertion value for PhaseSelection()
    pm.DFYsM = SizeFactor * pa_p->DFYs; // pure condenced phase and its DC

    return SizeFactor;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Internal memory allocation for IPM performance optimization
/// (since version 2.2.0)
//
void TMultiBase::Alloc_A_B( long int newN )
{
  if( AA && BB && (newN == sizeN) )
    return;
  Free_A_B();
  AA = new  double[newN*newN];
  BB = new  double[newN];
  sizeN = newN;
}

void TMultiBase::Free_A_B()
{
  if(AA)
    { delete[] AA; AA = 0; }
  if(BB)
    { delete[] BB; BB = 0; }
  sizeN = 0;
}

#define  a(j,i) ((*(pm.A+(i)+(j)*pm.N)))
/// Building an index list of non-zero elements of the matrix pm.A
void TMultiBase::Build_compressed_xAN()
{
 long int ii, jj, k;

 // Calculate number of non-zero elements in A matrix
 k = 0;
 for( jj=0; jj<pm.L; jj++ )
   for( ii=0; ii<pm.N; ii++ )
     if( fabs( a(jj,ii) ) > 1e-12 )
       k++;

   // Free old memory allocation
    Free_compressed_xAN();

   // Allocate memory
   arrL = new long int[pm.L+1];
   arrAN = new long int[k];

   // Set indexes in the index arrays
   k = 0;
   for( jj=0; jj<pm.L; jj++ )
   { arrL[jj] = k;
     for( ii=0; ii<pm.N; ii++ )
       if( fabs( a(jj,ii) ) > 1e-12 )
       {
        arrAN[k] = ii;
        k++;
       }
   }
   arrL[jj] = k;
}
#undef a

void TMultiBase::Free_compressed_xAN()
{
  if( arrL  )
    { delete[] arrL; arrL = 0;  }
  if( arrAN )
    { delete[] arrAN; arrAN = 0;  }
}

void TMultiBase::Free_internal()
{
  Free_compressed_xAN();
  Free_A_B();
 }

/// Internal memory allocation for IPM performance optimization
void TMultiBase::Alloc_internal()
{
// optimization 08/02/2007
 Alloc_A_B( pm.N );
 Build_compressed_xAN();
}

// add09
void TMultiBase::setErrorMessage( long int num, const char *code, const char * msg)
{
  size_t len_code, len_msg;
  pm.Ec  = num;
  len_code = strlen(code);
  if(len_code > 99)
      len_code = 99;
  memcpy( pm.errorCode, code, len_code );
  pm.errorCode[len_code] ='\0';
  len_msg = strlen(msg);
  if(len_msg > 1023)
      len_msg = 1023;
  memcpy( pm.errorBuf,  msg,  len_msg );
  pm.errorBuf[len_msg] ='\0';
}

void TMultiBase::addErrorMessage( const char * msg)
{
  auto len = strlen(pm.errorBuf);
  auto lenm = strlen( msg );
  if( len + lenm < 1023 )
  {
    memcpy(pm.errorBuf+len, msg, lenm);
    pm.errorBuf[len+lenm] ='\0';
  }
}

/// Added for implementation of divergence detection in dual solution 06.05.2011 DK
void TMultiBase::Alloc_uDD( long int newN )
{
    if( U_mean && U_M2 && U_CVo && U_CV && ICNud && (newN == nNu) )
      return;
    Free_uDD();
    U_mean = new  double[newN]; // w3 u mean values for r
    U_M2 = new  double[newN];   // w3 u mean values for r-1
    U_CVo = new  double[newN];  // w3 u mean difference for r-1
    U_CV = new  double[newN];   // w3 u mean difference for r
    ICNud = new long int[newN];
    nNu = newN;
}

void TMultiBase::Free_uDD()
{
    if( U_mean  )
      { delete[] U_mean; U_mean = 0; }
    if( U_M2  )
      { delete[] U_M2; U_M2 = 0; }
    if( U_CVo  )
      { delete[] U_CVo; U_CVo = 0; }
    if( U_CV  )
      { delete[] U_CV; U_CV = 0; }
    if( ICNud )
      { delete[] ICNud; ICNud = 0; }
    nNu = 0;
}

/// initializing data for u divergence detection
void TMultiBase::Reset_uDD( long int nr, bool trace )
{
    long int i;
    cnr = nr;
    for( i=0; i<nNu; i++)
    {
      U_mean[i] = 0.; U_M2[i] = 0.;
      U_CVo[i] = 0.; U_CV[i] = 0;
      ICNud[i] = -1L;
    }
    nCNud = 0;
    if ( trace )
    {
        ipm_logger->debug("UD3 trace: {}  SIA={}  Itr   C_D:  {}",
                          char_array_to_string(pm.stkey, EQ_RKLEN), pm.pNP, char_array_to_string(pm.SB1[0],MAXICNAME));
    }
    if( base_param()->PSM >= 3 )
    {
      TNode::ipmlog_file->debug(" UD3 trace: {}  SIA= {} Itr   C_D: {}",
                           char_array_to_string(pm.stkey, EQ_RKLEN), pm.pNP, char_array_to_string(pm.SB1[0],MAXICNAME));
    }
}

/// Incrementing mean u values for r-th (current) IPM iteration
void TMultiBase::Increment_uDD( long int r, bool trace )
{
    long int i;
    //double delta;
    cnr = r; // r+1;
    if( cnr == 0 )
        return;
    if( base_param()->PSM >= 3 )
    {
       TNode::ipmlog_file->debug("ncrement_uDD {}  {}", r, pm.PCI);
    }
    if( trace )
    {
        ipm_logger->debug("ncrement_uDD {}  {}", r, pm.PCI);
    }

    for( i=0; i<nNu; i++)
    {
// Calculating moving average of three u_i values
      switch( cnr )
      {
          case 1: U_mean[i] = pm.U[i];
                  U_M2[i] = U_mean[i];
                  U_CV[i] = 0.;
                  pm.Uc[i][0] = pm.U[i];
                  pm.Uc[i][1] = pm.U[i];
                  break;
          case 2: U_M2[i] = U_mean[i];
                  U_mean[i] = (pm.U[i] + pm.Uc[i][0] + pm.Uc[i][0] )/3.;
                  pm.Uc[i][1] = pm.Uc[i][0];
                  pm.Uc[i][0] = pm.U[i];
                  break;
          default:U_M2[i] = U_mean[i];
                  U_mean[i] = (pm.U[i] + pm.Uc[i][0] + pm.Uc[i][1])/3.;
                  pm.Uc[i][1] = pm.Uc[i][0];
                  pm.Uc[i][0] = pm.U[i];
                  break;
      }
      U_CVo[i] = U_CV[i];
      U_CV[i] = U_mean[i] - U_M2[i];
      //delta = fabs(U_CV[i] - U_CVo[i]);
      if( trace )
      {
        ipm_logger->debug("U={}  U_mean={} U_CV={}", pm.U[i],  U_mean[i], U_CV[i]);
      }
      if( base_param()->PSM >= 3 )
      {
         TNode::ipmlog_file->debug("U={}  U_mean={} U_CV={}", pm.U[i],  U_mean[i], U_CV[i]);

      }
//      delta = pm.U[i] - U_mean[i];
//      U_mean[i] += delta / cnr;
//      U_M2[i] += delta * ( pm.U[i] - U_mean[i] );
//      if( cnr > 2 )
//          U_CVo[i] = U_CV[i];
//      U_CV[i] = sqrt( U_M2[i]/cnr ) / fabs( U_mean[i] );
//      if( cnr < 2 )
//          U_CVo[i] = U_mean[i];
      // Copy of dual solution approximation
//      if( cnr < 2 )
//          pm.Uc[i][0] = pm.U[i];
//
    } // end for i
}

/// Checking for divergence in coef.variation of dual solution approximation.
/// Compares with CV value tolerance (mode = 0) or with CV increase
///          tolerance (mode = 1)
/// \return  0 if no divergence has been detected
///          >0 - number of diverging dual chemical potentials
///            (their IC names are collected in the ICNud list)
///
long int TMultiBase::Check_uDD( long int mode, double DivTol,  bool trace )
{
    long int i;
    double delta = 0., tol_gen=1., tolerance=1., log_bi=0.;
    bool FirstTime = true;

    tol_gen = DivTol;
    if( pm.PCI < 1 )
        tol_gen *= pm.PCI;
    if( cnr <= 1 )
        return 0;

    //  Check here that pm.PCI is reasonable (i.e. C_D < 1)?
    //
    for( i=0; i<nNu; i++)
    {
        if( DivTol >= 1e6 )
            continue;     // Disabling divergence checks for complete tracing
        // Checking absolute ranges of u[i] - to be checked for 'exotic' systems!
        if( i == nNu-1 && pm.E && pm.U[i] >= -50. && pm.U[i] <= 100.) // charge
            continue;
        else if( pm.U[i] >= -600. && pm.U[i] <= 400. ) // range for other ICs
        {
            tolerance = tol_gen;
            log_bi = log( pm.B[i] );  // Fixed 11.07.2011 DK
            if( log_bi > 0. )
                tolerance = tol_gen / log_bi;
            if( log_bi < 0. )
                tolerance = tol_gen * -log_bi;
            if( tolerance < 1. )
                tolerance = 1.;     // To prevent dangerous low tolerances
            if( tolerance > DivTol)
                tolerance = DivTol; // To prevent useless high tolerances

            if( !mode ) // Monitor difference between new and old mean3 u_i
            {
                // Calculation of abs.difference of moving averages at r and r-1
                delta = fabs(U_mean[i] - U_M2[i]);
                if( delta <= tolerance || cnr <= 2 )
                    continue;
            }
            if( mode ) // Monitor the difference between differences between new and old mean3 u_i
            {
                // Calculation of abs.difference of moving average differences at r and r-1
                delta = fabs(U_CV[i] - U_CVo[i]);
                if( delta <= tolerance || cnr <= 2 )
                    continue;
            }
        }
        // Divergence detected
        ICNud[nCNud++] = i;
        if( FirstTime )
        {
            if( trace )
            {
                ipm_logger->debug(" Tol = {} | uDD ITG = {}", tol_gen, pm.ITG);
            }
            if( base_param()->PSM >= 3 )
            {
                TNode::ipmlog_file->debug(" Tol = {} | uDD ITG = {}", tol_gen, pm.ITG);
            }
            FirstTime = false;
        }
        if( trace )
        {
            ipm_logger->debug("Divergent ICs: {} | ln_bi= {} | Tol= {} |",
                              char_array_to_string(pm.SB[i], MAXICNAME), log_bi, tolerance);
        }
        if( base_param()->PSM >= 3 )
        {
            TNode::ipmlog_file->debug("Divergent ICs: {} | ln_bi= {} | Tol= {} |",
                                    char_array_to_string(pm.SB[i], MAXICNAME), log_bi, tolerance);
        }
    } // for i
    return nCNud;
}

//--------------------- End of ipm_main.cpp ---------------------------
