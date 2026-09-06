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

#include "ms_multi.h"
#include "jama_lu.h"
#include "jama_cholesky.h"
#include "kinetics.h"
#include "v_service.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <chrono>

// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> TMultiBase::ipm_logger = spdlog::stdout_color_mt("ipm");

#define uDDtrace false

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Native-solver event trace - see the declaration in ms_multi.h for what it is
/// for and why it is env-gated rather than NDEBUG-gated.
///
/// Thread note: the open is a function-local static, so its initialisation is
/// thread-safe (C++11 magic statics); concurrent fprintf() calls on one FILE*
/// are individually locked by glibc, which is sufficient for a line-oriented
/// diagnostic. Interleaving between threads is possible and would show as
/// out-of-order events, not as corrupted lines.
FILE* native_trace_file()
{
    static FILE* fp = []() -> FILE*
    {
        const char* fn = std::getenv( "GEMS3K_NATIVE_TRACE_FILE" );
        return fn ? fopen( fn, "a" ) : nullptr;
    }();
    return fp;
}

/// Companion to native_trace_file() for the one thing an event-level trace
/// cannot show: the shape of the IPM descent, one line per iteration. See the
/// call site in InteriorPointsMethod() for what it measured and why.
FILE* ipm_probe_file()
{
    static FILE* fp = []() -> FILE*
    {
        const char* fn = std::getenv( "GEMS3K_IPM_PROBE" );
        return fn ? fopen( fn, "a" ) : nullptr;
    }();
    return fp;
}

/// Worst per-IC mass-balance residual of the CURRENT primal pm.Y, recomputed
/// locally rather than read out of pm.C[] - that array is MBR's own scratch and
/// is one update stale by the time MBR returns. Reports the same quantities
/// MBR's own convergence test compares (ipm_main.cpp, the pa_DT branches):
/// worst RELATIVE |C[i]|/(B[i]*DHBM) and worst ABSOLUTE |C[i]|, each with the
/// IC that carries it, so the "which test was binding" question the pa_DT /
/// `||`-vs-`&&` item turns on is answered at no extra instrumentation cost.
/// Scans the same [0, N - E) range MBR does, i.e. excluding the charge row.
/// Worst relative and worst absolute mass-balance residual of the amount vector `amt`,
/// over the ordinary IC range [0, Z) - the charge-balance IC in [Z, pm.N) is excluded,
/// exactly as MBR's own convergence loops exclude it. `rel` is normalised so that
/// rel > 1 means "this state fails the relative test MBR applies".
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

/// One MBRX line naming the exit path MassBalanceRefinement() actually took.
///
/// WHY THIS AND NOT JUST eRet. The return code does not distinguish "every IC
/// passed the test" from "gave up, and the strict check that would have reported
/// that did not apply". The latter is reached whenever the post-loop guard
///
///     if( pa_p->DW && ( WhereCalledFrom == 0L || pm.pNP ) )
///
/// is false - so a COLD call's second MBR (WhereCalledFrom = pm.K2 >= 1,
/// pm.pNP = 0) skips it whatever pa_DW is set to, and returns iRet = 0 on a
/// state its own per-IC test rejects. That is the mechanism behind native
/// returning an answer whose relative mass-balance residual exceeds pa_DHB, and
/// hence behind the eight projects whose own SIA cannot re-solve their own
/// converged answer: a cold call never re-checks what it returns, and a warm one
/// does. Observed directly on Resources/gems3k-fail/Al-species_G_sys_2_0_0_101,
/// where the cold call exits "lenient" at rel = 2.61x tolerance and the warm
/// call then correctly refuses the same state.
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
  // DATABR values as received, before any internal processing -- lets a caller-side
  // bug (e.g. a bulk composition silently zeroed before GEMS3K ever runs) be visually
  // distinguished from a solver-side one without needing a second historical build.
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

   // Exact-zero census of the answer this call returns. Pairs with the PSSC
   // lines above: an interior-point method on a box with a positive lower bound
   // cannot produce a single exact zero, so every zero counted here was written
   // by an explicit removal step, and the PSSC lines say which phase and why.
   // This is the observation behind the "native deletes the phase, it does not
   // solve it" reading of the psina failures.
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

   // ---- Mass-balance verdict on the ANSWER this call returns.
   //
   // WARN ONLY - deliberately, and this is a project-owner decision (2026-09-05),
   // not an oversight. Native's COLD path structurally never checks the state it
   // hands back: MBR's strict guard is gated on ( WhereCalledFrom == 0L || pm.pNP ),
   // so a cold call's SECOND MBR (K2 >= 1, pNP == 0) is exempt whatever pa_DW is and
   // returns iRet = 0 on a state its own per-IC test rejects. A warm call's FIRST MBR
   // has WhereCalledFrom == 0, the guard applies, and the same state is refused -
   // which is the whole mechanism behind the ten projects whose native SIA cannot
   // re-solve their own converged answer (plan v5 section 60.5).
   //
   // Dropping that clause was measured: it turns cold OK into FAIL on exactly the six
   // projects whose warm restart already fails. That is a user-facing behaviour change
   // on real projects, so the verdict is left alone and the FACT is surfaced instead.
   // A caller that wants to act on it has the message; one that does not is unaffected.
   // This is also the free signal an outcome-driven solver chooser needs (plan v5,
   // handoff item 8: "build the signal first").
   //
   // Unconditional, not NDEBUG-gated, for the same reason as the clamp warnings in
   // ipm_chemical.cpp: the event is rare and its whole point is to make a silently
   // accepted state visible in a production build. Suppressed when pm.MK/pm.PZ already
   // say the solution is bad, since testMulti() reports that case - the signal worth
   // having is the one on a run that would otherwise read as clean.
   //
   // Cost when it does not fire: one O(N*L) pass, against a solve that has just done
   // hundreds of them.
   if( !pm.MK && !pm.PZ )
   {
       long int iRel = -1, iAbs = -1; double rel = 0., absr = 0.;
       native_trace_mb_of( pm, pm.X, iRel, rel, iAbs, absr );
       const double dtExp  = fabs( (double)base_param()->DT );
       const double absCut = ( dtExp < 2. ) ? pm.DHBM : pow( 10., -dtExp );
       // Same per-IC rule the convergence branches apply: with DT == 0 only the
       // relative test exists; with DT != 0 an IC must exceed BOTH.
       const bool fails = !base_param()->DT ? ( rel > 1. )
                                            : ( rel > 1. && absr > absCut );
       if( fails && iRel >= 0 )
           gems_logger->warn(
               "GEM answer accepted with an unsatisfied mass balance: IC {} is {:.3e}x its own "
               "tolerance (|residual| {:.3e} mol against pa_DHB*b = {:.3e}); worst absolute "
               "|residual| {:.3e} mol at IC {}. Returned unchanged - native's cold path does "
               "not gate on this - but a warm (SIA) re-solve of the same state will reject it. "
               "Consider pa_DT (an absolute floor) or a tighter pa_DHB for this project.",
               char_array_to_string( pm.SB[iRel], MAXICNAME ),
               rel, fabs( pm.B[iRel] * pm.DHBM * rel ), pm.B[iRel] * pm.DHBM,
               absr, iAbs >= 0 ? char_array_to_string( pm.SB[iAbs], MAXICNAME )
                               : std::string( "-" ) );
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
                     pmbuf += char_array_to_string(pm.SF[k_miss],20);
                 std::string pubuf = std::to_string(k_unst)+ ": ";
                 if(k_unst >=0 )
                    pubuf += char_array_to_string(pm.SF[k_unst],20);

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
                       pmbuf += char_array_to_string(pm.SF[k_miss],20);
                  std::string pubuf = std::to_string(k_unst)+ ": ";
                  if(k_unst >=0 )
                      pubuf += char_array_to_string(pm.SF[k_unst],20);

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
                          buf += char_array_to_string(pm.SB[i],3)+" = "+std::to_string(pm.B[i]);
              setErrorMessage( 20, "W20IPM: IPM Main Descent:", buf.c_str());
    	   }
           else
           {
              addErrorMessage((std::string(", ")+char_array_to_string(pm.SB[i],3)+" = "+std::to_string(pm.B[i])).c_str());
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
        DC_RaiseZeroedOff( 0, pm.L );
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
                    buf += char_array_to_string( pm.SM[jK], MAXDCNAME);
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

    // Best-residual tracking. PORTED 2026-09-02 from branch ipm_contraints
    // commit d3ae685 (implemented and validated there 2026-08-21), which had
    // never been merged onto develop_optima - the same cross-branch gap as the
    // MBR Jacobi preconditioner, found by the same survey. The Tier A titrant
    // terms that version snapshots alongside Y do not exist on this branch, so
    // the metric below is written directly on pm.C[] - which is exactly what
    // that version's own helper reduces to when Tier A is inactive.
    //
    // WHY. Every early-exit path below - the LM-too-small break, the
    // stall-detector break, and simply exhausting pa_p->DP iterations - returns
    // whatever pm.Y happens to hold at that moment, which is the state AFTER
    // the last update was applied. Nothing guarantees that is the best state
    // MBR actually visited; in the oscillating case it is routinely worse than
    // an earlier iteration already reached. bestY snapshots pm.Y at whichever
    // iteration has the smallest ResidualMetric() seen so far, taken right
    // where the residual is measured and BEFORE that iteration's own update.
    //
    // UNCONDITIONAL, not behind a new BASE_PARAM flag. It runs only on paths
    // that are already non-ideal exits - the clean "balance residuals OK"
    // branch returns before reaching the restore - and it can only substitute
    // a state with a strictly SMALLER ResidualMetric() for the one about to be
    // returned, on the same metric the convergence and stall checks already
    // use. So it cannot make a converging run worse, and a toggle would only
    // add a way to keep returning a known-worse state. Unlike stall detection
    // itself - a heuristic that changes WHEN MBR gives up, and so needs its own
    // A/B - this changes only WHICH already-computed state is reported once MBR
    // has independently decided to give up.
    std::vector<double> bestY( pm.L );
    double bestResidual = std::numeric_limits<double>::max();
    bool haveBest = false, trRestored = false;

    // Worst normalized mass-balance residual across all N ICs (>1 means at
    // least one IC is outside its tolerance) - the same relative/absolute
    // combined tolerance logic as the per-IC convergence checks below, but
    // scanned over the FULL [0,N) range every time rather than from the
    // first-failing index onward, so that states from different iterations -
    // which can have different first-failing indices - stay comparable on one
    // consistent scale.
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
           // PER-IC-CLASS rule, opt-in, default OFF (pa_MbClassRule == 0 takes the
           // pre-existing branches below, byte for byte).
           //
           // This is what Kulik 2013 App. 2.2 actually prescribes: a RELATIVE
           // threshold for minor/trace ICs and an ABSOLUTE one for major ICs. The
           // branches below apply BOTH tests to EVERY IC when pa_DT != 0, which is
           // strictly stricter and never a per-class selection - see the field's
           // own comment in ms_multi.h for the measured consequences.
           //
           // Classification is by RATIO to the largest IC, not by an absolute
           // amount: pm.B[] is internally rescaled to pa_DG total moles, so an
           // absolute classification would not be scale-invariant.
           double maxB = 0.;
           for( I=0; I<Z; I++ )
               if( pm.B[I] > maxB ) maxB = pm.B[I];
           double AbsMbCutoff_cls;
           {
               const double e = fabs( (double)pa_p->DT );
               // |DT| >= 2 names the major absolute cutoff explicitly and is the
               // recommended usage; otherwise fall back so the switch works alone.
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
       else { // combined balance accuracy - an absolute FLOOR under the relative test
           // Each IC must exceed BOTH thresholds to count as not converged, i.e.
           // the effective bar is max(absolute, relative) per IC. That makes the
           // absolute cutoff a FLOOR: a trace IC, whose relative bar B[I]*DHBM is
           // vanishingly small, is judged on the absolute one instead, while a
           // major IC keeps its relative bar because that is already the larger
           // of the two. This is the per-IC-class behaviour pa_DT has always
           // documented and never had.
           //
           // CHANGED 2026-09-02 from `||` to `&&`. With `||` an IC failed if
           // EITHER threshold was exceeded, so pa_DT != 0 was strictly STRICTER
           // than pa_DT == 0 and could never relax anything - confirmed
           // empirically before the change (pa_DT = -9 on scratch copies of
           // Al-species, FeNaCl_FyGt_Precip and 07PSIna_G_iron changed nothing at
           // all, because their relative test was already the binding one).
           //
           // WHY IT MATTERS, measured with tools/trace_ladder on
           // Resources/gems3k-fail/07PSIna_G_iron (scale bIC[Fe] across a ladder,
           // everything else bit-identical, fresh TNode per rung): native's own
           // SIA cannot re-solve its own converged answer at any Fe below 3e-7,
           // and EVERY failing rung has an absolute residual between 1e-14 and
           // 1e-8 mol - 0.5 picomole at Fe = 3e-10, 35 femtomoles at 3e-12. No
           // physical criterion would reject those; they fail only because the
           // test is relative and Fe's own total is tiny.
           //
           // AND IT DOES NOT MAKE THE SYSTEM DETERMINATE - record it that way.
           // The same ladder, reporting SIGNED H and O residuals, shows the error
           // lying along the WATER direction (H/O ~ 2.0, same sign, in five of
           // eight rungs) - the same near-rank-1 H2O H:O = 2:1 dependence already
           // diagnosed for MBR's own Schur-complement matrix - and the worst rung
           // (Fe = 3e-9) is exactly where the redox constraint evaporates
           // (H2(aq) collapses to 8.7e-27 with O2(aq) still 0, neither couple
           // present). So this makes the TEST physical; it does not resolve the
           // degeneracy the test is detecting.
           //
           // Corpus-wide no-op at the time of the change: all 52 projects in
           // Resources/{gems3k,gems3k-fail,gems3k-psina} either set pa_DT = 0 or
           // omit it, so none of them takes this branch at all. A project must
           // opt in by setting pa_DT (|DT| >= 2 names the absolute cutoff as
           // 10^-|DT|, which is the recommended usage).
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
       // Best-residual snapshot. Taken here because pm.C[] is this iteration's
       // freshly measured residual for the PRE-update pm.Y (the LM-scaled
       // update below has not been applied yet) and AbsMbCutoff_stall has just
       // been set by the residual test above - see the declaration for why.
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
                      buf += ")): Too small LM step size - cannot converge (check Pa_DG?)";
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
              gems_logger->warn("MBR({}): stall at IT1={} stalledIter={} "
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

    // Best-residual restore. Every path that reaches this point is a non-ideal
    // exit - the LM-too-small break, the stall-detector break, or exhausting
    // pa_p->DP iterations - and each returns whatever pm.Y was left at, with no
    // guarantee it is the best state MBR actually visited. Recompute the
    // residual of the state about to be returned and put back the in-loop
    // snapshot if it is strictly better. See the bestY/ResidualMetric
    // declarations above for why this is unconditional.
    if( haveBest )
    {
        MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C );
        double curResidualFinal = ResidualMetric();
        if( bestResidual < curResidualFinal )
        {
            gems_logger->warn("MBR({}): restoring best-residual state (best={:.3e} vs current={:.3e})",
                               WhereCalledFrom, bestResidual, curResidualFinal);
            for( j=0; j<pm.L; j++ )
                pm.Y[j] = bestY[j];
            MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C );
            trRestored = true;
        }
    }

    //  Prescribed mass balance precision cannot be reached
                    // Temporary workaround for pathological systems 06.05.2010 DK
   if( pa_p->DW && ( WhereCalledFrom == 0L || pm.pNP ) )  // Now controlled by DW flag
   {  // Strict mode of mass balance control
       iRet = 2;
       std::string buf = "(MBR("+std::to_string(WhereCalledFrom);
                   buf += ")) Maximum allowed number of MBR iterations (";
                   buf += std::to_string(pa_p->DP) +") exceeded! ";
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
    // Noise-stall accept, gated on pa_IpmStallWindow (default 0 = off). See the
    // field's own comment in ms_multi.h for the mechanism and the measurement.
    //
    // FOUR SIGNALS, ALL REQUIRED, AND NONE OF THEM REFERENCES pa_DK. That is the
    // whole point: an earlier version guarded on "pm.PCI is already within 30x of
    // pm.DXM", which works but ties the accepted ACCURACY to the tolerance the
    // rule is meant to replace - and flipping its default then broke
    // proposed.aop's "native G is settings-independent" assertion. Each signal
    // below is individually insufficient and was individually measured unsafe:
    //   energy alone      -> 115 % error on f_/j_TestSUP98 (plan v5 74.9)
    //   criterion alone   -> 11 % error, fires on 25 of 42 (80.2)
    //   mass balance      -> unusable here: it is established by MBR and then
    //                        DEGRADES monotonically inside the IPM loop by design
    // Together they are safe, because a transient plateau in one is not a
    // simultaneous plateau in all.
    const long int kIpmStallMaxW = 200;
    const double kIpmStallFXTol   = 1.e-9;  // energy flat over the window
    const double kIpmStallCompTol = 1.e-6;  // sumX and maxX flat over the window
    const double kIpmStallIncLo   = 0.35;   // PCI increases on 35-65 % of steps,
    const double kIpmStallIncHi   = 0.65;   //   i.e. it is bouncing, not moving
    std::vector<double> stW_pci, stW_fx, stW_sum, stW_max;
    const BASE_PARAM *pa_p = base_param();

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
   " It is not possible to obtain a valid GEM IPM solution.\n"  );
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
            addErrorMessage((std::string(" ")+char_array_to_string(pm.SB[ICNud[Z]],6)).c_str());
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

        // Per-iteration descent record, env-gated (GEMS3K_IPM_PROBE=<path>).
        //
        // WHY IT EXISTS. The IPM loop terminates on pm.PCI <= pm.DXM, and on part
        // of this corpus PCI stops carrying signal long before that test is met:
        // once the composition has settled, PCI is computed from differences of
        // nearly-equal numbers and becomes noise, wandering in a band whose median
        // sits ABOVE DXM with only its lower tail below. Termination is then a
        // waiting time for a lucky draw, which is what makes the iteration count
        // irreproducible under a 1e-15 bIC nudge while the answer is not. That is
        // only visible per iteration - the event-level trace above cannot show it -
        // and it is measured from this record's PCI/FX columns:
        //
        //   f_Kaolinite, 9 nudges: energy final to 1e-11 by iteration 25 on EVERY
        //   run; PCI thereafter median 6.6e-6 against DXM 1e-6, rising on 49 % of
        //   steps (no trend), P(PCI <= DXM) = 0.0063 per iteration. Predicted mean
        //   25 + 1/0.0063 = 184.8 against an observed 183.8 over 62-363.
        //
        // Cost when unset: one null test on an already-resolved static pointer,
        // the same shape as native_trace_file(). L is small on the projects this
        // diagnoses; do not use it on the 1392-species giants without redirecting
        // to scratch.
        if( FILE* ipf = ipm_probe_file() )
        {
            double sumX = 0., minX = 1e300, maxX = 0.;
            for( long int jj = 0; jj < pm.L; jj++ )
            {
                sumX += pm.X[jj];
                if( pm.X[jj] > maxX ) maxX = pm.X[jj];
                if( pm.X[jj] > 0. && pm.X[jj] < minX ) minX = pm.X[jj];
            }
            // Mass-balance residual of the CURRENT primal pm.X, normalised so that
            // mbRel > 1 means "this state fails the per-IC test MBR applies". It is
            // measured against pa_DHB - a PHYSICAL tolerance - which is what makes it
            // usable as an independent second signal alongside PCI (a numerical one).
            // O(N*L) per iteration, so this probe is for diagnosis, not for the giants.
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
                if( spreadOK( stW_fx,  kIpmStallFXTol )
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
            "IPM convergence criterion tolerance (Pa_DK) could not be reached"
    		" (more than Pa_IIM iterations done);\n" );
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
            buf +=  char_array_to_string(pm.SB[i],3);
            setErrorMessage( 2, "E02IPM: PSSC(): ", buf.c_str());
        }
        else
        {
            addErrorMessage((std::string(", ")+char_array_to_string(pm.SB[i],3)).c_str());
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
// The following diagnostics (condition-number estimation + per-phase timing)
// add real overhead per linear solve — up to several hundred percent of the
// Cholesky/LU solve itself for small systems (see gems-benchmark/CLAUDE.md,
// 2026-07-28). Gated so normal/production use of GEMS3K compiles none of
// this in; only builds that explicitly opt in (GEMS3K's CMake option
// ENABLE_BENCHMARK_DIAGNOSTICS, used by gems-benchmark) pay the cost. Future
// benchmark/diagnostics-only instrumentation should reuse this same macro
// rather than introducing a new one per feature.

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
    auto solve_t0 = std::chrono::high_resolution_clock::now();
    double diag_ms = 0.;   // condition-number-diagnostics-only time, this call
    pm.SolveCallCount++;
#endif

    long int ii, i, jj, kk, k, Na = pm.N;

    // ------------------------------------------------------------------
    // Appendix A of Leal et al. (2017), Eqs. 130-136: a pivot/non-pivot
    // SPLIT of this reduction, gated on pa_MbPivotSplit (default 0 = off).
    // See BASE_PARAM::MbPivotSplit (ms_multi.h) for the derivation, the
    // classification rule and the honest bound on what it can achieve.
    //
    // Implemented as a self-contained early branch rather than by widening
    // AA/BB, so that with the field off this function is byte-identical to
    // what it was before - AA's stride is N everywhere below, and changing
    // it would touch every indexing site in the naive path.
    //
    // Falls through to that naive path when the split is empty (Eq. 136
    // then degenerates to Eq. 132 exactly, so there is nothing to gain and
    // the existing diagnostics are worth keeping) or when the non-pivot set
    // is larger than N (the augmented solve would then be more than twice
    // the size of the thing it replaces, against the paper's own claim that
    // |I_n| is "typically small and not greater than the number of
    // elements" - a system that violates that is not the case Appendix A
    // was written for).
    // ------------------------------------------------------------------
    if( initAppr && base_param()->MbPivotSplit )
    {
        // Eq. 133b, with D_jj = 1/W[j] and C_(i,j) = a(j,i):
        //     non-pivot  <=>  1/W[j] < max_i |a(j,i)|
        std::vector<long int> np;              // I_n, in species order
        std::vector<char> isNP( (size_t)pm.L, 0 );
        for( jj = 0; jj < pm.L; jj++ )
        {
            if( pm.Y[jj] <= min( pm.lowPosNum, pm.DcMinM ) )
                continue;                      // same filter as the assembly
            if( !( pm.W[jj] > 0. ) )
                continue;
            double colmax = 0.;
            for( i = arrL[jj]; i < arrL[jj+1]; i++ )
            {   ii = arrAN[i];
                if( ii >= N )
                    continue;
                double v = fabs( a(jj,ii) );
                if( v > colmax ) colmax = v;
            }
            if( colmax > 0. && pm.W[jj] * colmax > 1. )
            {   np.push_back( jj );  isNP[(size_t)jj] = 1;  }
        }

        const long int nn = (long int)np.size();
        if( nn > 0 && nn <= N )
        {
            const long int M = N + nn;
            std::vector<double> AM( (size_t)M * (size_t)M, 0. );
            std::vector<double> BM( (size_t)M, 0. );
            // AM is row-major: AM[r*M + c] is row r, column c. The matrix is
            // symmetric, so this agrees with the naive path's own (column-
            // major) convention regardless.

            // Top-left N x N: the same Gram matrix, over PIVOT species only.
            for( jj = 0; jj < pm.L; jj++ )
            {
                if( pm.Y[jj] <= min( pm.lowPosNum, pm.DcMinM ) )
                    continue;
                if( isNP[(size_t)jj] )
                    continue;
                for( k = arrL[jj]; k < arrL[jj+1]; k++)
                    for( i = arrL[jj]; i < arrL[jj+1]; i++ )
                    {   ii = arrAN[i];
                        kk = arrAN[k];
                        if( ii >= N || kk >= N )
                            continue;
                        AM[(size_t)ii*(size_t)M + (size_t)kk] += a(jj,ii) * a(jj,kk) * pm.W[jj];
                    }
            }

            // Coupling blocks C_n / B_n and the non-pivot diagonal D_n.
            // Row N+t is the retained equation  sum_k a(j,k) y_k - x_t/W[j] = 0
            // (negated from Eq. 134's own row so the whole matrix stays
            // symmetric; the right-hand side is zero either way).
            for( long int t = 0; t < nn; t++ )
            {
                const long int j = np[(size_t)t];
                const long int r = N + t;
                for( i = arrL[j]; i < arrL[j+1]; i++ )
                {   ii = arrAN[i];
                    if( ii >= N )
                        continue;
                    const double v = a(j,ii);
                    AM[(size_t)ii*(size_t)M + (size_t)r] = v;
                    AM[(size_t)r *(size_t)M + (size_t)ii] = v;
                }
                AM[(size_t)r*(size_t)M + (size_t)r] = -1. / pm.W[j];
            }

            for( ii = 0; ii < N; ii++ )
                BM[(size_t)ii] = pm.C[ii];     // BM[N..M) stay 0

            // Same symmetric Jacobi scaling as the naive path - see the block
            // below. Applied here too on purpose: comparing an unpreconditioned
            // Appendix A against a preconditioned baseline would measure the
            // loss of the preconditioner rather than the gain of the split.
            std::vector<double> Ds( (size_t)M, 1. );
            for( long int r = 0; r < M; r++ )
            {
                double d = fabs( AM[(size_t)r*(size_t)M + (size_t)r] );
                if( d > 1e-300 )
                    Ds[(size_t)r] = 1. / sqrt( d );
            }
            for( long int r = 0; r < M; r++ )
            {
                for( long int c = 0; c < M; c++ )
                    AM[(size_t)r*(size_t)M + (size_t)c] *= Ds[(size_t)r] * Ds[(size_t)c];
                BM[(size_t)r] *= Ds[(size_t)r];
            }

            // The augmented matrix is symmetric but INDEFINITE (the non-pivot
            // diagonal block is -1/W[j] < 0), so Cholesky cannot apply and is
            // not attempted; LU with partial pivoting is the whole point of
            // Appendix A in the first place.
            Array2D<double> AAm( M, M, AM.data() );
            Array1D<double> BBm( M, BM.data() );
            JAMA::LU<double> lum( AAm );
            if( !lum.isNonsingular() )
            {
                ipm_logger->warn("MakeAndSolveSystemOfLinearEquations (Appendix A "
                                 "pivot split, |I_n|={}): augmented matrix singular", nn);
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
                pm.SolveTimeMs += std::chrono::duration<double, std::milli>(
                                      std::chrono::high_resolution_clock::now() - solve_t0 ).count();
                pm.CondNumTimeMs += diag_ms;
#endif
                return 1;
            }
            BBm = lum.solve( BBm );
            for( ii = 0; ii < N; ii++ )
                pm.Uefd[ii] = BBm[(int)ii] * Ds[(size_t)ii];
            ipm_logger->trace("Appendix A pivot split: |I_n|={} of {} active species, "
                              "augmented system {} x {}", nn, pm.L, M, M);
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
            pm.SolveTimeMs += std::chrono::duration<double, std::milli>(
                                  std::chrono::high_resolution_clock::now() - solve_t0 ).count();
            pm.CondNumTimeMs += diag_ms;
#endif
            return 0;
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

    // Diagonal (Jacobi) preconditioning of the initAppr (MBR) Schur-complement
    // matrix. PORTED 2026-09-02 from branch ipm_contraints commit 26d9d57a,
    // where it was implemented and validated on 2026-08-21 - it had never been
    // merged onto develop_optima, which branched from master, so this branch's
    // native MBR ran unpreconditioned for the whole of the Optima work. Same
    // cross-branch trap already on record twice (PSTALL, and
    // ENABLE_BENCHMARK_DIAGNOSTICS).
    //
    // On at least one real project (Cu-Pourbaix) AA's diagonal spans up to ~17
    // orders of magnitude - diagMax stays flat ~4.4e5 while diagMin collapses
    // ~10x per iteration toward ~1e-12 - which is a SCALING problem, not rank
    // deficiency. Symmetric scaling A' = D*A*D, B' = D*B with
    // D = diag(1/sqrt(|A_ii|)) preserves symmetry and positive-definiteness
    // (A is a Gram matrix, a(j,i)*a(j,k)*W[j] summed over j, so A' is one too)
    // while normalising every nonzero diagonal entry to exactly 1. The solved
    // dual is unscaled back (U = D*U') at the unpack site below.
    //
    // Rows/columns with a structurally-zero diagonal (no species touches that
    // IC) get Dscale = 1 rather than dividing by zero - such a row already
    // makes the matrix singular regardless of preconditioning, and the existing
    // singular-matrix path handles it unchanged.
    //
    // Measured when first adopted: the real eigenvalue-based condition number
    // drops ~7 orders (~3.4e17 implied -> ~3.06e10 measured) with zero
    // regression across all 25 gems-benchmark Resources/gems3k projects. It
    // does NOT touch the remaining ~10 orders of NON-diagonal ill-conditioning,
    // which was traced to water's fixed H:O = 2:1 stoichiometry making the H
    // and O rows of the assembled matrix near-linearly-dependent - a genuine
    // near-singularity of A D^-1 A^T, which no rescaling can repair.
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
    auto diag_t0 = std::chrono::high_resolution_clock::now();
    pm.CondNumDiag = std::max( pm.CondNumDiag, DiagRatioConditionProxy( AA, N ) );
    double lambda_max = PowerIterationMaxEig( AA, N, 6 );
    diag_ms += std::chrono::duration<double, std::milli>(
                   std::chrono::high_resolution_clock::now() - diag_t0 ).count();
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
        auto inv_t0 = std::chrono::high_resolution_clock::now();
        double lambda_min = InverseIterationMinEig( chol, N, 6 );
        diag_ms += std::chrono::duration<double, std::milli>(
                       std::chrono::high_resolution_clock::now() - inv_t0 ).count();
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
                                  std::chrono::high_resolution_clock::now() - solve_t0 ).count();
            pm.CondNumTimeMs += diag_ms;
#endif
            return 1; // Singular matrix - too bad! No solution ...
        }

        B = lu.solve( B );
#ifdef GEMS3K_BENCHMARK_DIAGNOSTICS
        auto inv_t0 = std::chrono::high_resolution_clock::now();
        double lambda_min = InverseIterationMinEig( lu, N, 6 );
        diag_ms += std::chrono::duration<double, std::milli>(
                       std::chrono::high_resolution_clock::now() - inv_t0 ).count();
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
                          std::chrono::high_resolution_clock::now() - solve_t0 ).count();
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
        ipm_logger->info("UD3 trace: {}  SIA={}  Itr   C_D:  {}",
                          char_array_to_string(pm.stkey, EQ_RKLEN), pm.pNP, char_array_to_string(pm.SB1[0],MAXICNAME));
    }
    if( base_param()->PSM >= 3 )
    {
      TNode::ipmlog_file->info(" UD3 trace: {}  SIA= {} Itr   C_D: {}",
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
       TNode::ipmlog_file->info("ncrement_uDD {}  {}", r, pm.PCI);
    }
    if( trace )
    {
        ipm_logger->info("ncrement_uDD {}  {}", r, pm.PCI);
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
        ipm_logger->info("U={}  U_mean={} U_CV={}", pm.U[i],  U_mean[i], U_CV[i]);
      }
      if( base_param()->PSM >= 3 )
      {
         TNode::ipmlog_file->info("U={}  U_mean={} U_CV={}", pm.U[i],  U_mean[i], U_CV[i]);

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
                ipm_logger->info(" Tol = {} | uDD ITG = {}", tol_gen, pm.ITG);
            }
            if( base_param()->PSM >= 3 )
            {
                TNode::ipmlog_file->info(" Tol = {} | uDD ITG = {}", tol_gen, pm.ITG);
            }
            FirstTime = false;
        }
        if( trace )
        {
            ipm_logger->info("Divergent ICs: {} | ln_bi= {} | Tol= {} |",
                              char_array_to_string(pm.SB[i], MAXICNAME), log_bi, tolerance);
        }
        if( base_param()->PSM >= 3 )
        {
            TNode::ipmlog_file->info("Divergent ICs: {} | ln_bi= {} | Tol= {} |",
                                    char_array_to_string(pm.SB[i], MAXICNAME), log_bi, tolerance);
        }
    } // for i
    return nCNud;
}

//--------------------- End of ipm_main.cpp ---------------------------
