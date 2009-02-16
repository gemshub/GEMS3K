//-------------------------------------------------------------------
// $Id: ipm_main.cpp 705 2006-04-28 19:39:01Z gems $
//
// Copyright (C) 1992-2007 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
// Uses: JAMA/C++ Linear Algebra Package based on the Template
// Numerical Toolkit (TNT) - an interface for scientific computing in C++,
// (c) Roldan Pozo, NIST (USA), http://math.nist.gov/tnt/download.html
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry and of the
// standalone GEMIPM2K code.
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//

#include "m_param.h"
#include "jama_lu.h"
#include "jama_cholesky.h"
using namespace TNT;
using namespace JAMA;

#ifndef IPMGEMPLUGIN
#include "service.h"
#include "stepwise.h"
#endif
#include "node.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Main sequence of IPM calculations
//  Main place for implementation of diagnostics and setup
//  of IPM precision and convergence
//  rLoop is the index of the primal solution refinement loop (for tracing)
//   or -1 if this is main GEMIPM2 call
//
void TMulti::MultiCalcMain( long int rLoop )
{
    long int i, j, k, eRet, iB, status=0;
    SPP_SETTING *pa = &TProfil::pm->pa;

    pmp->W1=0; pmp->K2=0;               // internal counters
    pmp->Ec = pmp->MK = pmp->PZ = 0;    // Return codes
    setErrorMessage( 0, "" , "");  // empty error info
    if( TProfil::pm->pa.p.PLLG == 0 )  // SD 10/02/2009
    	TProfil::pm->pa.p.PLLG = 32000;


    if( pmp->pULR && pmp->PLIM )
        Set_DC_limits( DC_LIM_INIT );

    if( rLoop < 0 )
    {   // cleaning the f_alpha vector (phase stability criteria) - not done during refinement loops
        for( k=0; k<pmp->FI; k++ )
            pmp->Falp[k] = 0.0;
        for( j=0; j<pmp->L; j++ )
        {
        	pmp->EMU[j] = 0.0;
        	pmp->NMU[j] = 0.0;
        }
    }

    // testing the entry into feasible domain
mEFD:
	if( rLoop>=0 || pmp->W1)
     {
        for( i=0; i<pmp->N; i++ )
           pmp->Uc[i]=pmp->U[i];
     }

    if(pa->p.DW && pmp->W1)
     {
        for( j=0; j<pmp->L; j++ )
           pmp->Y[j]=pmp->X[j];
        TotalPhases( pmp->Y, pmp->YF, pmp->YFA );
     }

     eRet = EnterFeasibleDomain( ); // Here the IPM-2 EFD() algorithm is called

#ifdef GEMITERTRACE
to_text_file( "MultiDumpC.txt" );   // Debugging
#endif

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
    pVisor->Update(false);
#endif
// STEPWISE (2)  - stop point to examine output from EFD()
STEP_POINT("After FIA");
#endif

    switch( eRet )
    {
     case 0:  // OK
              break;
     case 5:  // Initial Lagrange multiplier for metastability broken for DC
     case 4:  // Initial mass balance broken for IC
     case 3:  // too small step length in descent algorithm
     case 2:  // max number of iterations has been exceeded in EnterFeasibleDomain()
     case 1: // degeneration in R matrix  for EnterFeasibleDomain()
   	         if( pmp->pNP )
   	         {   // bad PIA mode - trying the AIA mode
       	        pmp->MK = 2;   // Set to check in calcMulti() later on
#ifdef GEMITERTRACE
f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT << " ! PIA->AIA on E04IPM" << endl;
#endif
               goto FORCED_AIA;
   	         }
   	         else
   	        	 Error( pmp->errorCode ,pmp->errorBuf );
              break;
    }

   // call of the main IPM-2 minimization algorithm
   eRet = InteriorPointsMethod( status, rLoop );

#ifdef GEMITERTRACE
to_text_file( "MultiDumpD.txt" );   // Debugging
#endif

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
    pVisor->Update(false);
#endif
// STEPWISE (3)  - stop point to examine output from IPM()
   STEP_POINT("After IPM");
#endif

// Diagnostics of IPM results
   switch( eRet )
   {
     case 0:  // OK
              break;
     case 2:  // max number of iterations has been exceeded in InteriorPointsMethod()
     case 1: // degeneration in R matrix  for InteriorPointsMethod()
         if( pmp->pNP )
         {   // bad PIA mode - trying the AIA mode
	        pmp->MK = 2;   // Set to check in calcMulti() later on
#ifdef GEMITERTRACE
f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT << " ! PIA->AIA on E06IPM" << endl;
#endif
	        goto FORCED_AIA;
         }
         else
        	 if( pmp->DX < 1e-4 || pmp->DHBM < 1e-6 )
             {  // Attempt to get result using most tolerant thresholds for convergence and balance accuracy
        	    pmp->DX = 1e-4;
                pmp->DHBM = 1e-6;
#ifdef GEMITERTRACE
f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT << " AIA: DX->1e-4, DHBM->1e-6 on E06IPM" << endl;
#endif
                goto mEFD;
             }
       	 Error( pmp->errorCode ,pmp->errorBuf );
         break;
     case 3:  // bad GammaCalc() status in PIA mode
     case 4: // Mass balance broken after DualTh recover of DC amounts
         if( pmp->pNP )
         {   // bad PIA mode - trying the AIA mode
	        pmp->MK = 2;   // Set to check in calcMulti() later on
	        goto FORCED_AIA;
         }
       	 Error( pmp->errorCode ,pmp->errorBuf );
         break;
   }

    pmp->FI1 = 0;
    pmp->FI1s = 0;
    for( i=0; i<pmp->FI; i++ )
        if( pmp->YF[i] > 1e-18 )
        {
            pmp->FI1++;
            if( i < pmp->FIs )
                pmp->FI1s++;
        }

    if( !pa->p.PC || pmp->PZ )    //  No PhaseSelect() operation allowed or
    {   if( pmp->PD >= 2 )           //  pmp->PZ == 2 no impotant error
           for( i=0; i<pmp->L; i++)
             pmp->G[i] = pmp->G0[i];
        return;  // solved
    }

    if( pa->p.PC == 2 )
        XmaxSAT_IPM2();  // Install upper limits to xj of surface species

//========= calling Selekt2 algorithm =======
   long int ps_rcode, k_miss, k_unst;

 if( rLoop < 0 )
    ps_rcode = PhaseSelect( k_miss, k_unst, rLoop );
else
   ps_rcode = 1; // do not call Selekt2() in primal solution refinement loops (experimental!!!!!)

iB = CheckMassBalanceResiduals( pmp->X );
if( iB >= 0 )
{
    addErrorMessage( " (After PhaseSelect())" );
    if( pmp->pNP )
    {   // bad PIA mode - trying the AIA mode
       pmp->MK = 2;   // Set to check in calcMulti() later on
       goto FORCED_AIA;
    }
  	Error( pmp->errorCode ,pmp->errorBuf );
}
    //   MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->X, pmp->B, pmp->C);

   switch( ps_rcode )
   {
      case 1:   // IPM solution is final and consistent, no phases were inserted
                break;
      case 0:   // some phases were inserted and a new IPM loop is needed
#ifdef GEMITERTRACE
f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " K2=" << pmp->K2 <<
      " k_miss=" << k_miss << " k_unst=" << k_unst <<  " ! (new Selekt loop)" << endl;
#endif
    	        goto mEFD;
      default:
      case -1:  // the IPM solution is inconsistent after 3 Selekt2() loops
      {
          gstring pmbuf("");
          if(k_miss >=0 )
           pmbuf = gstring(pmp->SF[k_miss],0,20);
          gstring pubuf("");
          if(k_unst >=0 )
          pubuf = gstring(pmp->SF[k_unst],0,20);
          char buf[400];
          sprintf( buf,
   " The computed phase assemblage remains inconsistent after 3 PhaseSelect() loops.\n"
   " Problematic phase(s): %ld %s   %ld %s \n"
   " Probably, vector b is not balanced, or DC stoichiometries\n"
   " or standard g0 are wrong, or some relevant phases or DCs are missing.",
                   k_miss, pmbuf.c_str(), k_unst, pubuf.c_str() );
          setErrorMessage( 8, "W08IPM: PhaseSelect():", buf );
    	  if( pmp->pNP )
    	  {   // bad PIA mode - there are inconsistent phases after 3 attempts. Attempting AIA mode
        	  pmp->MK = 2;   // Set to check in calcMulti() later on
#ifdef GEMITERTRACE
f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT <<
      " k_miss=" << k_miss << " k_unst=" << k_unst << " ! PIA->AIA on E08IPM (Selekt)" << endl;
#endif
        	  goto FORCED_AIA;
    	  }
    	  else
    	  { pmp->PZ = 3; // IPM solution could not be improved in PhaseSelect()
    	                 //   some inconsistent phases remain
    	    return;
          }
      }
   }

#ifndef IPMGEMPLUGIN
   pmp->t_end = clock();
   pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);
// STEPWISE (4) Stop point after PhaseSelect()
   STEP_POINT("Before Refine()");
#ifndef Use_mt_mode
   pVisor->Update( false );
#endif
#endif
   if( pmp->MK == 2 )
       goto FORCED_AIA;

   if( pa->p.DW )
   {
     if( !pmp->W1 )
     {
       pmp->W1++;            // IPM-2 precision refinement - 1st run
       goto mEFD;
     }
     else
       if( pmp->W1 <  pa->p.DW )
       {
          for(i=0;i<pmp->N-pmp->E;i++)
          {
            if( fabs(pmp->C[i]) > pmp->DHBM // * pa->p.GAS
            || fabs(pmp->C[i]) > pmp->B[i] * pa->p.GAS )
            {
               if(pmp->W1 < pa->p.DW-1)
               {
                 pmp->W1++;  // IPM-2 precision enhancement - further runs
                 goto mEFD;
               }
           	   if( pmp->pNP )
               {   // bad PIA mode - precision cannot be reached. Attempting the AIA mode
                 	  pmp->MK = 2;   // Set to check in calcMulti() later on
#ifdef GEMITERTRACE
  f_log << " ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT <<
      " W1= " << pmp->W1 << " ! PIA-AEA on Prec.refin.errors E09IPM or E10IPM" << endl;
#endif
                     goto FORCED_AIA;
                }
             	char  buf[200];
                long int jj=0;
                for( j=i; j<pmp->N-pmp->E; j++ )
                   if( fabs(pmp->C[j]) > pmp->B[j] * pa->p.GAS )
                   {
          			  if(!jj )
          			  {
           				 sprintf(buf, "Prescribed mass balance accuracy cannot be reached\n"
           						      "for trace independent components:  %-2.2s  ", pmp->SB[j]);
          				 setErrorMessage( 9, "W09IPM: IPM-main(): ",buf);
                         jj = 1;  // Error state is activated
          			  }
          			  else
          			  {
           				 sprintf(buf,"  %-2.2s" ,  pmp->SB[j] );
          				 addErrorMessage(buf);
          			  }
                    }
                 if( !jj )
                 {
                	 setErrorMessage( 10, "E10IPM: IPM-main(): " ,
                       "Inconsistent GEM solution: Inaccurate mass balance\n"
                       "for some major independent components " );
               	  	Error( pmp->errorCode ,pmp->errorBuf );

                 }
                 pmp->PZ = 1; // no good solution
                 break;
            } // if
          } // end of i loop
      }
   }
#ifdef GEMITERTRACE
f_log << "ITF=" << pmp->ITF << " ITG=" << pmp->ITG << " IT=" << pmp->IT << " MBPRL="
   << pmp->W1 << " rLoop=" << rLoop;
    if( pmp->pNP )
	   f_log << " Ok after PIA" << endl;
    else
	   f_log << " Ok after AIA" << endl;
#endif

FORCED_AIA:
   for( i=0; i<pmp->L; i++)
      pmp->G[i] = pmp->G0[i];
   // Normal return after successful improvement of the mass balance precision
   pmp->t_end = clock();   // Fix pure runtime
   pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);
}

//Call for IPM iteration sequence
void TMulti::MultiCalcIterations( long int rLoop )
{
#ifdef GEMITERTRACE
to_text_file( "MultiDumpB.txt" );   // Debugging
#endif

     MultiCalcMain( rLoop );
if( !pmp->pNP )
     pmp->ITaia = pmp->IT;

#ifdef GEMITERTRACE
to_text_file( "MultiDumpE.txt" );   // Debugging
#endif

    // calculation of demo data for gases
    for( long int ii=0; ii<pmp->N; ii++ )
        pmp->U_r[ii] = pmp->U[ii]*pmp->RT;
    GasParcP();

}

// Finding whether the automatic initial approximation is necessary for
// launching the IPM algorithm.
//
// Uses a modified simplex method with two-side constraints (Karpov ea 1997)
//
// Return code:
// false - OK for IPM
// true  - OK solved
//
bool TMulti::AutoInitialApprox(  )
{
    long int i, j, NN, eCode=-1L;
    double minB, sfactor;
    char buf[300];
    SPP_SETTING *pa = &TProfil::pm->pa;

#ifdef GEMITERTRACE
to_text_file( "MultiDumpA.txt" );   // Debugging
#endif

// Scaling the IPM numerical controls for the system total amount and minimum b(IC)
    NN = pmp->N - pmp->E;
    minB = pmp->B[0]; // pa->p.DB;
    for(i=0;i<NN;i++)
    {
        if( pmp->B[i] < pa->p.DB )
        {
           if( eCode < 0  )
    	   {
              eCode = i;  // Error state is activated
              pmp->PZ = 3;
              sprintf(buf, "Too small input amount of independent component %-3.3s = %lg",
            	pmp->SB[i], pmp->B[i] );
    		  setErrorMessage( 20, "W20IPM: IPM-main():" ,buf);
    	   }
           else
           {
        	  sprintf(buf,"  %-3.3s = %lg" ,  pmp->SB[i], pmp->B[i] );
        	  addErrorMessage( buf );
           }
           pmp->B[i] = pa->p.DB;
        }
    	if( pmp->B[i] < minB )
           minB = pmp->B[i];      // Looking for the smallest IC amount
    }

    if( eCode >= 0 )
    {
       TProfil::pm->testMulti();
       pmp->PZ = 0;
       setErrorMessage( -1, "" , "");
    }

    //  check Ymin (cutoff)
   if(pmp->lowPosNum>minB*0.01)
      pmp->lowPosNum=minB*0.01;

   sfactor = calcSfactor();
   pmp->DHBM = sfactor * pa->p.DHB;  // 2.5 root
   pmp->DHBM *= (0.097+0.95/(1+exp(-(log10(minB)+11)/0.4)));

   pmp->DX = pa->p.DK;
   pmp->DX *= sfactor;
   pmp->DX *= (0.097+0.95/(1+exp(-(log10(minB)+6.1)/0.54)));
   if( pmp->DX < 0.01 * pa->p.DK )
       pmp->DX = 0.01 * pa->p.DK;
   pmp->DSM = pa->p.DS;  // Shall we add  * sfactor ?

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
   pVisor->Update(false);
#endif
#endif

// Analyzing if the Simplex approximation is necessary
    if( !pmp->pNP  )
    {   // Preparing to call Simplex method
        pmp->FitVar[4] = pa->p.AG;  //  initializing the smoothing parameter
        pmp->ITaia = 0;             // resetting the previous number of AIA iterations
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
//      pmp->IC = 0.0;  For reproducibility of simplex FIA?
        if( pmp->FIs )
            GammaCalc(LINK_FIA_MODE);
        if( pa->p.PC == 2 )
           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
        pmp->IT = 0; pmp->ITF += 1; // Assuming simplex() time equal to one iteration of EFD()
        pmp->PCI = 0;

     // Calling the simplex method here
        SimplexInitialApproximation( );

//  STEPWISE (0) - stop point for examining results from simplex IA
#ifndef IPMGEMPLUGIN
STEP_POINT( "End Simplex" );
#endif

        // no multi-component phases?
        if( !pmp->FIs )
            return true; // If so, the GEM problem is already solved !

        // Setting default trace amounts to DCs that were zeroed off
        RaiseZeroedOffDCs( 0, pmp->L, sfactor );
        // this operation greatly affects the accuracy of mass balance!

        TotalPhases( pmp->Y, pmp->YF, pmp->YFA );

        for( j=0; j< pmp->L; j++ )
            pmp->X[j] = pmp->Y[j];
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
    }
    else  // Taking previous GEMIPM result as an initial approximation
    {
		TotalPhases( pmp->Y, pmp->YF, pmp->YFA );
        for( j=0; j< pmp->L; j++ )
            pmp->X[j] = pmp->Y[j];
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
        if( pmp->PD==3 /* && pmp->Lads==0 */ )    // added for stability at PIA 06.03.2008 DK
            GammaCalc( LINK_UX_MODE);

        if( pmp->pNP <= -1 )
        {  // With raising species and phases zeroed off by simplex()
           // Setting default trace amounts of DCs that were zeroed off
           RaiseZeroedOffDCs( 0, pmp->L, sfactor );
        }
    }

// STEPWISE (1) - stop point to see IA from old solution or raised simplex
#ifndef IPMGEMPLUGIN
STEP_POINT("Before FIA");
#endif

    return false;
}

// ------------------- ------------------ ----------------
// Calculation of a feasible IPM initial approximation point
//
// Algorithm: see Karpov, Chudnenko, Kulik 1997 Amer.J.Sci. vol 297 p. 798-799
// (Appendix B)
//
// Returns: 0 -  OK,
//          1 -  no SLE colution at the specified precision pa.p.DHB
//          2  - used up more than pa.p.DP iterations
//          3  - too small step length (< pa.p.DG), no descent possible
//          4  - error in Initial mass balance residuals (debugging)
//          5  - error in MetastabilityLagrangeMultiplier() (debugging)
//
long int TMulti::EnterFeasibleDomain()
{
    long int IT1;
    long int I, J, Z,  N, sRet, iRet=0, jK, iB;
    double LM;
    SPP_SETTING *pa = &TProfil::pm->pa;

    ErrorIf( !pmp->MU || !pmp->W, "EnterFeasibleDomain()",
                              "Error of memory allocation for pmp->MU or pmp->W." );

    // Initial rough check of mass balance residuals  - disabled by DK 13.02.2009
//    iB = CheckMassBalanceResiduals( pmp->Y );
//    if( iB >= 0 )
//    {  // Experimental
//      	 return 4;
//    }

    // calculation of total mole amounts of phases
    TotalPhases( pmp->Y, pmp->YF, pmp->YFA );

    if( pmp->PLIM )
        Set_DC_limits(  DC_LIM_INIT );

    // Adjustment of primal approximation according to kinetic constraints
    // Now returns <0 (OK) or index of DC that caused a problem
    //
    jK = MetastabilityLagrangeMultiplier();
    if( jK >= 0 )
    {  // Experimental
        char buf[200];
        sprintf( buf, "Invalid initial Lagrange multiplier for metastability-constrained DC %16s ", pmp->SM[jK] );
		setErrorMessage( 17, "E17IPM: EnterFeasibleDomain(): ", buf);
      	return 5;
    }

//----------------------------------------------------------------------------
// BEGIN:  main loop
    for( IT1=0; IT1 < pa->p.DP; IT1++, pmp->ITF++ )
    {
        // get size of task
        pmp->NR=pmp->N;
        if( pmp->LO )
        {   if( pmp->YF[0] < pmp->DSM && pmp->YFA[0] < pmp->lowPosNum*100.)
                pmp->NR= pmp->N-1;
        }
        N=pmp->NR;
       // Calculation of mass-balance residuals in IPM
       MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C);

      // Testing mass balance residuals
       Z = pmp->N - pmp->E;
       for(I=0;I<Z;I++)
         if( fabs(pmp->C[I]) > pmp->DHBM  ||
             ( pmp->W1 && fabs(pmp->C[I]) > pmp->B[I] * pa->p.GAS ) )
           break;
       if( I == Z ) // OK
         return iRet;       // OK

       WeightMultipliers( true );

       // Assembling and solving the system of linearized equations
       sRet = SolverLinearEquations( N, true );
       if( sRet == 1 )  // error: no SLE solution!
       {
    	 iRet = 1;
 		 setErrorMessage( 5, "E05IPM: EnterFeasibleDomain(): " ,
           "Degeneration in R matrix (fault in the linearized system solver).\n"
           "Invalid initial approximation - cannot start IPM.");
    	 break;
       }

      // SOLVED: solution of linear matrix has been obtained
      pmp->PCI = calcDikin( N, true);  // calc of MU values and Dikin criterion

      LM = calcLM( true ); // Calculation of descent step size LM

      if( LM < min(pa->p.DG, 1e-5) )
      {  // Experimental
          iRet = 3;
    	  setErrorMessage( 3, "E03IPM: EnterFeasibleDomain():",
    			  "Too small LM step size - cannot converge (check Pa_DG?).");
          break;
       }
      if( LM > 1.)
         LM = 1.;
//      cout << "LM " << LM << endl;

      // calculation of new primal solution approximation
      // from step size LM and the descent vector MU
      for(J=0;J<pmp->L;J++)
            pmp->Y[J] += LM * pmp->MU[J];

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
  pVisor->Update( false );
#endif
// STEPWISE (5) Stop point at end of iteration of FIA()
STEP_POINT("FIA Iteration");
#endif
}  /* End loop on IT1 */
//----------------------------------------------------------------------------
    //  Prescribed mass balance precision cannot be reached
    //  Take a look at vector b or values of DHB and DS

    long int ScaleEFD = 1;  // added 11.02.2009 (TW)
   if( IT1 == ScaleEFD * pa->p.DP )
   {  // Experimental
       iRet = 2;
 	   setErrorMessage( 4, "E04IPM: EnterFeasibleDomain(): " ,
        "Maximum allowed number of EnterFeasibleDomain() iterations exceeded!\n ");
    }
   return iRet;   // no solution
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of chemical equilibrium using the Interior Points
//  Method algorithm (see Karpov et al., 1997, p. 785-786)
//  GEM IPM
// Returns: 0, if converged;
//          1, in the case of R matrix degeneration
//          2, (more than max iteration) - no convergence
//              or user's interruption
//          3, GammaCalc() returns bad (non-zero) status (in PIA mode only)
//          4, Mass balance broken  in DualTH (Mol_u)
long int TMulti::InteriorPointsMethod( long int &status, long int rLoop )
{
    long int N, IT1,J,Z,iRet,i;
    double LM=0., LM1=1., FX1;
    SPP_SETTING *pa = &TProfil::pm->pa;

    status = 0;
    if( pmp->FIs )
      for( J=0; J<pmp->Ls; J++ )
            pmp->lnGmo[J] = pmp->lnGam[J];

    for(J=0;J<pmp->N;J++)
    {   pmp->U[J]=0.;
        pmp->C[J]=0.;
    }

    pmp->FX=GX( LM  );  // calculation of G(x)

    if( pmp->FIs ) // multicomponent phases are present
      for(Z=0; Z<pmp->FIs; Z++)
        pmp->YFA[Z]=pmp->XFA[Z];

//----------------------------------------------------------------------------
//  Main loop of IPM iterations
    for( IT1 = 0; IT1 < pa->p.IIM; IT1++, pmp->IT++, pmp->ITG++ )
    {
        pmp->NR=pmp->N;
        if( pmp->LO ) // water-solvent is present
        {
            if( pmp->YF[0]<pmp->DSM && pmp->Y[pmp->LO]< pmp->lowPosNum *100.)
                pmp->NR=pmp->N-1;
        }
        N = pmp->NR;

        PrimalChemicalPotentials( pmp->F, pmp->Y, pmp->YF, pmp->YFA );

        // Setting weight multipliers for DC
        WeightMultipliers( false );

        // Making and solving the R matrix of IPM linearized equations
        iRet = SolverLinearEquations( N, false );
        if( iRet == 1 )
        {
        	setErrorMessage( 7, "E07IPM: IPM-main(): ",
   " Degeneration in R matrix (fault in the linearized system solver).\n"
   " It is not possible to obtain a valid GEMIPM solution.\n"  );
          return 1;
        }

        if( /*false &&*/ (pmp->W1 || rLoop>=0) ) // no first loop
        {
          for(J=0;J<pmp->N;J++)
          {
      		 char buf[200];
          	  if( fabs(pmp->U[J]-pmp->Uc[J] ) > (double)(TProfil::pm->pa.p.PLLG)/*500.*/ ) // broken dual solution
        	  {
        	     if( pmp->Ec != 14 )
        	     { sprintf( buf, "Dual solution (vector u) has changed too much "
        	    		 " in a refinement loop %ld:%ld for IC %3.3s ",
        	    		 pmp->W1, rLoop, pmp->SB[J] );
        		   setErrorMessage( 14, "W14IPM: IPM-main():", buf);
        	     }
    			  else
    			  {
     				 sprintf(buf,"  %-3.3s" ,  pmp->SB[J] );
    				 addErrorMessage(buf);
    			  }
        	  }
          }
          if( pmp->Ec == 14 )
          {  // swapping the dual solution with its copy from previous refinement loop
            double uu;
        	  for( i=0; i<pmp->N; i++ )
        	  {  uu = pmp->U[i]; pmp->U[i] = pmp->Uc[i]; pmp->Uc[i] = uu;}
             pmp->PZ = 2;
             goto CONVERGED;
          }
        }

//SOLVED: got the dual solution u vector - calculating the Dikin criterion
//    of GEM IPM convergence
//       f_alpha( );  commented out 30.10.2007  DK
       pmp->PCI=calcDikin( N, false );

       // Determination of the descent step size LM
	   LM = calcLM( false );

	   LM1=LMD( LM ); // Finding an optimal value of the descent step
       FX1=GX( LM1 ); // New G(X) value after the descent step
       pmp->PCI = sqrt(pmp->PCI); // Dikin criterion

       if( (IT1 > 3) && (FX1 - pmp->FX > 10)  &&
          ( IT1 > pa->p.IIM/2 || pmp->PD==3 ) )
        {                                 // Adjusting Dikin threshold
            if( pmp->LO && pmp->X[pmp->LO]>pmp->lowPosNum && pmp->PD==3 )
                pmp->PD=2;     //  what does this mean? - Do not GammaCalc in all loops
            else  pmp->DX= 0.5 * pmp->PCI;
        }

       MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );

       pmp->FX=FX1;
       // Main IPM iteration done
       // Calculation of activity coefficients
        if( pmp->PD==3 )
            status = GammaCalc( LINK_UX_MODE );

        if( pmp->PHC[0] == PH_AQUEL && pmp->XF[0] <= pa->p.XwMin &&
             pmp->X[pmp->LO] <= pmp->lowPosNum*1e3 )    // bugfix 29.11.05 KD
        {
            pmp->XF[0]=0.;  // elimination of aqueous phase
            pmp->XFA[0] = 0.;
        }
        // Restoring vectors Y and YF
        Restoring_Y_YF();

if( pmp->pNP && rLoop < 0 && status )
{
	setErrorMessage( 18, "E18IPM: IPM-main():", "Bad GammaCalc() status in SIA mode");
	return 3L;
}

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
  pVisor->Update( false );
#endif
// STEPWISE (6)  Stop point at IPM() main iteration
STEP_POINT( "IPM Iteration" );
#endif
        if( pmp->PCI < pmp->DX || pmp->PZ == 2 )  // Dikin criterion satisfied - converged!
            goto CONVERGED;
    } // end of main IPM cycle

    setErrorMessage( 6, "E06IPM: IPM-main(): " ,
            "IPM convergence criterion threshold (Pa_DK) could not be reached"
    		" (more than Pa_IIM iterations done);\n" );
    return 2L; // bad convergence - too many IPM iterations!
//----------------------------------------------------------------------------
CONVERGED:
 // cout << "LM " << LM << " LM1 "	 << LM1 << " PCI " << pmp->PCI << endl;
 // Final calculation of phase amounts and activity coefficients
  TotalPhases( pmp->X, pmp->XF, pmp->XFA );

  if(!pmp->PZ && pmp->W1)
  { iRet =  Mol_u( pmp->Y, pmp->X, pmp->XF, pmp->XFA );
    if(iRet )
 	   return 4L;
  }

  if( pmp->PD == 1 || pmp->PD == 2  || pmp->PD == 3  )
        GammaCalc( LINK_UX_MODE );
//   else
  ConCalc( pmp->X, pmp->XF, pmp->XFA );

   MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->X, pmp->B, pmp->C);
   return 0L;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Calculation of mass-balance residuals in GEM IPM CSD
//  Params: N - number of IC in IPM problem
//          L -   number of DC in IPM problem
//          A - DC stoichiometry matrix (LxN)
//          Y - moles  DC quantities in IPM solution (L)
//          B - Input bulk chem. compos. (N)
//          C - mass balance residuals (N)
void
TMulti::MassBalanceResiduals( long int N, long int L, double *A, double *Y, double *B,
         double *C )
{
    long int ii, jj, i;
    for(ii=0;ii<N;ii++)
        C[ii]=B[ii];
    for(jj=0;jj<L;jj++)
     for( i=arrL[jj]; i<arrL[jj+1]; i++ )
     {  ii = arrAN[i];
         C[ii]-=(*(A+jj*N+ii))*Y[jj];
     }
}

// Diagnostics for a severe break of mass balance in the initial approximation
// or after GEM IPM PhaseSelect() (when pmp->X is passed as parameter)
// Returns -1 (Ok) or index of the chemical element for which the balance is
// broken
long int
TMulti::CheckMassBalanceResiduals(double *Y )
{
	double cutoff;
    long int iRet = -1L;
    char buf[300];

	// cutoff = pmp->DHBM * 1e4;
	cutoff = min (pmp->DHBM*1.0e5, 1.0e-3 );	// changed, 28.08.2008 (TW,DK)
	if( cutoff > 1e-3 )
		cutoff = 1e-3;
	MassBalanceResiduals( pmp->N, pmp->L, pmp->A, Y, pmp->B, pmp->C);
    for(long int i=0; i<pmp->N; i++)
	{
	   if( pmp->B[i] > cutoff )
	   {	   // Major IC
		  if( fabs( pmp->C[i] ) < pmp->B[i] )
			  continue;
	   }
	   else {  // Trace IC - less than (cutoff) mol
	      if( fabs( pmp->C[i] ) < cutoff )
	    	  continue;
	   }
	   if( iRet < 0  )
	   {
          iRet = i;  // Error state is activated
		  sprintf(buf, "Mass balance broken on iteration %ld  for independent components %-3.3s",
					     pmp->ITG, pmp->SB[i] );
		  setErrorMessage( 2, "E02IPM: IPM-main():" ,buf);
		}
		else
		{
		 sprintf(buf," %-3.3s" ,  pmp->SB[i] );
		 addErrorMessage( buf );
		 }
	}
	return iRet;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Interior Points Method:
// subroutine for unconditional minimization of the descent step length
// on the interval 0 to LM
// uses the "Golden Section" algorithm
//
// Returns: optimal value of LM which provides the largest monotonous
// decrease in G(X)
//
double TMulti::LMD( double LM )
{
    double A,B,C,LM1,LM2,FX1,FX2;
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
    if( FX1>FX2)
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

// Cleaning the unstable phase with index k >= 0

void TMulti::ZeroDCsOff( long int jStart, long int jEnd, long int k )
{
  if( k >=0 )
     pmp->YF[k] = 0.;

  for(long int j=jStart; j<jEnd; j++ )
     pmp->Y[j] =  0.0;
}

// Inserting minor quantities of DC which were zeroed off after simplex()
// (important for the automatic initial approximation with solution phases
//  (k = -1)  or after Selekt2() algorithm (k >= 0)
//
void TMulti::RaiseZeroedOffDCs( long int jStart, long int jEnd, double scalingFactor, long int k )
{
  double sfactor = scalingFactor;
  SPP_SETTING *pa = &TProfil::pm->pa;

  if( fabs( sfactor ) > 1. )   // can reach 30 at total moles in system above 300000 (DK 11.03.2008)
	  sfactor = 1.;       // Workaround for very large systems (insertion breaks the EFD convergence)
  if( k >= 0 )
       pmp->YF[k] = 0.;

  for(long int j=jStart; j<jEnd; j++ )
  {
     switch( pmp->DCC[j] )
     {
       case DC_AQ_PROTON:
       case DC_AQ_ELECTRON:
       case DC_AQ_SPECIES:
          if( k >= 0 || pmp->Y[j] < pa->p.DFYaq * sfactor )
               pmp->Y[j] =  pa->p.DFYaq * sfactor;
           break;
       case DC_AQ_SOLVCOM:
       case DC_AQ_SOLVENT:
            if( k >= 0 || pmp->Y[j] < pa->p.DFYw * sfactor )
                pmp->Y[j] =  pa->p.DFYw * sfactor;
            break;
       case DC_GAS_H2O:
       case DC_GAS_CO2:
       case DC_GAS_H2:
       case DC_GAS_N2:
       case DC_GAS_COMP:
       case DC_SOL_IDEAL:
            if( k >= 0 || pmp->Y[j] < pa->p.DFYid*sfactor )
                  pmp->Y[j] = pa->p.DFYid * sfactor;
             break;
       case DC_SOL_MINOR:
            if( k >= 0 || pmp->Y[j] < pa->p.DFYh*sfactor )
                   pmp->Y[j] = pa->p.DFYh * sfactor;
             break;
       case DC_SOL_MAJOR:
            if( k >= 0 || pmp->Y[j] < pa->p.DFYr * sfactor )
                  pmp->Y[j] =  pa->p.DFYr * sfactor;
             break;
       case DC_SCP_CONDEN:
             if( k >= 0 )
             {                // Added 05.11.2007 DK
                 pmp->Y[j] =  pa->p.DFYs * sfactor;
                 break;
             }
             if( pmp->Y[j] < pa->p.DFYc * sfactor )
                  pmp->Y[j] =  pa->p.DFYc * sfactor;
              break;
                    // implementation for adsorption?
       default:
             if( k>=0 || pmp->Y[j] < pa->p.DFYaq *sfactor )
                   pmp->Y[j] =  pa->p.DFYaq * sfactor;
             break;
     }
     if( k >=0 )
     pmp->YF[k] += pmp->Y[j];
   } // i
}

// Adjustment of primal approximation according to kinetic constraints
long int TMulti::MetastabilityLagrangeMultiplier()
{
    double E = TProfil::pm->pa.p.DKIN; //1E-8;  Default min value of Lagrange multiplier p
//    E = 1E-30;

    for(long int J=0;J<pmp->L;J++)
    {
        if( pmp->Y[J] < 0. )   // negative number of moles!
        	return J;
    	if( pmp->Y[J] < pmp->lowPosNum )
            continue;

        switch( pmp->RLC[J] )
        {
        case NO_LIM:
        case LOWER_LIM:
            if( pmp->Y[J]<=pmp->DLL[J])
                pmp->Y[J]=pmp->DLL[J]+E;
            break;
        case BOTH_LIM:
            if( pmp->Y[J]<=pmp->DLL[J])
                pmp->Y[J]=pmp->DLL[J]+E;
            if( pmp->Y[J]>=pmp->DUL[J])     // SD 22/01/2009
            {
                if( pmp->DUL[J] == 1e6 )
                   return J;   // Broken initial approximation!
                pmp->Y[J]=pmp->DUL[J]-E;
                if( pmp->Y[J]<=pmp->DLL[J])
                  	pmp->Y[J]=(pmp->DUL[J]+pmp->DLL[J])/2.;
             }
             break;
        case UPPER_LIM:
            if( pmp->Y[J]>=pmp->DUL[J])
            {
               	if( pmp->DUL[J] == 1e6 )
               	    return J;   // Broken initial approximation!
            	pmp->Y[J]=pmp->DUL[J]-E;
                if( pmp->Y[J]<=0)         // SD 22/01/2009
                	pmp->Y[J]=(pmp->DUL[J])/2.;
            }
            break;
        }
    }   // J
    return -1L;
}

// Calculation of weight multipliers for DCs
void TMulti::WeightMultipliers( bool square )
{
  long int J;
  double  W1, W2;

  for( J=0; J<pmp->L; J++)
  {
    switch( pmp->RLC[J] )
    {
      case NO_LIM:
      case LOWER_LIM:
           W1=(pmp->Y[J]-pmp->DLL[J]);
           if( square )
             pmp->W[J]= W1 * W1;
             else
             pmp->W[J] = max( W1, 0. );
           break;
      case UPPER_LIM:
           W1=(pmp->DUL[J]-pmp->Y[J]);
           if( square )
             pmp->W[J]= W1 * W1;
           else
             pmp->W[J] = max( W1, 0.);
           break;
      case BOTH_LIM:
           W1=(pmp->Y[J]-pmp->DLL[J]);
           W2=(pmp->DUL[J]-pmp->Y[J]);
           if( square )
           {
             W1 = W1*W1;
             W2 = W2*W2;
           }
           pmp->W[J]=( W1 < W2 ) ? W1 : W2 ;
           if( !square && pmp->W[J] < 0. ) pmp->W[J]=0.;
           break;
      default: // error
    	  setErrorMessage( 16, "E16IPM: IPM-main():", "Error in codes of some DC metastability constraints" );
          Error( pmp->errorCode, pmp->errorBuf );
    }
  } // J
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Make and Solve a system of linear equations to find the dual vector
// approximation using a method of Cholesky Decomposition (good if a
// square matrix R happens to be symmetric and positive defined).
// If Cholesky Decomposition does not solve the problem, an attempt is done
// to solve the SLE using method of LU Decomposition
// (A = L*U , L is lower triangular ( has elements only on the diagonal and below )
//   U is is upper triangular ( has elements only on the diagonal and above))
// Parameters:
// bool initAppr - Inital approximation point(true) or iteration of IPM (false)
// int N - dimension of the matrix R (number of equations)
// pmp->U - dual solution vector (N).
// Return values: 0  - solved OK;
//                1  - no solution, degenerated or inconsistent system
//
#define  a(j,i) ((*(pmp->A+(i)+(j)*Na)))
//#define  a(j,i) ((*(pmp->A+(i)+(j)*N)))  obsolete
long int TMulti::SolverLinearEquations( long int N, bool initAppr )
{
  long int ii,i, jj, kk, k, Na = pmp->N;
//  double aa;
  Alloc_A_B( N );

  // Making the  matrix of IPM linear equations
  for( kk=0; kk<N; kk++)
   for( ii=0; ii<N; ii++ )
      (*(AA+(ii)+(kk)*N)) = 0.;

  for( jj=0; jj<pmp->L; jj++ )
   if( pmp->Y[jj] > pmp->lowPosNum   )
   {
      for( k=arrL[jj]; k<arrL[jj+1]; k++)
        for( i=arrL[jj]; i<arrL[jj+1]; i++ )
        { ii = arrAN[i];
          kk = arrAN[k];
          if( ii>= N || kk>= N )
           continue;
          (*(AA+(ii)+(kk)*N)) += a(jj,ii) * a(jj,kk) * pmp->W[jj];
        }
    }

   if( initAppr )
     for( ii=0; ii<N; ii++ )
         BB[ii] = pmp->C[ii];
   else
      {
         for( ii=0; ii<N; ii++ )
            BB[ii] = 0.;
          for( jj=0; jj<pmp->L; jj++)
           if( pmp->Y[jj] > pmp->lowPosNum  )
              for( i=arrL[jj]; i<arrL[jj+1]; i++ )
              {  ii = arrAN[i];
                 if( ii>= N )
                  continue;
                 BB[ii] += pmp->F[jj] * a(jj,ii) * pmp->W[jj];
              }
      }

#ifndef PGf90

  Array2D<double> A(N,N, AA);
  Array1D<double> B(N, BB);

#else

  Array2D<double> A(N,N);
  Array1D<double> B(N);

  for( kk=0; kk<N; kk++)
   for( ii=0; ii<N; ii++ )
      A[kk][ii] = (*(AA+(ii)+(kk)*N));

   for( ii=0; ii<N; ii++ )
     B[ii] = BB[ii];

#endif
// this routine constructs its Cholesky decomposition, A = L x LT .
  Cholesky<double>  chol(A);

  if( chol.is_spd() )  // is positive definite A.
  {
    B = chol.solve( B );
  }
  else
  {
// no solution by Cholesky decomposition; Trying the LU Decompositon
// The LU decompostion with pivoting always exists, even if the matrix is
// singular, so the constructor will never fail.

   LU<double>  lu(A);

// The primary use of the LU decomposition is in the solution
// of square systems of simultaneous linear equations.
// This will fail if isNonsingular() returns false.
   if( !lu.isNonsingular() )
     return 1; // Singular matrix - it's a pity.

  B = lu.solve( B );
  }

  for( ii=0; ii<N; ii++ )
   pmp->U[ii] = B[(int)ii];
  return 0;
}
#undef a

// Calculation of MU values (dual DC chemical potentials) and Dikin criterion
// Parameters:
// bool initAppr - Inital approximation (true) or main iteration of IPM (false)
// int N - dimension of the matrix R (number of equations)
double TMulti::calcDikin(  long int N, bool initAppr )
{
  long int  J;
  double Mu, PCI=0.;

  for(J=0;J<pmp->L;J++)
  {
    if( pmp->Y[J] > pmp->lowPosNum )
    {
      Mu = DualChemPot( pmp->U, pmp->A+J*pmp->N, N, J );
      if( !initAppr )
        Mu -= pmp->F[J];
      PCI += Mu*pmp->W[J]*Mu;
      pmp->MU[J] = Mu*pmp->W[J];
    }
    else
      pmp->MU[J]=0.; // initializing dual potentials
  }
  if( initAppr )
  {
     if( PCI > pmp->lowPosNum  )
          PCI=1/sqrt(PCI);
     else PCI=1.; // zero Psi value ?
  }

  return PCI;
}

// Calculation of the descent step length LM
// Parameters:
// bool initAppr - Inital approximation (true) or iteration of IPM (false)
double TMulti::calcLM(  bool initAppr )
{
   long int J, Z = -1;
   double LM=1., LM1=1., Mu;

   for(J=0;J<pmp->L;J++)
   {
     Mu = pmp->MU[J];
    if( pmp->RLC[J]==NO_LIM ||
        pmp->RLC[J]==LOWER_LIM || pmp->RLC[J]==BOTH_LIM )
    {
       if( Mu < 0 && fabs(Mu)>pmp->lowPosNum*100.)
       {
         if( Z == -1 )
         { Z = J;
           LM = (-1)*(pmp->Y[Z]-pmp->DLL[Z])/Mu;
         }
         else
         {
           LM1 = (-1)*(pmp->Y[J]-pmp->DLL[J])/Mu;
           if( LM > LM1)
             LM = LM1;
         }
       }
    }
    if( pmp->RLC[J]==UPPER_LIM || pmp->RLC[J]==BOTH_LIM )
    {
       if( Mu > pmp->lowPosNum*100.)
       {
         if( Z == -1 )
         { Z = J;
           LM = (pmp->DUL[Z]-pmp->Y[Z])/Mu;
         }
         else
         {
           LM1=(pmp->DUL[J]-pmp->Y[J])/Mu;
           if( LM>LM1)
               LM=LM1;
         }
      }
    }
  }

  if( initAppr )
  { if( Z == -1 )
     LM = pmp->PCI;
    else
     LM *= .95;     // Smoothing of final lambda value
  }
  else
  {  if( Z == -1 )
       LM = 1./sqrt(pmp->PCI);  // Might cause infinite loop in LMD() if PCI is too low?
//     LM = min( LM, 10./pmp->DX );
       LM = min( LM, 1.0e10 );  // Set an empirical upper limit for LM to prevent freezing
  }
  return LM;
}

// Restoring primal vectors Y and YF
void TMulti::Restoring_Y_YF()
{
 long int Z, I, JJ = 0;

 for( Z=0; Z<pmp->FI ; Z++ )
 {
   if( pmp->XF[Z] <= pmp->DSM ||
       ( pmp->PHC[Z] == PH_SORPTION &&
       ( pmp->XFA[Z] < TProfil::pm->pa.p.ScMin) ) )
   {
      pmp->YF[Z]= 0.;
      if( pmp->FIs && Z<pmp->FIs )
         pmp->YFA[Z] = 0.;
      for(I=JJ; I<JJ+pmp->L1[Z]; I++)
      {
        pmp->Y[I]=0.;
        pmp->lnGam[I] = 0.;
      }
   }
   else
   {
     pmp->YF[Z] = pmp->XF[Z];
     if( pmp->FIs && Z < pmp->FIs )
        pmp->YFA[Z] = pmp->XFA[Z];
     for(I = JJ; I < JJ+pmp->L1[Z]; I++)
        pmp->Y[I]=pmp->X[I];
   }
   JJ += pmp->L1[Z];
 } // Z

}

// Calculation of the scaling factor for the IPM-2 algorithm
double TMulti::calcSfactor()
{
    double molB=0.0, sfactor;
    long int i, NN = pmp->N - pmp->E;

    for(i=0;i<NN;i++)
      molB += pmp->B[i];

   sfactor = pow( molB, 0.4 )/7.7;
   return sfactor;
}

//===================================================================
// Checking Karpov phase stability criteria Fa for phases and DCs
//  using Selekt2() algorithm by Karpov & Chudnenko (1989)
//  modified by DK in 1995 and in 2007
//  Returns 0, if some phases were inserted and a new IPM loop is needed
//             (up to 3 loops possible);
//          1, if the IPM solution is final and consistent, no phases were inserted
//          -1, if the IPM solution is inconsistent after 3 Selekt2() loops
//  In this case, the index of most problematic phase is passed through kf or
//  ku parameter (parameter value -1 means that no problematic phases were found)
//
long int TMulti::PhaseSelect( long int &kfr, long int &kur, long int rLoop )
{
    long int k, j, jb;
    long int kf, ku;
    double F1, F2, *F0, sfactor;
    SPP_SETTING *pa = &TProfil::pm->pa;

    sfactor = calcSfactor();
    f_alpha( );  // calculation of Karpov phase stability criteria (in pmp->Falp)
    F0 = pmp->Falp;

    (pmp->K2)++;
    kf = -1; ku = -1;  // Index for phase diagnostics
    F1 = pa->p.DF;  // Fixed 29.10.2007  DK
    F2 = -pa->p.DFM;  // Meaning of DFM changed 02.11.2007

    for(k=0;k<pmp->FI;k++)
    {
        if( F0[k] > F1 && pmp->YF[k] < pa->p.DS )  // < pmp->lowPosNum?
        {            // stable phase not in mass balance - to be inserted
            F1=F0[k];
            kf=k;
        }
        if( F0[k] < F2 && pmp->YF[k] >= pa->p.DS )  // Fixed 2.11.2007
        {            // unstable phase in mass balance - to be excluded
            F2=F0[k];
            ku=k;
        }
    }
kfr = kf;
kur = ku;
    if( kfr < 0 && kur < 0 )
    {    // No phases to insert/exclude or no Fa distortions found
          // Successful end of iterations of SELEKT2()
        return 1L;
    }

    if( (F2 < -pa->p.DFM ) && ( ku >= 0 ) )
    {
        if( pmp->K2 > 4 ) // Three Selekt2() loops have already been done!
        	return -1L;   // Persistent presence of unstable phase(s) - bad system!

        // Excluding problematic phases
        do
        {  // excluding all phases with  F2 < DF*sfactor
            for( jb=0, k=0; k < ku; k++ )
                 jb += pmp->L1[k];

            ZeroDCsOff( jb, jb+pmp->L1[ku], ku ); // Zeroing the phase off
            pmp->FI1--;
            // find a new phase to exclude, if any exists
            F2= -pa->p.DFM;
            ku = -1;
            for( k=0; k<pmp->FI; k++ )
                if( F0[k] < F2 && pmp->YF[k] >= pa->p.DS )
                {
                    F2=F0[k];
                    ku=k;
                }
        }
        while( ( F2 <= -pa->p.DFM ) && ( ku >= 0 ) );
    } //if ku

    // Inserting problematic phases
    if( F1 > pa->p.DF && kf >= 0 )
    {
        if( pmp->K2 > 4 )
           return -1L;   // Persistent absence of stable phase(s) - bad system!

        // There is a phase for which DF*sfactor threshold is exceeded
        do
        {   // insert this phase and set Y[j] for its components
            // with account for asymmetry and non-ideality
            for( jb=0, k=0; k < kf; k++ )
                 jb += pmp->L1[k];

            RaiseZeroedOffDCs( jb, jb+pmp->L1[kf], sfactor, kf );

            pmp->FI1++;  // check phase rule

            if( pmp->FI1 >= pmp->NR+1 )
               break;   // No more phases can be inserted

            // find a new phase to insert, if any exists
            F1= pmp->lowPosNum*10000.; // 1e-16
            kf = -1;
            for( k=0; k<pmp->FI; k++ )
                if( F0[k] > F1 && pmp->YF[k] < pmp->lowPosNum )
                {
                    F1=F0[k];
                    kf=k;
                }
        }
        while( F1 > pa->p.DF && kf >= 0 );
        // end of insertion cycle
    } // if kf changed SD 03/02/2009
// Raise zeros in DC amounts in phases-solutions in the case if some phases
// were inserted or excluded                          - experimental option!
        double RaiseZeroVal = pmp->DHBM;  // Added 29.10.07  by DK
        jb=0;
        for(k=0;k<pmp->FIs;k++)
        {
            if( ( pmp->YF[k] >= pa->p.DS )  || (pmp->pNP && rLoop < 0 ) ) // Only in phase present in mass balance!
            {                            // (acc. to definition of L_S set) PIA only if initial!
                 pmp->YF[k]=0.;
                 for(j=jb;j<jb+pmp->L1[k];j++)
                 {
                    if( pmp->Y[j] < pmp->lowPosNum )
                        pmp->Y[j] = RaiseZeroVal; // bugfix 29.10.07
                    pmp->YF[k] += pmp->Y[j]; // calculate new amounts of phases
                 }
            }
            jb+=pmp->L1[k];
        }
    // } SD 03/02/2009
#ifndef IPMGEMPLUGIN
    STEP_POINT("Select2()");
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelection"
#endif
#endif
    if( pmp->K2 > 1 )
    { // more then the first step - but the IPM solution has not been improved
       double RaiseZeroVal = pmp->DHBM*0.1;   // experimental
       for(j=0;j<pmp->L;j++)
          if( fabs(pmp->Y[j]- pmp->XY[j]) > RaiseZeroVal ) //
               goto S6;
       // pmp->PZ=2; // No significant change has been done by Selekt2()
       return 1L;
    }
S6: // copy of X vector has been changed by Selekt2() algorithm - store
    for(j=0;j<pmp->L;j++)
        pmp->XY[j]=pmp->Y[j];

    return 0L;  // Another loop is needed
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Internal memory allocation for IPM performance optimization
// (since version 2.2.0)
//
void TMulti::Alloc_A_B( long int newN )
{
  if( AA && BB && (newN==sizeN) )
    return;
  Free_A_B();
  AA = new  double[newN*newN];
  BB = new  double[newN];
  sizeN = newN;
}

void TMulti::Free_A_B()
{
  if( AA  )
    { delete[] AA; AA = 0; }
  if( BB )
    { delete[] BB; BB = 0; }
  sizeN = 0;
}


// Building an index list of non-zero elements of the matrix pmp->A
#define  a(j,i) ((*(pmp->A+(i)+(j)*pmp->N)))
void TMulti::Build_compressed_xAN()
{
 long int ii, jj, k;

 //if( arrL && arrAN && (sizeL == pmp->L+1) && ( sizeAN == pmp->N ) )
 //  return; // The index arrays are intact - no need to remake (added by DK 27.05.08)  SD 26/11/2008

 // Calculate number of non-zero elements in A matrix
 k = 0;
 for( jj=0; jj<pmp->L; jj++ )
   for( ii=0; ii<pmp->N; ii++ )
     if( fabs( a(jj,ii) ) > 1e-12 )
       k++;

   // Free old memory allocation
    Free_compressed_xAN();

   // Allocate memory
   arrL = new long int[pmp->L+1]; // sizeL = pmp->L+1;
   arrAN = new long int[k];       // sizeAN = pmp->N; // sizeAN = k;

   // Set indexes in the index arrays
   k = 0;
   for( jj=0; jj<pmp->L; jj++ )
   { arrL[jj] = k;
     for( ii=0; ii<pmp->N; ii++ )
       if( fabs( a(jj,ii) ) > 1e-12 )
       {
        arrAN[k] = ii;
        k++;
       }
   }
   arrL[jj] = k;
}
#undef a

void TMulti::Free_compressed_xAN()
{
  if( arrL  )
    { delete[] arrL; arrL = 0;  }
  if( arrAN )
    { delete[] arrAN; arrAN = 0;  }
}

void TMulti::Free_internal()
{
  Free_compressed_xAN();
  Free_A_B();
}

void TMulti::Alloc_internal()
{
// optimization 08/02/2007
 Alloc_A_B( pmp->N );
 Build_compressed_xAN();
}

// added SD 03/02/2009
void TMulti::setErrorMessage( long int num, const char *code, const char * msg)
{
  pmp->Ec  = num;
  strncpy( pmp->errorCode, code, 99 );
  pmp->errorCode[99] ='\0';
  strncpy( pmp->errorBuf,  msg,  499 );
  pmp->errorBuf[499] ='\0';
}

void TMulti::addErrorMessage( const char * msg)
{
  long int ln = strlen(pmp->errorBuf);
  strncpy(pmp->errorBuf+ln, msg, 499-ln  );
  pmp->errorBuf[499] ='\0';
}

//--------------------- End of ipm_main.cpp ---------------------------
