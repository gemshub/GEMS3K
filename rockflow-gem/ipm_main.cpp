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
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
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
//  Returns: true   if IPM result is OK
//           false  if a good result could not be obtained
//
void TMulti::MultiCalcMain()
{
    int i, j, RepeatSel=0, eRet;
    SPP_SETTING *pa = &TProfil::pm->pa;

    pmp->W1=0;
    if( pmp->pULR && pmp->PLIM )
        Set_DC_limits( DC_LIM_INIT );

    // testing the entry into feasible domain
mEFD:
     if(pmp->PZ && pmp->W1)
     { for( i=0; i<pmp->L; i++ )
        pmp->Y[i]=pmp->X[i];
      TotalPhases( pmp->Y, pmp->YF, pmp->YFA );
     }

    eRet = EnterFeasibleDomain( );

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
     case 2:  // max number of iterations has been exceeded in EnterFeasibleDomain()
                 Error("E04IPM EFD(): " ,
     "Prescribed precision of mass balance could not be reached because the vector b or\n"
     "DC stoichiometries or standard-state thermodynamic data are inconsistent.\n");
              break;
     case 1: // degeneration in R matrix  for EnterFeasibleDomain()
           if( pmp->DHBM<1e-5 )
            {
               pmp->DHBM *= 10.;
               goto mEFD;
            }
           else
                 Error("E05IPM EFD(): " ,
          "Degeneration in R matrix (fault in the linearized system solver).\n"
          "Invalid initial approximation - further IPM calculations are not possible");
           break;
    }

   // call of main IPM minimization algorithm
   eRet = InteriorPointsMethod( );

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
               if( pmp->DX<1e-3 )
               {
                   pmp->DX *= 10.;
                   goto mEFD;
                }
               else
                 Error("E06IPM IPM-main(): " ,
     "Given IPM convergence criterion could not be reached;\n Perhaps, vector b is not balanced,\n"
     "or DC stoichiometries or standard-state thermodynamic data are inconsistent. \n");
     case 1: // degeneration in R matrix  for InteriorPointsMethod()
           if( pmp->DHBM<1e-5 )
            {
               pmp->DHBM *= 10.;
               goto mEFD;
            }
           else
               Error("E07IPM IPM-main(): ",
        "Degeneration in R matrix (fault in the linearized system solver).\n"
        "No valid IPM solution could be obtained. Probably, vector b is not balanced,\n"
        "or DC stoichiometries or standard-state thermodynamic data are inconsistent,\n"
        "or some relevant phases or DC are missing, or some kinetic constraints are too stiff.\n"
        );
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

    if( !pa->p.PC )    //  No PhaseSelect() operation allowed
    {  if( pmp->PD >= 2 )
           for( i=0; i<pmp->L; i++)
             pmp->G[i] = pmp->G0[i];
        return;  // solved
    }

//========= call Selekt2 algorithm =======
   PhaseSelect();
   if( pa->p.PC == 2 )
       XmaxSAT_IPM2();  // Install upper limits to xj of surface species

   if( !pmp->MK )
     if( RepeatSel<3 )
       { RepeatSel++;
         goto mEFD;
       }
     else
       Error( "E08IPM PhaseSelect(): "," Insertion of phases was incomplete!");

  MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->X, pmp->B, pmp->C);

#ifndef IPMGEMPLUGIN
// STEPWISE (4) Stop point after PhaseSelect()
STEP_POINT("PhaseSelect");
#ifndef Use_mt_mode
pVisor->Update( false );
#endif
#endif

   if(pmp->PZ )
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
               else
               {
                 gstring  buf,buf1;
                 vstr pl(5);
                 int jj=0;
                 for( j=0; j<pmp->N-pmp->E; j++ )
//  if( fabs(pmp->C[j]) > pmp->DHBM  || fabs(pmp->C[j]) > pmp->B[j] * pa->p.GAS )
                  if( fabs(pmp->C[j]) > pmp->B[j] * pa->p.GAS )
                  { sprintf( pl, " %-2.2s  ", pmp->SB[j] );
                    buf1 +=pl;
                    jj=1;
                  }
                 if( jj )
                 {
                    buf = "Prescribed balance precision cannot be reached\n";
                    buf += "for some trace independent components: ";
                    buf += buf1;
                    Error("E09IPM IPM-main(): ", buf.c_str() );
                 }
                else
                 Error("E10IPM IPM-main(): " ,
                 "Inconsistent GEM solution: imprecise mass balance\n for some major independent components: " );
              }
            }
          } // end of i loop
      }
   }
   for( i=0; i<pmp->L; i++)
      pmp->G[i] = pmp->G0[i];
   // Normal return after successfull improvement of mass balance precision
}

//Call for IPM iteration sequence
void TMulti::MultiCalcIterations()
{

    MultiCalcMain();

    // calculation of demo data for gases
    for( int ii=0; ii<pmp->N; ii++ )
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
bool TMulti::AutoInitialApprox( )
{
    int i, j, k, NN;
    double minB, sfactor;
    SPP_SETTING *pa = &TProfil::pm->pa;

    NN = pmp->N - pmp->E;
    minB=pmp->B[0];
    for(i=0;i<NN;i++)
      if( pmp->B[i] < minB )
          minB = pmp->B[i];
    if( minB < pa->p.DB )  // foolproof
       minB = pa->p.DB;

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
    // Analyzing if Simplex approximation is needed
    if( !pmp->pNP  )
    {   // Preparing to call Simplex method
        pmp->FitVar[4] = pa->p.AG;  // smoothing parameter initialization
        pmp->pRR1 = 0;
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
//      pmp->IC = 0.0;  For reproducibility of simplex FIA?
        if( pa->p.PSM && pmp->FIs )
            GammaCalc(LINK_FIA_MODE);
        if( pa->p.PC == 2 )
           XmaxSAT_IPM2_reset();  // Reset upper limits for surface species
        pmp->IT = 0;
        pmp->pNP = 0;
        pmp->K2 = 0;
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
    }
    else  // Taking previous GEMIPM result as initial approximation
    {
        int jb, je=0, jpb, jpe=0, jdb, jde=0, ipb, ipe=0;
        double LnGam, FitVar3;

        FitVar3 = pmp->FitVar[3];
        pmp->FitVar[3] = 1.0;
        TotalPhases( pmp->Y, pmp->YF, pmp->YFA );
        for( j=0; j< pmp->L; j++ )
            pmp->X[j] = pmp->Y[j];
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
        ConCalc( pmp->X, pmp->XF, pmp->XFA );
        if( pmp->E && pmp->LO )
            IS_EtaCalc();

        for( k=0; k<pmp->FI; k++ )
        {
            jb = je;
            je += pmp->L1[k];
// Indexes for extracting data from IPx, PMc and DMc arrays
    ipb = ipe;                  // added 07.12.2006 by KD
    ipe += pmp->LsMod[k*3]*pmp->LsMod[k*3+1];
            jpb = jpe;
            jpe += pmp->LsMod[k*3]*pmp->LsMod[k*3+2];  // Changed 07.12.2006  by KD
            jdb = jde;
            jde += pmp->LsMdc[k]*pmp->L1[k];

            if( pmp->PHC[k] == PH_SORPTION || pmp->PHC[k] == PH_POLYEL )
            {
               if( pmp->E && pmp->LO )
                   GouyChapman( jb, je, k );
//               SurfaceActivityTerm( jb, je, k );
                 SurfaceActivityCoeff(  jb, je, jpb, jdb, k );
            }
            for( j=jb; j<je; j++ )
            {
                LnGam = pmp->lnGam[j];
                if( fabs( LnGam ) > 1e-9 )
                    pmp->Gamma[j] = exp( LnGam );
                else pmp->Gamma[j] = 1.0;
                pmp->F0[j] = Ej_init_calc( 0.0, j, k  );
                pmp->G[j] = pmp->G0[j] + pmp->GEX[j] + pmp->F0[j];
                pmp->lnGmo[j] = LnGam;
            }  // j
        } // k
        pmp->FitVar[3] = FitVar3; // Restoring smoothing parameter

        if( pmp->pNP <= -1 )
        {  // With raising species and phases zeroed off by simplex()
           // Setting default trace amounts of DCs that were zeroed off
           RaiseZeroedOffDCs( 0, pmp->L, 1/1000. );
        }
    }  // else
    pmp->MK = 1;
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
// Returns: 0 - OK,
//          1 - no convergence at the specified precision pa.p.DHB
//          2  - used up more than pa.p.DP iterations
//
int TMulti::EnterFeasibleDomain()
{
    short IT1;
    int I, J, Z,  N, sRet, iRet=0;
    double DHB, LM;
    SPP_SETTING *pa = &TProfil::pm->pa;

    ErrorIf( !pmp->MU || !pmp->W, "EnterFeasibleDomain()",
                              "Error alloc pmp->MU or pmp->W." );
    DHB = pmp->DHBM;  // Convergence (balance precision) criterion
    pmp->Ec=0;  // Return code

    // calculation of total mole amounts of phases
    TotalPhases( pmp->Y, pmp->YF, pmp->YFA );

    if( pmp->PLIM )
        Set_DC_limits(  DC_LIM_INIT );

    // Adjustment of primal approximation according to kinetic constraints
    LagrangeMultiplier();

//----------------------------------------------------------------------------
// BEGIN:  main loop
    for( IT1=0; IT1 < pa->p.DP; IT1++, pmp->IT++ )
    {
        // get size of task
        pmp->NR=pmp->N;
        if( pmp->LO )
        {   if( pmp->YF[0] < pmp->DSM && pmp->YFA[0] < pmp->lowPosNum*100.)
                pmp->NR= (short)(pmp->N-1);
        }
        N=pmp->NR;

       // Calculation of mass-balance residuals in IPM
       MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C);

      // Testing mass balance residuals
       Z = pmp->N - pmp->E;
       for(I=0;I<Z;I++)
         if( fabs(pmp->C[I]) > DHB  ||
             ( pmp->W1 && fabs(pmp->C[I]) > pmp->B[I] * pa->p.GAS ) )
           break;
       if( I == Z ) // OK
       { pmp->IT -= IT1;
         return iRet;       // OK
       }

       WeightMultipliers( true );
       // Assembling and solving the system of linearized equations
       sRet = SolverLinearEquations( N, true );
       if( sRet == 1 )  // error: no SLE solution!
          break;

// SOLVED: solution of linear matrix has been obtained
      pmp->PCI = calcDikin( N, true);  // calc of MU values and Dikin criterion

      LM = calcLM( true ); // Calculation of descent step size LM
      if( LM < 1E-5 )
         Error( "E03IPM FIA-iteration:",
            "Too small step size - too slow convergence");
      if( LM > 1.)
         LM = 1.;

      // calculation of new prime solution approximation
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

   pmp->Ec=1;
   if( IT1 == pa->p.DP )
      iRet = 2;
   else
      iRet = 1;
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
//
int TMulti::InteriorPointsMethod( )
{
    int N, IT1,J,Z,iRet;
    double LM=0., LM1, FX1;
    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pmp->FIs )
      for( J=0; J<pmp->Ls; J++ )
            pmp->lnGmo[J] = pmp->lnGam[J];

    for(J=0;J<pmp->N;J++)
    {   pmp->U[J]=0.;
        pmp->C[J]=0.;
    }

    pmp->Ec=0;
    pmp->FX=GX( LM  );  // calculation of G(x)

    if( pmp->FIs ) // multicomponent phases are present
      for(Z=0; Z<pmp->FIs; Z++)
        pmp->YFA[Z]=pmp->XFA[Z];

//----------------------------------------------------------------------------
//  Main loop of IPM iterations
    for( IT1 = 0; IT1 < pa->p.IIM; IT1++, pmp->IT++ )
    {
        pmp->NR=pmp->N;
        if( pmp->LO ) // water-solvent is present
        {
            if( pmp->YF[0]<pmp->DSM && pmp->Y[pmp->LO]< pmp->lowPosNum *100.)
                pmp->NR=(short)(pmp->N-1);
        }
        N = pmp->NR;
//        memset( pmp->F, 0, pmp->L*sizeof(double));

        PrimalChemicalPotentials( pmp->F, pmp->Y, pmp->YF, pmp->YFA );

        // Setting weight multipliers for DC
        WeightMultipliers( false );
        // Making and solving the R matrix of IPM linearized equations
        iRet = SolverLinearEquations( N, false );
        if( iRet == 1 )
        { pmp->Ec=1;
          return 1;
        }

//SOLVED: got the dual solution u vector - calculating the Dikin criterion
//    of GEM IPM convergence
       f_alpha( );
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
                pmp->PD=2;     //  what does this mean?
            else  pmp->DX= 0.5 * pmp->PCI;
        }

       MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );

       pmp->FX=FX1;
       // Main IPM iteration done
       // Calculation of activity coefficients
        if( pmp->PD==3 )
            GammaCalc( LINK_UX_MODE );

        if( pmp->PHC[0] == PH_AQUEL && pmp->XF[0] <= pa->p.XwMin &&
             pmp->X[pmp->LO] <= pmp->lowPosNum*1e3 )    // bugfix 29.11.05 KD
        {
            pmp->XF[0]=0.;  // elimination of aqueous phase
            pmp->XFA[0] = 0.;
        }
        // Restoring vectors Y and YF
        Restoring_Y_YF();

#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
  pVisor->Update( false );
#endif
// STEPWISE (6)  Stop point at IPM() main iteration
STEP_POINT( "IPM Iteration" );
#endif
        if( pmp->PCI < pmp->DX )  // Dikin criterion satisfied - converged!
            break;
    } // end of main IPM cycle
//----------------------------------------------------------------------------

  if( IT1 >= pa->p.IIM)
   return 2; // bad convergence - too many IPM iterations!

 // Final calculation phase amounts and activity coefficients
  TotalPhases( pmp->X, pmp->XF, pmp->XFA );

  if(pmp->PZ && pmp->W1)
      Mol_u( pmp->Y, pmp->X, pmp->XF, pmp->XFA );

  if( pmp->PD==1 || pmp->PD == 2  /*|| pmp->PD == 3*/  )
        GammaCalc( LINK_UX_MODE);
   else
        ConCalc( pmp->X, pmp->XF, pmp->XFA );

   MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->X, pmp->B, pmp->C);
   return 0;
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
TMulti::MassBalanceResiduals( int N, int L, float *A, double *Y, double *B,
         double *C )
{
    int ii, jj, i;
//    for(J=0;J<N;J++)    obsolete - replaced by an optimized one 09.02.2007
//    {
//        C[J]=B[J];
//        for(I=0;I<L;I++)
//            C[J]-=(double)(*(A+I*N+J))*Y[I];
//    }

    for(ii=0;ii<N;ii++)
        C[ii]=B[ii];
    for(jj=0;jj<L;jj++)
     for( i=arrL[jj]; i<arrL[jj+1]; i++ )
     {  ii = arrAN[i];
         C[ii]-=(double)(*(A+jj*N+ii))*Y[jj];
     }
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
    if( LM<2.)
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

//Added by Sveta
//===================================================================

// Inserting minor quantities of DC which were zeroed off after simplex()
// (important for the automatic initial approximation with solution phases
//  or after Selekt2() algorithm
//
void TMulti::RaiseZeroedOffDCs( int iStart, int iEnd, double sfactor, int JJ )
{
  SPP_SETTING *pa = &TProfil::pm->pa;

  for(int i=iStart; i<iEnd; i++ )
  {
    if( JJ >=0 )
      pmp->YF[JJ] = 0.;

     switch( pmp->DCC[i] )
     {
       case DC_AQ_PROTON:
       case DC_AQ_ELECTRON:
       case DC_AQ_SPECIES:
          if( JJ>=0 || pmp->Y[i] < pa->p.DFYaq * sfactor )
               pmp->Y[i] =  pa->p.DFYaq * sfactor;
           break;
       case DC_AQ_SOLVCOM:
       case DC_AQ_SOLVENT:
            if( JJ>=0 || pmp->Y[i] < pa->p.DFYw * sfactor )
                pmp->Y[i] =  pa->p.DFYw * sfactor;
            break;
       case DC_GAS_H2O:
       case DC_GAS_CO2:
       case DC_GAS_H2:
       case DC_GAS_N2:
       case DC_GAS_COMP:
       case DC_SOL_IDEAL:
            if( JJ>=0 || pmp->Y[i] < pa->p.DFYid*sfactor )
                  pmp->Y[i] = pa->p.DFYid * sfactor;
             break;
       case DC_SOL_MINOR:
            if( JJ>=0 || pmp->Y[i] < pa->p.DFYh*sfactor )
                   pmp->Y[i] = pa->p.DFYh * sfactor;
             break;
       case DC_SOL_MAJOR:
            if( JJ>=0 || pmp->Y[i] < pa->p.DFYr * sfactor )
                  pmp->Y[i] =  pa->p.DFYr * sfactor;
             break;
       case DC_SCP_CONDEN:
             if( JJ>=0 || pmp->Y[i] < pa->p.DFYc * sfactor )
                  pmp->Y[i] =  pa->p.DFYc * sfactor;
              break;
                    // implementation for adsorption?
       default:
             if( JJ>=0 || pmp->Y[i] < pa->p.DFYaq *sfactor )
                   pmp->Y[i] =  pa->p.DFYaq * sfactor;
             break;
     }
    if( JJ >=0 )
     pmp->YF[JJ] += pmp->Y[i];
   } // i
}

// Adjustment of primal approximation according to kinetic constraints
void TMulti::LagrangeMultiplier()
{
    double E = 1E-8;  // pa.p.DKIN? Default min value of Lagrange multiplier p
//    E = 1E-30;

    for(int J=0;J<pmp->L;J++)
    {
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
                pmp->Y[J]=pmp->DLL[J]+E;  /* 1e-10: pa.DKIN ? */
        case UPPER_LIM:
            if( pmp->Y[J]>=pmp->DUL[J])
                pmp->Y[J]=pmp->DUL[J]-E;
            break;
        }
    }   // J
}

// Calculation of weight multipliers for DCs
void TMulti::WeightMultipliers( bool square )
{
  int J;
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
            Error( "E04IPM IPM-Iteration:",
                 "Error in codes of metastability constraints" );
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
#define  a(j,i) ((double)(*(pmp->A+(i)+(j)*Na)))
//#define  a(j,i) ((double)(*(pmp->A+(i)+(j)*N)))  obsolete
int TMulti::SolverLinearEquations( int N, bool initAppr )
{
  int ii,i, jj, kk, k, Na = pmp->N;
  double aa;
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

  Array2D<double> A(N,N, AA);
  Array1D<double> B(N, BB);

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
   pmp->U[ii] = B[ii];
  return 0;
}
#undef a

// Calculation of MU values (dual DC chemical potentials) and Dikin criterion
// Parameters:
// bool initAppr - Inital approximation (true) or main iteration of IPM (false)
// int N - dimension of the matrix R (number of equations)
double TMulti::calcDikin(  int N, bool initAppr )
{
  int  J;
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
     if( PCI > pmp->lowPosNum /* 1E-19*/ )
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
   int J;
   int Z = -1;
   double LM, LM1, Mu;

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
       LM = 1./sqrt(pmp->PCI);
  }
  return LM;
}

// Restoring primal vectors Y and YF
void TMulti::Restoring_Y_YF()
{
 int Z, I, JJ = 0;

 for( Z=0; Z<pmp->FI ; Z++ )
 {
   if( pmp->XF[Z] <= pmp->DSM ||
       pmp->PHC[Z] == PH_SORPTION &&
       ( pmp->XFA[Z] < TProfil::pm->pa.p.ScMin) )
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
    int i, NN = pmp->N - pmp->E;

    for(i=0;i<NN;i++)
      molB += pmp->B[i];

   sfactor = pow( molB, 0.4 )/7.7;
   return sfactor;
}

//===================================================================

// Checking Karpov phase stability criteria Fa for phases and DCs
//  using Selekt2() algorithm by Karpov & Chudnenko (1989)
//  modified by DK in 1995.
//  Returns 0, if a new IPM loop is needed;
//          1, if the solution is final.
//
void TMulti::PhaseSelect()
{
    short JJ,Z,I,J;
    double F1,F2,*F0;
    SPP_SETTING *pa = &TProfil::pm->pa;

    f_alpha( );
    F0 = pmp->Falp;

    (pmp->K2)++;
    pmp->MK=0;
    JJ= -1;
 //   II= -1;
    F1= F2= pmp->lowPosNum*10000.;  // 1E-16;

    for(Z=0;Z<pmp->FI;Z++)
    {
        if( F0[Z]>F1 && pmp->YF[Z]<pmp->lowPosNum )
        {
            F1=F0[Z];  // selection of max Fa and corresponding phase index
            JJ=Z;
        }
        if( F0[Z]>F2 && pmp->YF[Z]>pmp->lowPosNum )
        {
            F2=F0[Z];
           // II=Z;
        }
    }
    if( F1 > pa->p.DF && JJ >= 0 )
    {
        double sfactor;
        // There is a phase for which DF (0.01) is exceeded
        sfactor = calcSfactor();
        do
        {  // inserting all phases with  F1 > DFM
            // insert this phase and set Y[j] for its components
            // with account for asymmetry and non-ideality
            for( J=0, I=0; I<JJ; I++ )
                 J+=pmp->L1[I];

            RaiseZeroedOffDCs( J, J+pmp->L1[JJ], sfactor, JJ );

            pmp->FI1++;  // check phase rule
            if( pmp->FI1 >= pmp->NR+1 )
            { // No more phase can be inserted
                break;
            }
            // find a new phase to insert, if any exists
            F1= pmp->lowPosNum*10000.; // 1e-16
            JJ = -1;
            for( Z=0; Z<pmp->FI; Z++ )
                if( F0[Z] > F1 && pmp->YF[Z] < pmp->lowPosNum )
                {
                    F1=F0[Z];
                    JJ=Z;
                }
        }
        while( F1 > pa->p.DF && JJ >= 0 );
        // end of insertion cycle
        J=0;   // insert primeal GEM IPM solution
        for(Z=0;Z<pmp->FIs;Z++)
        {
            if( pmp->YF[Z ] > pmp->lowPosNum /* 1E-18 */ )
            {
                pmp->YF[Z]=0.;
                for(I=J;I<J+pmp->L1[Z];I++)
                {
                   if( pmp->Y[I] < pmp->lowPosNum ) // Check what to insert !
                       pmp->Y[I] = pa->p.DFM; // lowPosNum should be used ?
                    pmp->YF[Z]+=pmp->Y[I]; // calculate new amounts of phases
                }
            }
            J+=pmp->L1[Z];
        }
#ifndef IPMGEMPLUGIN
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelection"
#endif
#endif
        if( pmp->K2>1 )
        { // more then the first step - but the IPM solution has not been improved
            for(I=0;I<pmp->L;I++)
                if( fabs(pmp->Y[I]- pmp->XY[I]) > pa->p.DF*pmp->Y[I] ) // Check!
                    goto S6;
            pmp->PZ=2;
            goto S4;
        }
S6: // copy of X vector changed by Selekt2() algorithm
        for(I=0;I<pmp->L;I++)
            pmp->XY[I]=pmp->Y[I];
// Copy here also to pmp->X[I]?
        return;
    }
S4: // No phases to insert or no Fa distortions found
    // Successful end of iterations of SELEKT2
    pmp->MK=1;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Internal memory allocation for IPM performance optimization
// (since version 2.2.0)
//
void TMulti::Alloc_A_B( int newN )
{
  if( AA && BB && newN==sizeN )
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
#define  a(j,i) ((double)(*(pmp->A+(i)+(j)*pmp->N)))
void TMulti::Build_compressed_xAN()
{
 int ii, jj, k;

 // Free old memory allocation
 Free_compressed_xAN();

 // Calculate sizes
 k = 0;
 for( jj=0; jj<pmp->L; jj++ )
   for( ii=0; ii<pmp->N; ii++ )
     if( fabs( a(jj,ii) ) > 1e-12 )
       k++;

 // Allocate memory
 arrL = new int[pmp->L+1];
 arrAN = new int[k];

 // Set indexes
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
    { delete[] arrL; arrL = 0; }
  if( arrAN )
    { delete[] arrAN; arrAN = 0; }
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
//--------------------- End of ipm_main.cpp ---------------------------
