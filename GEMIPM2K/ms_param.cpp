//-------------------------------------------------------------------
// $Id: ms_param.cpp 1392 2009-08-10 13:39:26Z gems $
//
// Copyright  (C) 1992,2007 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation  of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry and
// of the GEMIPM2K standalone code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifdef __unix__
#include <unistd.h>
#endif

#include <math.h>
#include "m_param.h"
#include "num_methods.h"
#include "gdatastream.h"
#include "node.h"

TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481,
                             H2O_mol_to_kg = 55.50837344,
                               Min_phys_amount = 1.66e-24;

enum volume_code {  // Codes of volume parameter ???
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

SPP_SETTING pa_ = {
    "GEM-Selektor v3.0.0: Numerical flags, controls and thresholds",
    {   // Typical default set (interim, 08.12.2009)
        1,  /* PC */  2,     /* PD */   3,   /* PRD */
        1,  /* PSM  */ 300,  /* DP */   30,   /* DW */
        -2, /* DT */     200,   /* PLLG */   1,  /* PE */
        1000,   /* IIM */
        1e-5, /* DG */   1e-11,  /* DHB */  1e-13,  /* DS */
        7e-7,  /* DK */  0.01,  /* DF */  0.1,  /* DFM */
        1e-6,  /* DFYw */  1e-6,  /* DFYaq */    1e-6,  /* DFYid */
        1e-6,  /* DFYr,*/  1e-6,  /* DFYh,*/   1e-6,  /* DFYc,*/
        1e-7, /* DFYs, */  1e-17,  /* DB */   -1.0,   /* AG */
        -0.5,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
        0.001, /* GAS */   12.05,  /* DNS */   1e-9,  /* XwMin, */
        1e-7,  /* ScMin, */  1e-23, /* DcMin, */   1e-13, /* PhMin, */
        3e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
        1e-7,  /* DKIN  */ 0,  /* tprn */
  },
}; // SPP_SETTING


void BASE_PARAM::write(fstream& oss)
{
  short arr[10];

  arr[0] = PC;
  arr[1] = PD;
  arr[2] = PRD;
  arr[3] = PSM;
  arr[4] = DP;
  arr[5] = DW;
  arr[6] = DT;
  arr[7] = PLLG;
  arr[8] = PE;
  arr[9] = IIM;

	oss.write( (char*)arr, 10*sizeof(short) );
    oss.write( (char*)&DG, 28*sizeof(double) );
    oss.write( (char*)&tprn, sizeof(char*) );
}

void SPP_SETTING::write(fstream& oss)
{
    oss.write( ver, TDBVERSION );
    p.write( oss );
}

TProfil::TProfil( TMulti* amulti )
{
    pa= pa_;
    multi = amulti;
    pmp = multi->GetPM();
}

// test result GEM IPM calculation of equilibrium state in MULTI
long int TProfil::testMulti(  )
{
  if( pmp->MK || pmp->PZ )
  {
	if( pa.p.PSM == 2 )
	{
      fstream f_log("ipmlog.txt", ios::out|ios::app );
      f_log << "Warning " << pmp->stkey << ": " <<  pmp->errorCode << ":" << endl;
      f_log << pmp->errorBuf << endl;
	}
   return 1L;
  }
  return 0L	;
}

// GEM IPM calculation of equilibrium state in MULTI
// Modified on 10.09.2007 to return elapsed GEMIPM2 runtime in seconds
// Modified on 15.11.2007 to return more detailed info on FIA and IPM iterations
// and precision refinement loops
// Modified on 22.01.2008 to implement "smart PIA" mode
//
double TProfil::calcMulti( long int& RefinLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
    pmp = multi->GetPM();
    pmp->t_start = clock();
    pmp->t_end = pmp->t_start;

//    multi->MultiCalcInit( 0 );
    pmp->ITF = pmp->ITG = 0;
FORCED_AIA:
    multi->MultiCalcInit( "GEMIPM2K" );
    if( pmp->pNP )
    {
	   if( pmp->ITaia <=30 )       // Foolproof
		  pmp->IT = 30;
	   else
		  pmp->IT = pmp->ITaia;     // Setting number of iterations for the smoothing parameter
    }
    bool IAstatus;
    IAstatus = multi->AutoInitialApprox( );
    if( IAstatus == false )
    {
        multi->MultiCalcIterations(-1 );    // Calling main IPM2 sequence
    }
    RefinLoops_ = pmp->W1 + pmp->K2 - 1;
    NumIterFIA_ = pmp->ITF;
    NumIterIPM_ = pmp->ITG;

    if( IAstatus == true )
        goto FINISHED;  // Only pure phases - simplex solution is Ok

    if( pmp->MK || pmp->PZ ) // no good solution
        goto FINISHED;

    pmp->IT = pmp->ITG;   // This is to provide correct number of IPM iterations to upper levels

    if( pa.p.PRD < 0 && pa.p.PRD > -50 ) // max 50 loops
    {  // Test refinement loops for highly non-ideal systems Added here by KD on 15.11.2007
       long int pp, pNPo = pmp->pNP,  TotW1 = pmp->W1+pmp->K2-1,
ITstart = 10,        TotIT = pmp->IT;  //  ITold = pmp->IT,
       pmp->pNP = 1;
       for( pp=0; pp < abs(pa.p.PRD); pp++ )
       {
    	  pmp->IT = ITstart; // Value is important for refinement in highly non-ideal systems!
          if( multi->AutoInitialApprox( ) == false )
          {
             multi->MultiCalcIterations( pp );
          }
          TotIT += pmp->IT - ITstart;
          TotW1 += pmp->W1+pmp->K2-1;
          if( pmp->MK || pmp->PZ ) // no good solution
            break;	 // goto FINISHED;

       }  // end pp loop

       pmp->pNP = pNPo;
       pmp->IT = TotIT; // ITold;

       RefinLoops_ = TotW1;
       NumIterFIA_ = pmp->ITF;  //   TotITF;
       NumIterIPM_ = pmp->ITG;  //   TotITG;
    }

FINISHED:

	if( pmp->MK == 2 )
	{	if( pmp->pNP )
         {
    	    pmp->pNP = 0;
    	    pmp->MK = 0;
//cout << pmp->stkey << " return simplex W1 " <<  pmp->W1 << " pmp->Ec " << pmp->errorBuf << endl;
    	    goto FORCED_AIA;  // Trying again with AIA set after bad SIA
         }
    	else
    		Error( pmp->errorCode ,pmp->errorBuf );
	}
    if( pmp->MK || pmp->PZ ) // no good solution
    {
//cout << "Iter"  << " MK " << pmp->MK << " PZ " << pmp->PZ << " " << pmp->errorCode << endl;
    }
    else // only test 30/01/2009 SD
    {	int iB = multi->CheckMassBalanceResiduals( pmp->X );
    	if( iB >= 0 )
    	{
    	   	    Error( "Point1  -  After Finish calculation",pmp->errorBuf );
    	}

    }
    pmp->t_end = clock();
    pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);
return pmp->t_elap_sec;
}

void TProfil::outMulti( GemDataStream& ff, gstring& /*path*/  )
{
	 short arr[10];

	  arr[0] = pa.p.PC;
	  arr[1] = pa.p.PD;
	  arr[2] = pa.p.PRD;
	  arr[3] = pa.p.PSM;
	  arr[4] = pa.p.DP;
	  arr[5] = pa.p.DW;
	  arr[6] = pa.p.DT;
	  arr[7] = pa.p.PLLG;
	  arr[8] = pa.p.PE;
	  arr[9] = pa.p.IIM;

	ff.writeArray( arr, 10 );
    ff.writeArray( &pa.p.DG, 28 );
    multi->to_file( ff );
}

void TProfil::outMultiTxt( const char *path, bool append  )
{
    multi->to_text_file( path, append );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( GemDataStream& ff )
{
	 short arr[10];

	 ff.readArray( arr, 10 );
	  pa.p.PC = arr[0];
	  pa.p.PD = arr[1];
	  pa.p.PRD = arr[2];
	  pa.p.PSM = arr[3];
	  pa.p.DP = arr[4];
	  pa.p.DW = arr[5];
	  pa.p.DT = arr[6];
	  pa.p.PLLG = arr[7];
	  pa.p.PE = arr[8];
	  pa.p.IIM = arr[9];

	  ff.readArray( &pa.p.DG, 28 );
      multi->from_file( ff );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( const char* path )
{
      multi->from_text_file_gemipm( path);
}

 bool load = false;

 #include "io_arrays.h"
/*
 void TProfil::test_G0_V0_H0_Cp0_DD_arrays( long int nT, long int nP )
 {
   long int kk, jj, ii, l1, l2, xTP, lev = 11;
   double cT, cP, dT, dP;
   double *G0, *V0, *H0,* *Cp0, *S0, *A0, *U0,* *denW, *epsW, *denWg, *epsWg;
   DATACH  *CSD = TNode::na->pCSD();
   fstream ff("lagrange_T_11_3.out.txt", ios::out );
   ErrorIf( !ff.good() , "lagrange.out", "Fileopen error");

   G0 =  new double[CSD->nDC*nT*nP];
//   V0 =  new double[CSD->nDC*nT*nP];
//   H0 =  new double[TProfil::pm->mup->L*nT*nP];
   Cp0 = new double[CSD->nDC*nT*nP];
   S0 = new double[CSD->nDC*nT*nP];
//   A0 = new double[TProfil::pm->mup->L*nT*nP];
//   U0 = new double[TProfil::pm->mup->L*nT*nP];
   denW =  new double[5*nT*nP];
   epsW =  new double[5*nT*nP];

   cT = CSD->TCval[0];
   cP = CSD->Pval[0];
   if(nT > 1)
     dT = (CSD->TCval[CSD->nTp-1]-CSD->TCval[0])/(nT-1);
   else
	dT = 0;
   if(nP > 1 )
     dP = (CSD->Pval[CSD->nPp-1]-CSD->Pval[0])/(nP-1);
   else
	  dP = 0;
   for(  ii=0; ii<nT; ii++)
   {
//     cT += dT;
     for(  jj=0; jj<nP; jj++)
     {
//       cP += dP;
       xTP = TNode::na->check_grid_TP( cT, cP ); // changed to Kelvin
       for( kk=0; kk<5; kk++)
       {
         l1 = ( kk * nP + jj) * nT + ii;
         l2  =  kk * CSD->nPp * CSD->nTp;

         if( xTP >= 0 )
          {
           denW[l1] = CSD->denW[l2+xTP];
           epsW[l1] = CSD->epsW[l2+xTP];
          }
          else
          {
            denW[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->denW+l2,
                            cP, cT, CSD->nTp, CSD->nPp,lev );
            epsW[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->epsW+l2,
                            cP, cT, CSD->nTp, CSD->nPp,lev );
          }
       }
       // copy to arrays
       for( kk=0; kk<CSD->nDC; kk++)
       {
           l1 = ( kk * nP + jj) * nT + ii;
           l2  =  kk * CSD->nPp * CSD->nTp;
          if( xTP >= 0 )
           {
        	 G0[l1] = CSD->G0[ l2+xTP];
           	 S0[l1] = CSD->S0[ l2+xTP];
           	 Cp0[l1] = CSD->Cp0[ l2+xTP];
            }
            else
            {
              G0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->G0+l2,
                              cP, cT, CSD->nTp, CSD->nPp,lev );
              S0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->S0+l2,
                              cP, cT, CSD->nTp, CSD->nPp,lev );
              Cp0[l1] = LagranInterp( CSD->Pval, CSD->TCval, CSD->Cp0+l2,
                              cP, cT, CSD->nTp, CSD->nPp,lev );
            }
        }
       cP += dP;
     }
     cT += dT;
   }

//  fstream ff("lagrange_T_5_3.out.txt", ios::out );
//  ErrorIf( !ff.good() , "lagrange.out", "Fileopen error");
  TPrintArrays  prar(0, 0, ff);
  prar.writeArray(  "denW", denW, 5**nP*nT, nP*nT );
  prar.writeArray(  "epsW", epsW, 5**nP*nT, nP*nT );
  prar.writeArray(  "S0_Ca+2", S0, nP*nT,  nP*nT ); //Ca+2
  prar.writeArray(  "S0_H2O", S0+13*nP*nT, nP*nT, nP*nT); //H2O@
  prar.writeArray(  "S0_Brc", S0+18*nP*nT, nP*nT, nP*nT ); //Brc
  prar.writeArray(  "G0_Ca+2", G0, nP*nT, nP*nT ); //Ca+2
  prar.writeArray(  "G0_H2O", G0+13*nP*nT, nP*nT, nP*nT ); //H2O@
  prar.writeArray(  "G0_Brc", G0+18*nP*nT, nP*nT, nP*nT ); //Brc
  prar.writeArray(  "Cp0_Ca+2", Cp0, nP*nT, nP*nT ); //Ca+2
  prar.writeArray(  "Cp0_H2O", Cp0+13*nP*nT, nP*nT, nP*nT ); //H2O@
  prar.writeArray(  "Cp0_Brc", Cp0+18*nP*nT, nP*nT, nP*nT ); //Brc
//  prar.writeArray(  "V0", V0, CSD->nDC*nP*nT, nP*nT );
//  prar.writeArray(  "G0", G0, CSD->nDC*nP*nT, nP*nT );

   delete[] G0;
//   delete[] V0;
//   delete[] H0;
   delete[] Cp0;
   delete[] S0;
//   delete[] A0;
//   delete[] U0;
   delete[] denW;
   delete[] epsW;
 }
*/
// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
// (only used in standalone GEMIPM2K version)
 //
void TMulti::CompG0Load()
{
  long int j, jj, k, xTP, jb, je=0;
  double Go, Gg, Vv, h0=0., S0 = 0., Cp0= 0., a0 = 0., u0 = 0.;
  double TK, P, PPa;

  DATACH  *dCH = TNode::na->pCSD();
//  DATABR  *dBR = TNodeArray::na->pCNode();

  TK = TNode::na->cTK();
  PPa = TNode::na->cP();
  P = PPa/bar_to_Pa;
// if( dCH->nTp <=1 && dCH->nPp <=1 )
  if( dCH->nTp <1 || dCH->nPp <1 || TNode::na->check_TP( TK, PPa ) == false )
  {
	  char buff[256];
	  sprintf( buff, " Temperature %g or pressure %g out of range, or no T/D data are provided\n",
			  TK, PPa );
	  Error( "ECompG0Load: " , buff );
      return;
  }
  xTP = TNode::na->check_grid_TP( TK, PPa );

 if( load && fabs( pmp->Tc - TK ) < dCH->Ttol /*1.e-10*/ &&
            fabs( pmp->P - P ) < dCH->Ptol /*1.e-10*/ )
   return;    //T, P not changed - problematic for UnSpace!

 pmp->T = pmp->Tc = TK;
 pmp->TC = pmp->TCc = TK-C_to_K;
 pmp->P = pmp->Pc = P;

// if( dCH->ccPH[0] == PH_AQUEL )
// {
   for( k=0; k<5; k++ )
   {
     jj =  k * dCH->nPp * dCH->nTp;
     if( xTP >= 0 )
      {
       pmp->denW[k] = dCH->denW[jj+xTP]/1e3;
       pmp->epsW[k] = dCH->epsW[jj+xTP];
       pmp->denWg[k] = dCH->denWg[jj+xTP]/1e3;
       pmp->epsWg[k] = dCH->epsWg[jj+xTP];
      }
     else
     {
       pmp->denW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denW+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,6 )/1e3;// from test denW enough
       pmp->epsW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsW+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );// from test epsW enough
       pmp->denWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denWg+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 )/1e3;
       pmp->epsWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsWg+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
     }
  }
// }
// else
// {
//   pmp->denW[0] = 1.;
//   pmp->epsW[0] = 78.;
// }

 pmp->RT = R_CONSTANT * pmp->Tc;
 pmp->FRT = F_CONSTANT/pmp->RT;
 pmp->lnP = 0.;
 if( P != 1. ) // ???????
   pmp->lnP = log( P );

 for( k=0; k<pmp->FI; k++ )
 {
   jb = je;
   je += pmp->L1[k];
   // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
   // depending on the presence of these arrays in DATACH and Multi structures
    for( j=jb; j<je; j++ )
    {
      jj =  j * dCH->nPp * dCH->nTp;
      if( xTP >= 0 )
      {
        Go = dCH->G0[ jj+xTP];
        Vv = dCH->V0[ jj+xTP]*1e5;
        if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
        if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
        if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
        if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
        if( dCH->U0 ) h0 = dCH->U0[ jj+xTP];
      }
     else
     {
       Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
       Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 5 )*1e5;
       if( dCH->S0 ) S0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->S0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 4 ); // from test S0[Ca+2] enough
       if( dCH->H0 ) h0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->H0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
       if( dCH->Cp0 ) Cp0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->Cp0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp, 3 ); // from test Cp0[Ca+2] not more
       if( dCH->A0 ) a0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->A0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
       if( dCH->U0 ) u0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->U0+jj,
                          PPa, TK, dCH->nTp, dCH->nPp,5 );
     }
     if( pmp->tpp_G )
    	  pmp->tpp_G[j] = Go;
     if( pmp->Guns )
           Gg = pmp->Guns[j];
     else
           Gg = 0.;
     pmp->G0[j] = Cj_init_calc( Go+Gg, j, k ); // Inside this function, pmp->YOF[k] can be added!
     switch( pmp->PV )
     { // put mol volumes of components into A matrix or into the vector of molar volumes
       // to be checked!
       case VOL_CONSTR:
           if( pmp->Vuns )
              Vv += pmp->Vuns[jj];
           // pmp->A[j*pmp->N+xVol] = tpp->Vm[jj]+Vv;
             pmp->A[j*pmp->N] = Vv; // !!  error
       case VOL_CALC:
       case VOL_UNDEF:
    	     if( pmp->tpp_Vm )
    	    	  pmp->tpp_Vm[j] = Vv;
              if( pmp->Vuns )
                     Vv += pmp->Vuns[j];
 	          pmp->Vol[j] = Vv  * 10.;
              break;
     }
     if( pmp->S0 ) pmp->S0[j] = S0;
     if( pmp->H0 ) pmp->H0[j] = h0;
     if( pmp->Cp0 ) pmp->Cp0[j] = Cp0;
     if( pmp->A0 ) pmp->A0[j] = a0;
     if( pmp->U0 ) pmp->U0[j] = u0;
   }
 }
 load = true;
}

// GEM IPM calculation of equilibrium state in MULTI
void TMulti::MultiCalcInit( const char* /*key*/ )
{
  long int  j,k;

  SPP_SETTING *pa = &TProfil::pm->pa;
   // SD 09/03/09 strncpy( pmp->stkey, key, EQ_RKLEN );  // needed for the ipmlog error diagnostics
    pmp->Ec = pmp->K2 = pmp->MK = 0;
    pmp->PZ = pa->p.DW; // IPM-2 default
    pmp->W1 = 0;
    pmp->is = 0;
    pmp->js = 0;
    pmp->next  = 0;
    pmp->ln5551 = log( H2O_mol_to_kg ); //  4.016533882  ln(55.50837344)  4.0165339;
//    pmp->lowPosNum = pa->p.DcMin;  obsolete
    pmp->lowPosNum = Min_phys_amount;
    pmp->logXw = -16.;
    pmp->logYFk = -9.;
    pmp->FitVar[0] = 0.0640000030398369; // pa->aqPar[0]; setting T,P dependent b_gamma parameters
    pmp->DXM = pa->p.DK;
    RescaleToSize( true );  // Added to set default cutoffs/inserts 30.08.2009 DK
    pmp->T0 = 273.15;    // not used anywhere
    pmp->FX = 7777777.;
    pmp->YMET = 0;
//    pmp->PCI = 0.0;
    pmp->PCI = 1.0;

    // calculating mass of the system
    pmp->MBX = 0.0;
    for(long int i=0; i<pmp->N; i++ )
    {
        pmp->MBX += pmp->B[i] * pmp->Awt[i];
    }
    pmp->MBX /= 1000.;

    // optimization 08/02/2007 - allocation of A matrix index lists and IPM work arrays
    Alloc_internal(); // performance optimization 08/02/2007

    if(  pmp->pNP )     // Checking if this is SIA or AIA mode
    {   //  Smart Initial Approximation mode
        for( j=0; j< pmp->L; j++ )
          pmp->X[j] = pmp->Y[j];
 //       pmp->IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
        TotalPhases( pmp->X, pmp->XF, pmp->XFA );
        ConCalc( pmp->X, pmp->XF, pmp->XFA);
 //       pmp->IT = pmp->ITaia;
       if( pmp->FIs )
       {
           // Load activity coeffs for phases-solutions
          long int jb, je=0;
             for( k=0; k<pmp->FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += pmp->L1[k];
            // Load activity coeffs for phases-solutions
            for( j=jb; j< je; j++ )
            {
               pmp->lnGmo[j] = pmp->lnGam[j];
               if( fabs( pmp->lnGam[j] ) <= 84. )
      //                pmp->Gamma[j] = exp( pmp->lnGam[j] );
                      pmp->Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
               else pmp->Gamma[j] = 1.0;
             } // j
          }  // k
       }
    }
    else // Simplex initial approximation to be done (AIA mode)
    {
    	for( j=0; j<pmp->L; j++ )
    	{                           // cleaning work vectors
    		pmp->X[j] = pmp->Y[j] = pmp->lnGam[j] = pmp->lnGmo[j] = 0.0;
    		pmp->Gamma[j] = 1.0;
                pmp->MU[j] = 0.;// SD 06/12/2009
                pmp->XU[j] = 0.;// SD 06/12/2009
                pmp->EMU[j] = 0.;// SD 06/12/2009
                pmp->NMU[j] = 0.;// SD 06/12/2009
                pmp->W[j] = 0.;// SD 06/12/2009
                pmp->F[j] = 0.;// SD 06/12/2009
                pmp->F0[j] = 0.;// SD 06/12/2009
    	}
    }
    CompG0Load(); // Loading thermodynamic data into MULTI structure

 
     pmp->PD = abs(pa->p.PD);
	 //           SolModLoad() and Phase scripts cannot be used here!
    
    // recalculating kinetic restrictions for DC amounts
     if( pmp->pULR && pmp->PLIM )
          Set_DC_limits(  DC_LIM_INIT );

     bool AllPhasesPure = true;
     // checking if all phases are pure
     for( k=0; k < pmp->FI; k++ )
         if( pmp->L1[k] > 1 )
             AllPhasesPure = false;
     if( !AllPhasesPure )     // bugfix DK 09.03.2010   was if(!pmp->FIs)
    {
      long int jb, je=0;
     // Calc Eh, pe, pH,and other stuff
    if( pmp->E && pmp->LO && pmp->pNP )
    {
    	ConCalc( pmp->X, pmp->XF, pmp->XFA);
    	IS_EtaCalc();
        if( pmp->Lads )  // Calling this only when sorption models are present

        {
    	   for( k=0; k<pmp->FIs; k++ )
    	   { // loop on solution phases
    	      jb = je;
    	      je += pmp->L1[k];
    	      if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
    	      {
    		     if( pmp->PHC[0] == PH_AQUEL && pmp->XF[k] > pmp->DSM
    		       && (pmp->XFA[0] > pmp->ScMinM && pmp->XF[0] > pmp->XwMinM ))  // fixed 30.08.2009 DK
    		       GouyChapman( jb, je, k );  // getting PSIs - electrical potentials on surface planes
    	      }
    	   }
        }
    }
    GammaCalc( LINK_TP_MODE);   // Computing DQF, FugPure and G wherever necessary
                                      // Activity coeffs are restored from lnGmo
 }
   else {  // no multi-component phases
        pmp->PD = 0;  pmp->pNP = 0;
  //      pmp->pIPN = 1;
    }

}
//-------------------------------------------------------------------------
// internal functions

// read string as: "<characters>"
istream& f_getline(istream& is, gstring& str, char delim)
{
    char ch;
    is.get(ch);
    str="";

    while( is.good() && ( ch==' ' || ch=='\n' || ch== '\t') )
        is.get(ch);
    if(ch == '\"')
        is.get(ch);
    while( is.good() &&  ch!=delim && ch!= '\"' )
    {
        str += ch;
        is.get(ch);
    }
    while( is.good() &&  ch!=delim )
            is.get(ch);

   return is;
}

gstring u_makepath(const gstring& dir,
           const gstring& name, const gstring& ext)
{
    gstring Path(dir);
    if( dir != "")
      Path += "/";
    Path += name;
    Path += ".";
    Path += ext;

    return Path;
}

void u_splitpath(const gstring& Path, gstring& dir,
            gstring& name, gstring& ext)
{
    size_t pos = Path.rfind("/");
    if( pos != npos )
        dir = Path.substr(0, pos), pos++;
    else
        dir = "",    pos = 0;

    size_t pose = Path.rfind(".");
    if( pose != npos )
    {
        ext = Path.substr( pose+1, npos );
        name = Path.substr(pos, pose-pos);
    }
    else
    {
        ext = "";
        name = Path.substr(pos, npos);
    }
}

const long int bGRAN = 20;

// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path,
		long int& nElem, char delim ))[fileNameLength]
{
  long int ii, bSize = bGRAN;
  char  (*filesList)[fileNameLength];
  char  (*filesListNew)[fileNameLength];
  filesList = new char[bSize][fileNameLength];
  gstring name;

// Get path
   gstring path_;
   gstring flst_name = f_name;
   unsigned long int pos = flst_name.rfind("/");
   path_ = "";
   if( pos < npos )
      path_ = flst_name.substr(0, pos+1);
   strncpy( Path, path_.c_str(), 256-fileNameLength);
   Path[255] = '\0';

//  open file stream for the file names list file
   fstream f_lst( f_name/*flst_name.c_str()*/, ios::in );
   ErrorIf( !f_lst.good(), f_name, "Fileopen error");

// Reading list of names from file
  nElem = 0;
  while( !f_lst.eof() )
  {
	f_getline( f_lst, name, delim);
    if( nElem >= bSize )
    {    bSize = bSize+bGRAN;
         filesListNew = new char[bSize][fileNameLength];
         for( ii=0; ii<nElem-1; ii++ )
		   strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	     delete[] filesList;
		 filesList =  filesListNew;
	}
    strncpy( filesList[nElem], name.c_str(), fileNameLength);
	filesList[nElem][fileNameLength-1] = '\0';
    nElem++;
  }

  // Realloc memory for reading size
  if( nElem != bSize )
  {
    filesListNew = new char[nElem][fileNameLength];
    for(  ii=0; ii<nElem; ii++ )
	  strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	delete[] filesList;
	filesList =  filesListNew;
  }

  return filesList;
}


// ------------------ End of ms_param.cpp -----------------------




