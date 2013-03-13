//-------------------------------------------------------------------
// $Id$
//
// Stub implementation of TKinMet class for further development
//
// Copyright (C) 2013  D.Kulik, B.Thien, S.Dmitrieva
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

#include <cmath>
#include <cstdio>
#include "verror.h"
#include "s_kinmet.h"


//=============================================================================================
// TKinMet base class constructor (simulation of phase metastability and MWR kinetics)
// (c) DK March 2013
//
//=============================================================================================

/// generic constructor (new)
TKinMet::TKinMet( const KinMetData *kmd ):
    KinRateCode(kmd->KinRateCod_),    KinDirCode(kmd->KinDirCod_),   KinUptCode(kmd->KinUptCod_),
    KinLnkCode(kmd->KinLnkCod_),  KinSizedCode(kmd->KinSizedCod_),  KinRezdCode(kmd->KinRezdCod_),
    NComp(kmd->NComp_), nlPh(kmd->nlPh_), nlPc(kmd->nlPc_), nPRk(kmd->nPRk_), nSkr(kmd->nSkr_),
    nrpC(kmd->nrpC_), naptC(kmd->naptC_), nAscC(kmd->nAscC_), numpC(kmd->numpC_), iRes4(kmd->iRes4_),
    R_CONST(8.31451), T_k(kmd->T_k_), P_bar(kmd->P_bar_), kTau(kmd->kTau_), kdT(kmd->kdT_),
    IS(kmd->IS_), pH(kmd->pH_),  pe(kmd->pe_),  Eh(kmd->Eh_),  nPh(kmd->nPh_), mPh(kmd->mPh_),
    vPh(kmd->vPh_), sAPh(kmd->sAPh_), LaPh(kmd->LaPh_), OmPh(kmd->OmPh_),
    sSA(kmd->sSA_),  sgw(kmd->sgw_),  sgg(kmd->sgg_),  rX0(kmd->rX0_),  hX0(kmd->hX0_),
    sVp(kmd->sVp_), sGP(kmd->sGP_), nPul(kmd->nPul_), nPll(kmd->nPll_)
{


    // pointer assignments
    arfeSAr = kmd->arfeSAr_;   //  [nPRk]
    arAscp = kmd->arAscp_;     //  [nAscC]
    SM = kmd->SM_;             //  [NComp]
    arDCC = kmd->arDCC_;       //  [NComp]
    arlPhC = kmd->arlPhC_;     //  [nlPh]
    arocPRk = kmd->arocPRk_;   // [nPRk]
    arxSKr = kmd->arxSKr_;     // [nSKr]
    arxlPh = kmd->arxlPh_;     // [nlPh]
    arym = kmd->arym_;         // L
    arla = kmd->arla_;         // L
    arnx = kmd->arnx_;         // [NComp]
    arnxul = kmd->arnxul_;     // [NComp]  (DUL)
    arnxll = kmd->arnxll_;     // [NComp]  (DLL)
    arWx = kmd->arWx_;         // NComp
    arVol = kmd->arVol_;       // NComp

// 2D and 3D arrays
    arlPhc = NULL;
    arrpCon = NULL;
    arapCon = NULL;
    arUmpCon = NULL;
    alloc_kinrtabs( );
    init_kinrtabs( kmd->arlPhc_, kmd->arrpCon_,  kmd->arapCon_,  kmd->arUmpCon_);

/*
    // Work data and kinetic law calculation results

    TKinReact arPRt[]; /// work array of parameters and results for 'parallel reaction' terms [nPRk]

    double spcfu[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]
    double spcfl[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]

    double kTot;   /// Total rate constant (per m2 phase surface area)
    double rTot;   /// Current total MWR rate (mol/s)
    double vTot;   /// Total surface propagation velocity (nm/s)

    double sSAcor; /// Corrected specific surface area (m2/g)
    double sAph_c; /// Corrected surface area of the phase (m2/g)





    lnGamConf = new double[NComp];
    lnGamRecip = new double[NComp];
    lnGamEx = new double[NComp];  // Work arrays for lnGamma components - should we zero off?
    for (long int i=0; i<NComp; i++)
    {
       lnGamConf[i] = 0.0;
       lnGamRecip[i] = 0.0;
       lnGamEx[i] = 0.0;
    }
    // initialize phase properties
    Gex = 0.0;

*/
}

// Destructor
TKinMet::~TKinMet()
{
  free_kinrtabs();
}

/// allocates memory for multisite data
void TKinMet::alloc_kinrtabs()
{
   long int j, s, lp, pr;

   if( nlPh && nlPc )
   {
     arlPhc = new double *[nlPh];
     for( lp=0; lp<nlPh; lp++)
     {
          arlPhc[lp] = new double[nlPc];
     }
   }

   if( nPRk && nrpC )
   {
     arrpCon = new double *[nPRk];
     for( pr=0; pr<nPRk; pr++)
     {
          arrpCon[pr] = new double[nrpC];
     }
   }

   if( NComp && numpC )
   {
     arUmpCon = new double *[NComp];
     for( j=0; j<NComp; j++)
     {
          arUmpCon[j] = new double[numpC];
     }
   }

   if( nPRk && nSkr && naptC )
   {
        arapCon = new double **[nPRk];
        for(j=0; j<nPRk; j++)
        {
            arapCon[j]   = new double *[nSkr];
            for(s=0; s<nSkr; s++)
            {
                arapCon[j][s] = new double [naptC];
            }
        }
    }
}

/// returns 0 if o.k. or some arrays were not allocated.
/// returns -1 if error was encountered.
//
long int TKinMet::init_kinrtabs( double *p_arlPhc, double *p_arrpCon,  double *p_arapCon,  double *p_arUmpCon  )
{
    long int j, i, s, lp, pr;

    if( arlPhc ) {

        for( lp=0; lp<nlPh; lp++)
            for( i=0; i<nlPc; i++)
                arlPhc[lp][i] = p_arlPhc[nlPc*lp+i];
    }
    if( arrpCon ) {

        for( pr=0; pr<nPRk; pr++)
            for( i=0; i<nrpC; i++)
                arrpCon[pr][i] = p_arrpCon[nrpC*pr+i];
    }
    if( arUmpCon ) {

        for( j=0; j<NComp; j++)
            for( i=0; i<numpC; i++)
                arUmpCon[j][i] = p_arUmpCon[numpC*j+i];
    }
    if( arapCon ) {

        for( j=0; j<nPRk; j++)
            for( s=0; s<nSkr; s++)
                for( i=0; i<numpC; i++)
                    arapCon[j][s][i] = p_arapCon[numpC*nSkr*j+numpC*s+i];
    }
    return 0;
}

/// frees memory for TKinMet tables
//
void TKinMet::free_kinrtabs()
{
    long int j, s, lp, pr;

    if( arlPhc )
    {
      for( lp=0; lp<nlPh; lp++)
      {
           delete[]arlPhc[lp];
      }
      delete[]arlPhc;
    }
    if( arrpCon )
    {
      for( pr=0; pr<nPRk; pr++)
      {
           delete[]arrpCon[pr];
      }
      delete[]arrpCon;
    }
    if( arUmpCon )
    {
      for( j=0; j<NComp; j++)
      {
           delete[]arUmpCon[j];
      }
      delete[]arUmpCon;
    }
    if( arapCon )
    {
         for(j=0; j<nPRk; j++)
         {
             for(s=0; s<nSkr; s++)
             {
                 delete[]arapCon[j][s];
             }
         }
         for(j=0; j<nPRk; j++)
         {
             delete[]arapCon[j];
         }
         delete[]arapCon;
    }
}

// sets the specific surface area of the phase and 'parallel reactions' area fractions
long int TKinMet::UpdateFSA( const double *fSAf_p, const double As )
{

}

// returns modified specific surface area of the phase and 'parallel reactions' area fractions
double TKinMet::ModifiedFSA ( double *fSAf_p )
{

}

// sets new system TP state
long int TKinMet::UpdatePT ( const double T_k, const double P_bar )
{

}

// sets new time and time step
bool TKinMet::UpdateTime( const double Tau, const double dTau )
{

}

// Checks dimensions in order to re-allocate class instance, if necessary
bool TKinMet::testSizes( const KinMetData *kmd )
{

}













//--------------------- End of s_kinmet.cpp ---------------------------

