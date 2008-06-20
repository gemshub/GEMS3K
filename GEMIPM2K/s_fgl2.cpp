//-------------------------------------------------------------------
// $Id: s_fgl2.cpp 1091 2008-06-20 14:28:02Z wagner $
//
// Copyright (c) 2003-2007   S.Churakov, Th.Wagner,
//    D.Kulik, S.Dmitrieva
//
// Implementation of parts of TPRSVcalc and TCFGcalc classes
// called from m_dcomp.cpp
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------

#include <math.h>

#include "s_fgl.h"
#include "m_const.h"

#ifndef IPMGEMPLUGIN
#include "s_tpwork.h"
//--------------------------------------------------------------------//
//
int TPRSVcalc::CalcFugPure( void )
{
    double T, P, Fugcoeff = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
    float *Coeff;
    double Eos2parPT[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 } ;
    int retCode = 0;

    ErrorIf( !aW.twp, "PRSV EoS", "Undefined twp");

    P = aW.twp->P;    /* P in 10^5 Pa? */
    T = aW.twp->TC+273.15;   /* T?in K */

    Coeff = aW.twp->CPg;     /* pointer to coeffs of CG EOS */

// Calling PRSV EoS functions here

    if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1e-5 && P < 1e5 )
       retCode = PRFugacityPT( P, T, Coeff, Eos2parPT, Fugcoeff, Volume,
            DeltaH, DeltaS );
    else {
            Fugcoeff = 1.;
            Volume = 8.31451*T/P;
            aW.twp->V = Volume;
            aW.twp->Fug = Fugcoeff*P;
            return retCode;
          }

    aW.twp->G += 8.31451 * T * log( Fugcoeff );   // from fugacity coeff
    /* add enthalpy and enthropy increments */
    aW.twp->H +=  DeltaH;   // in J/mol - to be completed
    aW.twp->S +=  DeltaS;   // to be completed
    aW.twp->V = Volume; /* /10.  in J/bar */
    aW.twp->Fug = Fugcoeff * P;   /* fugacity at P */

//  passing corrected EoS coeffs to calculation of fluid mixtures (3 added, 31.05.2008 TW)
    aW.twp->wtW[6] = Eos2parPT[0];	// a
    aW.twp->wtW[7] = Eos2parPT[1];  // b
    aW.twp->wtW[8] = Eos2parPT[2];	// sqrAL
    aW.twp->wtW[9] = Eos2parPT[3];	// ac
    aW.twp->wtW[10] = Eos2parPT[4];	// dALdT

    return retCode;
}

#endif

#ifndef IPMGEMPLUGIN
//--------------------------------------------------------------------//
//
int TCGFcalc::CGcalcFug( void )
{
    double T, P, Fugacity = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
    double X[1]={1.};
    float *Coeff, Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
                  Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;
    int retCode = 0;

    ErrorIf( !aW.twp, "CG EoS", "Undefined twp");

    P = aW.twp->P;    /* P in 10^5 Pa? */
    T = aW.twp->TC+273.15;   /* T?in K */

    Coeff = aW.twp->Cemp;     /* pointer to coeffs of CG EOS */

// Calling CG EoS functions here

    if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1e-6 && P < 1e5 )
       retCode = CGFugacityPT( Coeff, Eos4parPT, Fugacity, Volume, P, T );
    else {
            Fugacity = P;
            Volume = 8.31451*T/P;
            aW.twp->V = Volume;
            aW.twp->Fug = Fugacity;
            aW.twp->wtW[6] = Coeff[0];
            if( aW.twp->wtW[6] < 1. || aW.twp->wtW[6] > 10. )
                aW.twp->wtW[6] = 1.;                 // foolproof temporary
            aW.twp->wtW[7] = Coeff[1];
            aW.twp->wtW[8] = Coeff[2];
            aW.twp->wtW[9] = Coeff[3];
            return retCode;
          }

//    if( retCode < 0 )
//    {  //  error - too low pressure
//       Fugacity = P;
//      Volume = 8.31451*T;
//    }

    aW.twp->G += 8.31451 * T * log( Fugacity / P );
    aW.twp->V = Volume /* /10.  in J/bar */;
//    aW.twp->U = ((aW.twp->H/4.184)-RP*fg.VLK*fg.P2)*4.184;
    aW.twp->Fug = Fugacity;   /* fugacity at P */
// For passing corrected EoS coeffs to calculation of fluid
// mixtures
    aW.twp->wtW[6] = Eos4parPT[0];
if( aW.twp->wtW[6] < 1. || aW.twp->wtW[6] > 10. )
  aW.twp->wtW[6] = 1.;                            // foolproof temporary
    aW.twp->wtW[7] = Eos4parPT[1];
    aW.twp->wtW[8] = Eos4parPT[2];
    aW.twp->wtW[9] = Eos4parPT[3];
    
    // add enthalpy and enthropy increments 
    retCode = CGFugacityPT( Coeff, Eos4parPT1, Fugacity, Volume, P, T+T*DELTA );   
    CGEnthalpy( X, Eos4parPT, Eos4parPT1, 1, P, T, DeltaH, DeltaS );
    aW.twp->H +=  DeltaH;   
    aW.twp->S +=  DeltaS;   
//
    return retCode;
}

#endif


// -----------------------------------------------------------------------------
// Implementation of the TSolMod class
// Started by Th.Wagner and D.Kulik on 07.03.2007



// Generic constructor for the TSolMod class
//
TSolMod::TSolMod( int NSpecies, int NParams, int NPcoefs, int MaxOrder,
       int NPperDC, double T_k, double P_bar, char Mod_Code,
       short* arIPx, float* arIPc, float* arDCc,
       double *arWx, double *arlnGam, double dW, double eW, double iS )
{
    R_CONST = 8.31451;
    NComp = NSpecies;
    NPar = NParams;
    NPcoef = NPcoefs;
    MaxOrd = MaxOrder;
    NP_DC = NPperDC;
    Tk = T_k;
    Pbar = P_bar;
    ModCode = Mod_Code;
    RhoW = dW;
    EpsW = eW;
    IonStr = iS;

    // pointer assignment
    aIPx = arIPx;   // Direct access to index list and parameter coeff arrays!
    aIPc = arIPc;
    aDCc = arDCc;
    x = arWx;
    lnGamma = arlnGam;
}



TSolMod::~TSolMod()
{
// In principle, the stuff below is not necessary if the memory is not
// allocated within the class
	aIPx = NULL;
	aIPc = NULL;
	aDCc = NULL;
	x = NULL;
	lnGamma = NULL;
}



// Van Laar model for solid solutions (c) TW March 2007
// Calculates T,P corrected binary interaction parameters
// Returns 0 if Ok or 1 if error
int
TSolMod::VanLaarPT()
{
// calculates P-T dependence of binary interaction parameters
	int ip;
	double Wij[4];

        if ( /* ModCode != SM_VANLAAR || */ NPcoef < 4 || NPar < 1 )
           return 1;  // foolproof!

	for (ip=0; ip<NPar; ip++)
	{
        Wij[0] = (double)aIPc[NPcoef*ip];
        Wij[1] = (double)aIPc[NPcoef*ip+1];
        Wij[2] = (double)aIPc[NPcoef*ip+2];
	    Wij[3] = Wij[0]+ Wij[1]*Tk + Wij[2]*Pbar;
	    aIPc[NPcoef*ip+3] = (float)Wij[3];
	}
	return 0;
}



// Van Laar model for solid solutions (c) TW March 2007
// References:  Holland & Powell (2003)
// Calculates activity coefficients and excess functions
// Returns 0 if Ok or not 0 if error
int
TSolMod::VanLaarMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2;
   double dj, dk;
   double sumPhi; // Sum of Phi terms
   double gEX, vEX, hEX, sEX, cpEX, uEX;
   double *Wu;
   double *Ws;
   double *Wv;
   double *Wpt;   // Interaction coeffs at P-T
   double *Phi;   // Mixing terms
   double *PsVol; // End member volume parameters

   if ( /* ModCode != SM_VANLAAR || */ NPcoef < 4 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   Wu = new double [NPar];
   Ws = new double [NPar];
   Wv = new double [NPar];
   Wpt = new double [NPar];
   Phi = new double [NComp];
   PsVol = new double [NComp];

   if( !Wpt || !Phi || !PsVol )
        return -1;  // memory alloc failure

	// read P-T corrected interaction parameters
   for (ip=0; ip<NPar; ip++)
   {
        Wu[ip] = (double)aIPc[NPcoef*ip];
	Ws[ip] = (double)aIPc[NPcoef*ip+1];
	Wv[ip] = (double)aIPc[NPcoef*ip+2];
	Wpt[ip] = (double)aIPc[NPcoef*ip+3]; // were stored in VanLaarPT()
   }

   // calculating Phi values
   sumPhi = 0.;
   for (j=0; j<NComp; j++)
   {
       PsVol[j] = (double)aDCc[NP_DC*j];  // reading pseudo-volumes
       sumPhi +=  x[j]*PsVol[j];
   }
   if( fabs(sumPhi) < 1e-30 )
       return 2;    // to prevent zerodivide!

   for (j=0; j<NComp; j++)
       Phi[j] = x[j]*PsVol[j]/sumPhi;

   // calculate activity coefficients
   for (j=0; j<NComp; j++)      // index end members with j
   {
	lnGamRT = 0.;
	for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
	{
        index1 = (int)aIPx[MaxOrd*ip];
	    index2 = (int)aIPx[MaxOrd*ip+1];

   	    if( j == index1 )
		dj = 1.;
	    else
		dj = 0.;
	    if( j == index2 )
		dk = 1.;
	    else
		dk = 0.;
	    lnGamRT -= (dj-Phi[index1])*(dk-Phi[index2])*Wpt[ip]
                         *2.*PsVol[j]/(PsVol[index1]+PsVol[index2]);
	}
        lnGam = lnGamRT/(R_CONST*Tk);
//	Gam = exp(lnGam);
	lnGamma[j] = lnGam;
	}

   // calculate bulk phase excess properties
   gEX = 0.0;
   vEX = 0.0;
   hEX = 0.0;
   sEX = 0.0;
   cpEX = 0.0;
   uEX = 0.0;

   for (ip=0; ip<NPar; ip++)
   {
      index1 = (int)aIPx[MaxOrd*ip];
      index2 = (int)aIPx[MaxOrd*ip+1];
      gEX += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wpt[ip];
      vEX += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wv[ip];
      uEX += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wu[ip];
      sEX -= Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Ws[ip];
   }
   
   hEX = uEX+vEX*Pbar;
   Gex_ = gEX;
   Vex_ = vEX;
   Hex_ = hEX;
   Sex_ = sEX;
   CPex_ = cpEX;

   delete[]Wu;
   delete[]Ws;
   delete[]Wv;
   delete[]Wpt;
   delete[]Phi;
   delete[]PsVol;
   return 0;
}



// Regular model for multicomponent solid solutions (c) TW March 2007
// Calculates T,P corrected binary interaction parameters
// Returns 0 if Ok or 1 if error
int
TSolMod::RegularPT()
{
// calculates P-T dependence of binary interaction parameters
	int ip;
	double Wij[4];

        if ( /* ModCode != SM_REGULAR || */ NPcoef < 4 || NPar < 1 )
           return 1;  // foolproof!

	for (ip=0; ip<NPar; ip++)
	{
        Wij[0] = (double)aIPc[NPcoef*ip];
        Wij[1] = (double)aIPc[NPcoef*ip+1];
        Wij[2] = (double)aIPc[NPcoef*ip+2];
	    Wij[3] = Wij[0]+ Wij[1]*Tk + Wij[2]*Pbar;
	    aIPc[NPcoef*ip+3] = (float)Wij[3];
	}
	return 0;
}



// Regular model for multicomponent solid solutions (c) TW March 2007
// References:  Holland & Powell (1993)
// Calculates activity coefficients and excess functions
// Returns 0 if Ok or not 0 if error
int
TSolMod::RegularMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2;
   double dj, dk;
   double gEX, vEX, hEX, sEX, cpEX, uEX;
   double *Wu;
   double *Ws;
   double *Wv;
   double *Wpt;   // Interaction coeffs at P-T

   if ( /* ModCode != SM_REGULAR || */ NPcoef < 4 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   Wu = new double [NPar];
   Ws = new double [NPar];
   Wv = new double [NPar];
   Wpt = new double [NPar];

   if( !Wpt || !Wu || !Ws || !Wv )
        return -1;  // memory alloc failure

	// read P-T corrected interaction parameters
   for (ip=0; ip<NPar; ip++)
   {
    Wu[ip] = (double)aIPc[NPcoef*ip];
	Ws[ip] = (double)aIPc[NPcoef*ip+1];
	Wv[ip] = (double)aIPc[NPcoef*ip+2];
	Wpt[ip] = (double)aIPc[NPcoef*ip+3]; // were stored in RegularPT()
   }

   // calculate activity coefficients
   for (j=0; j<NComp; j++)      // index end members with j
   {
	lnGamRT = 0.;
	for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
	{
        index1 = (int)aIPx[MaxOrd*ip];
	    index2 = (int)aIPx[MaxOrd*ip+1];

   	    if( j == index1 )
		dj = 1.;
	    else
		dj = 0.;
	    if( j == index2 )
		dk = 1.;
	    else
		dk = 0.;
	    lnGamRT -= (dj-x[index1])*(dk-x[index2])*Wpt[ip];
	}
        lnGam = lnGamRT/(R_CONST*Tk);
//	Gam = exp(lnGam);
	lnGamma[j] = lnGam;
	}

   // calculate bulk phase excess properties
   gEX = 0.0;
   vEX = 0.0;
   hEX = 0.0;
   sEX = 0.0;
   cpEX = 0.0;
   uEX = 0.0;

   for (ip=0; ip<NPar; ip++)
   {
      index1 = (int)aIPx[MaxOrd*ip];
      index2 = (int)aIPx[MaxOrd*ip+1];
      gEX += x[index1]*x[index2]*Wpt[ip];
      vEX += x[index1]*x[index2]*Wv[ip];
      uEX += x[index1]*x[index2]*Wu[ip];
      sEX -= x[index1]*x[index2]*Ws[ip];
   }
   
   hEX = uEX+vEX*Pbar;
   Gex_ = gEX;
   Vex_ = vEX;
   Hex_ = hEX;
   Sex_ = sEX;
   CPex_ = cpEX;
   
   delete[]Wu;
   delete[]Ws;
   delete[]Wv;
   delete[]Wpt;
   return 0;
}



// Redlich-Kister model for multicomponent solid solutions (c) TW March 2007
// Calculates T,P corrected binary interaction parameters (4 per interaction)
// Returns 0 if Ok or 1 if error
int
TSolMod::RedlichKisterPT()
{
// calculates P-T dependence of binary interaction parameters
	int ip;
	double L0[5];
	double L1[5];
	double L2[5];
	double L3[5];

        if ( /* ModCode != SM_GUGGENM || */ NPcoef < 20 || NPar < 1 )
           return 1;  // foolproof!

	for (ip=0; ip<NPar; ip++)
	{
        L0[0] = (double)aIPc[NPcoef*ip];
        L0[1] = (double)aIPc[NPcoef*ip+1];
        L0[2] = (double)aIPc[NPcoef*ip+2];
        L0[3] = (double)aIPc[NPcoef*ip+3];
        L0[4] = L0[0] + L0[1]*Tk + L0[2]*Tk*log(Tk) + L0[3]*Pbar;
        L1[0] = (double)aIPc[NPcoef*ip+4];
        L1[1] = (double)aIPc[NPcoef*ip+5];
        L1[2] = (double)aIPc[NPcoef*ip+6];
        L1[3] = (double)aIPc[NPcoef*ip+7];
        L1[4] = L1[0] + L1[1]*Tk + L1[2]*Tk*log(Tk) + L1[3]*Pbar;
        L2[0] = (double)aIPc[NPcoef*ip+8];
        L2[1] = (double)aIPc[NPcoef*ip+9];
        L2[2] = (double)aIPc[NPcoef*ip+10];
        L2[3] = (double)aIPc[NPcoef*ip+11];
        L2[4] = L2[0] + L2[1]*Tk + L2[2]*Tk*log(Tk) + L2[3]*Pbar;
        L3[0] = (double)aIPc[NPcoef*ip+12];
        L3[1] = (double)aIPc[NPcoef*ip+13];
        L3[2] = (double)aIPc[NPcoef*ip+14];
        L3[3] = (double)aIPc[NPcoef*ip+15];
        L3[4] = L3[0] + L3[1]*Tk + L3[2]*Tk*log(Tk) + L3[3]*Pbar;

	    aIPc[NPcoef*ip+16] = (float)L0[4];
	    aIPc[NPcoef*ip+17] = (float)L1[4];
	    aIPc[NPcoef*ip+18] = (float)L2[4];
	    aIPc[NPcoef*ip+19] = (float)L3[4];
	}
	return 0;
}



// Redlich-Kister model for multicomponent solid solutions (c) TW March 2007
// References: Hillert (1998)
// Calculates activity coefficients and excess functions
// Returns 0 if Ok or not 0 if error
int
TSolMod::RedlichKisterMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2, L, I, J;
   double LU, LS, LCP, LV, LPT;
   double L0, L1, L2, L3;
   double gEX, vEX, hEX, sEX, cpEX, uEX;
   
   double **Lu;
   double **Ls;
   double **Lcp;
   double **Lv;
   double **Lpt;

   if ( /* ModCode != SM_REGULAR || */ NPcoef < 20 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   Lu = new double *[NPar];
   Ls = new double *[NPar];
   Lcp = new double *[NPar];
   Lv = new double *[NPar];
   Lpt = new double *[NPar];
   
   for (ip=0; ip<NPar; ip++)
   {
	   Lu[ip] = new double [4];
	   Ls[ip] = new double [4];
	   Lcp[ip] = new double [4];
	   Lv[ip] = new double [4];
	   Lpt[ip] = new double [4];
   }

   if( !Lu || !Ls || !Lcp || !Lv || !Lpt )
        return -1;  // memory alloc failure

	// read in interaction parameters
   	for (ip=0; ip<NPar; ip++)
   	{
	   	Lu[ip][0] = (double)aIPc[NPcoef*ip+0];
	   	Ls[ip][0] = (double)aIPc[NPcoef*ip+1];
	   	Lcp[ip][0] = (double)aIPc[NPcoef*ip+2];
	   	Lv[ip][0] = (double)aIPc[NPcoef*ip+3];
	   	Lpt[ip][0] = (double)aIPc[NPcoef*ip+16];
	   	Lu[ip][1] = (double)aIPc[NPcoef*ip+4];
	   	Ls[ip][1] = (double)aIPc[NPcoef*ip+5];
	   	Lcp[ip][1] = (double)aIPc[NPcoef*ip+6];
	   	Lv[ip][1] = (double)aIPc[NPcoef*ip+7];
	   	Lpt[ip][1] = (double)aIPc[NPcoef*ip+17];
	   	Lu[ip][2] = (double)aIPc[NPcoef*ip+8];
	   	Ls[ip][2] = (double)aIPc[NPcoef*ip+9];
	   	Lcp[ip][2] = (double)aIPc[NPcoef*ip+10];
	   	Lv[ip][2] = (double)aIPc[NPcoef*ip+11];
	   	Lpt[ip][2] = (double)aIPc[NPcoef*ip+18];
	   	Lu[ip][3] = (double)aIPc[NPcoef*ip+12];
	   	Ls[ip][3] = (double)aIPc[NPcoef*ip+13];
	   	Lcp[ip][3] = (double)aIPc[NPcoef*ip+14];
	   	Lv[ip][3] = (double)aIPc[NPcoef*ip+15];
	   	Lpt[ip][3] = (double)aIPc[NPcoef*ip+19];
	}

	// calculate activity coefficients
	for (j=0; j<NComp; j++)      // index end members with j
	{
		lnGamRT = 0.;
		for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
		{
			index1 = (int)aIPx[MaxOrd*ip];
			index2 = (int)aIPx[MaxOrd*ip+1];

			if ( j == index1 || j == index2) // interaction terms with j
			{
				if ( index1 == j ) // check order of idexes
				{
					L = index1;
					I = index2;
					L0 = Lpt[ip][0];
					L1 = Lpt[ip][1];
					L2 = Lpt[ip][2];
					L3 = Lpt[ip][3];
				}
				else
				{
					L = index2;
					I = index1;
					L0 = Lpt[ip][0];
					L1 = -Lpt[ip][1];
					L2 = Lpt[ip][2];
					L3 = -Lpt[ip][3];
				}

				lnGamRT += L0*x[I]*(1.-x[L])
					+ L1*x[I]*(2.*(1.-x[L])*(x[L]-x[I])+x[I])
					+ L2*x[I]*(x[L]-x[I])*(3.*(1.-x[L])*(x[L]-x[I])+2.*x[I])
					+ L3*x[I]*pow((x[L]-x[I]),2.)*(4.*(1.-x[L])*(x[L]-x[I])+3.*x[I]);
			}

			else // interaction terms without j
			{
				I = index1;
				J = index2;
				L0 = Lpt[ip][0];
				L1 = Lpt[ip][1];
				L2 = Lpt[ip][2];
				L3 = Lpt[ip][3];

				lnGamRT -= x[I]*x[J]*( L0 + L1*2.*(x[I]-x[J])
					+ L2*3.*pow((x[I]-x[J]),2.)
					+ L3*4.*pow((x[I]-x[J]),3.) );
			}
		}

		lnGam = lnGamRT/(R_CONST*Tk);
		//	Gam = exp(lnGam);
		lnGamma[j] = lnGam;
	}

   	// calculate bulk phase excess properties
   	gEX = 0.0;
   	vEX = 0.0;
   	hEX = 0.0;
   	sEX = 0.0;
   	cpEX = 0.0;
   	uEX = 0.0;

   	for (ip=0; ip<NPar; ip++)
   	{
   	   	index1 = (int)aIPx[MaxOrd*ip];
   	   	index2 = (int)aIPx[MaxOrd*ip+1];

      	LPT = Lpt[ip][0] + Lpt[ip][1]*(x[index1]-x[index2])
      			+ Lpt[ip][2]*pow((x[index1]-x[index2]),2.)
      			+ Lpt[ip][3]*pow((x[index1]-x[index2]),3.);

      	LV = Lv[ip][0] + Lv[ip][1]*(x[index1]-x[index2])
      			+ Lv[ip][2]*pow((x[index1]-x[index2]),2.)
      			+ Lv[ip][3]*pow((x[index1]-x[index2]),3.);

   	   	LU = (Lu[ip][0]-Lcp[ip][0]*Tk)
   	  			+ (Lu[ip][1]-Lcp[ip][1]*Tk)*(x[index1]-x[index2])
      			+ (Lu[ip][2]-Lcp[ip][2]*Tk)*pow((x[index1]-x[index2]),2.)
      			+ (Lu[ip][3]-Lcp[ip][3]*Tk)*pow((x[index1]-x[index2]),3.);

   	   	LS = (-Ls[ip][0]-Lcp[ip][0]*(1.+log(Tk)))
      			+ (-Ls[ip][1]-Lcp[ip][1]*(1.+log(Tk)))*(x[index1]-x[index2])
      			+ (-Ls[ip][2]-Lcp[ip][2]*(1.+log(Tk)))*pow((x[index1]-x[index2]),2.)
      			+ (-Ls[ip][3]-Lcp[ip][3]*(1.+log(Tk)))*pow((x[index1]-x[index2]),3.);

   	   	LCP = (-Lcp[ip][0]) + (-Lcp[ip][1])*(x[index1]-x[index2])
      			+ (-Lcp[ip][2])*pow((x[index1]-x[index2]),2.)
      			+ (-Lcp[ip][3])*pow((x[index1]-x[index2]),3.);

      	gEX += x[index1]*x[index2]*LPT;
      	vEX += x[index1]*x[index2]*LV;
      	uEX += x[index1]*x[index2]*LU;
      	sEX += x[index1]*x[index2]*LS;
      	cpEX += x[index1]*x[index2]*LCP;
  	}

   	hEX = uEX+vEX*Pbar;
   	Gex_ = gEX;
   	Vex_ = vEX;
   	Hex_ = hEX;
   	Sex_ = sEX;
   	CPex_ = cpEX;
   	
   	for (ip=0; ip<NPar; ip++)
   	{
   		delete[]Lu[ip];
   		delete[]Ls[ip];
   		delete[]Lcp[ip];
   		delete[]Lv[ip];
   		delete[]Lpt[ip];
   	}
   	delete[]Lu;
   	delete[]Ls;
   	delete[]Lcp;
   	delete[]Lv;
   	delete[]Lpt;
 
   	return 0;
}



// NRTL model for liquid solutions (c) TW June 2008
// Calculates T-corrected interaction parameters
// Returns 0 if OK or 1 if error
int
TSolMod::NRTL_PT()
{
	// calculates T-dependence of binary interaction parameters
	int ip;
	double A, B, C, D, E, F;
	double tau, dtau, d2tau, alp, dalp, d2alp;

        if ( /* ModCode != SM_NRTL || */ NPcoef < 12 || NPar < 1 )
           return 1;  // foolproof!

	for (ip=0; ip<NPar; ip++)
	{
		A = (double)aIPc[NPcoef*ip+0];
		B = (double)aIPc[NPcoef*ip+1];
		C = (double)aIPc[NPcoef*ip+2];
		D = (double)aIPc[NPcoef*ip+3];
		E = (double)aIPc[NPcoef*ip+4];
		F = (double)aIPc[NPcoef*ip+5];
		tau = A + B/Tk + C*Tk + D*log(Tk);	// partial derivatives of tau and alp
		dtau = - B/pow(Tk,2.) + C + D/Tk;
		d2tau = 2.*B/pow(Tk,3.) - D/pow(Tk,2.);
		alp = E + F*(Tk-273.15);
		dalp = F;
		d2alp = 0.0;
		aIPc[NPcoef*ip+6] = (float)tau;
		aIPc[NPcoef*ip+7] = (float)dtau;
		aIPc[NPcoef*ip+8] = (float)d2tau;
		aIPc[NPcoef*ip+9] = (float)alp;
		aIPc[NPcoef*ip+10] = (float)dalp;
		aIPc[NPcoef*ip+11] = (float)d2alp;
	}
	return 0;
}



// NRTL model for liquid solutions (c) TW June 2008
// References: Renon and Prausnitz (1968), Prausnitz et al. (1997)
// Calculates activity coefficients and excess functions
// heat capacity calculation added, 06.06.2008 (TW)
// Returns 0 if OK or 1 if error
int
TSolMod::NRTL_MixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
	int ip, j, i, k;
	int i1, i2;
	double K, L, M, N, O;
	double U, dU, V, dV, d2U, d2V;
	double g, dg, d2g, lnGam;
	double gEX, vEX, hEX, sEX, cpEX;
	double **Tau;
	double **dTau;
	double **d2Tau;
	double **Alp;
	double **dAlp;
	double **d2Alp;
	double **G;
	double **dG;
	double **d2G;

	if ( NPcoef < 12 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;  // foolproof!

	Tau = new double *[NComp];
	dTau = new double *[NComp];
	d2Tau = new double *[NComp];
	Alp = new double *[NComp];
	dAlp = new double *[NComp];
	d2Alp = new double *[NComp];
	G = new double *[NComp];
	dG = new double *[NComp];
	d2G = new double *[NComp];

    for (j=0; j<NComp; j++)
    {
    	Tau[j] = new double [NComp];
    	dTau[j] = new double [NComp];
    	d2Tau[j] = new double [NComp];
    	Alp[j] = new double [NComp];
    	dAlp[j] = new double [NComp];
    	d2Alp[j] = new double [NComp];
		G[j] = new double [NComp];
		dG[j] = new double [NComp];
		d2G[j] = new double [NComp];
	}

	// fill internal arrays of interaction parameters with standard value
	for (j=0; j<NComp; j++)
	{
		for ( i=0; i<NComp; i++ )
		{
			Tau[j][i] = 0.0;
			dTau[j][i] = 0.0;
			d2Tau[j][i] = 0.0;
			Alp[j][i] = 0.0;
			dAlp[j][i] = 0.0;
			d2Alp[j][i] = 0.0;
			G[j][i] = 1.0;
			dG[j][i] = 0.0;
			d2G[j][i] = 0.0;
		}
	}

	// read and convert parameters that have non-standard value
	for (ip=0; ip<NPar; ip++)
	{
		i1 = (int)aIPx[MaxOrd*ip];
		i2 = (int)aIPx[MaxOrd*ip+1];
		Tau[i1][i2] = (double)aIPc[NPcoef*ip+6];
		dTau[i1][i2] = (double)aIPc[NPcoef*ip+7];
		d2Tau[i1][i2] = (double)aIPc[NPcoef*ip+8];
		Alp[i1][i2] = (double)aIPc[NPcoef*ip+9];
		dAlp[i1][i2] = (double)aIPc[NPcoef*ip+10];
		d2Alp[i1][i2] = (double)aIPc[NPcoef*ip+11];
		
		G[i1][i2] = exp(-Alp[i1][i2]*Tau[i1][i2]);
		dG[i1][i2] = - ( dAlp[i1][i2]*Tau[i1][i2] + Alp[i1][i2]*dTau[i1][i2] )
				* exp(-Alp[i1][i2]*Tau[i1][i2]);
		d2G[i1][i2] = - ( (d2Alp[i1][i2]*Tau[i1][i2] + 2.*dAlp[i1][i2]*dTau[i1][i2]
				+ Alp[i1][i2]*d2Tau[i1][i2])*G[i1][i2]
				+ (dAlp[i1][i2]*Tau[i1][i2] + Alp[i1][i2]*dTau[i1][i2])*dG[i1][i2] );
		
		// old version with constant Alp
		// dG[i1][i2] = -Alp[i1][i2] * exp( -Alp[i1][i2]*Tau[i1][i2] ) * dTau[i1][i2];
		// d2G[i1][i2] = -Alp[i1][i2]*(-Alp[i1][i2]*exp(-Alp[i1][i2]*Tau[i1][i2])*dTau[i1][i2]*dTau[i1][i2]
		//		+ exp(-Alp[i1][i2]*Tau[i1][i2])*d2Tau[i1][i2]);		                                                                                                                                  
	}

	// calculate activity coefficients
	for (j=0; j<NComp; j++)
	{
		lnGam = 0.0;
		K = 0.0;
		L = 0.0;
		M = 0.0;
		for (i=0; i<NComp; i++)
		{
			N = 0.0;
			O = 0.0;
			K += ( x[i]*Tau[i][j]*G[i][j] );
			L += ( x[i]*G[i][j] );
			for (k=0; k<NComp; k++)
			{
				N += ( x[k]*G[k][i] );
				O += ( x[k]*Tau[k][i]*G[k][i] );
			}
			M += ( x[i]*G[j][i]/N * ( Tau[j][i] - O/N ) );
		}
		lnGam = K/L + M;
		lnGamma[j] = lnGam;
	}

	// calculate bulk phase excess properties
   	gEX = 0.0;
   	vEX = 0.0;
   	hEX = 0.0;
   	sEX = 0.0;
   	cpEX = 0.0;
   	g = 0.0;
   	dg = 0.0;
   	d2g = 0.0;

   	for (j=0; j<NComp; j++)
   	{
		U = 0.0;
		V = 0.0;
		dU = 0.0;
		dV = 0.0;
		d2U = 0.0;
		d2V = 0.0;
		for (i=0; i<NComp; i++)
		{
			U += x[i]*Tau[i][j]*G[i][j];
			V += x[i]*G[i][j];
			dU += x[i] * ( dTau[i][j]*G[i][j] + Tau[i][j]*dG[i][j] );
			dV += x[i]*dG[i][j];
			d2U += x[i] * ( d2Tau[i][j]*G[i][j] + 2.*dTau[i][j]*dG[i][j]
					+ Tau[i][j]*d2G[i][j] );
			d2V += x[i]*d2G[i][j];
		}
		g += x[j]*U/V;
		dg += x[j] * (dU*V-U*dV)/pow(V,2.);
		d2g += x[j] * ( (d2U*V+dU*dV)*pow(V,2.)/pow(V,4.) - (dU*V)*(2.*V*dV)/pow(V,4.)
				- (dU*dV+U*d2V)*pow(V,2.)/pow(V,4.) + (U*dV)*(2.*V*dV)/pow(V,4.) );
	}

	gEX = g*R_CONST*Tk;
	hEX = -R_CONST*pow(Tk,2.)*dg;
	sEX = (hEX-gEX)/Tk;
	cpEX = -R_CONST * ( 2.*Tk*dg + pow(Tk,2.)*d2g );

	Gex_ = gEX;
	Vex_ = vEX;
	Hex_ = hEX;
	Sex_ = sEX;
	CPex_ = cpEX;

   	// cleaning memory
   	for (j=0; j<NComp; j++)
   	{
   		delete[]Tau[j];
   		delete[]dTau[j];
   		delete[]d2Tau[j];
   		delete[]Alp[j];
   		delete[]dAlp[j];
   		delete[]d2Alp[j];
		delete[]G[j];
		delete[]dG[j];
		delete[]d2G[j];
	}
   	delete[]Tau;
   	delete[]dTau;
   	delete[]d2Tau;
   	delete[]Alp;
   	delete[]dAlp;
   	delete[]d2Alp;
	delete[]G;
	delete[]dG;
	delete[]d2G;
	return 0;
}



// Wilson model for liquid solutions (c) TW June 2008
// Calculates T-corrected interaction parameters
// Returns 0 if OK or 1 if error
int
TSolMod::Wilson_PT()
{
	// calculates T-dependence of binary interaction parameters
	int ip;
	double A, B, C, D;
	double lam, dlam, d2lam;

        if ( /* ModCode != SM_NRTL || */ NPcoef < 7 || NPar < 1 )
           return 1;  // foolproof!

	for (ip=0; ip<NPar; ip++)
	{
		A = (double)aIPc[NPcoef*ip+0];
		B = (double)aIPc[NPcoef*ip+1];
		C = (double)aIPc[NPcoef*ip+2];
		D = (double)aIPc[NPcoef*ip+3];
		lam = exp( A + B/Tk + C*Tk + D*log(Tk) );
		dlam = lam*( - B/pow(Tk,2.) + C + D/Tk );
		d2lam = dlam*( - B/pow(Tk,2.) + C + D/Tk ) + lam*( 2.*B/pow(Tk,3.) - D/pow(Tk,2.) );
		aIPc[NPcoef*ip+4] = (float)lam;
		aIPc[NPcoef*ip+5] = (float)dlam;
		aIPc[NPcoef*ip+6] = (float)d2lam;
	}
	return 0;
}



// Wilson model for liquid solutions (c) TW June 2008
// References: Prausnitz et al. (1997)
// Calculates activity coefficients and excess functions
// heat capacity calculation added, 06.06.2008 (TW)
// Returns 0 if OK or 1 if error
int
TSolMod::Wilson_MixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
	int ip, j, i, k;
	int i1, i2;
	double K, L, M;
	double U, dU, d2U;
	double g, dg, d2g, lnGam;
	double gEX, vEX, hEX, sEX, cpEX;
	double **Lam;
	double **dLam;
	double **d2Lam;

	if ( NPcoef < 7 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;  // foolproof!

	Lam = new double *[NComp];
	dLam = new double *[NComp];
	d2Lam = new double *[NComp];

    for (j=0; j<NComp; j++)
    {
		Lam[j] = new double [NComp];
		dLam[j] = new double [NComp];
		d2Lam[j] = new double [NComp];
	}

	// fill internal arrays of interaction parameters with standard value
	for (j=0; j<NComp; j++)
	{
		for ( i=0; i<NComp; i++ )
		{
			Lam[j][i] = 1.0;
			dLam[j][i] = 0.0;
			d2Lam[j][i] = 0.0;
		}
	}

	// read and convert parameters that have non-standard value
	for (ip=0; ip<NPar; ip++)
	{
		i1 = (int)aIPx[MaxOrd*ip];
		i2 = (int)aIPx[MaxOrd*ip+1];
		Lam[i1][i2] = (double)aIPc[NPcoef*ip+4];
		dLam[i1][i2] = (double)aIPc[NPcoef*ip+5];
		d2Lam[i1][i2] = (double)aIPc[NPcoef*ip+6];
	}

	// calculate activity coefficients (Wilson)
	for (j=0; j<NComp; j++)
	{
		lnGam = 0.0;
		K = 0.0;
		L = 0.0;
		for (i=0; i<NComp; i++)
		{
			M = 0.0;
			K += x[i]*Lam[j][i];
			for (k=0; k<NComp; k++)
			{
				M += x[k]*Lam[i][k];
			}
			L += x[i]*Lam[i][j]/M;
		}
		lnGam = 1.-log(K)-L;
		lnGamma[j] = lnGam;
	}

	// calculate bulk phase excess properties
	gEX = 0.0;
	vEX = 0.0;
	hEX = 0.0;
	sEX = 0.0;
	cpEX = 0.0;
	g = 0.0;
	dg = 0.0;
	d2g = 0.0;

	for (j=0; j<NComp; j++)
	{
		U = 0.0;
		dU = 0.0;
		d2U = 0.0;
		for (i=0; i<NComp; i++)
		{
			U += x[i]*Lam[j][i];
			dU += x[i]*dLam[j][i];
			d2U += x[i]*d2Lam[j][i];
		}
		g -= x[j]*log(U);
		dg -= x[j]*(1./U)*dU;
		d2g -= x[j] * ( (-1./pow(U,2.))*dU*dU + (1./U)*d2U );  // fixed, 11.06.2008 (TW)
	}

	// final calculations and assignments
	gEX = g*R_CONST*Tk;
	hEX = -R_CONST*pow(Tk,2.)*dg;
	sEX = (hEX-gEX)/Tk;
	cpEX = -R_CONST * ( 2.*Tk*dg + pow(Tk,2.)*d2g );

	Gex_ = gEX;
	Vex_ = vEX;
	Hex_ = hEX;
	Sex_ = sEX;
	CPex_ = cpEX;

   	// cleaning memory
   	for (j=0; j<NComp; j++)
   	{
		delete[]Lam[j];
		delete[]dLam[j];
		delete[]d2Lam[j];
	}
	delete[]Lam;
	delete[]dLam;
	delete[]d2Lam;
	return 0;
}



// add other solution models here





//--------------------- End of s_fgl2.cpp ---------------------------
