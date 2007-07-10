//-------------------------------------------------------------------
// $Id: s_fgl2.cpp 901 2007-04-03 12:41:09Z gems $
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
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
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

//  passing corrected EoS coeffs to calculation of fluid mixtures
    aW.twp->wtW[6] = Eos2parPT[0];      // a
    aW.twp->wtW[7] = Eos2parPT[1];      // b
// three more to add !!!
    return retCode;
}

#endif

#ifndef IPMGEMPLUGIN
//--------------------------------------------------------------------//
//
int TCGFcalc::CGcalcFug( void )
{
    double T, P, Fugacity = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
    float *Coeff, Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 } ;
    int retCode = 0;

    ErrorIf( !aW.twp, "CG EoS", "Undefined twp");

    P = aW.twp->P;    /* P in 10^5 Pa? */
    T = aW.twp->TC+273.15;   /* T?in K */

    Coeff = aW.twp->Cemp;     /* pointer to coeffs of CG EOS */

// Calling CG EoS functions here

    if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1. && P < 1e5 )
       retCode = CGFugacityPT( Coeff, Eos4parPT, Fugacity, Volume,
            DeltaH, DeltaS, P, T );
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
    /* add enthalpy and enthropy increments */
    aW.twp->H +=  DeltaH;   // in J/mol - to be completed
    aW.twp->S +=  DeltaS;   // to be completed
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
       double *arWx, double *arlnGam )
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
//
int
TSolMod::VanLaarMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2;
   double dj, dk;
   double sumPhi; // Sum of Phi terms
   double *Wh;
   double *Ws;
   double *Wv;
   double *Wpt;   // Interaction coeffs at P-T
   double *Phi;   // Mixing terms
   double *PsVol; // End member volume parameters

   if ( /* ModCode != SM_VANLAAR || */ NPcoef < 4 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   Wh = new double [NPar];
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
        Wh[ip] = (double)aIPc[NPcoef*ip];
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
   Gex = 0.;
   Vex = 0.;
   Hex = 0.;
   Sex = 0.;
   CPex = 0.;

   for (ip=0; ip<NPar; ip++)
   {
      index1 = (int)aIPx[MaxOrd*ip];
      index2 = (int)aIPx[MaxOrd*ip+1];
      Gex += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wpt[ip];
      Vex += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wv[ip];
      Hex += Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wh[ip];
      Sex -= Phi[index1]*Phi[index2]*2.*sumPhi/(PsVol[index1]+PsVol[index2])*Wv[ip];
   }

   Gex_ = Gex;
   Vex_ = Vex;
   Hex_ = Hex;
   Sex_ = Sex;
   CPex_ = CPex;

   delete[]Wh;
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
//
int
TSolMod::RegularMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2;
   double dj, dk;
   double *Wh;
   double *Ws;
   double *Wv;
   double *Wpt;   // Interaction coeffs at P-T

   if ( /* ModCode != SM_REGULAR || */ NPcoef < 4 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   Wh = new double [NPar];
   Ws = new double [NPar];
   Wv = new double [NPar];
   Wpt = new double [NPar];

   if( !Wpt || !Wh || !Ws || !Wv )
        return -1;  // memory alloc failure

	// read P-T corrected interaction parameters
   for (ip=0; ip<NPar; ip++)
   {
    Wh[ip] = (double)aIPc[NPcoef*ip];
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
   Gex = 0.;
   Vex = 0.;
   Hex = 0.;
   Sex = 0.;
   CPex = 0.;

   for (ip=0; ip<NPar; ip++)
   {
      index1 = (int)aIPx[MaxOrd*ip];
      index2 = (int)aIPx[MaxOrd*ip+1];
      Gex += x[index1]*x[index2]*Wpt[ip];
      Vex += x[index1]*x[index2]*Wv[ip];
      Hex += x[index1]*x[index2]*Wh[ip];
      Sex -= x[index1]*x[index2]*Wv[ip];
   }

   Gex_ = Gex;
   Vex_ = Vex;
   Hex_ = Hex;
   Sex_ = Sex;
   CPex_ = CPex;

   delete[]Wh;
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
// References:  Hillert (1998)
// Calculates activity coefficients and excess functions
// Returns 0 if Ok or not 0 if error
//
int
TSolMod::RedlichKisterMixMod( double &Gex_, double &Vex_, double &Hex_, double &Sex_,
         double &CPex_ )
{
   int ip, j;
   int index1, index2, L, I, J;
   double Lh, Ls, Lcp, Lv, Lpt;
   double L0, L1, L2, L3;
   // double dj, dk;

   double *L0h;
   double *L0s;
   double *L0cp;
   double *L0v;
   double *L0pt;

   double *L1h;
   double *L1s;
   double *L1cp;
   double *L1v;
   double *L1pt;

   double *L2h;
   double *L2s;
   double *L2cp;
   double *L2v;
   double *L2pt;

   double *L3h;
   double *L3s;
   double *L3cp;
   double *L3v;
   double *L3pt;

   if ( /* ModCode != SM_REGULAR || */ NPcoef < 20 || NPar < 1 || NComp < 2
         || MaxOrd < 2 || !x || !lnGamma )
           return 1;  // foolproof!

   L0h = new double [NPar];
   L0s = new double [NPar];
   L0cp = new double [NPar];
   L0v = new double [NPar];
   L0pt = new double [NPar];

   L1h = new double [NPar];
   L1s = new double [NPar];
   L1cp = new double [NPar];
   L1v = new double [NPar];
   L1pt = new double [NPar];

   L2h = new double [NPar];
   L2s = new double [NPar];
   L2cp = new double [NPar];
   L2v = new double [NPar];
   L2pt = new double [NPar];

   L3h = new double [NPar];
   L3s = new double [NPar];
   L3cp = new double [NPar];
   L3v = new double [NPar];
   L3pt = new double [NPar];

   if( !L0h || !L0s || !L0cp || !L0v || !L0pt || !L1h || !L1s || !L1cp
   		|| !L1v || !L1pt || !L2h || !L2s || !L2cp || !L2v || !L2pt
   			|| !L3h || !L3s || !L3cp || !L3v || !L3pt )
        return -1;  // memory alloc failure

	// read in interaction parameters
   	for (ip=0; ip<NPar; ip++)
   	{
	   	L0h[ip] = (double)aIPc[NPcoef*ip+0];
	   	L0s[ip] = (double)aIPc[NPcoef*ip+1];
	   	L0cp[ip] = (double)aIPc[NPcoef*ip+2];
	   	L0v[ip] = (double)aIPc[NPcoef*ip+3];
	   	L0pt[ip] = (double)aIPc[NPcoef*ip+16];

	   	L1h[ip] = (double)aIPc[NPcoef*ip+4];
	   	L1s[ip] = (double)aIPc[NPcoef*ip+5];
	   	L1cp[ip] = (double)aIPc[NPcoef*ip+6];
	   	L1v[ip] = (double)aIPc[NPcoef*ip+7];
	   	L1pt[ip] = (double)aIPc[NPcoef*ip+17];

	   	L2h[ip] = (double)aIPc[NPcoef*ip+8];
	   	L2s[ip] = (double)aIPc[NPcoef*ip+9];
	   	L2cp[ip] = (double)aIPc[NPcoef*ip+10];
	   	L2v[ip] = (double)aIPc[NPcoef*ip+11];
	   	L2pt[ip] = (double)aIPc[NPcoef*ip+18];

	   	L3h[ip] = (double)aIPc[NPcoef*ip+12];
	   	L3s[ip] = (double)aIPc[NPcoef*ip+13];
	   	L3cp[ip] = (double)aIPc[NPcoef*ip+14];
	   	L3v[ip] = (double)aIPc[NPcoef*ip+15];
	   	L3pt[ip] = (double)aIPc[NPcoef*ip+19];
	}


	// calculate activity coefficients (under construction)
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
					L0 = L0pt[ip];
					L1 = L1pt[ip];
					L2 = L2pt[ip];
					L3 = L3pt[ip];
				}

				else
				{
					L = index2;
					I = index1;
					L0 = L0pt[ip];
					L1 = -L1pt[ip];
					L2 = L2pt[ip];
					L3 = -L3pt[ip];
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
				L0 = L0pt[ip];
				L1 = L1pt[ip];
				L2 = L2pt[ip];
				L3 = L3pt[ip];

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
   	Gex = 0.;
   	Vex = 0.;
   	Hex = 0.;
   	Sex = 0.;
   	CPex = 0.;

   	for (ip=0; ip<NPar; ip++)
   	{
   	   	index1 = (int)aIPx[MaxOrd*ip];
   	   	index2 = (int)aIPx[MaxOrd*ip+1];

      	Lpt = L0pt[ip] + L1pt[ip]*(x[index1]-x[index2])
      			+ L2pt[ip]*pow((x[index1]-x[index2]),2.)
      			+ L3pt[ip]*pow((x[index1]-x[index2]),3.);

      	Lv = L0v[ip] + L1v[ip]*(x[index1]-x[index2])
      			+ L2v[ip]*pow((x[index1]-x[index2]),2.)
      			+ L3v[ip]*pow((x[index1]-x[index2]),3.);

   	   	Lh = (L0h[ip]-L0cp[ip]*Tk)
   	  			+ (L1h[ip]-L1cp[ip]*Tk)*(x[index1]-x[index2])
      			+ (L2h[ip]-L2cp[ip]*Tk)*pow((x[index1]-x[index2]),2.)
      			+ (L3h[ip]-L3cp[ip]*Tk)*pow((x[index1]-x[index2]),3.);

   	   	Ls = (-L0s[ip]-L0cp[ip]*(1.+log(Tk)))
      			+ (-L1s[ip]-L1cp[ip]*(1.+log(Tk)))*(x[index1]-x[index2])
      			+ (-L2s[ip]-L2cp[ip]*(1.+log(Tk)))*pow((x[index1]-x[index2]),2.)
      			+ (-L3s[ip]-L3cp[ip]*(1.+log(Tk)))*pow((x[index1]-x[index2]),3.);

   	   	Lcp = (-L0cp[ip]) + (-L1cp[ip])*(x[index1]-x[index2])
      			+ (-L2cp[ip])*pow((x[index1]-x[index2]),2.)
      			+ (-L3cp[ip])*pow((x[index1]-x[index2]),3.);

      	Gex += x[index1]*x[index2]*Lpt;
      	Vex += x[index1]*x[index2]*Lv;
      	Hex += x[index1]*x[index2]*Lh;
      	Sex += x[index1]*x[index2]*Ls;
      	CPex += x[index1]*x[index2]*Lcp;
  	}

   	Gex_ = Gex;
   	Vex_ = Vex;
   	Hex_ = Hex;
   	Sex_ = Sex;
   	CPex_ = CPex;

   	delete[]L0h;
   	delete[]L0s;
   	delete[]L0cp;
   	delete[]L0v;
   	delete[]L0pt;

   	delete[]L1h;
   	delete[]L1s;
   	delete[]L1cp;
   	delete[]L1v;
   	delete[]L1pt;

   	delete[]L2h;
   	delete[]L2s;
   	delete[]L2cp;
   	delete[]L2v;
   	delete[]L2pt;

   	delete[]L3h;
   	delete[]L3s;
   	delete[]L3cp;
   	delete[]L3v;
   	delete[]L3pt;

   	return 0;
}






// add other solution models here





//--------------------- End of s_fgl2.cpp ---------------------------
