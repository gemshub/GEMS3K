//-------------------------------------------------------------------
// $Id: s_fgl2.cpp 1143 2008-12-10 14:40:41Z gems $
//
// Copyright (c) 2007,2008  Th.Wagner, D.Kulik, S.Dmitrieva
//
// Implementation of the TSolMod class
// and TVanLaar, TRegular, TRedlichKister, TNRTL and TWilson classes
// Started by Th.Wagner and D.Kulik on 07.03.2007
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

//--------------------------------------------------------------------------------------------------------------
// Generic constructor for the TSolMod class
//
TSolMod::TSolMod( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, long int NPTPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
    ModCode(Mod_Code), NComp(NSpecies),  NPar(NParams), NPcoef(NPcoefs),
    MaxOrd(MaxOrder),  NP_DC(NPperDC), NPTP_DC(NPTPperDC), R_CONST(8.31451), RhoW(dW),
    EpsW(eW),  Tk(T_k), Pbar(P_bar)

{
    // pointer assignment
    aIPx = arIPx;   // Direct access to index list and parameter coeff arrays!
    aIPc = arIPc;
    aIP = new double[ NPar ];
    aDCc = arDCc;
    if( NPTP_DC )
    {	aDC = new double *[NComp];
       for (long int i=0; i<NComp; i++)
           aDC[i] = new double[NPTP_DC];
    }
    else aDC = 0;
    x = arWx;
    phVOL = aphVOL;
//    aZ = arZ;
//    aM =	arM;
    lnGamma = arlnGam;
}


TSolMod::~TSolMod()
{
   delete[] aIP;
   if( aDC )
   {
	for ( long i=0; i<NComp; i++)
      delete[] aDC[i];
    delete[] aDC;
   }
}



//=============================================================================================
// Van Laar model for solid solutions (c) TW March 2007
// References:  Holland & Powell (2003)
//=============================================================================================


// Generic constructor for the TVanLaar class
TVanLaar::TVanLaar( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0, 
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  alloc_internal();
  // PTparam();
}


TVanLaar::~TVanLaar()
{
  free_internal();
}


void TVanLaar::alloc_internal()
{
	   Wu = new double [NPar];
	   Ws = new double [NPar];
	   Wv = new double [NPar];
	   Wpt = new double [NPar];
	   Phi = new double [NComp];
	   PsVol = new double [NComp];
}


void TVanLaar::free_internal()
{
	 if(Wu)  delete[]Wu;
	 if(Ws)  delete[]Ws;
	 if(Wv)  delete[]Wv;
	 if(Wpt)  delete[]Wpt;
	 if(Phi)  delete[]Phi;
	 if(PsVol)  delete[]PsVol;
}


// Calculates T,P corrected binary interaction parameters
long int TVanLaar::PTparam()
{
	long int ip;
    if ( NPcoef < 3 || NPar < 1 )
       return 1;

// read P-T corrected interaction parameters
	   for (ip=0; ip<NPar; ip++)
	   {
	     Wu[ip] = aIPc[NPcoef*ip];
		 Ws[ip] = aIPc[NPcoef*ip+1];
		 Wv[ip] = aIPc[NPcoef*ip+2];
		 Wpt[ip] = Wu[ip]+ Ws[ip]*Tk + Wv[ip]*Pbar;
	     aIP[ip] = Wpt[ip];
		 // aIPc[NPcoef*ip+3] = Wpt[ip]; // obsolete
	   }
	   return 0;
}


// Calculates activity coefficients and excess functions
long int TVanLaar::MixMod()
{
   long int ip, j, i1, i2;
   double dj, dk;
   double sumPhi; // Sum of Phi terms
   double gEX, vEX, hEX, sEX, cpEX, uEX;

   if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
           return 1;

   // calculating Phi values
   sumPhi = 0.;
   for (j=0; j<NComp; j++)
   {
       PsVol[j] = aDCc[NP_DC*j];  // reading pseudo-volumes
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
        i1 = aIPx[MaxOrd*ip];
	    i2 = aIPx[MaxOrd*ip+1];

   	    if( j == i1 )
		dj = 1.;
	    else
		dj = 0.;
	    if( j == i2 )
		dk = 1.;
	    else
		dk = 0.;
	    lnGamRT -= (dj-Phi[i1])*(dk-Phi[i2])*Wpt[ip]
                         *2.*PsVol[j]/(PsVol[i1]+PsVol[i2]);
	}

    lnGam = lnGamRT/(R_CONST*Tk);
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
      i1 = aIPx[MaxOrd*ip];
      i2 = aIPx[MaxOrd*ip+1];
      gEX += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wpt[ip];
      vEX += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wv[ip];
      uEX += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wu[ip];
      sEX -= Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Ws[ip];
   }

   hEX = uEX+vEX*Pbar;

   return 0;
}



//=============================================================================================
// Regular model for multicomponent solid solutions (c) TW March 2007
// References:  Holland & Powell (1993)
//=============================================================================================


// Generic constructor for the TRegular class
TRegular::TRegular( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0, 
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  alloc_internal();
  // PTparam();
}


TRegular::~TRegular()
{
  free_internal();
}


void TRegular::alloc_internal()
{
	   Wu = new double [NPar];
	   Ws = new double [NPar];
	   Wv = new double [NPar];
	   Wpt = new double [NPar];
}


void TRegular::free_internal()
{
	 if(Wu)  delete[]Wu;
	 if(Ws)  delete[]Ws;
	 if(Wv)  delete[]Wv;
	 if(Wpt)  delete[]Wpt;
}


//   Calculates T,P corrected binary interaction parameters
long int TRegular::PTparam()
{
	long int ip;
   if ( NPcoef < 3 || NPar < 1 )
	           return 1;

// read interaction parameters and correct to T,P
	   for (ip=0; ip<NPar; ip++)
	   {
	     Wu[ip] = aIPc[NPcoef*ip];
		 Ws[ip] = aIPc[NPcoef*ip+1];
		 Wv[ip] = aIPc[NPcoef*ip+2];
		 Wpt[ip] = Wu[ip]+ Ws[ip]*Tk + Wv[ip]*Pbar;
	     aIP[ip] = Wpt[ip];
		 // aIPc[NPcoef*ip+3] = Wpt[ip]; // obsolete
	   }
	   return 0;
}


// Calculates activity coefficients and excess functions
long int
TRegular::MixMod()
{
   long int ip, j, i1, i2;
   double dj, dk;
   double gEX, vEX, hEX, sEX, cpEX, uEX;

   if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
           return 1;

   // calculate activity coefficients
   for (j=0; j<NComp; j++)      // index end members with j
   {
	lnGamRT = 0.;
	for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
	{
        i1 = aIPx[MaxOrd*ip];
	    i2 = aIPx[MaxOrd*ip+1];

   	    if( j == i1 )
		dj = 1.;
	    else
		dj = 0.;
	    if( j == i2 )
		dk = 1.;
	    else
		dk = 0.;
	    lnGamRT -= (dj-x[i1])*(dk-x[i2])*Wpt[ip];
	 }

    lnGam = lnGamRT/(R_CONST*Tk);
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
      i1 = aIPx[MaxOrd*ip];
      i2 = aIPx[MaxOrd*ip+1];
      gEX += x[i1]*x[i2]*Wpt[ip];
      vEX += x[i1]*x[i2]*Wv[ip];
      uEX += x[i1]*x[i2]*Wu[ip];
      sEX -= x[i1]*x[i2]*Ws[ip];
   }

   hEX = uEX+vEX*Pbar;

   return 0;
}



//=============================================================================================
// Redlich-Kister model for multicomponent solid solutions (c) TW March 2007
// References: Hillert (1998)
//=============================================================================================


// Generic constructor for the TRedlichKister class
TRedlichKister::TRedlichKister( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  alloc_internal();
  // PTparam();
}


TRedlichKister::~TRedlichKister()
{
  free_internal();
}


void TRedlichKister::alloc_internal()
{
	   Lu = new double [NPar][4];
	   Ls = new double [NPar][4];
	   Lcp = new double [NPar][4];
	   Lv = new double [NPar][4];
	   Lpt = new double [NPar][4];
}


void TRedlichKister::free_internal()
{
	 if(Lu)  delete[]Lu;
	 if(Ls)  delete[]Ls;
	 if(Lv)  delete[]Lv;
	 if(Lpt)  delete[]Lpt;
	 if(Lcp)  delete[]Lcp;
}


//   Calculates T,P corrected binary interaction parameters
long int TRedlichKister::PTparam()
{
   long int ip;

   if ( NPcoef < 16 || NPar < 1 )
      return 1;

   // read in interaction parameters
  	for (ip=0; ip<NPar; ip++)
  	{
	   	Lu[ip][0] = aIPc[NPcoef*ip+0];
	   	Ls[ip][0] = aIPc[NPcoef*ip+1];
	   	Lcp[ip][0] = aIPc[NPcoef*ip+2];
	   	Lv[ip][0] = aIPc[NPcoef*ip+3];
	   	Lpt[ip][0] = Lu[ip][0] + Ls[ip][0]*Tk + Lcp[ip][0]*Tk*log(Tk) + Lv[ip][0]*Pbar;
	    aIP[ip] = Lpt[ip][0];
	// aIPc[NPcoef*ip+16] = Lpt[ip][0]; // obsolete
	   	Lu[ip][1] = aIPc[NPcoef*ip+4];
	   	Ls[ip][1] = aIPc[NPcoef*ip+5];
	   	Lcp[ip][1] = aIPc[NPcoef*ip+6];
	   	Lv[ip][1] = aIPc[NPcoef*ip+7];
	   	Lpt[ip][1] = Lu[ip][1] + Ls[ip][1]*Tk + Lcp[ip][1]*Tk*log(Tk) + Lv[ip][1]*Pbar;
	// aIPc[NPcoef*ip+17] = Lpt[ip][1]; // obsolete
	   	Lu[ip][2] = aIPc[NPcoef*ip+8];
	   	Ls[ip][2] = aIPc[NPcoef*ip+9];
	   	Lcp[ip][2] = aIPc[NPcoef*ip+10];
	   	Lv[ip][2] = aIPc[NPcoef*ip+11];
	   	Lpt[ip][2] = Lu[ip][2] + Ls[ip][2]*Tk + Lcp[ip][2]*Tk*log(Tk) + Lv[ip][2]*Pbar;
	// aIPc[NPcoef*ip+18] = Lpt[ip][2]; // obsolete
	   	Lu[ip][3] = aIPc[NPcoef*ip+12];
	   	Ls[ip][3] = aIPc[NPcoef*ip+13];
	   	Lcp[ip][3] = aIPc[NPcoef*ip+14];
	   	Lv[ip][3] = aIPc[NPcoef*ip+15];
	   	Lpt[ip][3] = Lu[ip][3] + Ls[ip][3]*Tk + Lcp[ip][3]*Tk*log(Tk) + Lv[ip][3]*Pbar;
	// aIPc[NPcoef*ip+19] = Lpt[ip][3]; // obsolete
	}
   return 0;
}


// Calculates activity coefficients and excess functions
long int
TRedlichKister::MixMod()
{
   long int ip, j;
   long int i1, i2, L, I, J;
   double LU, LS, LCP, LV, LPT;
   double L0, L1, L2, L3;
   double gEX, vEX, hEX, sEX, cpEX, uEX;

   if ( NPcoef < 16 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
           return 1;

	// calculate activity coefficients
	for (j=0; j<NComp; j++)      // index end members with j
	{
		lnGamRT = 0.;
		for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];

			if ( j == i1 || j == i2) // interaction terms with j
			{
				if ( i1 == j ) // check order of idexes
				{
					L = i1;
					I = i2;
					L0 = Lpt[ip][0];
					L1 = Lpt[ip][1];
					L2 = Lpt[ip][2];
					L3 = Lpt[ip][3];
				}
				else
				{
					L = i2;
					I = i1;
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
				I = i1;
				J = i2;
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
   	   	i1 = aIPx[MaxOrd*ip];
   	   	i2 = aIPx[MaxOrd*ip+1];

      	LPT = Lpt[ip][0] + Lpt[ip][1]*(x[i1]-x[i2])
      			+ Lpt[ip][2]*pow((x[i1]-x[i2]),2.)
      			+ Lpt[ip][3]*pow((x[i1]-x[i2]),3.);

      	LV = Lv[ip][0] + Lv[ip][1]*(x[i1]-x[i2])
      			+ Lv[ip][2]*pow((x[i1]-x[i2]),2.)
      			+ Lv[ip][3]*pow((x[i1]-x[i2]),3.);

   	   	LU = (Lu[ip][0]-Lcp[ip][0]*Tk)
   	  			+ (Lu[ip][1]-Lcp[ip][1]*Tk)*(x[i1]-x[i2])
      			+ (Lu[ip][2]-Lcp[ip][2]*Tk)*pow((x[i1]-x[i2]),2.)
      			+ (Lu[ip][3]-Lcp[ip][3]*Tk)*pow((x[i1]-x[i2]),3.);

   	   	LS = (-Ls[ip][0]-Lcp[ip][0]*(1.+log(Tk)))
      			+ (-Ls[ip][1]-Lcp[ip][1]*(1.+log(Tk)))*(x[i1]-x[i2])
      			+ (-Ls[ip][2]-Lcp[ip][2]*(1.+log(Tk)))*pow((x[i1]-x[i2]),2.)
      			+ (-Ls[ip][3]-Lcp[ip][3]*(1.+log(Tk)))*pow((x[i1]-x[i2]),3.);

   	   	LCP = (-Lcp[ip][0]) + (-Lcp[ip][1])*(x[i1]-x[i2])
      			+ (-Lcp[ip][2])*pow((x[i1]-x[i2]),2.)
      			+ (-Lcp[ip][3])*pow((x[i1]-x[i2]),3.);

      	gEX += x[i1]*x[i2]*LPT;
      	vEX += x[i1]*x[i2]*LV;
      	uEX += x[i1]*x[i2]*LU;
      	sEX += x[i1]*x[i2]*LS;
      	cpEX += x[i1]*x[i2]*LCP;
  	}

   	hEX = uEX+vEX*Pbar;
   	return 0;
}



//=============================================================================================
// NRTL model for liquid solutions (c) TW June 2008
// References: Renon and Prausnitz (1968), Prausnitz et al. (1997)
//=============================================================================================


// Generic constructor for the TNRTL class
TNRTL::TNRTL( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  alloc_internal();
  // PTparam();
}


TNRTL::~TNRTL()
{
  free_internal();
}


void TNRTL::alloc_internal()
{
	Tau = new double *[NComp];
	dTau = new double *[NComp];
	d2Tau = new double *[NComp];
	Alp = new double *[NComp];
	dAlp = new double *[NComp];
	d2Alp = new double *[NComp];
	G = new double *[NComp];
	dG = new double *[NComp];
	d2G = new double *[NComp];
    for (long int j=0; j<NComp; j++)
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

}


void TNRTL::free_internal()
{
  	// cleaning memory
	   	for (long int j=0; j<NComp; j++)
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
}


//   Calculates T,P corrected binary interaction parameters
long int TNRTL::PTparam()
{
	long int ip, i, j, i1, i2;
	double A, B, C, D, E, F;
	double tau, dtau, d2tau, alp, dalp, d2alp;

    if ( NPcoef < 6 || NPar < 1 )
       return 1;

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
		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		A = aIPc[NPcoef*ip+0];
		B = aIPc[NPcoef*ip+1];
		C = aIPc[NPcoef*ip+2];
		D = aIPc[NPcoef*ip+3];
		E = aIPc[NPcoef*ip+4];
		F = aIPc[NPcoef*ip+5];
		tau = A + B/Tk + C*Tk + D*log(Tk);	// partial derivatives of tau and alp
		dtau = - B/pow(Tk,2.) + C + D/Tk;
		d2tau = 2.*B/pow(Tk,3.) - D/pow(Tk,2.);
		alp = E + F*(Tk-273.15);
		dalp = F;
		d2alp = 0.0;

		Tau[i1][i2] = tau;
		dTau[i1][i2] =  dtau;
		d2Tau[i1][i2] = d2tau;
		Alp[i1][i2] = alp;
		dAlp[i1][i2] = dalp;
		d2Alp[i1][i2] =  d2alp;
		//aIPc[NPcoef*ip+6] = tau;          //obsolete
		//aIPc[NPcoef*ip+7] = dtau;
		//aIPc[NPcoef*ip+8] = d2tau;
		//aIPc[NPcoef*ip+9] = alp;
		//aIPc[NPcoef*ip+10] = dalp;
		//aIPc[NPcoef*ip+11] = d2alp;

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
   return 0;
}


// Calculates activity coefficients and excess functions
// heat capacity calculation added, 06.06.2008 (TW)
long int
TNRTL::MixMod()
{
	long int  j, i, k;
	double K, L, M, N, O;
	double U, dU, V, dV, d2U, d2V;
	double g, dg, d2g, lnGam;
	double gEX, vEX, hEX, sEX, cpEX;

	if ( NPcoef < 6 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

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

	return 0;
}



//=============================================================================================
// Wilson model for liquid solutions (c) TW June 2008
// References: Prausnitz et al. (1997)
//=============================================================================================


// Generic constructor for the TWilson class
TWilson::TWilson( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  alloc_internal();
  // PTparam();
}


TWilson::~TWilson()
{
  free_internal();
}


void TWilson::alloc_internal()
{
	Lam = new double *[NComp];
	dLam = new double *[NComp];
	d2Lam = new double *[NComp];

    for (long int j=0; j<NComp; j++)
    {
		Lam[j] = new double [NComp];
		dLam[j] = new double [NComp];
		d2Lam[j] = new double [NComp];
	}
}


void TWilson::free_internal()
{
   	// cleaning memory
   	for (long int j=0; j<NComp; j++)
   	{
		delete[]Lam[j];
		delete[]dLam[j];
		delete[]d2Lam[j];
	}
	delete[]Lam;
	delete[]dLam;
	delete[]d2Lam;
}


// Calculates T-corrected interaction parameters
long int TWilson::PTparam()
{
	long int ip, i, j, i1, i2;
	double A, B, C, D;
	double lam, dlam, d2lam;

    if ( NPcoef < 4 || NPar < 1 )
           return 1;

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
		A = aIPc[NPcoef*ip+0];
		B = aIPc[NPcoef*ip+1];
		C = aIPc[NPcoef*ip+2];
		D = aIPc[NPcoef*ip+3];
		lam = exp( A + B/Tk + C*Tk + D*log(Tk) );
		dlam = lam*( - B/pow(Tk,2.) + C + D/Tk );
		d2lam = dlam*( - B/pow(Tk,2.) + C + D/Tk ) + lam*( 2.*B/pow(Tk,3.) - D/pow(Tk,2.) );

		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		Lam[i1][i2] = lam;
		dLam[i1][i2] = dlam;
		d2Lam[i1][i2] = d2lam;
		//aIPc[NPcoef*ip+4] = lam;     //obsolete
		//aIPc[NPcoef*ip+5] = dlam;
		//aIPc[NPcoef*ip+6] = d2lam;
	}
	return 0;
}


// Calculates activity coefficients and excess functions
// heat capacity calculation added, 06.06.2008 (TW)
long int
TWilson::MixMod( )
{
	long int  j, i, k;
	double K, L, M;
	double U, dU, d2U;
	double g, dg, d2g, lnGam;
	double gEX, vEX, hEX, sEX, cpEX;

	if ( NPcoef < 4 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

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

	return 0;
}


//--------------------- End of s_fgl2.cpp ---------------------------
