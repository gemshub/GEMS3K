//-------------------------------------------------------------------
// $Id: s_fgl.cpp 1209 2009-01-31 18:13:53Z wagner $
//
// Copyright (C) 2004-2009  T.Wagner, S.Churakov, D.Kulik
//
// Implementation of TPRSVcalc, TCGFcalc and TSRKcalc classes
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// and part of the GEMIPM2K standalone code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------

#include <math.h>
#include <stdio.h>

#include "s_fgl.h"
#include "verror.h"
#ifndef IPMGEMPLUGIN
  #include "m_const.h"
#endif



//=======================================================================================================
// Peng-Robinson-Stryjek-Vera (PRSV) model for fluid mixtures
// References: Stryjek and Vera (1986), Proust and Vera (1989)
// Implementation of the TPRSVcalc class
//=======================================================================================================

// Constructor
TPRSVcalc::TPRSVcalc( long int NCmp, double Pp, double Tkp ):
	TSolMod( NCmp, 0, 0, 0, 0, 4, 'P',
         0, 0, 0, 0, 0, 0, Tkp, Pp, 0, 0  )

{
	aGEX = 0;
    aVol = 0;
	Pparc = 0;
	alloc_internal();
}


TPRSVcalc::TPRSVcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int *arIPx, double *arIPc, double *arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arPparc,
        double *arGEX, double *arVol, double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 4,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
{
	aGEX = arGEX;
	aVol = arVol;
	Pparc = arPparc;
	alloc_internal();
}


TPRSVcalc::~TPRSVcalc()
{
	free_internal();
}


// allocate work arrays for pure fluid and fluid mixture properties
void TPRSVcalc::alloc_internal()
{
	Eosparm = new double [NComp][6];
	Pureparm = new double [NComp][4];
	Fugpure = new double [NComp][6];
	Fugci = new double [NComp][4];
	DepPh = new double [7];
	KK = new double *[NComp];
	dKK = new double *[NComp];
	d2KK = new double *[NComp];
	AA = new double *[NComp];

    for (long int i=0; i<NComp; i++)
    {
    	KK[i] = new double[NComp];
    	dKK[i] = new double[NComp];
    	d2KK[i] = new double[NComp];
    	AA[i] = new double[NComp];
    }
}


void TPRSVcalc::free_internal()
{
	long int i;

	for (i=0; i<NComp; i++)
	{
		delete[]KK[i];
		delete[]dKK[i];
		delete[]d2KK[i];
		delete[]AA[i];
	}

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
	delete[]DepPh;
	delete[]KK;
	delete[]dKK;
	delete[]d2KK;
	delete[]AA;
}


// High-level method to retrieve pure fluid fugacities
long int TPRSVcalc::PureSpecies()
{
	double Fugcoeff, Volume, DeltaH, DeltaS;
	long int j, retCode = 0;

	for( j=0; j<NComp; j++)
    {
		// Calling PRSV EoS for pure fugacity
        retCode =  PRFugacityPT( j, Pbar, Tk, aDCc+j*NP_DC,
						aDC[j], Fugcoeff, Volume, DeltaH, DeltaS );

        aGEX[j] = log( Fugcoeff );  // now here (since 26.02.2008) DK
        Pparc[j] = Fugcoeff * Pbar;  // Necessary only for performance
        aVol[j] = Volume * 10.;  // molar volume of pure fluid component, J/bar to cm3
	} // j

	if ( retCode )
	{
		char buf[150];
		sprintf(buf, "PRSV Fluid(): calculation of pure fugacity failed");
		Error( "E71IPM IPMgamma: ",  buf );
	}
	return 0;
}


// High-level method to calculate T,P corrected binary interaction parameters
long int TPRSVcalc::PTparam()
{
	long int j, i, ip;
	long int i1, i2;
	double p0, p1, p2;
	double k, dk, d2k;

	PureSpecies();

	if( NPcoef > 0 )
	{
		// fill internal array of interaction parameters with standard value
		for( j=0; j<NComp; j++ )
		{
			for( i=0; i<NComp; i++ )
			{
				KK[j][i] = 0.;
				dKK[j][i] = 0.;
				d2KK[j][i] = 0.;
			}
		}

		// transfer those interaction parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			p0 = aIPc[NPcoef*ip];
			k = p0;
			KK[i1][i2] = k;
			KK[i2][i1] = k;   // symmetric case
		}
	}

	return 0;
}


// High-level method to retrieve activity coefficients of the fluid mixture
long int TPRSVcalc::MixMod()
{
	long int j, iRet;

    iRet = FugacitySpec( Pparc );

    phVOL[0] = PhVol * 10.;

    for(j=0; j<NComp; j++)
    {
    	if( Fugci[j][3] > 1e-23 )
    		lnGamma[j] = log( Fugci[j][3] );
        else
        	lnGamma[j] = 0;
    }
    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PRSV Fluid(): calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }
    return iRet;
}


// High-level method to retrieve departure functions of the fluid mixture
long int TPRSVcalc::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	long int iRet;

	iRet = DepartureFunct( Pparc );

    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PRSV Fluid(): calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }

	// assignments
	Gex_ = DepPh[0];
	Sex_ = DepPh[1];
	Hex_ = DepPh[2];
	CPex_ = DepPh[3];
	Vex_ = DepPh[4];

	return iRet;
}


// retrieve pure fluid properties
long int TPRSVcalc::PRFugacityPT( long int i, double P, double Tk, double *EoSparam, double *Eos2parPT,
        double &Fugacity, double &Volume, double &DeltaH, double &DeltaS )
{

	long int iRet = 0;
    double Tcrit, Pcrit, omg, k1, k2, k3;
    double apure, bpure, da, d2a;

    // reads EoS parameters from database into work array
    if( !EoSparam )
    	return -1;  // Memory alloc error

    Eosparm[i][0] = EoSparam[0];   // critical temperature in K
    Eosparm[i][1] = EoSparam[1];   // critical pressure in bar
    Eosparm[i][2] = EoSparam[2];   // Pitzer acentric factor omega
    Eosparm[i][3] = EoSparam[3];   // empirical EoS parameter k1
    Eosparm[i][4] = EoSparam[4];   // empirical EoS parameter k2
    Eosparm[i][5] = EoSparam[5];   // empirical EoS parameter k3
    Tcrit = Eosparm[i][0];
    Pcrit = Eosparm[i][1];
    omg = Eosparm[i][2];
    k1 = Eosparm[i][3];
	k2 = Eosparm[i][4];
	k3 = Eosparm[i][5];

	AB( Tcrit, Pcrit, omg, k1, k2, k3, apure, bpure, da, d2a );

	Pureparm[i][0] = apure;
	Pureparm[i][1] = bpure;
	Pureparm[i][2] = da;
	Pureparm[i][3] = d2a;
	Eos2parPT[0] = apure;
	Eos2parPT[1] = bpure;
	Eos2parPT[2] = da;
	Eos2parPT[3] = d2a;

	iRet = FugacityPure( i );
	if( iRet)
		return iRet;

	Fugacity = Fugpure[i][0];  // Fugacity coefficient
	DeltaH = Fugpure[i][2];  // H departure function
	DeltaS = Fugpure[i][3];  // S departure function
	Volume = Fugpure[i][4];  //  J/bar

	return iRet;
}


// Calculates attractive (a) and repulsive (b) parameter of PRSV equation of state
// and partial derivatives of alpha function
long int TPRSVcalc::AB( double Tcrit, double Pcrit, double omg, double k1, double k2, double k3,
		double &apure, double &bpure, double &da, double &d2a )
{
	double Tred, k0, k, alph, ac, sqa, dsqa, d2sqa;

	Tred = Tk/Tcrit;
	k0 = 0.378893 + 1.4897153*omg - 0.17131848*pow(omg,2.) + 0.0196554*pow(omg,3.);
	if(Tk >= Tcrit)
	{
		k1 = 0.0;
		k2 = 0.0;
		k3 = 0.0;
	}
	k = k0 + (k1 + k2*(k3-Tred)*(1.-sqrt(Tred))) * (1.+sqrt(Tred)) * (0.7-Tred);
	alph = pow(1. + k*(1.-sqrt(Tred)), 2.);
	ac = 0.457235*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit;
	apure = alph*ac;
	bpure = 0.077796*R_CONST*Tcrit/Pcrit;
	sqa = 1.+k*(1.-sqrt(Tred));
	// dsqa = (-1.)*k0/(2.*sqrt(Tk*Tcrit)) - 1.7*k1/Tcrit + 2.*k1*Tk/(pow(Tcrit,2.));  // extend dA/dT for k2, k3
	dsqa = ( k1*(0.7-Tred)/(2.*sqrt(Tred)*Tcrit) - k1*(1.+sqrt(Tred))/Tcrit ) * (1.-sqrt(Tred))
				- (k0 + k1*(1.+sqrt(Tred))*(0.7-Tred))/(2*sqrt(Tred)*Tcrit);
	da = 2.*ac*(sqa*dsqa);
	d2sqa = ( - (k1*(0.7-Tred))/(4.*pow(Tred,1.5)*pow(Tcrit,2.)) - k1/(sqrt(Tred)*pow(Tcrit,2.)) ) * (1.-sqrt(Tred))
				+ ( - k1*(0.7-Tred)/(2.*sqrt(Tred)*Tcrit) - k1*(1.+sqrt(Tred))/Tcrit ) / (sqrt(Tred)*Tcrit)
				+ ( k0 + k1*(1.+sqrt(Tred))*(0.7-Tred) ) / (4.*pow(Tred,1.5)*pow(Tcrit,2.));
	d2a = 2.*ac*(dsqa*dsqa + sqa*d2sqa);

	return 0;
}


// Calculates fugacities and departure functions of pure fluid species
long int TPRSVcalc::FugacityPure( long int i )
{
	double Tcrit, Pcrit, Tred, aprsv, bprsv, alph, da, d2a;
	double k, A, B, a2, a1, a0, z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, cpig, fugpure;
	double gdep, hdep, sdep, cpdep;
	double cv, dPdT, dPdV, dVdT;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONST*log(Pbar);
	gig = hig - Tk*sig;
	cpig = 0.;

	// retrieve a and b terms of cubic EoS
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	Tred = Tk/Tcrit;
	aprsv = Pureparm[i][0];
	bprsv = Pureparm[i][1];
	da = Pureparm[i][2];
	d2a = Pureparm[i][3];

	// solve cubic equation
	A = aprsv*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bprsv*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		z = z2; vol = vol2; lnf = lnf2;
	}
	else
	{
		z = z1; vol = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		z = z3; vol = vol3; lnf = lnf3;
	}
	else
	{
		z = z; vol = vol; lnf = lnf;
	}

	// calculate thermodynamic properties
	alph = aprsv/(0.457235*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit);
	k = (sqrt(alph)-1.)/(1.-sqrt(Tred));
	gdep = R_CONST*Tk*(z-1.-log(z-B)-A/(B*sqrt(8.))
				*log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B)));
	hdep = R_CONST*Tk*(z-1.-log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B))
				*A/(B*sqrt(8.))*(1+k*sqrt(Tred)/sqrt(alph)));
	sdep = (hdep-gdep)/Tk;

	// heat capacity part
	cv = Tk*d2a/(bprsv*sqrt(8.))
			 * log( (z+B*(1.+sqrt(2.)))/(z+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vol-bprsv) - da/( vol*(vol+bprsv) + bprsv*(vol-bprsv) );
	dPdV = - R_CONST*Tk/pow((vol-bprsv),2.) + 2*aprsv*(vol+bprsv)/pow((vol*(vol+bprsv)+bprsv*(vol-bprsv)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	cpdep = cv + Tk*dPdT*dVdT - R_CONST;

	// increment thermodynamic properties
	fugpure = exp(lnf);
	Fugpure[i][0] = fugpure;
	Fugpure[i][1] = gdep;  // changed to departure functions, 31.05.2008 (TW)
	Fugpure[i][2] = hdep;
	Fugpure[i][3] = sdep;
    Fugpure[i][4] = vol;
    Fugpure[i][5] = cpdep;

    return 0;
}


// Cubic equation root solver based on Cardanos method
long int TPRSVcalc::Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 )
{
	double q, rc, q3, rc2, theta, ac, bc;

	q = (pow(a2,2.) - 3.*a1)/9.;
	rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
	q3 = pow(q,3.);
	rc2 = pow(rc,2.);
	if (rc2 < q3)  // three real roots
	{
		theta = acos(rc/sqrt(q3));
		z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
		z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
		z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
	}
	else  // one real root
	{
		ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
		if (ac != 0.)
			bc = q/ac;
		else
			bc = 0.;
		z1 = ac+bc-a2/3.;
		z2 = ac+bc-a2/3.;
		z3 = ac+bc-a2/3.;
	}
	return 0;
}


// Calculates mixing properties of the fluid mixture
long int TPRSVcalc::MixParam( double &amix, double &bmix )
{
	long int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
            K = KK[i][j];
			AA[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}
	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + x[i]*x[j]*AA[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + x[i]*Pureparm[i][1];
	}
	return 0;
}


// Calculates fugacity of the bulk fluid mixture
long int TPRSVcalc::FugacityMix( double amix, double bmix, double &fugmix, double &zmix,
		double &vmix )
{
	double A, B, a2, a1, a0, z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	A = amix*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bmix*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano( a2, a1, a0, z1, z2, z3 );

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}
	fugmix = exp(lnf);
        PhVol = vmix;
	return 0;
}


// Calculates fugacities and activities of fluid species in the mixture,
long int TPRSVcalc::FugacitySpec( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double A, B, lnfci, fci;

    // Reload params to Pureparm
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + x[j]*AA[i][j];
		}
		lnfci = Pureparm[i][1]/bmix*(zmix-1.) - log(zmix-B)
		      + A/(sqrt(8.)*B)*(2.*sum/amix-Pureparm[i][1]/bmix)
                      * log((zmix+B*(1.-sqrt(2.)))/(zmix+B*(1.+sqrt(2.))));
		fci = exp(lnfci);
		Fugci[i][0] = fci;  // fugacity coefficient using engineering convention
		Fugci[i][1] = x[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (x[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/x[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;
	}

	return iRet;
}


// calculates departure functions in the mixture bla
long int TPRSVcalc::DepartureFunct( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0.;
	double A, B;
	double Gig, Hig, Sig, CPig, Gdep, Hdep, Sdep, CPdep;
	double K, dK, d2K, Q, dQ, d2Q;
	double damix, d2amix, ai, aj, dai, daj, d2ai, d2aj;
	double cv, dPdT, dPdV, dVdT;

    // Reload params to Pureparm (probably now obsolete?)
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// ideal gas changes from 1 bar to P (at T of interest)
	Hig = 0.;
	Sig = (-1.)*R_CONST*log(Pbar);
	Gig = Hig - Tk*Sig;
	CPig = 0.;

	// calculate total state functions of the mixture
	damix = 0.;
	d2amix = 0.;
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// pull parameters
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];
			K = KK[i][j];
			dK = dKK[i][j];
			d2K = d2KK[i][j];

			// increments to derivatives
			Q = sqrt(ai*aj);
			dQ = 0.5*( sqrt(aj/ai)*dai + sqrt(ai/aj)*daj );
			d2Q = 0.5*( dai*daj/sqrt(ai*aj) + d2ai*sqrt(aj)/sqrt(ai) + d2aj*sqrt(ai)/sqrt(aj)
					- 0.5*( pow(dai,2.)*sqrt(aj)/sqrt(pow(ai,3.))
					+ pow(daj,2.)*sqrt(ai)/sqrt(pow(aj,3.)) ) );
			damix = damix + x[i]*x[j] * ( dQ*(1.-K) - Q*dK );
			d2amix = d2amix + x[i]*x[j] * ( d2Q*(1.-K) - 2.*dQ*dK - Q*d2K );
		}
	}

	// calculate thermodynamic properties
	Gdep = (amix/(R_CONST*Tk*sqrt(8.)*bmix) * log((vmix+(1.-sqrt(2.))*bmix)
		/ (vmix+(1.+sqrt(2.))*bmix))-log(zmix*(1.-bmix/vmix))+zmix-1.)*R_CONST*Tk;
	Hdep = ((amix-Tk*damix)/(R_CONST*Tk*sqrt(8.)*bmix)*log((vmix+(1.-sqrt(2.))
		*bmix)/(vmix+(1.+sqrt(2))*bmix))+zmix-1.)*R_CONST*Tk;
	Sdep = (Hdep - Gdep)/Tk;

	// heat capacity part
	cv = Tk*d2amix/(bmix*sqrt(8.))
			 * log( (zmix+B*(1.+sqrt(2.)))/(zmix+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vmix-bmix) - damix/( vmix*(vmix+bmix) + bmix*(vmix-bmix) );
	dPdV = - R_CONST*Tk/pow((vmix-bmix),2.) + 2*amix*(vmix+bmix)/pow((vmix*(vmix+bmix)+bmix*(vmix-bmix)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	CPdep = cv + Tk*dPdT*dVdT - R_CONST;

	// assignments
	DepPh[0] = Gdep + Gig;
	DepPh[1] = Sdep + Sig;
	DepPh[2] = Hdep;
	DepPh[3] = CPdep;
	DepPh[4] = vmix;

	return iRet;
}


#ifndef IPMGEMPLUGIN
#include "s_tpwork.h"

// Calculates properties of pure fluids when called from RTParm
long int TPRSVcalc::PRCalcFugPure( void )
{
    double T, P, Fugcoeff = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
		// float*Coeff;
    double Coeff[7];  // MAXCRITPARAM
    double Eos2parPT[4] = { 0.0, 0.0, 0.0, 0.0 } ;
    long int retCode = 0;

    ErrorIf( !aW.twp, "PRSV EoS", "Undefined twp");

    P = aW.twp->P;    /* P in 10^5 Pa? */
    T = aW.twp->TC+273.15;   /* T?in K */

    for(long int ii=0; ii<7; ii++ )
      Coeff[ii] = aW.twp->CPg[ii];
		// Coeff = aW.twp->CPg;


    // Calling PRSV EoS functions here
    if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1e-5 && P < 1e5 )
       retCode = PRFugacityPT( 0, P, T, Coeff, Eos2parPT, Fugcoeff, Volume,
            DeltaH, DeltaS );
    else {
            Fugcoeff = 1.;
            Volume = 8.31451*T/P;
            aW.twp->V = Volume;
            aW.twp->Fug = Fugcoeff*P;
            return retCode;
          }

    // increment thermodynamic properties
    aW.twp->G += 8.31451 * T * log( Fugcoeff );  // from fugacity coeff
    aW.twp->H +=  DeltaH;  // to be completed
    aW.twp->S +=  DeltaS;  // to be completed
    aW.twp->V = Volume;  // in J/bar
    aW.twp->Fug = Fugcoeff * P;  // fugacity at P

    // passing corrected EoS coeffs to calculation of fluid mixtures
    aW.twp->wtW[6] = Eos2parPT[0];  // a
    aW.twp->wtW[7] = Eos2parPT[1];  // b
    aW.twp->wtW[8] = Eos2parPT[2];  // da
    aW.twp->wtW[9] = Eos2parPT[3];  // d2a

    return retCode;
}

#endif



//=======================================================================================================
// Churakov-Gottschalk (CG) model for fluid mixtures
// References: Churakov and Gottschalk (2003a, 2003b)
// Implementation of the TCGFcalc class
//=======================================================================================================

// Constructor
TCGFcalc::TCGFcalc( long int NCmp, double Pp, double Tkp ):
    TSolMod( NCmp, 0, 0, 0, 0, 0, 'F',
         0, 0, 0, 0, 0, 0, Tkp, Pp, 0., 0. )
{
	Pparc = 0;
	phWGT = 0;
	aX = 0;
    aGEX = 0;
	aVol = 0;

    set_internal();
	alloc_internal();
}


TCGFcalc::TCGFcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int *arIPx, double *arIPc, double *arDCc,
        double *arWx, double *arlnGam, double *aphVOL,
        double * arPparc, double *arphWGT,double *arX,
        double *arGEX, double *arVol, double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 8,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
{
	Pparc = arPparc;
	phWGT = arphWGT;
	aX = arX;
	aGEX = arGEX;
	aVol = arVol;
	set_internal();
	alloc_internal();
}


// Destructor
TCGFcalc::~TCGFcalc()
{
	free_internal();
}


// set internally used parameters
void TCGFcalc::set_internal()
{
	PI_1 = 3.141592653589793120;  // pi
	TWOPI = 6.283185307179586230;  // 2.*pi
	PISIX = 0.523598775598298927;  // pi/6.
	TWOPOW1SIX = 1.12246204830937302;  // 2^ = 1/6)
	DELTA  = 0.00001;
	DELTAMOLLIM  = 0.0000001;
	R = 8.31439;  // R constant
	NA = 0.6023;
	P1 = 1.186892378996;
	PP2 = -0.4721963005527;
	P3 = 3.259515855283;
	P4 = 3.055229342609;
	P5 = 1.095409321023;
	P6 = 1.282306659774E-2;
	P7 = 9.55712461425E-2;
	P8 = 13.67807693107;
	P9 = 35.75464856619;
	P10 = 16.04724381643;
	AA1 = -0.120078459237;
	AA2 = -.808712488307;
	AA3 = .321543801337;
	A4 = 1.16965477132;
	A5 = -.410564939543;
	A6 = -.516834310691;
	BB1 = -2.18839961483;
	BB2 = 1.59897428009;
	BB3 = -.392578806128;
	B4 = -.189396607904;
	B5 = -.576898496254;
	B6 = -0.0185167641359;
	A00 = .9985937977069455;
	A01 = .5079834224407451;
	A10 = 1.021887697885469;
	A11 = -5.136619463333883;
	A12 = -5.196188074016755;
	A21 = -6.049240839050804;
	A22 = 18.67848155616692;
	A23 = 20.10652684217768;
	A31 = 9.896491419756988;
	A32 = 14.6738380473899;
	A33 = -77.44825116542995;
	A34 = -4.82871082941229;
}


void TCGFcalc::alloc_internal()
{
	paar = 0;
	paar1 = 0;
    FugCoefs =  0;
    EoSparam =  0;
    EoSparam1 = 0;
    DepPh = new double [7];
}


void TCGFcalc::free_internal()
{
	if( paar )   delete paar;
	paar = 0;
	if( paar1 )  delete paar1;
    paar1 = 0;
	if( FugCoefs )  delete[]FugCoefs;
	if( EoSparam )  delete[]EoSparam;
	if( EoSparam1 ) delete[]EoSparam1;
	if (DepPh)	delete[]DepPh;
}


// High-level method to retrieve pure fluid fugacities
long int TCGFcalc::PureSpecies()
{
	double Fugacity = 0.1, Volume = 0.0;
	double X[1] = {1.};
	double Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
            Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;
	double roro;  // added, 21.06.2008 (TW)
	long int j, retCode = 0;

	for( j=0; j<NComp; j++)
	{
		// Calling CG EoS for pure fugacity
        if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
        	retCode = CGFugacityPT( aDCc+j*NP_DC, Eos4parPT, Fugacity, Volume, Pbar, Tk, roro );
        else {
            Fugacity = Pbar;
            Volume = 8.31451*Tk/Pbar;
            aDC[j][0] = aDCc[j*NP_DC];
            if( aDC[j][0] < 1. || aDC[j][0] > 10. )
            	aDC[j][0] = 1.;                 // foolproof temporary
            aDC[j][1] = aDCc[j*NP_DC+1];
            aDC[j][2] = aDCc[j*NP_DC+2];
            aDC[j][3] = aDCc[j*NP_DC+3];
            aDC[j][4] = 0.;
            aDC[j][5] = 0.;
            aDC[j][6] = 0.;
            aDC[j][7] = 0.;
            continue;
        }

        aGEX[j] = log( Fugacity / Pbar );  // now here (since 26.02.2008)  DK
        Pparc[j] = Fugacity;  // Necessary only for performance
        aVol[j] = Volume * 10.;  // molar volume of pure fluid component, J/bar to cm3

        // passing corrected EoS coeffs to calculation of fluid mixtures
        aDC[j][0] = Eos4parPT[0];
        if( aDC[j][0] < 1. || aDC[j][0] > 10. )
        	aDC[j][0] = 1.;                            // foolproof temporary
        aDC[j][1] = Eos4parPT[1];
        aDC[j][2] = Eos4parPT[2];
        aDC[j][3] = Eos4parPT[3];

        CGFugacityPT( aDCc+j*NP_DC, Eos4parPT1, Fugacity, Volume, Pbar, Tk+Tk*GetDELTA(), roro );

        // passing corrected EoS coeffs for T+T*DELTA
        aDC[j][4] =Eos4parPT1[0];
        if( aDC[j][4] < 1. || aDC[j][4] > 10. )
        	aDC[j][4] = 1.;                            // foolproof temporary
        aDC[j][5] = Eos4parPT1[1];
        aDC[j][6] = Eos4parPT1[2];
        aDC[j][7] = Eos4parPT1[3];

        // Calculation of departure functions
        CGDepartureFunct( X, Eos4parPT, Eos4parPT1, 1, roro, Tk );  // changed, 21.06.2008 (TW)
    }  // j

	if ( retCode )
	{
		char buf[150];
		sprintf(buf, "CG2004Fluid(): calculation of pure fugacity failed");
		Error( "E71IPM IPMgamma: ",  buf );
	}
	return 0;
}


// Calculates T,P corrected binary interaction parameters
long int TCGFcalc::PTparam()
{
	long int i,j;

	if( FugCoefs )  delete[]FugCoefs;
	if( EoSparam )  delete[]EoSparam;
	if( EoSparam1 ) delete[]EoSparam1;

    FugCoefs = new double[ NComp ];
    EoSparam = new double[ NComp*4 ];
    EoSparam1 = new double[ NComp*4 ];

    PureSpecies();

    // Copying T,P corrected coefficients
    for( j=0; j<NComp; j++)
    {
    	for( i=0; i<4; i++)
    		EoSparam[j*4+i] = aDC[j][i];
    	for( i=0; i<4; i++)
    		EoSparam1[j*4+i] = aDC[j][i+4];
    }
    return 0;
}


// High-level method to retrieve activity coefficients in the fluid mixture
long int TCGFcalc::MixMod()
{
	long int j;
	double roro; // changed, 21.06.2008 (TW)

	if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
	{
		CGActivCoefPT( aX, EoSparam, FugCoefs, NComp, Pbar, Tk, roro );  // changed, 21.06.2008 (TW)
		if (roro <= 0. )
		{
			char buf[150];
			sprintf(buf, "CGFluid(): bad calculation of density ro= %lg", roro);
			Error( "E71IPM IPMgamma: ",  buf );
		}

		// Phase volume of the fluid in cm3 (not needed any more?)
		phVOL[0] = phWGT[0] / roro;

	}

	else  // Setting Fugcoefs to 0 outside TP interval
		for( j=0; j<NComp; j++ )
			FugCoefs[ j ] = 0.0;

		for( j=0; j<NComp; j++  )
		{
			if( FugCoefs[j] > 1e-23 )
				lnGamma[j] = log(FugCoefs[j]/Pparc[j]);
			else
				lnGamma[j] = 0;
		}  // j
	return 0;
}


long int TCGFcalc::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	double roro; // changed, 21.06.2008 (TW)
	double Gig, Sig, Hig, CPig;

	// ideal gas changes from 1 bar to P (at T of interest)
	Hig = 0.;
	Sig = (-1.)*R_CONST*log(Pbar);
	Gig = Hig - Tk*Sig;
	CPig = 0.;

	if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
	{
		CGActivCoefPT( aX, EoSparam, FugCoefs, NComp, Pbar, Tk, roro );  // changed, 21.06.2008 (TW)
		if (roro <= 0. )
		{
			char buf[150];
			sprintf(buf, "CGFluid(): bad calculation of density ro= %lg", roro);
			Error( "E71IPM IPMgamma: ",  buf );
		}

		// calculate departure functions
		CGDepartureFunct( aX, EoSparam, EoSparam1, NComp, roro, Tk );

	}

	else  // setting departure functions to 0 outside TP interval
	{
		DepPh[0] = 0.;
		DepPh[1] = 0.;
		DepPh[2] = 0.;
		DepPh[3] = 0.;
		DepPh[4] = 0.;
	}

	// assignments
	Gex_ = DepPh[0] + Gig;
	Sex_ = DepPh[1] + Sig;
	Hex_ = DepPh[2];
	CPex_ = DepPh[3];
	Vex_ = DepPh[4];

	return 0;
}


// High-level method to retrieve pure fluid properties
long int TCGFcalc::CGFugacityPT( double *EoSparam, double *EoSparPT, double &Fugacity,
        double &Volume, double P, double T, double &roro )
{
	long int iRet = 0;
	// double ro;
	double X[1] = {1.};
	double FugPure[1];

	// modification to simplify CG database structure, 20.03.2007 (TW)
	EoSparPT[0] = EoSparam[0]+EoSparam[4]*exp(T*EoSparam[5]);
	EoSparPT[1] = EoSparam[1]+EoSparam[6]*exp(T*EoSparam[7]);
	EoSparPT[2] = EoSparam[2]+EoSparam[8]/(T+EoSparam[9]);
	EoSparPT[3] = EoSparam[3]+EoSparam[10]/(T+EoSparam[11]);

	// returns density
	CGActivCoefPT( X, EoSparPT, FugPure, 1, P, T, roro );  // changed, 21.06.2008 (TW)
	if( roro < 0.  )
	{
		return -1;
	};

	Fugacity = FugPure[0];
	roro = DENSITY( X, EoSparPT, 1, P, T );

	if( roro < 0 )
	{  // error - density could not be calculated
		iRet = -2;
		roro = 1.0;
	}

	Volume = 0.1/roro;  // in J/bar
	// roro = ro;  // added, 21.06.2008 (TW)

	return iRet;
}


long int TCGFcalc::CGActivCoefPT( double *X,double *param, double *act,
		   unsigned long int NN,   double Pbar, double T, double &roro )
{
	double *xtmp,*Fx;
	double P = Pbar/10.;
	xtmp = new double [NN];
	Fx = new double [NN];

	if(!paar)
		paar = new  EOSPARAM(X, param, NN);
	else
		paar->init( X, param, NN );

	double F0,Z,F1,fideal;
	double ro,delta = DELTA,ax,dx /*,tmp*/;
	long int i;

	norm(paar->XX0,paar->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());

	paar->ParamMix(xtmp);

	ro = ROTOTALMIX(P,T,paar);

	if( ro < 0.0 )  // Too low pressure, no corrections will be done
		return ( -1 );

	Z = P/(R*T*ro);
	F0 = FTOTALMIX(T,ro,paar);

	// fideal=log(R*T*ro/BARMPA);
	fideal = log(R*T*ro/0.1);
	ax = Z - 1.+fideal;

	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			copy(paar->XX0,xtmp,paar->NCmp());
			dx = xtmp[i]*delta;
			xtmp[i] += dx;
			norm(xtmp,paar->NCmp());
			paar->ParamMix(xtmp);
			F1 = FTOTALMIX(T,ro,paar)*(1.+dx);
			Fx[i] = (F1-F0)/(dx);
		}
		else Fx[i] = 0.;
	};

	// GMix=0.;
	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. && Fx[i]< 100. )
		{
			// tmp=log(paar.XX0[i]);
			// GMix+=tmp*paar.XX0[i];
			act[i] = exp(ax+Fx[i]);
		}
		else
		{
			act[i] = 0.;
		}
	};

	// GMix+=F0 + ax;
	// MLPutRealList(stdlink,act,paar.NCmp());
	delete[]xtmp;
	delete[]Fx;
	roro = ro;  // added, 21.06.2008 (TW)

	return 0;  // changed, 21.06.2008 (TW)
}


// Calculate departure functions through numerical derivative
long int TCGFcalc::CGDepartureFunct( double *X, double *param, double *param1, unsigned long int NN,
		double ro, double T )
{
	double F0, Z, F1;
	double delta = DELTA;
	double *xtmp = new double [NN];
	double Gdep, Sdep, Hdep, CPdep, vmix;

	if(!paar)
		paar = new  EOSPARAM(X, param, NN);
	else
		paar->init( X, param, NN );

	if(!paar1)
		paar1 = new  EOSPARAM(X, param1, NN);
	else
		paar1->init( X, param1, NN );

	norm(paar->XX0,paar->NCmp());
	norm(paar1->XX0,paar1->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());
	paar->ParamMix(xtmp);
	paar1->ParamMix(xtmp);
	Z = ZTOTALMIX(T,ro,paar);

	F0 = FTOTALMIX(T,ro,paar);
	// recalculate param1 for T+T*delta
	F1 = FTOTALMIX(T+T*delta,ro,paar1);
	// F1 = FTOTALMIX(T+T*delta,ro,paar);

	Sdep = - ( (F1-F0)/(delta*Tk)*Tk + F0 ) * R_CONST;	// corrected, 20.06.2008 (TW)
	Hdep = (F0*Tk*R_CONST + Tk*Sdep) + Z*R_CONST*Tk;
	Gdep = Hdep - Tk*Sdep;
	CPdep = 0.;
	vmix = Z*R_CONST*Tk/Pbar;

	// assignments
	DepPh[0] = Gdep;
	DepPh[1] = Sdep;
	DepPh[2] = Hdep;
	DepPh[3] = CPdep;
	DepPh[4] = vmix;

    delete [] xtmp;
    return 0;

}


// void ACTDENS(double *data,long nn, double *act )
long int TCGFcalc::CGActivCoefRhoT( double *X, double *param, double *act,
		unsigned long int NN, double ro, double T )
{
	double   F0,Z,F1,GMix,fideal;
	double delta = DELTA,ax,dx,tmp;
	long int i;
	double *Fx,*xtmp;
	xtmp = new double [NN];
	Fx = new double [NN];

	if(!paar)
		paar = new EOSPARAM(X, param, NN);
		else
			paar->init( X, param, NN );

	norm(paar->XX0,paar->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());
	paar->ParamMix(xtmp);
	Z = ZTOTALMIX(T,ro,paar);
	F0 = FTOTALMIX(T,ro,paar);
	fideal = log(R*T*ro/0.1);
	ax = Z - 1.+fideal;

	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			copy(paar->XX0,xtmp,NN);
			if ( xtmp[i]>DELTAMOLLIM )
			{
				dx = xtmp[i]*delta;
			}
			else
			{
				dx = DELTAMOLLIM*delta;
			}

			xtmp[i] += dx;
			norm(xtmp,paar->NCmp());
			paar->ParamMix(xtmp);
			F1 = FTOTALMIX(T,ro,paar)*(1.+dx);
			Fx[i] = (F1-F0)/(dx);

		}
		else Fx[i] = 0.;
	};

	GMix = 0.;
	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			tmp = log(paar->XX0[i]);
			GMix += tmp*paar->XX0[i];
			act[i] = exp(ax+Fx[i]);
		}
		else
		{
			act[i] = 0.;
		}
	};

	delete[]xtmp;
	delete[]Fx;
	return 0;
    // MLPutRealList(stdlink,act,paar.NCmp());
   };


double TCGFcalc::DIntegral( double T, double ro, unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
		{{-0.257431, 0.439229,  0.414783,  -0.457019, -0.145520,  0.299666},
		{-0.396724, 0.690721,  0.628935,  -0.652622, -0.201462, -0.23163 },
		{-0.488498, 0.863195,  0.761344,  -0.750086, -0.218562, -0.538463},
		{-0.556600, 0.995172,  0.852903,  -0.804710, -0.214736, -0.761700},
		{-0.611295, 1.103390,  0.921359,  -0.838804, -0.197999, -0.940714},
		{-0.657866, 1.196189,  0.975721,  -0.862346, -0.172526, -1.091678},
		{-0.698790, 1.278054,  1.020604,  -0.880027, -0.140749, -1.222733},
		{-0.735855, 1.351533,  1.058986,  -0.894024, -0.104174, -1.338626},
		{-0.769504, 1.418223,  1.092052,  -0.905347, -0.063730, -1.442391},
		{-0.800934, 1.479538,  1.121453,  -0.914864, -0.020150, -1.536070},
		{-0.829779, 1.535822,  1.147161,  -0.922381, 0.026157 , -1.621183},
		{-0.856655, 1.587957,  1.169885,  -0.928269, 0.074849 , -1.698853},
		{-0.881757, 1.636402,  1.190082,  -0.932668, 0.125590 , -1.769898},
		{-0.904998, 1.681421,  1.207610,  -0.935419, 0.178283 , -1.835070},
		{-0.926828, 1.723393,  1.223088,  -0.936667, 0.232649 , -1.894899},
		{-0.946773, 1.762571,  1.236007,  -0.936403, 0.288687 , -1.949858},
		{-0.965248, 1.799170,  1.246887,  -0.934650, 0.346207 , -2.000344}};

		// static double dt12[]=
		// {-2.139734,1.971553, 0.945513, -1.901492,-0.588630,-5.390941};
		// {-0.637684, 0.708107,  0.222086,  -0.481116, -0.332141, -3.492213};

	unsigned long int n;
	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		e = log(T);
		b = ro*ro;
		d = ro;
		c = ro*e;
		a = b*e;
	}

	// special case
	/*
	if ( IType==12 )
	{
		rez=(dt12[0]*T + dt12[1])*b +
		(dt12[2]*T + dt12[3])*ro + dt12[4]*T + dt12[5];
		return exp(rez);
	}
	*/

	n = IType-4;
	dtmp = data[n];
	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return exp(rez);
}


double TCGFcalc::LIntegral( double T, double ro,unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
	{{ -1.010391, 1.628552,  2.077476,  -2.30162 , -0.689931, -2.688117},
	{ -1.228611, 2.060090,  2.463396,  -2.453303, -0.573894, -3.350638},
	{ -1.354004, 2.402034,  2.718124,  -2.462814, -0.412252, -4.018632}};

	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}

	switch ( IType )
	{
		case 662:
			dtmp = data[0];
			break;
		case 1262:
			dtmp = data[1];
			break;
		case 12122:
			dtmp = data[2];
			break;
		default:
			return 0;
	}
	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return -exp(rez);
}


double TCGFcalc::KIntegral( double T, double ro,unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
	{{ -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720},
	{ -1.309550, 2.249120,  2.135877,  -2.278530, -0.773166, -3.704690},
	{ -1.490116, 2.619997,  2.404319,  -2.420706, -0.829466, -3.930928},
	{ -1.616385, 2.881007,  2.577600,  -2.484990, -0.828596, -4.175589},
	{ -1.940503, 3.552034,  2.940925,  -2.593808, -0.724353, -4.899975}};

	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}

	switch ( IType )
	{
		case 222333:
			dtmp = data[0];
			break;
		case 233344:
			dtmp = data[1];
			break;
		case 334445:
			dtmp = data[2];
			rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
			return -exp(rez);
		case 444555:
			dtmp = data[3];
			break;
		case 666777:
			dtmp = data[4];
			break;
		default:
			return 0;

	}

	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return exp(rez);
}


double TCGFcalc::K23_13( double T, double ro )
{
	static double TOld,roOld,KOLD;
	static double a,b,c,d,e;
	static double dtmp[]=
	{ -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720};

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}
	else return KOLD;

	KOLD = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	KOLD = exp(KOLD/3.);
	return KOLD;
}


double TCGFcalc::DENSITY( double *X,double *param, unsigned long NN ,double Pbar, double T )
{
	double P = Pbar * 0.1;
	double *xtmp;
	double ro;

	xtmp = new double [NN];
	if( !paar1 )
		paar1 = new EOSPARAM(X,param,NN);
	else
		paar1->init( X, param, NN );

	norm(paar1->XX0,paar1->NCmp());
	copy(paar1->XX0,xtmp,paar1->NCmp());
	paar1->ParamMix(xtmp);
	ro = ROTOTALMIX(P,T,paar1);

	delete [] xtmp;
	if( ro < 0. )
		Error( ""," Error - density cannot be found at this T,P" );
	return ro;
};


double TCGFcalc::PRESSURE( double *X,double *param,
		unsigned long int NN,double ro, double T )
{
	double *xtmp;
	xtmp = new double [NN];

	if( !paar1 )
		paar1 = new EOSPARAM(X,param,NN);
	else
		paar1->init( X, param, NN );

	norm(paar1->XX0,paar1->NCmp());
	copy(paar1->XX0,xtmp,paar1->NCmp());
	paar1->ParamMix(xtmp);
	double P = PTOTALMIX(T,ro,paar1);
	delete [] xtmp;
	return P*10.;
};


void TCGFcalc::copy( double* sours,double *dest,unsigned long int num )
{
	unsigned long int i;
	for ( i=0; i<num; i++)
	{
		dest[i]=sours[i];
	};
}


void TCGFcalc::norm( double *X,unsigned long int mNum )
{
	double tmp=0.;
	unsigned long int i;
	for ( i=0; i<mNum; i++ )
	{
		tmp += X[i];
	}
	tmp = 1./tmp;
	for ( i=0; i<mNum; i++ )
	{
		X[i] *= tmp;
	}
}


double TCGFcalc::RPA( double beta,double nuw )
{
	double fi1,fi2;
	fi1 = (1.20110+(0.064890+(-76.860+(562.686+(-2280.090+(6266.840+(-11753.40+(14053.8
			+(-9491.490 +2731.030*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw;
	fi2 = (0.588890+(-7.455360+(40.57590+(-104.8970+(60.25470+(390.6310+(-1193.080
			+(1576.350+(-1045.910+283.7580*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw*nuw;
	return  (-12.*fi1 + 192.*fi2*beta)*beta*beta/PI_1;
}


double TCGFcalc::dHS( double beta,double ro )
{
	// service constants
	double DV112 = 1./12.;
	double DV712 = 7./12.;
	// local variables
	double T12, T112, T712, B13, dB, delta, d;
	double a0, a1, a6, a3, a4, a7, a9, a12;
	double p0, p2, p6, p3, p5, p8, p11;
	double dbdl, ri6ro, ri6ro2, d3, d2, dnew, F0, F1;
	unsigned long int i;

	T12 = sqrt(beta);
	T112 = exp(DV112*log(beta));
	T712 = exp(DV712*log(beta));
	B13 = (1+beta);
	B13 = B13*B13*B13;

	dB = (P1*T112+PP2*T712+(P3+(P4+P5*beta)*beta)*beta)/B13;
	delta = (P6+P7*T12)/(1.+(P8+(P9+P10*T12)*T12)*T12);

	dbdl = dB*delta;
	ri6ro = PISIX*ro;
	ri6ro2 = ri6ro*ri6ro;

	a0 = dB+dbdl;
	a1 = -1.;
	a3 = (-1.5*dB -3.75*dbdl)*ri6ro;
	a4 = (1.5*ri6ro);
	a6 = (2.*dB + dbdl)*0.25*ri6ro2;
	a7 = -0.5*ri6ro2;
	a9 = -2.89325*ri6ro2*ri6ro*dbdl;
	a12 = -0.755*ri6ro2*ri6ro2*dbdl;

	p0 = -1.;
	p2 = a3*3.;
	p3 = a4*4.;
	p5 = a6*6.;
	p6 = a7*7.;
	p8 = a9*9.;
	p11 = a12*12.;

	d = dB;
	i = 0;

	while ( i++<21 )
	{
		d2 = d*d;
		d3 = d*d*d;
		F0 = a0+(a1+(a3+(a4+(a6+(a7+(a9+a12*d3)*d2)*d)*d2)*d)*d2)*d;
		F1 = p0+(p2+(p3+(p5+(p6+(p8+p11*d3)*d2)*d)*d2)*d)*d2;
		dnew = d-F0/F1;
		if ( fabs(dnew-d)<1.E-7 )
		{
			return dnew;
		}
		d = dnew;
	}

	if ( i>=20 )
	{
		return dB;
	}
	return dnew;
};


double TCGFcalc::FWCA( double T,double ro )
{
	static double TOld,roOld,F;
	double d,beta,nu,nuw;
	double nu1w1,nu1w2,nu1w3,nu1w4,nu1w5;
	double a0,a1,a2,a3;
	double I2;
	double I1_6,I1_12;
	double dW,dW12,dW6;
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	double F0,F1,FA;
	double rm,rmdw1,rmdw2,rmdw3,rmdw4,rmdw5;

	if ((T==TOld) && (ro==roOld))
	{
		return F;
	}
	else
	{
		TOld = T;
		roOld = ro;
	}

	rm = TWOPOW1SIX;
	beta = 1./T;
	d = dHS( beta, ro );
	tmp2 = PISIX*d*d*d;
	nu = tmp2*ro;
	tmp1 = (1. - nu/16.);
	nuw = nu*tmp1;
	dW = d*exp(1./3.*log(tmp1));

	nu1w1 = (1.-nuw);
	nu1w2 = nu1w1*nu1w1;
	nu1w3 = nu1w2*nu1w1;
	nu1w4 = nu1w2*nu1w2;
	nu1w5 = nu1w2*nu1w3;

	tmp1 = (1-nu);
	tmp1 = tmp1*tmp1;
	F0 = ((4.-3.*nu)*nu)/tmp1;

	a0 = fa0( nuw , nu1w2);
	a1 = fa1( nuw , nu1w3);
	a2 = fa2( nuw , nu1w4);
	a3 = fa3( nuw , nu1w5);

	I1_6 = fI1_6( nuw );
	I1_12 = fI1_12( nuw );

	rmdw1 = rm/dW;
	rmdw2 = rmdw1*rmdw1;
	rmdw3 = rmdw1*rmdw2;
	rmdw4 = rmdw2*rmdw2;
	rmdw5 = rmdw3*rmdw2;

	dW6 = dW*dW*dW;
	dW6 = 1./(dW6*dW6);
	dW12 = dW6*dW6;

	tmp1 = (a0/4.+ a1/12. + a2/24. + a3/24.)*dW6;
	tmp2 = (a0/10.+ a1/90. + a2/720. + a3/5040.)*(-dW12);
	tmp3 = (a0 - a1/3. + a2/12 - a3/60)/8.;
	tmp4 = (a0 - a1 + a2/2. - a3/6.)*rmdw2*(-9.)/40.;
	tmp5 = (a1 - a2 + a3/2)*rmdw3*(-2.)/9.;
	tmp6 = (a2 - a3)*rmdw4*(-9.)/64.;
	tmp7 = a3*(-3.)/35.*rmdw5;

	I2 = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

	F1 = 48.*nuw*(I1_12*dW12-I1_6*dW6 + I2)*beta;
	FA = RPA(beta,nuw);

	F = F0+F1+FA;

	return F;
}


double TCGFcalc::ZWCANum( double T,double ro )
{
	double delta = DELTA;
	double a0,a1;
	a1 = FWCA(T,ro*(1.+delta));
	a0 = FWCA(T,ro);
	return 1.+(a1-a0)/delta;
}


double TCGFcalc::UWCANum( double T,double ro )
{
	double delta = DELTA;
	double a0,a1,beta0,beta1;
	beta0 = 1./T;
	beta1 = beta0*(1.+delta);
	a1 = FWCA(1./beta1,ro);
	a0 = FWCA(T,ro);
	return (a1-a0)/(beta1-beta0);
}


double TCGFcalc::FDipPair( double T,double ro,double m2 )
{
	double kappa,Z,U,beta,F;
	kappa = m2*m2/(24.*T);
	beta = 1./T;
	Z = ZWCANum(T,ro);
	U = UWCANum(T,ro);
	F = kappa*(4.*beta*U-Z+1.);
	return F;
}


double TCGFcalc::J6LJ( double T,double ro )
{
	double kappa,Z,U,beta,F;
	beta = 1./T;
	Z = ZWCANum(T,ro);
	kappa = -16.*PI_1*ro*beta;
	U = UWCANum(T,ro);
	F = (4.*beta*U-Z+1.)/kappa;
	return F;
}


double TCGFcalc::FTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
{
	double FF,A0,A2,A3,AP,A1;
	// unsigned iall,inopol;
	double emix,s3mix,rotmp,T2R;
	double Jind,Jdp;
    long int /*itmp,jtmp,ktmp,*/ i,j,k;
    double s3tmp,mtmp,IK /*,atmp*/;
    double imtmp,jmtmp,iatmp,jatmp;
    double m2i,m2j,m2k;
    double s3tmpij,s3tmpik,s3tmpjk;
    double IKtmpij,IKtmpik,IKtmpjk;

    // iall=param.inonzero();
    // inopol=param.inonpolar();
    emix = param->EMIX();
    s3mix = param->S3MIX();

      rotmp = NA*ro_Real;
      T2R = T_Real*T_Real;

      A0 = FWCA(T_Real/emix,s3mix*rotmp);
     // if ( inopol< iall )
      {
        // dipole part
        A2 = 0.;
        for ( i=0; i<param->NCmp()-1; i++ )
        {
          for ( j=i+1; j<param->NCmp(); j++ )
          {
            s3tmp = param->MIXS3(i,j);
            Jdp = J6LJ(T_Real*s3tmp/param->MIXES3(i,j) , s3tmp*rotmp);
            A2 += param->M2R(i)*param->M2R(j)*Jdp*
                           param->X(i)*param->X(j)/s3tmp;
          }
        }
          A2 *= 2.;
          for ( i=0; i<param->NCmp(); i++ )
          {
            // itmp=param->ind(i);
            mtmp = param->M2R(i);
            s3tmp = param->SIG3(i);
            Jdp = J6LJ(T_Real/param->EPS(i),s3tmp*rotmp);
            A2 += mtmp*mtmp*Jdp*param->X(i)*param->X(i)/s3tmp;
          }

         A2 = -A2*TWOPI*rotmp/(3.*T2R);
         // A2 done

         if ( A2!=0. )
         {

          A3 = 0.;
          for ( i=0; i<param->NCmp(); i++ )
          {
            // itmp=param->ind(i);
            m2i = param->M2R(i);

            for ( j=0; j<param->NCmp(); j++  )
            {
             // jtmp=param->ind(j);
              m2j = param->M2R(j);

              s3tmpij = param->MIXS3(i,j);
              IKtmpij = K23_13(T_Real*s3tmpij/param->MIXES3(i,j),
                                                    s3tmpij*rotmp);
              for ( k=0; k<param->NCmp(); k++  )
              {
               // ktmp=param->ind(k);
               m2k = param->M2R(k);

               s3tmpik = param->MIXS3(i,k);
               s3tmpjk = param->MIXS3(j,k);

               IKtmpik = K23_13(T_Real*s3tmpik/param->MIXES3(i,k),
                                                    s3tmpik*rotmp);
               IKtmpjk = K23_13(T_Real*s3tmpjk/param->MIXES3(j,k),
                                                    s3tmpjk*rotmp);

               IK = IKtmpij*IKtmpik*IKtmpjk;
               A3 += m2i*m2j*m2k*IK*pow(s3tmpij*s3tmpik*s3tmpjk,-1./3.)*
               param->X(i)*param->X(j)*param->X(k);
              }
            }
          }
            A3 = A3*32.*sqrt(14.*PI_1/5.)*
                  rotmp*rotmp*PI_1*PI_1*PI_1/(135.*T_Real*T2R);
            AP = A2/(1. - A3/A2);
         }
         else AP = 0.;

        // induced interaction
        A1 = 0.;
        for ( i=0; i<param->NCmp(); i++ )
        {
         // itmp=param->ind(i);
          iatmp = param->A(i);
          imtmp = param->M2R(i);
          for ( j=0; j<param->NCmp(); j++ )
          {
            // jtmp=param->ind(j);
            jatmp = param->A(j);
            jmtmp = param->M2R(j);

            s3tmp = param->MIXS3(i,j);
            Jind = J6LJ(T_Real*s3tmp/param->MIXES3(i,j),s3tmp*rotmp);

           A1 += (iatmp*jmtmp + jatmp*imtmp)
                  *Jind*param->X(i)*param->X(j)/s3tmp;
          }
        }
        A1 = -A1*TWOPI*rotmp/T_Real;
// A1=-A1*FOURPI*rotmp/T_Real;
// A1=0.;

      }  // end of polar contribution

     FF = A0 + A1 + AP;
     //printf("%g %g %g %g %g",A0,A1,A2,A3,AP);
     //exit(1) ;

    return FF;
  }


double TCGFcalc::UTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
{
  double T /*,ro,s3 */;
  double delta = DELTA;
  double a0,a1,beta0,beta1,eps;
  eps = param->EMIX();
  T = T_Real/eps;

  beta0 = 1./T;
  beta1 = beta0*(1.+delta);
  a1 = FTOTALMIX((1./beta1)*eps,ro_Real,param);
  a0 = FTOTALMIX(T_Real,ro_Real,param);
  return (a1-a0)/(beta1-beta0);
 }


double TCGFcalc::ZTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
 {
  double delta = DELTA;
  double a0,a1;
  a1 = FTOTALMIX(T_Real,ro_Real*(1.+delta),param);
  a0 = FTOTALMIX(T_Real,ro_Real,param);

  return 1.+(a1-a0)/delta;
 }


double TCGFcalc::PTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
 {
  double Z;
    Z = ZTOTALMIX(T_Real,ro_Real,param);
    return Z*R*T_Real*ro_Real;
 }


// melting density
double TCGFcalc::Melt( double T )
 {

  return T*0.+.9;

 };


double TCGFcalc::Melt2(double T)
 {
  return T*0.+3.;

 };


 #define FIRSTSEED (15)
 #define ROMIN (1.E-2)
 #define NPOINT (5)


void  TCGFcalc::choose( double *pres, double P,unsigned long int &x1,unsigned long int &x2 )
 {
  unsigned long int i;
  double deltam = -10000000.,tmp;
  double deltap = 10000000.;

  for ( i=0; i<NPOINT; i++ )
  {
    tmp = P-pres[i];

    if ( tmp>0. )
    {
       if ( tmp<deltap )
       {
        deltap=tmp;
        x1 = i;
       }
    }
    else
    {
       if ( tmp>deltam )
       {
        deltam=tmp;
        x2 = i;
       }
    }
  }
     return ;

 }


double TCGFcalc::ROTOTALMIX( double P,double TT,EOSPARAM* param )
 {
     unsigned long int i;
     double T /*,ro*/;
     double fact, fact0, romax, dro, roarr[FIRSTSEED];
     double Ptmp[FIRSTSEED], ro0, ro1, rotest, PP0, PP1 /* ,Ptest */;
     double a,b;
     double inttofloat;
     double f[4],x[4],ff,dens[5],pres[5];
     unsigned long int x1,x2;
// double ptmp;

     T = TT/param->EMIX();
     fact0 = 1./(param->S3MIX()*NA);
     fact = R*TT*fact0;

     romax = Melt(T);
     inttofloat = FIRSTSEED-1;
     dro = (romax-ROMIN)/inttofloat;
     roarr[0] = ROMIN;
     roarr[1] = 2.*ROMIN;

     for ( i=2; i<FIRSTSEED; i++)
     {
       inttofloat = i;
       roarr[i] = ROMIN+inttofloat*dro;
     }

     for ( i=0; i<FIRSTSEED; i++)
     {
      Ptmp[i] = ZTOTALMIX(TT,roarr[i]*fact0,param);
      Ptmp[i] *= roarr[i] * fact;
      if ( Ptmp[i] > P )
      {
        break;
      }
     }

     if ( i==0 )  // Uses aproximation of ideal gas
     {
            return P/(R*TT);
     }

     // additional high pressure inteval
     if ( i==FIRSTSEED )
     {

     // roarr[0]=romax-0.0001;
     roarr[0] = roarr[FIRSTSEED-1];
     Ptmp[0] = Ptmp[FIRSTSEED-1];

     romax = Melt2(T);
     inttofloat = FIRSTSEED-1;
     dro = (romax-ROMIN)/inttofloat;
     for ( i=1; i<FIRSTSEED; i++)
     {
       inttofloat = i;
       roarr[i] = ROMIN+inttofloat*dro;
     }

     for ( i=1; i<FIRSTSEED; i++)
     {
      Ptmp[i] = ZTOTALMIX(TT,roarr[i]*fact0,param)*roarr[i]*fact;
      if ( Ptmp[i]>P )
      {
        break;
      }
     }

     if ( i==FIRSTSEED || i==0 )
     {
         printf("Input pressure is too high!\n");
//         exit(1);
         return (-1.0);
     }
     }

     ro0 = roarr[i-1];
     ro1 = roarr[i];
     PP0 = Ptmp[i-1];
     PP1 = Ptmp[i];
     i = 0;

   while ( i++<20 )
   {
     // Start interp
     ff = ro0;
     dens[0] = ro0;
     dens[1] = ro1;
     pres[0] = PP0;
     pres[1] = PP1;

     // first order
     x[0] = P-pres[0];
     f[0] = (dens[1]-dens[0])/(pres[1]-pres[0]);
     ff += f[0]*x[0];

     // second order
     dens[2] = ff;
     pres[2] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;

     if ( fabs(pres[2]-P)<1E-5 )
     {
       return ff*fact0;
     }

     x[1] = x[0]*(P-pres[1]);
     f[1] = (dens[2]-dens[1])/(pres[2]-pres[1]);

     f[0] = (f[1]-f[0])/(pres[2]-pres[0]);
     ff += f[0]*x[1];

     // third order
     dens[3] = ff;
     pres[3] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[3]-P)<1E-6 )
     {
      return ff*fact0;
     }
     x[2] = x[1]*(P-pres[2]);
     f[2] = (dens[3]-dens[2])/(pres[3]-pres[2]);
     f[1] = (f[2]-f[1])/(pres[3]-pres[1]);
     f[0] = (f[1]-f[0])/(pres[3]-pres[0]);
     ff += f[0]*x[2];
     dens[4] = ff;
     pres[4] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[4]-P)<1e-6 )
     {
      return ff*fact0;
     }

     choose(pres,P,x1,x2);

     ro0 = dens[x1];
     ro1 = dens[x2];
     PP0 = pres[x1];
     PP1 = pres[x2];

      if ( fabs((ro1-ro0))<0.001 )
      {
          a = (PP1-PP0)/(ro1-ro0);
          b = PP1-a*ro1;
          rotest = (P-b)/a;
          return rotest*(fact0);
      }
   }
        //return 10.;
         /// bad result

          a = (PP1-PP0)/(ro1-ro0);
          b = PP1-a*ro1;
          rotest = (P-b)/a;
          return rotest*(fact0);

 }


#ifndef IPMGEMPLUGIN
//--------------------------------------------------------------------//
// Calculates properties of pure fluids when called from RTParm
long int TCGFcalc::CGcalcFug( void )
{
	double T, P, Fugacity = 0.1, Volume = 0.0;
	double X[1] = {1.};
	double roro;  // added 21.06.2008 (TW)
	double Coeff[12];  // MAXEOSPARAM = 20;
	double Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
		Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;
	long int retCode = 0;

	ErrorIf( !aW.twp, "CG EoS", "Undefined twp");

	P = aW.twp->P;
	T = aW.twp->TC+273.15;

	for(long int ii=0; ii<12; ii++ )
		Coeff[ii] = aW.twp->Cemp[ii];

	// Calling CG EoS functions here
	if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1e-6 && P < 1e5 )
		retCode = CGFugacityPT( Coeff, Eos4parPT, Fugacity, Volume, P, T, roro );

	else
	{
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

	// increments to thermodynamic properties
	aW.twp->G += 8.31451 * T * log( Fugacity / P );
	aW.twp->V = Volume /* /10.  in J/bar */;
	// aW.twp->U = ((aW.twp->H/4.184)-RP*fg.VLK*fg.P2)*4.184;
	aW.twp->Fug = Fugacity;   /* fugacity at P */

	// passing corrected EoS coeffs to calculation of fluid mixtures
	// (possibly not required any more)
	aW.twp->wtW[6] = Eos4parPT[0];
	if( aW.twp->wtW[6] < 1. || aW.twp->wtW[6] > 10. )
		aW.twp->wtW[6] = 1.;  // foolproof temporary
	aW.twp->wtW[7] = Eos4parPT[1];
	aW.twp->wtW[8] = Eos4parPT[2];
	aW.twp->wtW[9] = Eos4parPT[3];

	// calculate departure functions and increment thermodynamic properties
	retCode = CGFugacityPT( Coeff, Eos4parPT1, Fugacity, Volume, P, T+T*DELTA, roro );
	CGDepartureFunct( X, Eos4parPT, Eos4parPT1, 1, roro, T );

	aW.twp->S +=  DepPh[1];
	aW.twp->H +=  DepPh[2];

    return retCode;
}

#endif



//=======================================================================================================
// Implementation of EOSPARAM class (used by TCGFcalc class)
//=======================================================================================================

void EOSPARAM::free()
{
	long int i;

	if ( NComp > 0)
	{
		for ( i=0;i<NComp;i++ )
			delete[]mixpar[i];
		delete[]mixpar;

		delete[]epspar;
		delete[]sig3par;
		delete[]XX;
		delete[]eps;
		delete[]eps05;
		delete[]sigpar;
		delete[]mpar;
		delete[]apar;
		delete[]aredpar;
		delete[]m2par;
		delete[]XX0;
		NComp = 0;
	}
}


void EOSPARAM::allocate()
{
	long int i;

	mixpar = new double*[NComp];
	for ( i=0; i<NComp; i++ )
		mixpar[i] = new double[NComp];

	epspar = new double[NComp];
	sig3par = new double[NComp];
	XX = new double[NComp];
	eps = new double[NComp];
	eps05 = new double[NComp];
	sigpar = new double[NComp];
	mpar = new double[NComp];
  	apar = new double[NComp];
	aredpar = new double[NComp];
	m2par = new double[NComp];
	XX0 = new double[NComp];
}


void EOSPARAM::init( double *Xinp, double * data, long int nn )
{
	long int i,j;
	double tmp;

	if( nn != NComp )
	{ // or error message
		free();
		NComp = nn;
		allocate();
	}

	for ( i=0;i<NComp;i++ )
	{
		XX0[i] = Xinp[i];

		sigpar[i] = data[i*4 ];
		eps[i] = data[i*4 + 1];
		mpar[i] = data[i*4 + 2];
		apar[i] = data[i*4 + 3];
	}
	for ( i=0; i<NComp; i++ )
	{
		tmp = sigpar[i];
		tmp = tmp*tmp*tmp;
		sig3par[i] = tmp;
		eps05[i] = sqrt(eps[i]);
		epspar[i] = tmp*eps[i];
		m2par[i] = mpar[i]*mpar[i]/(1.38048E-4);
		aredpar[i] = apar[i]/tmp;
	}

	// calculation of mixing properties
	for ( i=0; i<NComp-1; i++ )
	{
		for ( j=i+1; j<NComp; j++ )
		{
			tmp = (sigpar[i]+sigpar[j])*0.5;
			tmp = tmp*tmp*tmp;
			mixpar[i][j] = tmp;
			mixpar[j][i] = tmp*eps05[i]*eps05[j];
		}
	}
};


long int EOSPARAM::ParamMix( double *Xin )
  {
    long int j,i;
    double tmp,tmp1,tmp2;

    for ( i=0; i<NComp; i++ )
    	XX[i] = Xin[i];

    emix = 0.;
    s3mix = 0.;
    for ( i=0; i<NComp-1; i++ )
    {
      for ( j=i+1; j<NComp; j++ )
      {
          tmp = XX[i]*XX[j];
          tmp2 = mixpar[j][i];  //eps
          tmp1 = mixpar[i][j];  //signa
          s3mix += tmp1*tmp;
          emix += tmp2*tmp;
      }
    }
    s3mix *= 2.;
    emix *= 2.;
    for ( i=0; i<NComp; i++ )
    {
          tmp = XX[i]*XX[i];

          s3mix += sig3par[i]*tmp;
          emix += epspar[i]*tmp;
    }
    emix = emix/s3mix;
    return NComp;
  }



//=======================================================================================================
// Soave-Redlich-Kwong (SRK) model for fluid mixtures
// References: Soave (1972), Soave (1993)
// Implementation of the TSRKcalc class
//=======================================================================================================

// Constructor
TSRKcalc::TSRKcalc( long int NCmp, double Pp, double Tkp ):
    TSolMod( NCmp, 0, 0, 0, 0, 4, 'E',
         0, 0, 0, 0, 0, 0, Tkp, Pp, 0., 0. )
{
	aGEX = 0;
	aVol = 0;
	Pparc = 0;
	alloc_internal();
}


TSRKcalc::TSRKcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int *arIPx, double *arIPc, double *arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arPparc,
        double *arGEX, double *arVol, double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 4,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
{
	Pparc = arPparc;
	aGEX = arGEX;
	aVol = arVol;
	alloc_internal();
}


TSRKcalc::~TSRKcalc()
{
  free_internal();
}


// allocate work arrays for pure fluid and fluid mixture properties
void TSRKcalc::alloc_internal()
{
	Eosparm = new double [NComp][4];
	Pureparm = new double [NComp][4];
	Fugpure = new double [NComp][6];
	Fugci = new double [NComp][4];
	DepPh = new double [7];
	KK = new double *[NComp];
	dKK = new double *[NComp];
	d2KK = new double *[NComp];
	AA = new double *[NComp];

	for (long int i=0; i<NComp; i++)
	{
		KK[i] = new double[NComp];
		dKK[i] = new double[NComp];
		d2KK[i] = new double[NComp];
		AA[i] = new double[NComp];
	}
}


void TSRKcalc::free_internal()
{
	long int i;

	for (i=0; i<NComp; i++)
	{
		delete[]KK[i];
		delete[]dKK[i];
		delete[]d2KK[i];
		delete[]AA[i];
	}

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
	delete[]DepPh;
	delete[]KK;
	delete[]dKK;
	delete[]d2KK;
	delete[]AA;

}


// High-level method to retrieve pure fluid fugacities
long int TSRKcalc::PureSpecies()
{
	double Fugcoeff, Volume, DeltaH, DeltaS;
	long int j, retCode = 0;

	for( j=0; j<NComp; j++)
	{
		// Calling SRK EoS for pure fugacity
		retCode =  SRFugacityPT( j,  Pbar, Tk, aDCc+j*NP_DC,
				aDC[j], Fugcoeff, Volume, DeltaH, DeltaS );
		aGEX[j] = log( Fugcoeff );
		Pparc[j] = Fugcoeff * Pbar;  // Necessary only for performance
		aVol[j] = Volume * 10.;  // molar volume of pure fluid component, J/bar to cm3
	}  // j

	if ( retCode )
	{
		char buf[150];
		sprintf(buf, "SRK Fluid(): calculation of pure fugacity failed");
			Error( "E71IPM IPMgamma: ",  buf );
	}
	return 0;
}


// High-level method to calculate T,P corrected binary interaction parameters
long int TSRKcalc::PTparam()
{
	long int j, i, ip;
	long int i1, i2;
	double p0, p1, p2;
	double k, dk, d2k;

	PureSpecies();

	if( NPcoef > 0 )
	{
		// fill internal array of interaction parameters with standard value
		for( j=0; j<NComp; j++ )
		{
			for( i=0; i<NComp; i++ )
			{
				KK[j][i] = 0.;
				dKK[j][i] = 0.;
				d2KK[j][i] = 0.;
			}
		}

		// transfer those interaction parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			p0 = aIPc[NPcoef*ip];
			k = p0;
			KK[i1][i2] = k;
			KK[i2][i1] = k;   // symmetric case
		}
	}
	return 0;
}


// High-level method to retrieve activity coefficients of the fluid mixture
long int TSRKcalc::MixMod()
{
	long int j, iRet;

	iRet = FugacitySpec( Pparc );
	phVOL[0] = PhVol * 10.;

    for(j=0; j<NComp; j++)
    {
        if( Fugci[j][3] > 1e-23 )
        	lnGamma[j] = log( Fugci[j][3] );
        else
        	lnGamma[j] = 0;
    }
    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "SRK fluid(): calculation failed");
			Error( "E71IPM IPMgamma: ",  buf );
    }
    return iRet;
}


// High-level method to retrieve departure functions of the fluid mixture
long int TSRKcalc::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	// add excess property calculations
	long int iRet;

	iRet = DepartureFunct( Pparc );

    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "SRK fluid(): calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }

	// assignments
	Gex_ = DepPh[0];
	Sex_ = DepPh[1];
	Hex_ = DepPh[2];
	CPex_ = DepPh[3];
	Vex_ = DepPh[4];

	return iRet;

}


// High-level method to retrieve pure fluid properties
long int TSRKcalc::SRFugacityPT( long int i, double P, double Tk, double *EoSparam, double *Eos2parPT,
        double &Fugacity, double &Volume, double &DeltaH, double &DeltaS )
{
	long int iRet = 0;
	double Tcrit, Pcrit, omg, N;
	double apure, bpure, da, d2a;

	// reads EoS parameters from database into work array
	if( !EoSparam )
		return -1;  // Memory alloc error

	Eosparm[i][0] = EoSparam[0];  // critical temperature in K
	Eosparm[i][1] = EoSparam[1];  // critical pressure in bar
	Eosparm[i][2] = EoSparam[2];  // Pitzer acentric factor omega
	Eosparm[i][3] = EoSparam[3];  // empirical EoS parameter N
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	omg = Eosparm[i][2];
	N = Eosparm[i][3];

	AB(Tcrit, Pcrit, omg, N, apure, bpure, da, d2a);

	Pureparm[i][0] = apure;
	Pureparm[i][1] = bpure;
	Pureparm[i][2] = da;
	Pureparm[i][3] = d2a;
	Eos2parPT[0] = apure;
	Eos2parPT[1] = bpure;
	Eos2parPT[2] = da;
	Eos2parPT[3] = d2a;

	iRet = FugacityPure( i );
    if( iRet)
    	return iRet;

    Fugacity = Fugpure[i][0];  // Fugacity coefficient
    DeltaH = Fugpure[i][2];  // H departure function
    DeltaS = Fugpure[i][3];  // S departure function
    Volume = Fugpure[i][4];  //  J/bar

    return iRet;
}


// Calculates attractive (a) and repulsive (b) parameter of SRK equation of state
// and partial derivatives of alpha function
long int TSRKcalc::AB( double Tcrit, double Pcrit, double omg, double N,
		double &apure, double &bpure, double &da, double &d2a )
{
	double Tred, m, alph, ac, sqa, dsqa, d2sqa;

	Tred = Tk/Tcrit;
	m = 0.48 + 1.574*omg - 0.176*pow(omg,2.);
	alph = pow(1. + m*(1-sqrt(Tred)), 2.);
	ac = 0.42747*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit;
	apure = alph*ac;
	bpure = 0.08664*R_CONST*Tcrit/Pcrit;
	sqa = 1. + m*(1-sqrt(Tred));
	dsqa = -0.5*m/(sqrt(Tred)*Tcrit);
	da = 2.*ac*(sqa*dsqa);
	d2sqa = 0.25*m/(pow(Tred,1.5)*pow(Tcrit,2.));
	d2a = 2.*ac*(dsqa*dsqa + sqa*d2sqa);

	return 0;
}


// Calculates fugacities and departure functions of pure fluid species
long int TSRKcalc::FugacityPure( long int i )
{
	double Tcrit, Pcrit, Tred, asrk, bsrk, alph, da, d2a;
	double A, B, a2, a1, a0, z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, cpig;
	double fugpure, gdep, hdep, sdep, cpdep;
	double cv, dPdT, dPdV, dVdT;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONST*log(Pbar);
	gig = hig - Tk*sig;
	cpig = 0.;

	// retrieve a and b terms of cubic EoS
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	Tred = Tk/Tcrit;
	asrk = Pureparm[i][0];
	bsrk = Pureparm[i][1];
	da = Pureparm[i][2];
	d2a = Pureparm[i][3];

	// solve cubic equation
	A = asrk*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bsrk*Pbar/(R_CONST*Tk);
	a2 = (-1.);
	a1 = A-B-pow(B,2.);
	a0 = -A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;

	if (z1 > B)
		lnf1 = z1 - 1 - log(z1-B) - A/B*log(1.+B/z1);
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = z2 - 1 - log(z2-B) - A/B*log(1.+B/z2);
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = z3 - 1 - log(z3-B) - A/B*log(1.+B/z3);
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		z = z2; vol = vol2; lnf = lnf2;
	}
	else
	{
		z = z1; vol = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		z = z3; vol = vol3; lnf = lnf3;
	}
	else
	{
		z = z; vol = vol; lnf = lnf;
	}

	// calculate thermodynamic properties
	hdep = - ( 1 - z + 1./(bsrk*R_CONST*Tk) * (asrk-Tk*da) * log(1.+ bsrk/vol) )*R_CONST*Tk;
	sdep = ( log(z*(1.-bsrk/vol)) + 1./(bsrk*R_CONST) * da * log(1.+bsrk/vol) )*R_CONST;
	gdep = hdep - Tk*sdep;

	// heat capacity part
	cv = Tk*d2a/bsrk * log(1.+B/z);
	dPdT = R_CONST/(vol-bsrk) - da/(vol*(vol+bsrk));
	dPdV = - R_CONST*Tk/pow((vol-bsrk),2.) + asrk*(2.*vol+bsrk)/pow((vol*(vol+bsrk)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	cpdep = cv + Tk*dPdT*dVdT - R_CONST;

	// increment thermodynamic properties
	fugpure = exp(lnf);
	Fugpure[i][0] = fugpure;
	Fugpure[i][1] = gdep;
	Fugpure[i][2] = hdep;
	Fugpure[i][3] = sdep;
	Fugpure[i][4] = vol;
	Fugpure[i][5] = 0.;

	return 0;
}


// Cubic equation root solver based on Cardanos method
long int TSRKcalc::Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 )
{
	double q, rc, q3, rc2, theta, ac, bc;

	q = (pow(a2,2.) - 3.*a1)/9.;
	rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
	q3 = pow(q,3.);
	rc2 = pow(rc,2.);

	if (rc2 < q3)  // three real roots
	{
		theta = acos(rc/sqrt(q3));
		z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
		z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
		z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
	}

	else  // one real root
	{
		ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
		if (ac != 0.)
			bc = q/ac;
		else
			bc = 0.;
		z1 = ac+bc-a2/3.;
		z2 = ac+bc-a2/3.;
		z3 = ac+bc-a2/3.;
	}
	return 0;
}


// Calculates mixing properties of the fluid mixture
long int TSRKcalc::MixParam( double &amix, double &bmix )
{
	long int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			K = KK[i][j];
			AA[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}

	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + x[i]*x[j]*AA[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + x[i]*Pureparm[i][1];
	}
	return 0;
}


// Calculates fugacity of the bulk fluid mixture
long int TSRKcalc::FugacityMix( double amix, double bmix,
    double &fugmix, double &zmix, double &vmix )
{
	double A, B, a2, a1, a0, z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	A = amix*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bmix*Pbar/(R_CONST*Tk);
	a2 = (-1.);
	a1 = A-B-pow(B,2.);
	a0 = -A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;

	if (z1 > B)
		lnf1 = z1 - 1 - log(z1-B) - A/B*log(1.+B/z1);
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = z2 - 1 - log(z2-B) - A/B*log(1.+B/z2);
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = z3 - 1 - log(z3-B) - A/B*log(1.+B/z3);
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}

	fugmix = exp(lnf);
	PhVol = vmix;
	return 0;
}


// Calculates fugacities and activities of fluid species in the mixture,
long int TSRKcalc::FugacitySpec( double *fugpure )
{
	long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double A, B, bi, Bi, lnfci, fci;

	// Reload params to Pureparm (possibly not required any more)
	for( j=0; j<NComp; j++ )
	{
		Fugpure[j][0] = fugpure[j]/Pbar;
	}

	// calculate properties of the mixture
	iRet = MixParam( amix, bmix);
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix);
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		bi = Pureparm[i][1];
		Bi = bi*Pbar/(R_CONST*Tk);

		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + x[j]*AA[i][j];
		}

		lnfci = Bi/B*(zmix-1.) - log(zmix-B)
			+ A/B * ( Bi/B - 2./amix*sum ) * log(1.+B/zmix);
		fci = exp(lnfci);
		Fugci[i][0] = fci;  // fugacity coefficient using engineering convention
		Fugci[i][1] = x[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (x[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/x[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;
	}

	return iRet;
}


// calculates departure functions in the mixture
long int TSRKcalc::DepartureFunct( double *fugpure )
{
	long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0.;
	double A, B;
	double Gig, Hig, Sig, CPig, Gdep, Hdep, Sdep, CPdep;
	double K, dK, d2K, Q, dQ, d2Q;
	double damix, d2amix, ai, aj, dai, daj, d2ai, d2aj;
	double cv, dPdT, dPdV, dVdT;

	// Reload params to Pureparm (possibly not required any more)
	for( j=0; j<NComp; j++ )
	{
		Fugpure[j][0] = fugpure[j]/Pbar;
	}

	// calculate properties of the mixture
	iRet = MixParam( amix, bmix);
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix);
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// ideal gas changes from 1 bar to P (at T of interest)
	Hig = 0.;
	Sig = (-1.)*R_CONST*log(Pbar);
	Gig = Hig - Tk*Sig;
	CPig = 0.;

	// calculate total state functions of the mixture
	damix = 0.;
	d2amix = 0.;
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// pull parameters
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];
			K = KK[i][j];
			dK = dKK[i][j];
			d2K = d2KK[i][j];

			// increments to derivatives
			Q = sqrt(ai*aj);
			dQ = 0.5*( sqrt(aj/ai)*dai + sqrt(ai/aj)*daj );
			d2Q = 0.5*( dai*daj/sqrt(ai*aj) + d2ai*sqrt(aj)/sqrt(ai) + d2aj*sqrt(ai)/sqrt(aj)
					- 0.5*( pow(dai,2.)*sqrt(aj)/sqrt(pow(ai,3.))
					+ pow(daj,2.)*sqrt(ai)/sqrt(pow(aj,3.)) ) );
			damix = damix + x[i]*x[j] * ( dQ*(1.-K) - Q*dK );
			d2amix = d2amix + x[i]*x[j] * ( d2Q*(1.-K) - 2.*dQ*dK - Q*d2K );
		}
	}

	// calculate thermodynamic properties
	Hdep = - ( 1. - zmix + 1./(bmix*R_CONST*Tk) * (amix-Tk*damix )
			* log(1.+bmix/vmix) )*R_CONST*Tk;
	Sdep = ( log(zmix*(1.-bmix/vmix)) + 1./(bmix*R_CONST)*damix
			* log(1.+bmix/vmix) )*R_CONST;
	Gdep = Hdep - Tk*Sdep;

	// heat capacity part
	cv = Tk*d2amix/bmix * log(1.+B/zmix);
	dPdT = R_CONST/(vmix-bmix) - damix/(vmix*(vmix+bmix));
	dPdV = - R_CONST*Tk/pow((vmix-bmix),2.) + amix*(2.*vmix+bmix)/pow((vmix*(vmix+bmix)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	CPdep = cv + Tk*dPdT*dVdT - R_CONST;

	// assignments
	DepPh[0] = Gdep + Gig;
	DepPh[1] = Sdep + Sig;
	DepPh[2] = Hdep;
	DepPh[3] = CPdep;
	DepPh[4] = vmix;

	return iRet;
}


#ifndef IPMGEMPLUGIN
#include "s_tpwork.h"

// Calculates properties of pure fluids when called from RTParm
long int TSRKcalc::SRCalcFugPure( void )
{
	double T, P, Fugcoeff = 0.1, Volume = 0.0, DeltaH=0, DeltaS=0;
	double Coeff[7];  // MAXCRITPARAM
	double Eos2parPT[4] = { 0.0, 0.0, 0.0, 0.0 } ;
	long int retCode = 0;

	ErrorIf( !aW.twp, "SRK EoS", "Undefined twp");

	P = aW.twp->P;
	T = aW.twp->TC+273.15;
	for(long int ii=0; ii<7; ii++ )
		Coeff[ii] = aW.twp->CPg[ii];

	// Calling SRK EoS functions here
	if( T >= aW.twp->TClow +273.15 && T < 1e4 && P >= 1e-5 && P < 1e5 )
		retCode = SRFugacityPT( 0, P, T, Coeff, Eos2parPT, Fugcoeff, Volume, DeltaH, DeltaS );

	else
	{
		Fugcoeff = 1.;
		Volume = 8.31451*T/P;
		aW.twp->V = Volume;
		aW.twp->Fug = Fugcoeff*P;
		return retCode;
	}

	// increments to thermodynamic properties
	aW.twp->G += 8.31451 * T * log( Fugcoeff );  // from fugacity coefficient
	aW.twp->H +=  DeltaH;  // to be completed
	aW.twp->S +=  DeltaS;  // to be completed
	aW.twp->V = Volume;
	aW.twp->Fug = Fugcoeff * P;   // fugacity at P

	// passing corrected EoS coeffs to calculation of fluid mixtures
	aW.twp->wtW[6] = Eos2parPT[0];  // a
	aW.twp->wtW[7] = Eos2parPT[1];  // b
	aW.twp->wtW[8] = Eos2parPT[2];  // da
	aW.twp->wtW[9] = Eos2parPT[3];  // d2a

	return retCode;
}

#endif



//--------------------- End of s_fgl.cpp ---------------------------


