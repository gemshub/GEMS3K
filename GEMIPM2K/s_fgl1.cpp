//-------------------------------------------------------------------
// $Id: s_fgl1.cpp 1140 2008-12-08 19:07:05Z wagner $
//
// Copyright (C) 2008-2009  S.Dmitrieva, F.Hingerl, T.Wagner, D.Kulik
//
// Implementation of TSIT, TPitzer and TEUNIQUAC classes
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//

#include <math.h>
#include <stdio.h>
#include <iostream>
#include  <fstream>
using namespace std;
#include "verror.h"
#include "s_fgl.h"



//=============================================================================================
// SIT model (NEA version) reimplementation for aqueous electrolyte solutions
// References:
//=============================================================================================

// Generic constructor for the TSIT class
TSIT::TSIT( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arM, double *arZ,
        double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
{
  aZ = arZ;
  aM =	arM;
}


// Calculates activity coefficients and excess functions
long int TSIT::MixMod()
{
    long int j, icat, ian, /*ic, ia,*/  index1, index2, ip;
    double T, A, B, sqI, bgi=0, Z2, lgGam, SumSIT;
//    double nPolicy;

    I= IonicStr();
    if( I <  1e-6 /*TProfil::pm->pa.p.ICmin*/ )
    {
      for( j=0; j<NComp; j++)
    	  lnGamma[j] = 0.;
	  return 0;
    }
    T = Tk;
    A = 1.82483e6 * sqrt( RhoW ) / pow( T*EpsW, 1.5 );
    B = 50.2916 * sqrt( RhoW ) / sqrt( T*EpsW );

    sqI = sqrt( I );
    ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "SIT",
        "Error: A,B were not calculated - no values of RoW and EpsW !" );

    // Calculation of EDH equation
//  bgi = bg;
    ian= -1;
    icat = -1;
    for( j=0; j<NComp; j++ )
    {
// Determining the index of cation or anion
      if( aZ[j] < 0 )
          ian = j; // ian++;
      else if( aZ[j] > 0 )
          icat = j; // icat++;
//      else ;
      if( aZ[j] )
      {    // Charged species : calculation of the DH part
           Z2 = aZ[j]*aZ[j];
           lgGam = ( -A * sqI * Z2 ) / ( 1. + 1.5 * sqI );  // B * 4.562 = 1.5 at 25 C

// Calculation of SIT sum - new variant
           SumSIT = 0.;
           if( aZ[j] > 0 )  // cation
           {
              for( ip=0; ip<NPar; ip++ )
              {
                 index1 = aIPx[ip*MaxOrd];
                 if( index1 != icat )
                    continue;
                 index2 = aIPx[ip*MaxOrd+1];
                 SumSIT += aIPc[ip*NPcoef]   // epsilon
                        * aM[index2];
              }
           }
           else {   // anion
              for( ip=0; ip<NPar; ip++ )
              {
                 index2 = aIPx[ip*MaxOrd+1];
                 if( index2 != ian )
                    continue;
                 index1 = aIPx[ip*MaxOrd];  // index of cation
                 SumSIT += aIPc[ip*NPcoef]  // epsilon
                         * aM[index1];
              }
           }
           lgGam += SumSIT;
      }
      else { // Neutral species
         if( j != NComp-1 /*pmp->DCC[j] != DC_AQ_SOLVENT*/ ) // common salting-out coefficient ??
               lgGam = bgi * I;
            else // water-solvent - a0 - osmotic coefficient
               lgGam = 0.;
      }
      lnGamma[j] = lgGam * 2.302585093/*lg_to_ln*/;
    } // j

    return 0;
}


long int TSIT::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	// add excess property calculations

	// final assigments
	Gex_ = Gex;
	Vex_ = Vex;
	Hex_ = Hex;
	Sex_ = Sex;
	CPex_ = CPex;
	return 0;
}



//=============================================================================================
// Pitzer model for aqueous electrolyte solutions, Harvie-Moller-Weare (HMW) version
// References: Zhang et al. (2006)  Pitzer-Toughreact report
// Implemented by F.Hingerl as Matlab script,
//    converted by s.Dmytrieva into C++ program
//    in December 2008 for GEOTHERM CCES project
//=============================================================================================

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Generic constructor for the TPitzer class
TPitzer::TPitzer( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arM, double *arZ,
        double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW )
{
	aZ = arZ;
	aM = arM;

   // calculate sizes Nc, Na, Nn, Ns
   calcSizes();

  // reallocate internal arrays and set zeros
   alloc_internal();
}


TPitzer::~TPitzer()
{
  free_internal();
}


long int TPitzer::PTparam( )
{
	  // calculate vector of interaction parameters corrected to T,P of interest
		PTcalc( Tk );

	  // build conversion of species indexes between aqueous phase and Pitzer parameter tables
		setIndexes();

	  // put data from arIPx, arIPc to internal structure
		setValues();
    return 0;
}


// Calculation of activity coefficients
long int TPitzer::Pitzer_calc_Gamma( )
{
	long int M, N, X;

	//------------ Computing A- Factor
	  Aphi = A_Factor( Tk );
	//----------- Ionic Strength
	  Is = IonicStr( I );
	//----------- ------- F-Factor_______ Pitzer-Toughreact Report 2006 equation (A6)
	  Ffac = F_Factor( Aphi, I, Is );
	//----------- Z- Term________________ Pitzer-Toughreact Report 2006 equation (A8)
	  Zfac = Z_Term();

	lnGamma[Ns] = lnGammaH2O();

	for( M=0; M<Nc; M++ )
	{	lnGamma[xcx[M]] = lnGammaM( M );
//cout << "indexC = " <<	xcx[M] << " NComp " << NComp << endl;
	}

    for( X=0; X<Na; X++ )
    {	lnGamma[xax[X]] = lnGammaX( X );
//    cout << "indexA = " <<	xax[X] << " NComp " << lnGamma[xax[X]] << endl;
    	}
if( Nn > 0 )
    for( N=0; N<Nn; N++ )
    	lnGamma[xnx[N]] = lnGammaN( N );

// Pitzer_test_out( "test111.dat ");
    return 0;
}

// Moved macros here from s_fgl.h to restrict their visibility
// in other files (DK)
#define IPc( ii, jj )  ( aIPc[ (ii) * NPcoef + (jj) ])
#define IPx( ii, jj )  ( aIPx[ (ii) * MaxOrd + (jj) ])

#define mc( ii ) (aM[ xcx[(ii)] ])
#define ma( ii ) (aM[ xax[(ii)] ])
#define mn( ii ) (aM[ xnx[(ii)] ])
#define zc( ii ) (aZ[ xcx[(ii)] ])
#define za( ii ) (aZ[ xax[(ii)] ])

#define bet0( c,a ) ( abet0[ ((c)*Na+(a)) ])
#define bet1( c,a ) ( abet1[ ((c)*Na+(a)) ])
#define bet2( c,a ) ( abet2[ ((c)*Na+(a)) ])
#define Cphi( c,a ) ( aCphi[ ((c)*Na+(a)) ])
#define Lam( n,c )  ( aLam[ ((n)*Nc+(c)) ])
#define Lam1( n,a )  ( aLam1[ ((n)*Na+(a)) ])     // Bug fixed by SD 26.12.2008
#define Theta( c,c1 )  ( aTheta[ ((c)*Nc+(c1)) ])
#define Theta1( a,a1 ) ( aTheta1[ ((a)*Na+(a1)) ])

#define Psi( c,c1,a )  ( aPsi[(( (c) * Nc + (c1)  ) * Na + (a)) ])
#define Psi1( a,a1,c ) ( aPsi1[(( (a) * Na + (a1) ) * Nc + (c)) ])
#define Zeta( n,c,a )  ( aZeta[(( (n) * Nc + (c)  ) * Na + (a)) ])
//

// Output of test results into text file (standalone variant only)
void TPitzer::Pitzer_test_out( const char *path )
{

	long int ii, c, a, n;

	fstream ff(path, ios::out );
	ErrorIf( !ff.good() , path, "Fileopen error");

	ff << "Vector of interaction parameters corrected to T,P of interest" << endl;
	for( ii=0; ii<NPar; ii++ )
		ff << aIP[ii] << "  ";

	ff << endl << "list of indexes of Nc cations in aqueous phase" << endl;
	for( ii=0; ii<Nc; ii++ )
		ff << xcx[ii] << "  ";

	ff << endl << "list of indexes of Na anions in aq phase" << endl;
	for( ii=0; ii<Na; ii++ )
		ff << xax[ii] << "  ";

	ff << endl << "list of indexes of Nn neutral species in aq phase" << endl;
	for( ii=0; ii<Nn; ii++ )
		ff << xnx[ii] << "  ";

	ff << endl << "abet0" << endl;
	for( c=0; c<Nc; c++ )
	{	for( a=0; a<Na; a++ )
			ff << abet0[c*Na+a] << "  ";
		ff << endl;
	}

	ff << endl << "Theta" << endl;
	for( c=0; c<Nc; c++ )
	{	for( a=0; a<Na; a++ )
			ff << aTheta[c*Nc+a] << "  ";
		ff << endl;
	}

	ff << endl << "Lam1" << endl;
	for( n=0; n<Nn; n++ )
	{	for( a=0; a<Na; a++ )
			ff << Lam1( n,a ) << "  ";
		ff << endl;
	}

	ff << "\nAphi = " << Aphi << " I = " << I << " Is = " << Is <<
	       " Ffac = " << Ffac << " Zfac = " << Zfac  << endl;

	ff << endl << "ln activity coefficients of end members" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << lnGamma[ii] << "  ";
	ff << endl;

	ff << endl << "Activity coefficients of end members" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << exp(lnGamma[ii]) << "  ";

}


void TPitzer::alloc_internal()
{
	int i;
	// Input parameter arrays
	xcx = new long int[Nc];
	xax = new long int[Na];

	abet0 = new double[Nc*Na];
	abet1 = new double[Nc*Na];
	abet2 = new double[Nc*Na];
	aCphi = new double[Nc*Na];
	aTheta = new double[Nc*Nc];
	aTheta1 = new double[Na*Na];
    aPsi = new double[Nc*Nc*Na];
	aPsi1 = new double[Na*Na*Nc];
//	Cleaning allocated arrays
	for( i=0; i < Nc*Na; i++ )
	{
	   abet0[i] = 0.0;
	   abet1[i] = 0.0;
	   abet2[i] = 0.0;
	   aCphi[i] = 0.0;
	}
	for( i=0; i < Nc*Nc; i++ )
	  aTheta[i] = 0.0;
	for( i=0; i < Na*Na; i++ )
	  aTheta1[i] = 0.0;
	for( i=0; i < Nc*Nc*Na; i++ )
	  aPsi[i] = 0.0;
	for( i=0; i < Na*Na*Nc; i++ )
	  aPsi1[i] = 0.0;

	if( Nn > 0 )
	{
		xnx = new long int[Nn];
		aLam = new double[Nn*Nc];
		aLam1 = new double[Nn*Na];
		aZeta = new double[Nn*Nc*Na];
		for( i=0; i < Nn*Nc; i++ )
		  aLam[i] = 0.0;
		for( i=0; i < Nn*Na; i++ )
		  aLam1[i] = 0.0;
		for( i=0; i < Nn*Nc*Na; i++ )
		  aZeta[i] = 0.0;
	}
	else
	{
		xnx = 0;
		aLam = 0;
		aLam1 = 0;
		aZeta = 0;
	}
}

void TPitzer::free_internal()
{
	// Input parameter arrays
	if( xcx ) delete[] xcx;
	if( xax ) delete[] xax;

	if( abet0 ) delete[] abet0;
	if( abet1 ) delete[] abet1;
	if( abet2 ) delete[] abet2;
	if( aCphi ) delete[] aCphi;
	if( aTheta ) delete[] aTheta;
	if( aTheta1 ) delete[] aTheta1;
	if( aPsi ) delete[] aPsi;
	if( aPsi1 ) delete[] aPsi1;

	if( xnx ) delete[] xnx;
	if( aLam ) delete[] aLam;
	if( aLam1 ) delete[] aLam1;
	if( aZeta ) delete[] aZeta;

}


// Calculate temperature dependence of the interaction parameters.
// Two different versions: with 5 (4) or with 8 coefficients
void TPitzer::PTcalc( double T )
{
   long int ii;
   double Tr = 298.15;

   if( NPcoef == 5 ) // PHREEQPITZ and TOUGHREACT version
   {
	  for( ii=0; ii<NPar; ii++ )
       	aIP[ii] = IPc(ii,0) + IPc(ii,1)*(1./T-1./Tr) + IPc(ii,2)*log(T/Tr) +
       	          IPc(ii,3)*(T-Tr) + IPc(ii,4)*(T*T-Tr*Tr);
   }
   else
	if( NPcoef == 8 ) // original HMW version, also Felmy (GMIN) version
	{
	  for( ii=0; ii<NPar; ii++ )
	     aIP[ii] = IPc(ii,0) + IPc(ii,1)*T + IPc(ii,2)/T + IPc(ii,3)*log(T) + IPc(ii,4)/(T-263.) +
	       IPc(ii,5)*T*T + IPc(ii,6)/(680.-T) + IPc(ii,7)/(T-227.);
	}
	else Error( "", "PitzerHMW: Invalid number of coefficients to describe T dependence");
}


void TPitzer::calcSizes()
{
  long int jj;

  Nc = Na = Nn = 0;
  for( jj=0; jj<NComp-1; jj++ ) // -1 = no check H2O
  {
	 if( aZ[jj] > 0)
	   Nc++;
	 else
	  if( aZ[jj] < 0 )
		Na++;
	  else
		 Nn++;
  }
  Ns = NComp-1;  // index of water-solvent
}


// build conversion of species indexes between aq phase and Pitzer parameter tables
// list of indexes of Nc cations in aqueous phase
// list of indexes of Na anions in aq phase
// list of indexes of Nn neutral species in aq phase
void TPitzer::setIndexes()
{
  long int jj, ic, ia, in;

  ic = ia = in = 0;
  for( jj=0; jj<NComp-1; jj++ ) // -1 = no check H2O
  {
	 if( aZ[jj] > 0)
	 { xcx[ic] = jj;
	   ic++;
	 }
	 else
	  if( aZ[jj] < 0 )
	  {	  xax[ia] = jj;
		  ia++;
	  }
	  else
	  {  xnx[in] = jj;
		 in++;
	  }
  }
}


// put data from arIPx, arIPc, arDCc to internal structures
void TPitzer::	setValues()
{
  long int ii, ic, ia,in, i;

  for( ii=0; ii<NPar; ii++ )
  {
	switch( IPx(ii,3) )  // type of table
	{
	  case bet0_:
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Cation and anion index needed here"  );
	      bet0( ic, ia ) = aIP[ii];
//	      cout << "indexC = " << ic << "indexA = " << ia << " ind " << ((ic)*Na+(ia)) << endl;
	      break;
	  case bet1_:
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Cation and anion index needed here"  );
	      bet1( ic, ia ) = aIP[ii];
//	      cout << "indexC = " << ic << "indexA = " << ia << " ind " << ((ic)*Na+(ia)) << endl;
	      break;
	  case bet2_:
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Cation and anion indexes needed here"  );
	      bet2( ic, ia ) = aIP[ii];
//	      cout << "indexC = " << ic << "indexA = " << ia << " ind " << ((ic)*Na+(ia)) << endl;
	      break;
	  case Cphi_:
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Cation and anion indexes needed here"  );
	      Cphi( ic, ia ) = aIP[ii];
//	      cout << "indexC = " << ic << "indexA = " << ia << " ind " << ((ic)*Na+(ia)) << endl;
	      break;
	  case Lam_:
		  in = getIn( IPx(ii,0) );
	      if( in < 0 )
	      {   in = getIn( IPx(ii,1) );
		      ic = getIc( IPx(ii,0) );
	      }
	      else
	         ic = getIc( IPx(ii,1) );
	      ErrorIf( in<0||ic<0, "", "Cation and neutral species indexes needed here"  );
	      Lam( in, ic ) = aIP[ii];
//	      cout << "indexN = " << in << "indexC = " << ic << " ind " << ((in)*Nc+(ic)) << endl;
	      break;
	  case Lam1_:
		  in = getIn( IPx(ii,0) );
	      if( in < 0 )
	      {   in = getIn( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( in<0||ia<0, "", "Parameters must be anion and neutral species index"  );
	      Lam1( in, ia ) = aIP[ii];
//	      cout << "indexN = " << in << "indexA = " << ia << " ind " << ((in)*Na+(ia)) << endl;
	      break;
	  case Theta_:
	      ic = getIc( IPx(ii,0) );
          i = getIc( IPx(ii,1) );
	      ErrorIf( i<0||ic<0, "", "Only indexes of cations needed here"  );
	      Theta( ic, i ) = aIP[ii];
//	      cout << "indexC = " << ic << "indexC = " << i << " ind " << ((ic)*Nc+(i)) << endl;
	      break;
	  case Theta1_:
	      ia = getIa( IPx(ii,0) );
          i = getIa( IPx(ii,1) );
	      ErrorIf( i<0||ia<0, "", "Only indexes of anions needed here"  );
	      Theta1( ia, i ) = aIP[ii];
//	      cout << "indexA = " << ia << "indexA = " << i << " ind " << ((ia)*Na+(i)) << endl;
	      break;
	  case Psi_:
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );
		      ia = getIa( IPx(ii,0) );
		      i =  getIc( IPx(ii,2) );
	      }
	      else
	      {	 i =  getIc( IPx(ii,1) );
	         if( i<0 )
	         {
	        	 ia = getIa( IPx(ii,1) );
			     i =  getIc( IPx(ii,2) );
	         }
	         else
	          ia = getIa( IPx(ii,2) );
	      }
          ErrorIf( ic<0||ia<0||i<0, "", "Index of anion and 2 indexes of cations needed here"  );
	      Psi( ic, i, ia ) = aIP[ii];
//	      cout << "Psi = " << Psi( ic, i, ia ) << "index = " << (( (ic) * Nc + (i)  ) * Na + (ia)) << endl;
	      break;
	  case Psi1_:
		  ia = getIa( IPx(ii,0) );
	      if( ia < 0 )
	      {   ia = getIa( IPx(ii,1) );
		      ic = getIc( IPx(ii,0) );
		      i =  getIa( IPx(ii,2) );
	      }
	      else
	      {	 i =  getIa( IPx(ii,1) );
	         if( i<0 )
	         {
	        	 ic = getIc( IPx(ii,1) );
			     i =  getIa( IPx(ii,2) );
	         }
	         else
	          ic = getIc( IPx(ii,2) );
	      }
          ErrorIf( ic<0||ia<0||i<0, "", "Indexes of 2 anions and one cation needed here"  );
	      Psi1( ia, i, ic ) = aIP[ii];
//	      cout << "Psi1 = " << Psi1( ia, i, ic ) << "index = " << (( (ia) * Na + (i) ) * Nc + (ic)) << endl;
	      break;
	  case Zeta_:
		  in = getIn( IPx(ii,0) );
	      if( in < 0 )
	      {  in = getIn( IPx(ii,1) );
	         if( in < 0 )
	        	 in = getIn( IPx(ii,2) );
	      }
	      ic = getIc( IPx(ii,1) );
	      if( ic < 0 )
	      {  ic = getIc( IPx(ii,2) );
	         if( ic < 0 )
	        	 ic = getIc( IPx(ii,0) );
	      }
	      ia = getIa( IPx(ii,2) );
	      if( ia < 0 )
	      {  ia = getIa( IPx(ii,0) );
	         if( ia < 0 )
	        	 ia = getIa( IPx(ii,1) );
	      }
	      ErrorIf( ic<0||ia<0||in<0, "",
	    		  "Index of neutral species, index of cation and index of anion needed here"  );
	      Zeta( in, ic, ia ) = aIP[ii];
//	      cout << "Zeta = " << Zeta( in, ic, ia ) << "index = " << (( (in) * Nc + (ic)  ) * Na + (ia)) << endl;
	      break;
	}
  }
}


// Calculate Eta and Theta factors
// Reference: Anderson (2005), p. 610
void TPitzer::Ecalc( double z, double z1, double I, double Aphi,
		double& Etheta, double& Ethetap)
{
  double xMN, xMM, xNN,  x;
  double zet=0., dzdx=0.;
  double bk[23], dk[23];
  double JMN=0., JpMN=0., JMM=0., JpMM=0., JNN=0., JpNN=0.;

  // parameters for ak1 and ak2 values from Pitzer 1991 (p. 125, Table B1)
  static double ak1[21] = {  1.925154014814667, -0.060076477753119, -0.029779077456514,
                      -0.007299499690937,  0.000388260636404,  0.000636874599598,
                       0.000036583601823, -0.000045036975204, -0.000004537895710,
                       0.000002937706971,  0.000000396566462, -0.000000202099617,
                      -0.000000025267769,  0.000000013522610,  0.000000001229405,
                      -0.000000000821969, -0.000000000050847,  0.000000000046333,
                       0.000000000001943, -0.000000000002563, -0.000000000010991 };
					   // Prescribed constants
  static double ak2[23] = {  0.628023320520852,  0.462762985338493,  0.150044637187895,
                      -0.028796057604906, -0.036552745910311, -0.001668087945272,
                       0.006519840398744,  0.001130378079086, -0.000887171310131,
                      -0.000242107641309,  0.000087294451594,  0.000034682122751,
                      -0.000004583768938, -0.000003548684306, -0.000000250453880,
                       0.000000216991779,  0.000000080779570,  0.000000004558555,
                      -0.000000006944757, -0.000000002849257,  0.000000000237816,
                       0, 0 };  // Prescribed constants

  long int k, m;

  xMN= 6. * z*z1 * Aphi * pow(I,0.5);
  xMM= 6. * z* z * Aphi * pow(I,0.5);
  xNN= 6. *z1*z1 * Aphi * pow(I,0.5);

 for( k=1; k<=3; k++ )
 {
   if( k==1)
     x=xMN;
   else if( k==2 )
          x=xMM;
        else
          x=xNN;

  if( x<=1 )
  {   zet=4.0 * pow(x, 0.2) - 2.0;
      dzdx=0.8 * pow(x,(-0.8));

      bk[22]=0.; bk[21]=0.;
      dk[21]=0.; dk[22]=0.;
      for( m=20; m>=0; m--)
      {         bk[m]= zet * bk[m+1] - bk[m+2] + ak1[m];
                dk[m]= bk[m+1] + zet * dk[m+1]- dk[m+2];
      }
  }

   else
	  if( x>1)
	  {
          zet=-22.0/9.0 + (40.0/9.0) * pow(x,(-0.1));
          dzdx= (-40.0/90.) * pow(x,(-11./10.));

          bk[22]=0.; bk[21]=0.;
          dk[21]=0.; dk[22]=0.;
          for( m=20; m>=0; m--)
          {   bk[m] = zet *bk[m+1] - bk[m+2] + ak2[m];
              dk[m]=  bk[m+1] + zet *dk[m+1] - dk[m+2];
          }
	  }

   if( k==1 )
   {
    JMN=0.25*x -1. + 0.5* (bk[0]-bk[2]);
    JpMN=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
   }
  else
	if( k==2 )
	{
     JMM=0.25*x -1. + 0.5*(bk[0]-bk[2]);
     JpMM=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
	}
    else
    {
     JNN=0.25*x -1. +0.5*(bk[0]-bk[2]);
     JpNN=0.25 +0.5*dzdx*(dk[0]-dk[2]);
    }

  } //k
  Etheta=((z*z1) /(4.0*I)) * (JMN - 0.5*JMM - 0.5*JNN);
  Ethetap= - (Etheta/I) +((z*z1)/(8.0*I*I)) *(xMN*JpMN - 0.5*JpMM - 0.5*xNN*JpNN);
}


// Calculate Z-Term, Pitzer-Toughreact Report 2006, equation (A8)
double TPitzer::Z_Term()
{
    double Zan=0., Zca=0., Z;
	long int a, c;

	for( a=0; a<Na; a++)
	 Zan += za(a)*ma(a);

	for( c=0; c<Nc; c++)
	 Zca +=zc(c)*mc(c);

	Z = fabs(Zan) + Zca;
    return Z;
}


// Compute A-Factor
double TPitzer::A_Factor( double T )
{
	// Fix parameters
	double dens = 1.0;             // Density water
	double N0 = 6.0221415e23;      // Avogadro
	double eps = 78.38;            // Dielectricity constant solvent
	double k = 1.3806505e-23;      // Boltzmann
	double el = 1.60217635e-19;    // Coulomb charge unit
	double eps0 = 8.854187817e-12; // Dielectricity constant vacuum
	double pi = 3.141592654;
    double Aphi;

	//------------ Computing A- Factor
	Aphi = (1./3.) * pow((2.*pi*N0*dens*1000.),0.5) * pow((el*el)/(eps*4.*pi*eps0*k*T),1.5);

	return Aphi;
}


// Calculate Ionic Strength
double TPitzer::IonicStr( double& I )
{
    double Ia=0., Ic=0., IS;
	long int a, c;

	for( a=0; a<Na; a++ )
	  Ia += za(a)* za(a)* ma(a);

	for( c=0; c<Nc; c++ )
	  Ic += zc(c)* zc(c)* mc(c);

	IS =0.5*(Ia+Ic);
    I=IS;
	return pow(IS,0.5);
}



// Calculate osmotic coefficient, activity, and activity coefficient of water-solvent
double TPitzer::lnGammaH2O( )
{
    double Etheta=0., Ethetap=0.;
	long int a, c, n, c1, a1;
// Term OC1, Pitzer-Toughreact Report 2006, equation (A2)
	double OC1 = 2. * ( (-(Aphi*pow(I,1.5)) / (1.+1.2*Is) ));
// Term OC2
	double OC2=0., alp=0., alp1=0., C=0., h1=0., h2=0., B3=0.;
	for( c=0; c<Nc; c++)
	  for( a=0; a<Na; a++)
	  {
		 getAlp(  c,  a, alp, alp1 );
	     C = Cphi(c,a) / (2.*sqrt(fabs(za(a)*zc(c))));	// Pitzer-Toughreact Report 2006, equation (A7)
         h1=alp*Is;
	     h2=alp1*Is;
	     B3 = bet0(c,a)+ bet1(c,a)*exp(-h1)+(bet2(c,a)*exp(-h2)); //Pitzer-Toughreact Report 2006 equation (A9)
	     OC2 +=(mc(c)*ma(a)*(B3+Zfac*(C/(2.*sqrt(fabs(zc(c)*za(a)))))));
	  }
// Term OC3
	double OC3=0., z, z1, Phiphi;
	for( c=0; c<Nc; c++ )
	  for( c1=c+1; c1<Nc; c1++ )
	     for( a=0; a<Na; a++)
	     {
	    	 z=zc(c);
             z1=zc(c1);
	         Ecalc( z, z1, I, Aphi, Etheta,Ethetap);
	         Phiphi = Theta(c,c1) + Etheta + Ethetap * sqrt(I);	// Pitzer-Toughreact Report 2006, equation (A14)
	         OC3 += (mc(c)*mc(c1)*(Phiphi + (ma(a)*Psi(c,c1,a))));
	     }
// Term OC4
	double OC4=0., Phiphi1;
	for( a=0; a<Na; a++)
	  for( a1=a+1; a1<Na; a1++)
	    for( c=0; c<Nc; c++)
	    {  z=za(a);
	       z1=za(a1);
	       Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
	       Phiphi1 = Theta1(a,a1) + Etheta + Ethetap * sqrt(I);	// Pitzer-Toughreact Report, 2006 equation (A14)
	       OC4 += (ma(a)*ma(a1)*(Phiphi1+(mc(c)*Psi1(a,a1,c))));
	    }
// Term OC5
	double OC5, OC5a=0., OC5b=0.;
	for(  n=0; n<Nn; n++)
	  for( c=0; c<Nc; c++)
	     OC5a +=(mn(n)*mc(c)*Lam(n,c));

	for(  n=0; n<Nn; n++)
	  for( a=0; a<Na; a++)
	        OC5b +=(mn(n)*ma(a)*Lam1(n,a));
	OC5=OC5a+OC5b;
// Term OC6
	double OC6=0.;
	for(  n=0; n<Nn; n++)
	 for( c=0; c<Nc; c++)
	   for( a=0; a<Na; a++)
	        OC6 +=(mn(n)*mc(c)*ma(a)*Zeta(n,c,a));
// Addition of all sums
	double OCges=OC1+OC2+OC3+OC4+OC5+OC6;
// Summation of Molalities
	double   OCmol= sum(aM, xcx, Nc)+ sum(aM, xax, Na)+ sum(aM, xnx, Nn);
// Osmotic coefficient (OC)
	double OC = (1.+OCges) / OCmol;
// Activity of Water, Pitzer-Toughreact Report 2006, equation (A1)
	double Lna =(-18.1/1000.)*OC*OCmol;

	double activityH2O = exp(Lna);

//  lnGamma[Ns] = activityH2O/molefractionH2O;
	return Lna-log(x[Ns]);
}


void TPitzer::getAlp( long int c, long int a, double& alp, double& alp1 )
{
    if( zc(c) == 1. || za(a) == -1. )
    {    alp=2;
        alp1=12.;
    }
    else
      if( zc(c) == 2. && za(a) == -2. )
      {   alp=1.4;
          alp1=12.;
      }
      else
        if( zc(c) > 2. && za(a) <= -2. )
        { alp=2.0;
          alp1=50.;
        }
        else
          Error( " ", "alpha not defined");
}


// Calculate F-Factor, Pitzer-Toughreact Report 2006, equation (A6)
double TPitzer::F_Factor( double Aphi, double I, double Is )
{

  long int c, c1, a, a1;
  double z=0., z1=0., Etheta=0., Ethetap=0.;
// Term F1
  double F1=-Aphi*( (Is/(1.+1.2*Is)) + 2.*log(1.+1.2*Is)/1.2);
// Term F2
  double F2=0., Phip;
   for( c=0; c<Nc; c++ )
	 for( c1=c+1; c1<Nc; c1++ )
	 {
		 z=zc(c);
	     z1=zc(c1);
         Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
         Phip = Ethetap;					//Pitzer-Toughreact Report 2006, equation (A16)
         F2 +=(mc(c)*mc(c1)*(Phip));
	 }
// Term F3
  double F3=0., Phip1;
  for( a=0; a<Na; a++)
	 for( a1=a+1; a1<Na; a1++)
	 {  z=za(a);
        z1=za(a1);
        Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
        Phip1=Ethetap;      				//Pitzer-Toughreact Report 2006, equation (A16)
        F3 +=(ma(a)*ma(a1)*(Phip1));
	 }
// Term F4
  double F4=0., alp=0., alp1=0., h1, h2, g3, g4, B1;
  for( c=0; c<Nc; c++)
	 for( a=0; a<Na; a++)
	 {
		getAlp(  c,  a, alp, alp1 );
        h1=alp*Is;
        h2=alp1*Is;
        g3=(-2.*(1.-(1.+h1+(h1*h1)/2. )*exp(-h1)))/(h1*h1);	// Pitzer-Toughreact Report, 2006 equation (A13)
        g4=(-2.*(1.-(1.+h2+(h2*h2)/2 )*exp(-h2)))/(h2*h2);	// Pitzer-Toughreact Report, 2006 equation (A13)
        B1= (bet1(c,a)*g3)/I+ (bet2(c,a)*g4)/I;			// Pitzer-Toughreact Report 2006, equation (A12)
        F4 = F4+ (mc(c)*ma(a)*B1);
	 }
// Term F-Factor
   return F1+F2+F3+F4;
}


// Calculate lnGammaM - activity coefficient of a cation with index X
//
double TPitzer::lnGammaM(  long int M )
{
  double Etheta=0., Ethetap=0.;
  long int a, n, c1, a1;

// Calculate GM1, Pitzer-Toughreact Report 2006, equation (A3)

 double GM1=(zc(M)*zc(M))*Ffac;
// Term GM2
 double GM2=0., alp=0., alp1=0., h1, h2, g1,g2, B2, C;
 for( a=0; a<Na; a++)
 {
	 getAlp(  M,  a, alp, alp1 );
     C= Cphi(M,a)/(2.*sqrt(fabs(za(a)*zc(M))));	// Pitzer-Toughreact Report 2006, equation (A7)
     h1=alp*Is;
     h2=alp1*Is;
     g1=(2.*(1.-(1.+h1)*exp(-h1)))/(h1*h1);		  // Pitzer-Toughreact Report 2006, equation (A11)
     g2=(2.*(1.-(1.+h2)*exp(-h2)))/(h2*h1);		  // Pitzer-Toughreact Report 2006, equation (A11)
     B2= bet0(M,a)+(bet1(M,a)*g1)+ (bet2(M,a)*g2); // Pitzer-Toughreact Report 2006, equation (A10)
     GM2=GM2+(ma(a)*(2.*B2+Zfac*(C/(2.*sqrt(fabs(zc(M)*za(a)))))));
 }
// Term GM3
  double GM3=0., Phi, z, z1;
  for( c1=0; c1<Nc; c1++)
	 for( a=0; a<Na; a++)
	 {
	    if( M == c1)
	    {
	      Phi = 0.;
	      // Psi(M,c1,a) = 0.;
	    }
	    else
	    {  z=zc(M);
	       z1=zc(c1);
           Ecalc(z,z1,I,Aphi,Etheta,Ethetap);
           Phi=Theta(M,c1)+Etheta;  					// Pitzer-Toughreact Report 2006, equation (A15)

        GM3=GM3+mc(c1)*(2.*Phi+ ma(a)*Psi(M,c1,a)  );
	    } // Sveta 19/12/2008
	 }
// Term GM4
    double GM4=0.;
    for( a=0; a<Na; a++)
        for( a1=a+1; a1<Na; a1++)
            GM4=GM4+(ma(a)*ma(a1)*Psi1(a,a1,M));
// Term GM5
    double GM5a=0.;
    for( c1=0; c1<Nc; c1++)
     for( a=0; a<Na; a++)
     { C = Cphi(c1,a)/(2.*sqrt(fabs(za(a)*zc(c1))));			// Pitzer-Toughreact Report 2006, equation (A7)
       GM5a = GM5a+(mc(c1)*ma(a)* (C/(2.*sqrt(fabs(zc(M)*za(a)))))  );
     }
    double GM5=zc(M)*GM5a;
// Term GM6
   double GM6a=0;
   for( n=0; n<Nn; n++)
       GM6a += mn(n)*Lam(n,M);

   double GM6 = 2*GM6a;
// Term GM
   double GM=GM1+GM2+GM3+GM4+GM5+GM6;
   double actcoeffM = exp(GM);
   return GM;
}


// Calculate lnGammaX - activity coefficient of an anion with index X
double TPitzer::lnGammaX(  long int X )
{
  double Etheta=0., Ethetap=0.;
  long int c, n, c1, a1;

// Term GX1 (Pitzer-Toughreact Report 2006, equation A4)
  double   GX1=(za(X)*za(X))*Ffac;
// Term GX2
  double GX2=0., C=0., h1, h2, g1, g2, B2, alp, alp1;
  for( c=0; c<Nc; c++)
  {
 	 getAlp(  c,  X, alp, alp1 );
     C=Cphi(c,X)/(2.*sqrt(fabs(za(X)*zc(c))));    // Pitzer-Toughreact Report 2006, equation (A7)
     h1=alp*Is;
     h2=alp1*Is;
     g1=(2.*(1.-(1.+h1)*exp(-h1)))/(h1*h1);			 // Pitzer-Toughreact Report 2006, equation (A11)
     g2=(2.*(1.-(1.+h2)*exp(-h2)))/(h2*h2);			 // Pitzer-Toughreact Report 2006, equation (A11)
     B2= bet0(c,X)+bet1(c,X)*g1+ (bet2(c,X)*g2); 	// Pitzer-Toughreact Report 2006, equation (A10)
     GX2=GX2+(mc(c)*(2.*B2+Zfac*(C/(2.*sqrt(fabs(zc(c)*za(X)))))));
  }
// Term GX3
  double  GX3=0., z, z1, Phi1 ;
  for( a1=0; a1<Na; a1++)
	 for( c=0; c<Nc; c++)
	 {
	    if( X == a1)
	    {
	      Phi1 = 0.;
	      // Psi1(X,a1,c) = 0.; // Sveta 19/12/2008
	    }else
	      {
           z=za(X);
           z1=za(a1);
           Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
           Phi1=Theta1(X,a1)+Etheta; 			 // Pitzer-Toughreact Report 2006, equation (A15)

          GX3=GX3+ma(a1)*(2.*Phi1+mc(c)*Psi1(X,a1,c));  // Sveta 19/12/2008
	      }
	 }
// Term GX4
    double  GX4=0.;
    for( c=0; c<Nc; c++)
	   for( c1=c+1; c1<Nc; c1++)
              GX4=GX4+(mc(c)*mc(c1)*Psi(c,c1,X));
// Term GX5
     double GX5a=0.;
     for( c=0; c<Nc; c++)
  	   for( a1=0; a1<Na; a1++)
  	   {
          C =Cphi(c,a1)/(2.*sqrt(fabs(za(a1)*zc(c))));	 // Pitzer-Toughreact Report 2006, equation (A7)
          GX5a =GX5a+(mc(c)*ma(a1)* (C/(2.*sqrt(fabs(zc(c)*za(X))))) );
  	   }
      double GX5=fabs(za(X))*GX5a;
// Term GX6
     double GX6a=0.;
     for( n=0; n<Nn; n++)
         GX6a += mn(n)*Lam1(n,X);
     double GX6=2.*GX6a;
// Term GX
    double GX=GX1+GX2+GX3+GX4+GX5+GX6;
    double actcoeffX = exp(GX);
    return GX;
}


// Calculate lngammaN - activity coefficient of a neutral species with index N
double TPitzer::lnGammaN(  long int N )
{
  long int c, a;
// Term GN1, Pitzer-Toughreact Report 2006, equation (A5)
  double GN1=0.;
  for( a=0; a<Na; a++)
      GN1=GN1+(ma(a)*2.*Lam1(N,a));
// Term GN2
  double GN2=0.;
  for( c=0; c<Nc; c++)
      GN2=GN2+(mc(c)*2.*Lam(N,c));
// Term GN3
  double GN3=0.;
  for( c=0; c<Nc; c++)
  for( a=0; a<Na; a++)
      GN3=GN3+(mc(c)*ma(a)*Zeta(N,c,a));
// Term GN
  double GN=GN1+GN2+GN3;
  double actcoeffN = exp(GN);

  return GN;
}


long int TPitzer::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	// add excess property calculations

	// final assigments
	Gex_ = Gex;
	Vex_ = Vex;
	Hex_ = Hex;
	Sex_ = Sex;
	CPex_ = CPex;
	return 0;
}



//=============================================================================================
// Extended universal quasi-chemical (EUNIQUAC) model for aqueous electrolyte solutions
// References: Nicolaisen et al. (1993), Thomsen et al. (1996), Thomsen (2005)
//=============================================================================================


// Generic constructor for the TEUNIQUAC class
TEUNIQUAC::TEUNIQUAC( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, char Mod_Code,
        long int *arIPx, double *arIPc, double *arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arM, double *arZ,
        double T_k, double P_bar, double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, T_k, P_bar, dW, eW)
{
	alloc_internal();
	Z = arZ;
	M = arM;
}


TEUNIQUAC::~TEUNIQUAC()
{
	free_internal();
}


void TEUNIQUAC::alloc_internal()
{
	R = new double [NComp];
	Q = new double [NComp];
	Phi = new double [NComp];
	Theta = new double [NComp];
	U = new double *[NComp];
	dU = new double *[NComp];
	d2U = new double *[NComp];
	Psi = new double *[NComp];
	dPsi = new double *[NComp];
	d2Psi = new double *[NComp];
	Xi = new double *[NComp];

	for (long int j=0; j<NComp; j++)
	{
		U[j] = new double [NComp];
		dU[j] = new double [NComp];
		d2U[j] = new double [NComp];
		Psi[j] = new double [NComp];
		dPsi[j] = new double [NComp];
		d2Psi[j] = new double [NComp];
		Xi[j] = new double [NComp];
	}
}


void TEUNIQUAC::free_internal()
{
  	// cleaning memory
	for (long int j=0; j<NComp; j++)
	{
		delete[]U[j];
		delete[]dU[j];
		delete[]d2U[j];
		delete[]Psi[j];
		delete[]dPsi[j];
		delete[]d2Psi[j];
		delete[]Xi[j];
	}
	delete[]R;
	delete[]Q;
	delete[]Phi;
	delete[]Theta;
	delete[]U;
	delete[]dU;
	delete[]d2U;
	delete[]Psi;
	delete[]dPsi;
	delete[]d2Psi;
	delete[]Xi;
}


// Calculates T,P corrected binary interaction parameters
long int TEUNIQUAC::PTparam()
{
	long int j, i, ip, i1, i2;
	double u0, u1, u, du, d2u;
	double psi, dpsi, d2psi, v, dv;

    if ( NPcoef < 2 || NPar < 1 || NP_DC < 2 )
       return 1;

	// read and transfer species-dependent parameters
	for (j=0; j<NComp; j++)
	{
		R[j] = aDCc[NP_DC*j];   // volume parameter r
		Q[j] = aDCc[NP_DC*j+1];   // surface parameter q
	}

	// fill internal arrays of interaction parameters with standard value
	for (j=0; j<NComp; j++)
	{
		for (i=0; i<NComp; i++)
		{
			U[j][i] = 0.0;
			dU[j][i] = 0.0;
			d2U[j][i] = 0.0;
			Psi[j][i] = 1.0;
			dPsi[j][i] = 0.0;
			d2Psi[j][i] = 0.0;
		}
	}

	// read and convert interaction energies (uji) that have non-standard value
	for (ip=0; ip<NPar; ip++)
	{
		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		u0 = aIPc[NPcoef*ip+0];
		u1 = aIPc[NPcoef*ip+1];
		u = u0 + u1*(Tk-298.15);
		du = u1;
		d2u = 0.0;
		U[i1][i2] = u;
		dU[i1][i2] = du;
		d2U[i1][i2] = d2u;
		U[i2][i1] = u;
		dU[i2][i1] = du;
		d2U[i2][i1] = d2u;   // uij identical to uji
	}

	// calculate Psi and its partial derivatives
	for (j=0; j<NComp; j++)
	{
		for (i=0; i<NComp; i++)
		{
			psi = exp( -(U[j][i]-U[i][i])/Tk );
			v = (U[j][i]-U[i][i])/pow(Tk,2.) - (dU[j][i]-dU[i][i])/Tk;
			dv = (-2.)*(U[j][i]-U[i][i])/pow(Tk,3.) + 2.*(dU[j][i]-dU[i][i])/pow(Tk,2.)
					- (d2U[j][i]-d2U[i][i])/Tk;
			dpsi = psi * v;
			d2psi = dpsi*v + psi*dv;
			Psi[j][i] = psi;
			dPsi[j][i] = dpsi;
			d2Psi[j][i] = d2psi;
			Xi[j][i] = v;
		}
	}

	return 0;
}


// Calculates activity coefficients and excess functions
long int TEUNIQUAC::MixMod()
{
	int j, i, l, k, w;
	double Mw, Xw, IS, b, c;
	double A, RR, QQ, K, L, M;
	double gamDH, gamC, gamR, lnGam, Gam;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculation of DH parameters
	b = 1.5;
	c = 1.3287e-5;
	// A = c*sqrt(RhoW)/pow((EpsW*Tk),1.5);
	A = 1.131 + (1.335e-3)*(Tk-273.15) + (1.164e-5)*pow( (Tk-273.15), 2.);  // valid only for temperatures below 200 deg. C and Psat

	// calculation of ionic strength
	IS = 0.0;
	Mw = 0.018015;
	Xw = x[w];
	for (j=0; j<NComp; j++)
	{
		IS += 0.5*x[j]*pow(Z[j],2.)/(Xw*Mw);
	}

	// calculation of Phi and Theta terms
	for (j=0; j<NComp; j++)
	{
		RR = 0.0;
		QQ = 0.0;
		for (i=0; i<NComp; i++)
		{
			RR += x[i]*R[i];
			QQ += x[i]*Q[i];
		}
		Phi[j] = x[j]*R[j]/RR;
		Theta[j] = x[j]*Q[j]/QQ;
	}

	// loop over all species
	for (j=0; j<NComp; j++)
	{
		// species other than water solvent
		if (j < w)
		{
			K = 0.0;
			L = 0.0;
			for (k=0; k<NComp; k++)
			{
				M = 0.0;
				for (l=0; l<NComp; l++)
				{
					M += Theta[l]*Psi[l][k];
				}

				K += Theta[k]*Psi[k][j];
				L += Theta[k]*Psi[j][k]/M;
			}

			gamDH = - pow(Z[j],2.)*A*sqrt(IS)/(1.+b*sqrt(IS));
			gamC = log(Phi[j]/x[j]) - Phi[j]/x[j] - log(R[j]/R[w]) + R[j]/R[w]
					- 5.0*Q[j] * ( log(Phi[j]/Theta[j]) - Phi[j]/Theta[j]
					- log(R[j]*Q[w]/(R[w]*Q[j])) + R[j]*Q[w]/(R[w]*Q[j]) );
			gamR = Q[j] * ( - log(K) - L + log(Psi[w][j]) + Psi[j][w] );

			lnGam = gamDH + gamC + gamR;

			// convert activity coefficient to molality scale
			lnGam = lnGam + log(x[w]);
			lnGamma[j] = lnGam;

			// write debug results
			Gam = exp(lnGam);
			gammaDH[j] = gamDH;
			gammaC[j] = gamC;
			gammaR[j] = gamR;

		}

		// water solvent
		else
		{
			K = 0.0;
			L = 0.0;
			for (k=0; k<NComp; k++)
			{
				M = 0.0;
				for (l=0; l<NComp; l++)
				{
					M += Theta[l]*Psi[l][k];
				}

				K += Theta[k]*Psi[k][j];
				L += Theta[k]*Psi[j][k]/M;
			}

			gamDH = Mw*2.*A/pow(b,3.) * ( 1. + b*sqrt(IS) - 1./(1.+b*sqrt(IS)) - 2*log(1.+b*sqrt(IS)) );

			gamC = log(Phi[j]/x[j]) + 1. - Phi[j]/x[j] - 5.0*Q[j] * ( log(Phi[j]/Theta[j]) + 1. - Phi[j]/Theta[j] );

			gamR = Q[j] * (1. - log(K) - L );

			lnGam = gamDH + gamC + gamR;
			lnGamma[j] = lnGam;

			// write debug results
			Gam = exp(lnGam);
			gammaDH[j]=gamDH;
			gammaC[j] = gamC;
			gammaR[j] = gamR;
		}
	}

	return 0;
}


long int TEUNIQUAC::ExcessProp( double &Gex_, double &Vex_, double &Hex_, double &Sex_, double &CPex_ )
{
	// add excess property calculations
	int j, i, w;
	double Mw, Xw, IS, b, c;
	double A, dAdT, dAdP, d2AdT2;
	double phiti, phthi, RR, QQ, N, TPI, TPX, DTPX, CON;
	double gI, sI, gi, si;
	double gE, hE, sE, cpE, vE;
	double gDH, gC, gR, hR, cpR, gCI, gRI, gCX, gRX;   // DH, C and R contributions to properties
	double dg, d2g, dgRI, d2gRI, dgRX, d2gRX, dgDH, d2gDH;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculation of DH parameters
	b = 1.5;
	c = 1.3287e-5;
	// A = c*sqrt(RhoW)/pow((EpsW*Tk),1.5);
	// approximation valid only for temperatures below 200 deg. C and Psat
	A = 1.131 + (1.335e-3)*(Tk-273.15) + (1.164e-5)*pow( (Tk-273.15), 2.);
	dAdT = (1.335e-3) + 2.*(1.164e-5)*(Tk-273.15);
	d2AdT2 = 2.*(1.164e-5);

	// calculation of ionic strength
	IS = 0.0;
	Mw = 0.018015;
	Xw = x[w];
	for (j=0; j<NComp; j++)
	{
		IS += 0.5*x[j]*pow(Z[j],2.)/(Xw*Mw);
	}

	// calculation of Phi and Theta terms
	for (j=0; j<NComp; j++)
	{
		RR = 0.0;
		QQ = 0.0;
		for (i=0; i<NComp; i++)
		{
			RR += x[i]*R[i];
			QQ += x[i]*Q[i];
		}
		Phi[j] = x[j]*R[j]/RR;
		Theta[j] = x[j]*Q[j]/QQ;
	}

	// calculating bulk phase ideal mixing contributions
	gi = 0.0;
	si = 0.0;

	for (j=0; j<NComp; j++)
	{
		gi += x[j]*log(x[j]);
		si += x[j]*log(x[j]);
	}
	gI = R_CONST*Tk*gi;
	sI = - R_CONST*si;

	// calculation of bulk phase excess properties
	gE = 0.0;
	hE = 0.0;
	sE = 0.0;
	cpE = 0.0;
	vE = 0.0;
	gC = 0.0;
	gR = 0.0;
	hR = 0.0;
	cpR = 0.0;

	// infinite dilution part
	gCI = 0.;
	gRI = 0.;
	dgRI = 0.;
	d2gRI = 0.;

	for (j=0; j<NComp; j++)
	{
		phiti = R[j]/R[w];
		phthi = Q[j]/Q[w]*phiti;
		gCI += 1.0 + log(phiti) - phiti - 5.*Q[j]*(1.+log(phthi)-phthi);
		gRI += Q[j]*(1.-log(Psi[w][j])-Psi[j][w]);
		dgRI += - x[j]*Q[j] * (Xi[w][j] + dPsi[j][w]);
		d2gRI += - x[j]*Q[j] * (2.*(1./Tk)*(Xi[w][j]+dPsi[j][w])-Xi[j][w]*dPsi[j][w]);
	}

	// combinatorial and residual part
	gCX = 0.;
	gRX = 0.;
	dgRX = 0.;
	d2gRX = 0.;

	for (j=0; j<NComp; j++)
	{
		N = 0.0;
		TPI = 0.0;
		TPX = 0.0;
		DTPX = 0.0;
		for (i=0; i<NComp; i++)
		{
			N += Theta[i]*Psi[i][j];
			TPI += 1./N;
			TPX += Theta[i]*dPsi[i][j]*TPI;
			DTPX += - TPX*TPX + Theta[i]*d2Psi[i][j]*TPI;
		}
		gCX += x[j]*log(Phi[j]/x[j]) - 5.0*(Q[j]*x[j]*log(Phi[j]/Theta[j]));
		gRX += ( - x[j]*Q[j]*log(N) );
		dgRX += x[j]*Q[j]*TPX;
		d2gRX += x[j]*Q[j]*DTPX;
	}

	// DH part
	CON = - x[w]*Mw*4./pow(b,3.) * ( log(1.+b*sqrt(IS)) - b*sqrt(IS) + 0.5*pow(b,2.)*IS );
	gDH = CON*A;
	dgDH = - CON*dAdT;
	d2gDH = - CON*d2AdT2;
	// increment thermodynamic properties
	dg = ( dgDH + dgRX + dgRI );
	d2g = ( d2gDH + d2gRX + d2gRI );
	gE = ( gDH + gRX + gCX + gRI + gCI ) * R_CONST * Tk;
	hE = dg * pow(Tk,2.) * R_CONST;
	cpE = ( 2.*Tk*dg + pow(Tk,2.)*d2g ) * R_CONST;
	sE = (hE-gE)/Tk;

	// final assigments
	Gex_ = gE + gI;
	Vex_ = vE;
	Hex_ = hE;
	Sex_ = sE + sI;
	CPex_ = cpE;
	return 0;
}



// Output of test results into text file (standalone variant only)
void TEUNIQUAC::Euniquac_test_out( const char *path )
{

	long int ii, c, a, n;

//	const ios::open_mode OFSMODE = ios::out � ios::app;
	ofstream ff(path, ios::app );
	ErrorIf( !ff.good() , path, "Fileopen error");

	ff << endl << "Vector of interaction parameters corrected to T,P of interest" << endl;
	for( ii=0; ii<NPar; ii++ )
		ff << aIP[ii] << "  ";

	ff << endl << "Debye-H�ckel contribution to Activity Coefficients" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << gammaDH[ii] << "  ";

	ff << endl << "Contribution of Combinatorial Term to Activity Coefficients" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << gammaC[ii] << "  ";

	ff << endl << "Contribution of Residual Term to Activity Coefficients" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << gammaR[ii] << "  ";

/*	ff << endl << "abet0" << endl;
	for( c=0; c<Nc; c++ )
	{	for( a=0; a<Na; a++ )
			ff << abet0[c*Na+a] << "  ";
		ff << endl;
	}
*/
	ff << endl << "ln activity coefficients of end members" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << lnGamma[ii] << "  ";

	ff << endl << "Activity coefficients of end members" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << exp(lnGamma[ii]) << "  ";
	ff << endl;


}




//--------------------- End of s_fgl1.cpp ---------------------------

