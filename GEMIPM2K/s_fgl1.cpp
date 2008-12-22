//-------------------------------------------------------------------
// $Id: s_fgl1.cpp 1140 2008-12-08 19:07:05Z wagner $
//
// Copyright (C) 2008    S.Dmitrieva, F.Hingerl, D.Kulik, Th.Wagner
//
// Implementation of TSIT and TPitzer class
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
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arM, double *arZ,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
{
  aZ = arZ;
  aM =	arM;
//  PTparam();
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
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *aphVOL, double *arM, double *arZ,
        double dW, double eW ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 0,
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx,
        			 arlnGam, aphVOL, dW, eW )
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

    return 0;
}


// Output of test results into text file (standalone variant only)
void TPitzer::Pitzer_test_out( const char *path )
{

	long int ii, c, a;

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
	{	for( a=0; a<Nc; a++ )
			ff << aTheta[c*Nc+a] << "  ";
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
#define Lam1( n,a )  ( aLam[ ((n)*Na+(a)) ])
#define Theta( c,c1 )  ( aTheta[ ((c)*Nc+(c1)) ])
#define Theta1( a,a1 ) ( aTheta1[ ((a)*Na+(a1)) ])

#define Psi( c,c1,a )  ( aPsi[(( (c) * Nc + (c1)  ) * Na + (a)) ])
#define Psi1( a,a1,c ) ( aPsi1[(( (a) * Na + (a1) ) * Nc + (c)) ])
#define Zeta( n,c,a )  ( aZeta[(( (n) * Nc + (c)  ) * Na + (a)) ])
//

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


// Down here comes the Extended Uniquac model
//


//--------------------- End of s_fgl1.cpp ---------------------------

