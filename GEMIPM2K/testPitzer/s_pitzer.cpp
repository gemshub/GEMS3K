
#include <math.h>
#include <stdio.h>
#include <iostream>
#include  <fstream>
using namespace std;
#include "verror.h"
#include "s_pitzer.h"

    //--- ak1 and ak2 values from Pitzer 1991, p125, Table B1
	static double ak1[21] = {  1.925154014814667, -0.060076477753119, -0.029779077456514,
	                    -0.007299499690937,  0.000388260636404,  0.000636874599598,
	                     0.000036583601823, -0.000045036975204, -0.000004537895710,
	                     0.000002937706971,  0.000000396566462, -0.000000202099617,
	                    -0.000000025267769,  0.000000013522610,  0.000000001229405,
	                    -0.000000000821969, -0.000000000050847,  0.000000000046333,
	                     0.000000000001943, -0.000000000002563, -0.000000000010991 };  // Prescribed constants
	static double ak2[23] = {  0.628023320520852,  0.462762985338493,  0.150044637187895,
	                    -0.028796057604906, -0.036552745910311, -0.001668087945272,
	                     0.006519840398744,  0.001130378079086, -0.000887171310131,
	                    -0.000242107641309,  0.000087294451594,  0.000034682122751,
	                    -0.000004583768938, -0.000003548684306, -0.000000250453880,
	                     0.000000216991779,  0.000000080779570,  0.000000004558555,
	                    -0.000000006944757, -0.000000002849257,  0.000000000237816,
	                     0, 0 };  // Prescribed constants

//--------------------------------------------------------------------------------------------------------------
// Generic constructor for the TSolMod class
//
TSolMod::TSolMod( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *arM, double *arZ, 
        double dW, double eW, double iS ):
    ModCode(Mod_Code), NComp(NSpecies),  NPar(NParams), NPcoef(NPcoefs),
    MaxOrd(MaxOrder),  NP_DC(NPperDC), R_CONST(8.31451), RhoW(dW),
    EpsW(eW),  IonStr(iS), Tk(T_k), Pbar(P_bar)
        	
{
    // pointer assignment
    aIPx = arIPx;   // Direct access to index list and parameter coeff arrays!
    aIPc = arIPc;
    aIP = new double[ NPar ];
    aDCc = arDCc;
    x = arWx;
    aZ = arZ;
    aM =	arM;
    lnGamma = arlnGam;
}


TSolMod::~TSolMod()
{
   delete[] aIP;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Generic constructor for the TPitzer class
TPitzer::TPitzer( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
        long int NPperDC, double T_k, double P_bar, char Mod_Code,
        long int* arIPx, double* arIPc, double* arDCc,
        double *arWx, double *arlnGam, double *arM, double *arZ, 
        double dW, double eW, double iS ):
        	TSolMod( NSpecies, NParams, NPcoefs, MaxOrder, NPperDC, 
        			 T_k, P_bar, Mod_Code, arIPx, arIPc, arDCc, arWx, 
        			 arlnGam, arM, arZ, dW, eW, iS )    	
{
  // calculate sizes Nc, Na, Nn, Ns 
   calcSizes();	
	
  // realloc internal arrays and set zeros
   alloc_internal();

  // calc vector of interaction parameters corrected to T,P of interest
	P_T( Tk );

  // build conversion of species indexes between aq phase and Pitzer parameter tables
	setupIndexes();
  
  // put data from arIPx, arIPc to internal structure
	setupValues();	
}

TPitzer::~TPitzer()
{
  free_internal();
}

// Calculation of activity coefficients
long int TPitzer::Pitzer_calc_Gamma( )
{
	long int M, N, X; 

	//------------ Computing A- Factor
	  Aphi = A_Factor( Tk );
	//----------- Ionic Strength
	  Is = IonicStr( I );
	//----------- ------- F-Factor________________ Pitzer-Toughreact Report 2006 equation (A6)
	  Ffac = F_Factor( Aphi, I, Is );
	//----------- Z- Term________________ Pitzer-Toughreact Report 2006 equation (A8)
	  Zfac = Z_Term();
	
	lnGamma[Ns] = lnGammaH2O();

	for( M=0; M<Nc; M++ )
    	lnGamma[xcx[M]] = lnGammaM( M );	
	
    for( X=0; X<Na; X++ )
    	lnGamma[xax[X]] = lnGammaX( X );	
	
    for( N=0; N<Nn; N++ )
    	lnGamma[xnx[N]] = lnGammaN( N );
    
    return 0;
}

// Out test results for file
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
	{	for( a=0; a<Nc; a++ )
			ff << aTheta[c*Nc+a] << "  ";
		ff << endl;
	} 

	ff << "\nAphi = " << Aphi << " I = " << I << " Is = " << Is <<
	       " Ffac = " << Ffac << " Zfac = " << Zfac  << endl;
	
	ff << endl << "Pointer to ln activity coefficients of end members" << endl;
	for( ii=0; ii<NComp; ii++ )
		ff << lnGamma[ii] << "  ";
	
	ff << endl;

}


void TPitzer::alloc_internal()
{
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
	
	if( Nn > 0 )
	{
		xnx = new long int[Nn]; 
		aLam = new double[Nn*Nc];
		aLam1 = new double[Nn*Na];
		aZeta = new double[Nn*Nc*Na];
	}
	else
	{
		xnx = 0; 
		aLam = 0;
		aLam1 = 0;
		aZeta = 0;
	}
		
	/* Work parameter arrays
	B1 = new double[Nc*Na];
    B2 = new double[Nc*Na];
	B3 = new double[Nc*Na];
	Phi1 = new double[Na*Na];
	Phi2 = new double[Nc*Nc];
	Phi3 = new double[Na*Na];
	C = new double[Nc*Na];
	Etheta = new double[Nc*Nc];
	Ethetap = new double[Na*Na];
	*/
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
	
	/* Work parameter arrays
	if( B1 ) delete[] B1;
	if( B2 ) delete[] B2;
	if( B3 ) delete[] B3;
	if( Phi1 ) delete[] Phi1;
	if( Phi2 ) delete[] Phi2;
	if( Phi3 ) delete[] Phi3;
	if( C ) delete[] C;
	if( Etheta ) delete[] Etheta;
	if( Ethetap ) delete[] Ethetap;
	*/
}

// Calculate temperature dependence of the interaction parameters.
// There are different versions
void TPitzer::P_T( double T )
{
   long int ii;
   double Tr = 298.15; 
      
   if( NPcoef == 5 ) // PHREEQPITZ version
   {
	  for( ii=0; ii<NPar; ii++ )
       	aIP[ii] = IPc(ii,0) + IPc(ii,1)*(1/T-1/Tr) + IPc(ii,2)*log(T/Tr) + 
       	          IPc(ii,3)*(T-Tr) + IPc(ii,4)*(T*T-Tr*Tr);
   } 
   else
	if( NPcoef == 8 ) // original HMW version
	{
	  for( ii=0; ii<NPar; ii++ )
	     aIP[ii] = IPc(ii,0) + IPc(ii,1)*T + IPc(ii,2)/T + IPc(ii,3)*log(T) + IPc(ii,4)/(T-263) +
	       IPc(ii,5)*T*T + IPc(ii,6)/(680-T) + IPc(ii,7)/(T-227);
	}	    
	else Error( "", "Illegal equations to discribe T dependence");
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
  Ns = NComp-1;
}

// build conversion of species indexes between aq phase and Pitzer parameter tables
// list of indexes of Nc cations in aqueous phase
// list of indexes of Na anions in aq phase
// list of indexes of Nn neutral species in aq phase
void TPitzer::setupIndexes()
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

// put data from arIPx, arIPc, arDCc to internal structure
void TPitzer::	setupValues()
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
	      ErrorIf( ia<0||ic<0, "", "Parameters must be cation and anion index"  );
	      bet0( ic, ia ) = aIP[ii];
	      break;
	  case bet1_: 
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );  
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Parameters must be cation and anion index"  );
	      bet1( ic, ia ) = aIP[ii];
	      break;
	  case bet2_: 
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );  
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Parameters must be cation and anion index"  );
	      bet2( ic, ia ) = aIP[ii];
	      break;
	  case Cphi_: 
		  ic = getIc( IPx(ii,0) );
	      if( ic < 0 )
	      {   ic = getIc( IPx(ii,1) );  
		      ia = getIa( IPx(ii,0) );
	      }
	      else
	         ia = getIa( IPx(ii,1) );
	      ErrorIf( ia<0||ic<0, "", "Parameters must be cation and anion index"  );
	      Cphi( ic, ia ) = aIP[ii];
	      break;
	  case Lam_: 
		  in = getIn( IPx(ii,0) );
	      if( in < 0 )
	      {   in = getIn( IPx(ii,1) );  
		      ic = getIc( IPx(ii,0) );
	      }
	      else
	         ic = getIc( IPx(ii,1) );
	      ErrorIf( in<0||ic<0, "", "Parameters must be cation and neutral species index"  );
	      Lam( in, ic ) = aIP[ii];
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
	      break;
	  case Theta_: 
	      ic = getIc( IPx(ii,0) );
          i = getIc( IPx(ii,1) );
	      ErrorIf( i<0||ic<0, "", "Parameters must be cations"  );
	      Theta( ic, i ) = aIP[ii];
	      break;
	  case Theta1_: 
	      ia = getIa( IPx(ii,0) );
          i = getIa( IPx(ii,1) );
	      ErrorIf( i<0||ia<0, "", "Parameters must be anions"  );
	      Theta1( ia, i ) = aIP[ii];
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
          ErrorIf( ic<0||ia<0||i<0, "", "Parameters must be anion and 2 cations index"  );
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
         ErrorIf( ic<0||ia<0||i<0, "", "Parameters must be 2 anions and cation index"  );
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
	      ErrorIf( ic<0||ia<0||in<0, "", "Parameters must be  neutral species, cation  and anion index"  );
	      Zeta( in, ic, ia ) = aIP[ii];
	      break;
	}
  } 
	
}

//-------------------------------------------------------------------------
//----------------------Computation Eta and Theta Factor-------------------
//Formulation from G.ANDERSON:THERMODYNAMICS OF NATURAL SYSTEMS,2005; pp. 610
//-------------------------------------------------------------------------
void TPitzer::Ecalc( double z, double z1, double I, double Aphi, 
		double& Etheta, double& Ethetap)
{
  double xMN, xMM, xNN,  x;
  double zet, dzdx;
  double bk[23], dk[23];
  double JMN, JpMN, JMM, JpMM, JNN, JpNN;
  long int k, m;
  
  xMN= 6 * z*z1 * Aphi * pow(I,1/2);
  xMM= 6 * z* z * Aphi * pow(I,1/2);
  xNN= 6 *z1*z1 * Aphi * pow(I,1/2);

//-----------------------------------------------
 for( k=1; k<=3; k++ )
 {  
   if( k==1)
     x=xMN; 
   else if( k==2 )
          x=xMM;
        else 
          x=xNN;
//---------------------------------     
         
  if( x<=1 )
  {   zet=4.0 * pow(x, 0.2) - 2.0;
      dzdx=0.8 * pow(x,(-0.8));  

      bk[22]=0; bk[21]=0;
      dk[21]=0; dk[22]=0;
      for( m=20; m>=0; m--)     
      {         bk[m]= zet * bk[m+1] - bk[m+2] + ak1[m];
                dk[m]= bk[m+1] + zet * dk[m+1]- dk[m+2];
      }
  }   
//-------------
   else 
	  if( x>1)
	  {
          zet=-22.0/9.0 + (40.0/9.0) * pow(x,(-0.1)); 
          dzdx= (-40.0/90.) * pow(x,(-11./10.));

          bk[22]=0; bk[21]=0;
          dk[21]=0; dk[22]=0;
          for( m=20; m>=0; m--)     
          {   bk[m] = zet *bk[m+1] - bk[m+2] + ak2[m];
              dk[m]=  bk[m+1] + zet *dk[m+1] - dk[m+2];
          }    
	  }
//-----------------------------------------------
   if( k==1 )
   {
    JMN=0.25*x -1. + 0.5* (bk[0]-bk[2]);
    JpMN=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
   } 
  else
	if( k==2 )
	{
     JMM=0.25*x -1. + 0.5*(bk[0]-dk[2]);
     JpMM=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
	} 
    else
    {
     JNN=0.25*x -1. +0.5*(bk[0]-bk[2]);
     JpNN=0.25 +0.5*dzdx*(dk[0]-dk[2]);
    }
//---------------------------------------------------
  } //k 
  Etheta=(2.0 /(4.0*I)) * (JMN - 0.5*JMM - 0.5*JNN);
  Ethetap= - (Etheta/I) +(2.0/(8.0*I*I)) *(xMN*JpMN - 0.5*JpMM - 0.5*xNN*JpNN);
}

//----------- Z- Term________________ Pitzer-Toughreact Report 2006 equation (A8)
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

//------------ Computing A- Factor
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
    double A, Aphi;

	//------------ Computing A- Factor
	A = (1./3.) * pow((2*pi*N0*dens),0.5) * pow((el*el)/(eps*4*pi*eps0*k*T),1.5);

	// Agamma = ln(10)A; Aphi = Agamma/3
	Aphi=(log(10.)*(A))/3.;  
    return Aphi;
}    
	
//----------- Ionic Strength
double TPitzer::IonicStr( double& I )
{
    double Ia=0., Ic=0.;
	long int a, c;

	for( a=0; a<Na; a++ )
	  Ia += za(a)*za(a)*ma(a);

	for( c=0; c<Nc;c++ )
	  Ic += zc(c)*zc(c)*mc(c);

	I=0.5*(Ia+Ic);

	return pow(I,0.5);
}

//-----------------------------------------------------------------------
//------------------------CALCULATIONS-----------------------------------
//-------------Osmotic Coefficient and activity of water-----------------
//------------------------------------------------------------------------
double TPitzer::lnGammaH2O( )
{
    double Etheta, Ethetap;
	long int a, c, n, c1, a1;
//----------- OC1________________ Pitzer-Toughreact Report 2006 equation (A2)
	double OC1 = 2 * ( (-(Aphi*pow(I,1.5)) / (1.+1.2*Is) ));
//----------- OC2
	double OC2=0., alp, alp1, C, h1, h2, B3;
	for( c=0; c<Nc; c++)
	  for( a=0; a<Na; a++)
	  {    
		 getAlp(  c,  a, alp, alp1 );  
	     C = Cphi(c,a) / (2.*sqrt(fabs(za(a)*zc(c))));	// Pitzer-Toughreact Report 2006 equation (A7)
         h1=alp*Is;
	     h2=alp1*Is;
	     B3 = bet0(c,a)+ bet1(c,a)*exp(-h1)+(bet2(c,a)*exp(-h2)); //Pitzer-Toughreact Report 2006 equation (A9)
	     OC2 +=(mc(c)*ma(a)*(B3+Zfac*(C/(2.*sqrt(fabs(zc(c)*za(a)))))));
	  } 
//----------- OC3
	double OC3=0., z, z1, Phiphi;
	for( c=0; c<Nc; c++ )
	  for( c1=c+1; c1<Nc; c1++ )
	     for( a=0; a<Na; a++)
	     { 
	    	 z=zc(c);
             z1=zc(c1);
	         Ecalc( z, z1, I, Aphi, Etheta,Ethetap);
	         Phiphi = Theta(c,c1) + Etheta + Ethetap * sqrt(I);	//Pitzer-Toughreact Report 2006 equation (A14)
	         OC3 += (mc(c)*mc(c1)*(Phiphi + (ma(a)*Psi(c,c1,a))));
	     }
//----------- OC4
	double OC4=0., Phiphi1;
	for( a=0; a<Na; a++)
	  for( a1=a+1; a1<Na; a1++)
	    for( c=0; c<Nc; c++)
	    {  z=za(a);
	       z1=za(a1);
	       Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
	       Phiphi1 = Theta1(a,a1) + Etheta + Ethetap * sqrt(I);	//Pitzer-Toughreact Report 2006 equation (A14)
	       OC4 += (ma(a)*ma(a1)*(Phiphi1+(mc(c)*Psi1(a,a1,c))));    
	    }
//----------- OC5
	double OC5, OC5a=0., OC5b=0.; 
	for(  n=0; n<Nn; n++)
	  for( c=0; c<Nc; c++)
	     OC5a +=(mn(n)*mc(c)*Lam(n,c));

	for(  n=0; n<Nn; n++)
	  for( a=0; a<Na; a++)
	        OC5b +=(mn(n)*ma(a)*Lam1(n,a));
	OC5=OC5a+OC5b;
//----------- OC6
	double OC6=0.;
	for(  n=0; n<Nn; n++)
	 for( c=0; c<Nc; c++)
	   for( a=0; a<Na; a++)
	        OC6 +=(mn(n)*mc(c)*ma(a)*Zeta(n,c,a));
//----------- Addition of all sums:
	double OCges=OC1+OC2+OC3+OC4+OC5+OC6;
//----------- Summation of Molalities
	double   OCmol= sum(aM, xcx, Nc)+ sum(aM, xax, Na)+ sum(aM, xnx, Nn);
//----------- Osmotic coefficient (OC) = (1+Oges)/(OCmol) 
	double OC = (1.+OCges) / OCmol;
//-----------Activity of Water	Pitzer-Toughreact Report 2006 equation (A1)
	double Lna =(-18.1/1000.)*OC*OCmol;

	double activityH2O=exp(Lna);
	
//  lnGamma[Ns] = activityH2O/molefractionH2O;
	return Lna-log(x[Ns]);;
}

void TPitzer::getAlp( long int c, long int a, double& alp, double& alp1 )
{
    if( zc(c) || fabs(za(a)) ==1 )
    {    alp=2;
        alp1=12.;
    }
    else
      if( zc(c) && fabs(za(a)) ==2 )
      {   alp=1.4;
          alp1=12.;
      }    
      else
        if( zc(c) && fabs(za(a)) >=2 )
        { alp=2.0;
          alp1=50.;
        } 
        else
          Error( "", "alpha not defined");
}

//----------- ------- F-Factor________________ Pitzer-Toughreact Report 2006 equation (A6)
double TPitzer::F_Factor( double Aphi, double I, double Is )
{
	
  long int c, c1, a, a1;
  double z, z1, Etheta, Ethetap;
//---------- F1
  double F1=-Aphi*( (Is/(1.+1.2*Is)) + 2.*log(1.+1.2*Is)/1.2);
//----------- F2
  double F2=0., Phip;
   for( c=0; c<Nc; c++ )
	 for( c1=c+1; c1<Nc; c1++ )
	 {
		 z=zc(c);
	     z1=zc(c1);
         Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
         Phip = Ethetap;					//Pitzer-Toughreact Report 2006 equation (A16)
         F2 +=(mc(c)*mc(c1)*(Phip));
	 }
//---------- F3
  double F3=0., Phip1;
  for( a=0; a<Na; a++)
	 for( a1=a+1; a1<Na; a1++)
	 {  z=za(a);
        z1=za(a1);
        Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
        Phip1=Ethetap;      				//Pitzer-Toughreact Report 2006 equation (A16)
        F3 +=(ma(a)*ma(a1)*(Phip1));
	 }
//----------- F4
  double F4=0., alp, alp1, h1, h2, g3, g4, B1;
  for( c=0; c<Nc; c++)
	 for( a=0; a<Na; a++)
	 {   
		getAlp(  c,  a, alp, alp1 );  
        h1=alp*Is;
        h2=alp1*Is;
        g3=(-2.*(1.-(1.+h1+(h1*h1)/2. )*exp(-h1)))/(h1*h1);	//Pitzer-Toughreact Report 2006 equation (A13)
        g4=(-2.*(1.-(1.+h2+(h2*h2)/2 )*exp(-h2)))/(h2*h2);	//Pitzer-Toughreact Report 2006 equation (A13)
        B1= (bet1(c,a)*g3)/I+ (bet2(c,a)*g4)/I;			//Pitzer-Toughreact Report 2006 equation (A12)
        F4 = F4+ (mc(c)*ma(a)*B1);       
	 }
// ----------- F-Factor
   return F1+F2+F3+F4;
}   

//--------------------------lnGammaM-------------------------------------
//------------------------CALCULATIONS-----------------------------------
// Act coeff will be calculated for cation M(here cation number 1)
double TPitzer::lnGammaM(  long int M )
{
  double Etheta, Ethetap;
  long int a, n, c1, a1;

//----------- GM1________________ Pitzer-Toughreact Report 2006 equation (A3)
 
 double GM1=(zc(M)*zc(M))*Ffac;
//----------- GM2
 double GM2=0., alp, alp1, h1, h2, g1,g2, B2, C;
 for( a=0; a<Na; a++)
 {
	 getAlp(  M,  a, alp, alp1 );  
     C= Cphi(M,a)/(2.*sqrt(fabs(za(a)*zc(M))));	//Pitzer-Toughreact Report 2006 equation (A7)
     h1=alp*Is;
     h2=alp1*Is;
     g1=(2.*(1.-(1.+h1)*exp(-h1)))/(h1*h1);		  //Pitzer-Toughreact Report 2006 equation (A11)
     g2=(2.*(1.-(1.+h2)*exp(-h2)))/(h2*h1);		  //Pitzer-Toughreact Report 2006 equation (A11)
     B2= bet0(M,a)+(bet1(M,a)*g1)+ (bet2(M,a)*g2); //Pitzer-Toughreact Report 2006 equation (A10)
     GM2=GM2+(ma(a)*(2.*B2+Zfac*(C/(2.*sqrt(fabs(zc(M)*za(a)))))));   
 }
//----------- GM3
  double GM3=0., Phi, z, z1;
  for( c1=0; c1<Nc; c1++)
	 for( a=0; a<Na; a++)
	 {  
	    if( M == c1)
	    {
	      Phi = 0.;
	      Psi(M,c1,a) = 0.;
	    }
	    else
	    {  z=zc(M);
	       z1=zc(c1);
           Ecalc(z,z1,I,Aphi,Etheta,Ethetap);
           Phi=Theta(M,c1)+Etheta;  					//Pitzer-Toughreact Report 2006 equation (A15)
	    }   
        GM3=GM3+mc(c1)*(2.*Phi+ ma(a)*Psi(M,c1,a)  );
	 }
//----------- GM4
    double GM4=0.;
    for( a=0; a<Na; a++)
        for( a1=a+1; a1<Na; a1++)
            GM4=GM4+(ma(a)*ma(a1)*Psi1(a,a1,M));
//----------- GM5
    double GM5a=0.;
    for( c1=0; c1<Nc; c1++) 
     for( a=0; a<Na; a++)
     { C = Cphi(c1,a)/(2.*sqrt(fabs(za(a)*zc(c1))));			//Pitzer-Toughreact Report 2006 equation (A7)
       GM5a = GM5a+(mc(c1)*ma(a)* (C/(2.*sqrt(fabs(zc(M)*za(a)))))  );
     }
    double GM5=zc(M)*GM5a;
//----------- GM6
   double GM6a=0;
   for( n=0; n<Nn; n++)
       GM6a += mn(n)*Lam(n,M);
  
   double GM6 = 2*GM6a;
//----------- GM
   double GM=GM1+GM2+GM3+GM4+GM5+GM6;
   double actcoeffM=exp(GM);
   return GM;
}

//-------------------------lnGammaX--------------------------------------
//------------------------CALCULATIONS-----------------------------------
double TPitzer::lnGammaX(  long int X )
{
  double Etheta, Ethetap;
  long int c, n, c1, a1;

//----------- GX1________________ Pitzer-Toughreact Report 2006 equation (A4)
  double   GX1=(za(X)*za(X))*Ffac;
//----------- GX2
  double GX2=0., C, h1, h2, g1, g2, B2, alp, alp1;
  for( c=0; c<Nc; c++)
  {
 	 getAlp(  c,  X, alp, alp1 );  
     C=Cphi(c,X)/(2.*sqrt(fabs(za(X)*zc(c))));    // Pitzer-Toughreact Report 2006 equation (A7)
     h1=alp*Is;
     h2=alp1*Is;
     g1=(2.*(1.-(1.+h1)*exp(-h1)))/(h1*h1);			 // Pitzer-Toughreact Report 2006 equation (A11)
     g2=(2.*(1.-(1.+h2)*exp(-h2)))/(h2*h2);			 // Pitzer-Toughreact Report 2006 equation (A11)
     B2= bet0(c,X)+bet1(c,X)*g1+ (bet2(c,X)*g2); // Pitzer-Toughreact Report 2006 equation (A10)
     GX2=GX2+(mc(c)*(2.*B2+Zfac*(C/(2.*sqrt(fabs(zc(c)*za(X)))))));
  }           
//----------- GX3
  double  GX3=0., z, z1, Phi1 ;
  for( a1=0; a1<Na; a1++)
	 for( c=0; c<Nc; c++)
	 {  
	    if( X == a1)
	    {
	      Phi1 = 0.;
	      Psi1(X,a1,c) = 0.;
	    }else
	      {	
           z=zc(X);
           z1=zc(a1);
           Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
           Phi1=Theta(X,a1)+Etheta; 			 // Pitzer-Toughreact Report 2006 equation (A15)
	      }    
          GX3=GX3+ma(a1)*(2*Phi1+mc(c)*Psi(X,a1,c));
	 }
//----------- GX4
    double  GX4=0.;
    for( c=0; c<Nc; c++)
	   for( c1=c+1; c1<Nc; c1++)
              GX4=GX4+(mc(c)*mc(c1)*Psi(c,c1,X));
//----------- GX5
     double GX5a=0.;
     for( c=0; c<Nc; c++)
  	   for( a1=0; a1<Na; a1++)
  	   { 	
          C =Cphi(c,a1)/(2*sqrt(fabs(za(a1)*zc(c))));	 // Pitzer-Toughreact Report 2006 equation (A7)
          GX5a =GX5a+(mc(c)*ma(a1)* (C/(2*sqrt(fabs(zc(c)*za(X))))) );
  	   }   
      double GX5=za(X)*GX5a;
//----------- GX6
     double GX6a=0.;
     for( n=0; n<Nn; n++)
         GX6a += mn(n)*Lam(n,X);
     double GX6=2*GX6a;
//----------- GX
    double GX=GX1+GX2+GX3+GX4+GX5+GX6;
    double actcoeffX=exp(GX);
    return GX;
}


//--------------------------lngammaN-------------------------------------
//------------------------CALCULATIONS-----------------------------------
//   Act coeff will be calculated for this neutral species 
double TPitzer::lnGammaN(  long int N )
{
  long int c, a;
//----------- GN1________________ Pitzer-Toughreact Report 2006 equations (A5)
  double GN1=0.;
  for( a=0; a<Na; a++)
      GN1=GN1+(ma(a)*2.*Lam1(N,a));
//----------- GN2
  double GN2=0.;
  for( c=0; c<Nc; c++)
      GN2=GN2+(mc(c)*2.*Lam(N,c));
//----------- GN3
  double GN3=0.;
  for( c=0; c<Nc; c++)
  for( a=0; a<Na; a++)
      GN3=GN3+(mc(c)*ma(a)*Zeta(N,c,a));
//----------- GN
  double GN=GN1+GN2+GN3;
  double actcoeffN=exp(GN);
  
  return GN;
}


//--------------------- End of s_pitzer.cpp ---------------------------
