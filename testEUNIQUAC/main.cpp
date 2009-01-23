//--------------------------------------------------------------------
// $Id: main.cpp 182 2008-05-27 08:06:21Z gems $
//
// Demo test of usage of the TPitzer class for implementing a simple
//
// Copyright (C) 2008 T.Wagner, D.Kulik, F.F.Hingerl
//
// This file is part of GEMIPM2K code for thermodynamic modelling
// by Gibbs energy minimization
//
// This file may be distributed under the licence terms defined
// in GEMIPM2K.QAL
//
// See also http://les.web.psi.ch/Software/GEMS-PSI
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <iostream>
using namespace std;


#include "s_fgl.h"

int main( int argc, char* argv[] )
 {
   // Default arguments
    double phVOL[1] = {1.};
    long int NComp  = 6; //!!!!
    double alnGam[6]; // Na+, Cl-, SO42-, OH-, H+, H2O@
                       // 0    1    2	   3   4    5
    double aM[6] = { 2., 1., 0.5, 1e-7, 1e-7, 55.5083};
	double aWx[6] = { 0.033894, 0.016947, 0.00847, 2e-9, 2e-9, 0.940686} ;
    double aZ[6] = { 1., -1., -2., -1., 1., 0};

    long int MaxOrd = 2; // number of columns of aIPx
    long int NPcoef = 2; // r(i) and q(i)

    long int NPar = 21;  // u0 and ut number of rows
    long int aIPx[21*2] = {
    		0L, 0L,
    		0L, 1L,
    		0L, 2L,
    		0L, 3L,
    		0L, 4L,
    		0L, 5L,
    		1L, 1L,
    		1L, 2L,
    		1L, 3L,
    		1L, 4L,
    		1L, 5L,
    		2L, 2L,
    		2L, 3L,
    		2L, 4L,
    		2L, 5L,
    		3L, 3L,
    		3L, 4L,
    		3L, 5L,
    		4L, 4L,
    		4L, 5L,
    		5L, 5L

    };

    double aIPc[21*2] = {
    		0., 0.,
			1443.23, 15.635,
			845.135, 11.681,
    		1398.14, 20.278,
			1e10, 0.,
			733.286, 0.4872,
    		2214.81, 14.436,
			2036.06, 12.407,
			1895.52, 13.628,
            1e10, 0.,
			1523.39, 14.631,
			1265.83, 8.3194,
			1225.67, 8.5902,
			1e10, 0.,
            752.879, 9.4905,
			1562.88, 5.6169,
			1e10, 0.,
			600.495, 8.5455,
			0.,	0.,
			100000., 0.,
			0.,	0.
    };

    long int NP_DC = 2; // number of parameter per component
    double aDCc[6*2] = { 	// volume and surface area parameters
    		1.4034, 1.199,
    		10.386, 10.197,
    		12.794, 12.444,
			9.3973, 9.8171,
			0.13779, 1e-015,
			0.92, 1.4
    };

    double RhoW = 1., EpsW = 78.38, IS = 10. ;
    double T_k = 298.15;
    double P_bar = 1.;

    TEUNIQUAC *aEU = new TEUNIQUAC( NComp, NPar, NPcoef, MaxOrd,
	         NP_DC, 'Q',aIPx, aIPc, aDCc,
	         aWx, alnGam, phVOL, aM, aZ);

cout << "TEUNIQUAC class instance initialized" << endl;

    aEU->PTparam(T_k, P_bar, RhoW, EpsW);
//cout << "Pitzer parameters corrected to T,P" << endl;


	double Y=0;
	for (Y=0; Y<6; Y += 0.1)	{
		aM[0]=1+Y; aM[1]=Y;
		double sum = aM[0]+aM[1]+aM[2]+aM[3]+aM[4]+aM[5];
		aWx[0] = aM[0]/sum;
		aWx[1] = aM[1]/sum;
		aWx[2] = aM[2]/sum;
		aWx[3] = aM[3]/sum;
		aWx[4] = aM[4]/sum;
		aWx[5] = aM[5]/sum;
		aEU->MixMod();     //Euniquac_calc_Gamma();
		cout << "Euniquac activity coefficients calculated" << " m(NaCl) = " << Y << endl;
		aEU->Euniquac_test_out( "Euniquac_test.out" );
	}
		cout << "See file 'Euniquac_test.out' for results..." << endl;

	delete aEU;
    return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for gemnode - the TNode class usage example
