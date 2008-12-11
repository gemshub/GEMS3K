//--------------------------------------------------------------------
// $Id: main.cpp 182 2008-05-27 08:06:21Z gems $
//
// gemnode
// Demo test of usage of the TPitzer class for implementing a simple
//
// Copyright (C) 2008 S.Dmytrieva, D.Kulik
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
    long int NComp  = 3; //!!!!
    double aWx[3] = { 0.,0.,0.8222 } ;
    double alnGam[6]; // Na+, Cl-, H2O@
                       // 0    1    2
    double aM[6] = { 6., 6., 0.,};
    double aZ[6] = { 1., -1., 0.,};

    long int MaxOrd = 4;
    long int NPcoef = 5; // or 8

    long int NPar = 3; //!!!!
    long int aIPx[3*4] = {
    		0, 1, -1, -10,
    		0, 1, -1, -11,
    		0, 1, -1, -20,
    };

    double aIPc[3*5] = {
    		0.0765, 0., 0., 0., 0.,
    		0.2664, 0., 0., 0., 0.,
    		0.00127, 0., 0., 0., 0.,
    };

    long int NP_DC = 0;
    double *aDCc = NULL;

    double RhoW = 1., EpsW = 78.38, IS = 10. ;
    double T_k = 298.15;
    double P_bar = 1.;
    char letter;

    TPitzer *aPT = new TPitzer( NComp, NPar, NPcoef, MaxOrd,
	         NP_DC/**/, T_k, P_bar, 'Z',
	         aIPx, aIPc, aDCc/**/,
	         aWx /**/, alnGam, phVOL, aM, aZ,
	         RhoW, EpsW );
cout << "TPitzer class instance initialized: enter any character to proceed" << endl;
cin >> letter;
    aPT->Pitzer_calc_Gamma();
cout << "Pitzer activity coefficients calculated: enter any character to finish" << endl;
cin >> letter;
    aPT->Pitzer_test_out( "Pitzer_test.out" );
cout << "See file 'Pitzer_test.out' for results. Bye!" << endl;
    delete aPT;

    return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for gemnode - the TNode class usage example
