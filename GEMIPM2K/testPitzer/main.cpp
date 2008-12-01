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


#include "s_pitzer.h"

int main( int argc, char* argv[] )
 {
   // Default arguments
    long int NComp  = 6; //!!!!
    double aWx[6] = { 0.,0.,0.,0.,0.,0.9 } ;
    double alnGam[6]; // Cl-, Na+, NaCl@, H+, OH-, H2O@ 
                       // 0    1    2      3   4    5
    double aM[6] = { 5., 5., 0., 1e-7, 1e-7, 0., };
    double aZ[6] = { -1., 1., 0., 1., -1., 0.};

    long int MaxOrd = 4; 
    long int NPcoef = 5; // or 8

    long int NPar = 13; //!!!!
    long int aIPx[13*4] = { 
    		1, 0, -1, -10,
    		3, 0, -1, -10,
    		1, 4, -1, -10,
    		1, 0, -1, -11,
    		3, 0, -1, -11,
    		1, 4, -1, -11,
    		1, 0, -1, -20,
    		3, 0, -1, -20,
    		1, 4, -1, -20,
    		
    		1, 3, -1, -40,
    		1, 3, 0, -50,
    	    0, 4, -1, -41,
    		0, 4, 1, -51 };
    
    double aIPc[13*5] = { 
    		0.0765, 0., 0., 0., 0.,
    		0.1775, 0., 0., 0., 0.,
    		0.0864, 0., 0., 0., 0.,
    		0.2664, 0., 0., 0., 0.,
    		0.2945, 0., 0., 0., 0.,
    		0.253, 0., 0., 0., 0.,
    		0.00127, 0., 0., 0., 0.,
    		0.00080, 0., 0., 0., 0.,
    		0.0044, 0., 0., 0., 0.,
    		0.036, 0., 0., 0., 0.,
    		-0.004, 0., 0., 0., 0.,
    		-0.050, 0., 0., 0., 0.,
    		-0.006, 0., 0., 0., 0.};
    
    long int NP_DC = 0;
    double *aDCc = NULL;
    
    double RhoW = 1., EpsW = 2., IS = 1. ;
    double T_k = 298.15;
    double P_bar = 1.;

    TPitzer *aPT = new TPitzer( NComp, NPar, NPcoef, MaxOrd,
	         NP_DC/**/, T_k, P_bar, 'Z',
	         aIPx, aIPc, aDCc/**/,
	         aWx /**/, alnGam, aM, aZ, 
	         RhoW, EpsW, IS );

    aPT->Pitzer_calc_Gamma();
    aPT->Pitzer_test_out( "Pitzer_test.out" );

    delete aPT;  
    
    return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for gemnode - the TNode class usage example
