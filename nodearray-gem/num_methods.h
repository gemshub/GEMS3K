//-------------------------------------------------------------------
// $Id:  $
//
// C/C++ Some functions from NUMERICAL RECIPES IN C
//
// Copyright (C) 2006 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _num_methods_h_
#define _num_methods_h_

#include <math.h>

double enorm( int n, double *x );
int CholeskyDecomposition( int N, double* R, double* X, double* R1  );
int LUDecomposition( int N, double* A, double* X  );
//void InverseofMatrix( int N, double* A, double* Y );
//double DeterminantofMatrix( int N, double* A );

// Random numbers ==========================================================
double randuni(double& x); // uniform point
double randnorm(double& x); // normal point
// Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham
// shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0
float ran2(long& idum);
// According to Knuth, any large MBIG, and any smaller (but still large) MSEED
// can be substituted for the above values.
// Returns a uniform random deviate between 0.0 and 1.0.
float ran3(long& idum);

#endif   // _num_methods_h_

//-----------------------End of num_methods.h--------------------------

