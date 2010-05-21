//-------------------------------------------------------------------
// $Id: num_methods.h 705 2006-04-28 19:39:01Z gems $
//
// C/C++ Numerical Methods used in GEMS-PSI and GEMIPM2K
// (c) 2006,2008 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of the GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _num_methods_h_
#define _num_methods_h_

#include <math.h>


// Lagrangian interpolation functions
double LagranInterp(float *y, float *x, double *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(float *y, float *x, float *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(double *y, double *x, double *d, double yoi,
		double xoi, long int M, long int N, long int pp );


// generic functions for calculating partial derivatives
double quot( double u, double v, double du, double dv );

double quot( double u, double v, double du, double dv, double d2u, double d2v );


double prod2( double u, double v, double du, double dv );

double prod2 ( double u, double v, double du, double dv, double d2u, double d2v );


double prod3 ( double u, double v, double w, double du, double dv, double dw );

double prod3 ( double u, double v, double w, double du, double dv, double dw,
		double d2u, double d2v, double d2w );


double GoldenSection( double param[3], double funct[2], double (f_proc)(double val));


#endif   // _num_methods_h_

//-----------------------End of num_methods.h--------------------------

