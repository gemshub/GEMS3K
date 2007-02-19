//-------------------------------------------------------------------
// $Id: num_methods.h 705 2006-04-28 19:39:01Z gems $
//
// C/C++ Numerical Methods (Linear Algebra) used in GEMS-PSI and GEMIPM2K
// (c) 2006-2007 S.Dmytriyeva, D.Kulik
//
// Uses: JAMA/C++ Linear Algebra Package based on the Template
// Numerical Toolkit (TNT) - an interface for scientific computing in C++,
// (c) Roldan Pozo, NIST (USA), http://math.nist.gov/tnt/download.html
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

double LagranInterp(float *y, float *x, double *d, float yoi,
                    float xoi, int M, int N, int pp );
double LagranInterp(float *y, float *x, float *d, float yoi,
                    float xoi, int M, int N, int pp );

#endif   // _num_methods_h_

//-----------------------End of num_methods.h--------------------------

