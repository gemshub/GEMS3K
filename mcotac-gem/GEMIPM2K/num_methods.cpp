//-------------------------------------------------------------------
// $Id: num_methods.cpp 705 2006-04-28 19:39:01Z gems $
//
// C/C++ Numerical Methods (Lagrange interpolation)
// used in GEMS-PSI and GEMIPM2K
//
// (c) 2006-2007 S.Dmytriyeva, D.Kulik
//
// Uses: JAMA/C++ Linear Algebra Package based on the Template
// Numerical Toolkit (TNT) - an interface for scientific computing in C++,
// (c) Roldan Pozo, NIST (USA), http://math.nist.gov/tnt/download.html
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and the GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#include "v_user.h"
#include "num_methods.h"
#include "verror.h"

/*-----------------------------------------------------------------*/
// Interpolation over tabulated values (2D array d) using Lagrange method
//  y[N] - discrete values of argument over rows (ascending order)
//  x[M] - discrete values of arguments over columns (ascending order)
//  d[N][M] - discrete values of a function of x and y arguments
//  yoi - row (y) argument of interest ( y[0] <= yi <= y[N-1] )
//  xoi - column (x) argument of interest ( x[0] <= xi <= x[M-1] )
//  M - number of rows in y array
//  N - number of columns in y array
//  pp  -level of interpolation ( default 1)
//  Function returns an interpolated value of d(yoi,xoi) or error if
//  yoi or xoi are out of range
//
double LagranInterp(float *y, float *x, double *d, float yoi,
                    float xoi, int M, int N, int pp )
{
    double s,z,s1[21];
    int ppy, ppx, py, px, i, j, k, jx, jy, jy1;

    py = N-1;
    px = M-1;

   if (yoi < y[0] || yoi > y[py] )
     Error( "LagranInterp",
       "E34RErun: yoi < y[0] or yoi > y[py] ( row argument outside the range )");
   if(xoi < x[0] || xoi > x[px] )
   Error( "LagranInterp",
    "E34RErun: xoi < x[0] or xoi > x[px] ( column argument outside the range )");

   if( N==1 && M==1 ) // one dimension interpolation
      return d[0];

// find point in row
   ppy = min( N-1, pp );
   for(jy=0;jy<N;jy++)
     if ( yoi >= y[jy] && yoi <= y[jy+1])
        break;
   if( jy >= N-ppy)
       jy=N-ppy-1;
   jy1=jy;

// find point in column
   ppx = min( M-1, pp );
   for(jx=0;jx<M;jx++)
        if(xoi >= x[jx] && xoi <= x[jx+1])
            break;
   if(jx >= M-ppx)
     jx = M-ppx-1;

   for(j=0;j <= ppy; j++)
    {
        s=0.;
        for(i=0;i<=ppx;i++)
        {
            z=1; //z1=1;
            for(k=0;k<=ppx;k++)
                if(k!=i)
                    z*=(xoi-x[k+jx])/(x[i+jx]-x[k+jx]);
            s+=d[i+jx+(jy)*M]*z;
        }

       s1[j]=s;
       jy++;
    }
    s=0.;
    for(i=0;i<=ppy;i++)
    {
        z=1;
        for(k=0;k<=ppy;k++)
            if(k!=i)
                z*=(yoi-y[k+jy1])/(y[i+jy1]-y[k+jy1]);
        s+=s1[i]*z;
    }
    return(s);
}

double LagranInterp(float *y, float *x, float *d, float yoi,
                    float xoi, int M, int N, int pp )
{
  double *dd, res;

  dd = new double[N*M];
  for(int ii=0; ii<N*M; ii++ )
   dd[ii] = (double)d[ii];

  res = LagranInterp( y, x, dd, yoi, xoi, M, N, pp );

 delete[] dd;
 return res;
}

//-----------------------End of num_methods.cpp--------------------------


