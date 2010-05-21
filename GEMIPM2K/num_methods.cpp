//-------------------------------------------------------------------
// $Id: num_methods.cpp 705 2006-04-28 19:39:01Z gems $
//
// C/C++ Numerical Methods (Lagrange interpolation)
// used in GEMS-PSI and GEMIPM2K
//
// (c) 2006,20078 S.Dmytriyeva, D.Kulik
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
double LagranInterp(double *y, double *x, double *d, double yoi,
                    double xoi, long int M, long int N, long int pp )
{
    double s,z,s1[21];
    long int ppy, ppx, py, px, i, j, k, jx, jy, jy1;

    py = N-1;
    px = M-1;

   if (yoi < y[0] || yoi > y[py] )
     Error( "LagranInterp",
       "E34RErun: yoi < y[0] or yoi > y[py] ( row argument outside the range )");
   if(xoi < x[0] || xoi > x[px] )
   Error( "LagranInterp",
    "E34RErun: xoi < x[0] or xoi > x[px] ( column argument outside the range )");

   if( N==1 && M==1 ) // zero dimension interpolation
      return d[0];

// find the point in the row
   ppy = min( N-1, pp );
   for(jy=0;jy<N;jy++)
     if ( yoi >= y[jy] && yoi <= y[jy+1])
        break;
   if( jy >= N-ppy)
       jy=N-ppy-1;
   jy1=jy;

// find the point in the column
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
  double *dd, *yy, *xx, res;

  dd = new double[N*M];
  for(int ii=0; ii<N*M; ii++ )
    dd[ii] = (double)d[ii];
  yy = new double[N];
  for(int ii=0; ii<N; ii++ )
    yy[ii] = (double)y[ii];
  xx = new double[M];
  for(int ii=0; ii<M; ii++ )
    xx[ii] = (double)x[ii];

  res = LagranInterp( yy, xx, dd, (double)yoi, (double)xoi, M, N, pp );

 delete[] dd;
 delete[] yy;
 delete[] xx;
 return res;
}

double LagranInterp(float *y, float *x, double *d, float yoi,
                    float xoi, int M, int N, int pp )
{
  double *yy, *xx, res;

  yy = new double[N];
  for(int ii=0; ii<N; ii++ )
      yy[ii] = (double)y[ii];
  xx = new double[M];
  for(int ii=0; ii<M; ii++ )
      xx[ii] = (double)x[ii];

  res = LagranInterp( yy, xx, d, (double)yoi, (double)xoi, M, N, pp );

 delete[] yy;
 delete[] xx;
 return res;
}


// 1st partial derivative of quotient of two functions
double quot( double u, double v, double du, double dv )
{
	double derivative;
	derivative = ( du*v - u*dv ) / pow (v,2.);

	return derivative;
}


// 2nd partial derivative of quotient of two functions
double quot( double u, double v, double du, double dv, double d2u, double d2v )
{
	double derivative;
	derivative = (d2u*v + du*dv)/pow(v,2.) - (du*v)*(2.*dv)/pow(v,3.)
				- (du*dv + u*d2v)/pow(v,2.) + (u*dv)*(2.*dv)/pow(v,3.);

	return derivative;
}


// 1st partial derivative of product of two functions
double prod2( double u, double v, double du, double dv )
{
	double derivative;
	derivative = ( du*v + u*dv );

	return derivative;
}


// 2nd partial derivative of product of two functions
double prod2( double u, double v, double du, double dv, double d2u, double d2v )
{
	double derivative;
	derivative = ( d2u*v + 2*du*dv + u*d2v );

	return derivative;
}


// 1st partial derivative of product of three functions
double prod3( double u, double v, double w, double du, double dv, double dw )
{
	double derivative;
	derivative = ( du*v*w + u*dv*w + u*v*dw );

	return derivative;
}


// 2nd partial derivative of product of three functions
double prod3( double u, double v, double w, double du, double dv, double dw,
		double d2u, double d2v, double d2w )
{
	double derivative;
	derivative = ( d2u*v*w + du*dv*w + du*v*dw ) + ( du*dv*w + u*d2v*w + u*dv*dw )
				+ ( du*v*w + u*dv*w + u*v*d2w );

	return derivative;
}


// Method of Gold Section
//   param[3],  parameter x    - start, end, tolerance
//   funct[2],  function  f(x) - value, tolerance
//
double GoldenSection( double param[3], double funct[2], double (f_proc)(double val))
{
    double Fa, Fb, Fx1, Fx2, a, b, x1, x2;

    x1 = param[0];
    x2 = param[1];
    a = min( x1, x2 );
    b = max( x1, x2 );
    if( (b-a) < param[2]) goto DONE;
    x1 = a + .382*(b-a);
    x2 = a + .618*(b-a);
    Fa = f_proc( a);
    if( fabs(Fa) < funct[1] )
    {
        b = a;
        goto DONE;
    }
    Fb = f_proc( b);
    if(  fabs(Fb) < funct[1] )
    {
        a = b;
        goto DONE;
    }
    if( (Fa*Fb) > 0)
        Error( "GoldenSection",
  "W01PEexec: No result in specified interval! Change interval!");
    Fx1 = f_proc( x1);
    Fx2 = f_proc( x2);
    do
    {
        if( fabs( Fx1 ) < funct[1] )
        {
            a = b = x1;
            goto DONE;
        }
        if( fabs( Fx2 ) < funct[1] )
        {
            a = b = x2;
            goto DONE;
        }
        if( fabs( Fx1) > fabs( Fx2) )
        {
            a = x1;
            if( (b-a) < param[2])
                goto DONE;
            x1 = x2;
            Fx1 = Fx2;
            x2 = a + .618*(b-a);
            Fx2 = f_proc( x2);
        }
        else
        {
            b = x2;
            if( (b-a) < param[2])
                goto DONE;
            x2 = x1;
            Fx2 = Fx1;
            x1 = a + .382*(b-a);
            Fx1 = f_proc( x1);
        }
    }
    while( (b-a) > param[2] );
DONE:
    x1 = (a+b)/2;
    return x1;
}





//-----------------------End of num_methods.cpp--------------------------


