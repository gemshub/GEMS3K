//-------------------------------------------------------------------
// $Id:  $
//
// C/C++ Some functions from TNT-JAMA
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

#include "v_user.h"
#include "num_methods.h"
#include "verror.h"

/*-----------------------------------------------------------------*/
// Interpolation over tabulated values (array y) using the Lagrange method
//  y[N] - discrete values of argument over rows (ascending order)
//  x[M] - discrete values of arguments over columns (ascending order)
//  d[N][M] - discrete values of a function of x and y arguments
//  xoi - column (x) argument of interest ( x[0] <= xi <= x[M-1] )
//  yoi - row (y) argument of interest  ( y[0] <= yi <= y[N-1] )
//  N - number of rows in y array;
//  M - number of columns in y array.
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

// the squares of the following constants shall not under/overflow:
// these values seem good for an x86:
#define LM_SQRT_DWARF 1.e-160
#define LM_SQRT_GIANT 1.e150
// the following values should work on any machine:
// #define LM_SQRT_DWARF 3.834e-20
// #define LM_SQRT_GIANT 1.304e19
#define SQR(x)   (x)*(x)

double enorm( int n, double *x )
{
/*     given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     the euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. the sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. non-destructive underflows are permitted. underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     the definitions of small, intermediate and large components
 *     depend on two constants, LM_SQRT_DWARF and LM_SQRT_GIANT. the main
 *     restrictions on these constants are that LM_SQRT_DWARF**2 not
 *     underflow and LM_SQRT_GIANT**2 not overflow.
 *
 *     parameters
 *
 *	n is a positive integer input variable.
 *
 *	x is an input array of length n.
 */
    int i;
    double agiant, s1, s2, s3, xabs, x1max, x3max, temp;

    if( n <= 0 )
      return 0.0;
    s1 = 0;
    s2 = 0;
    s3 = 0;
    x1max = 0;
    x3max = 0;
    agiant = LM_SQRT_GIANT/( (double) n);

    for ( i=0; i<n; i++ )
    {
        xabs = fabs(x[i]);
        if ( xabs > LM_SQRT_DWARF && xabs < agiant )
        {
// **  sum for intermediate components.
            s2 += xabs*xabs;
            continue;
        }

        if ( xabs >  LM_SQRT_DWARF )
        {
// **  sum for large components.
            if (xabs > x1max)
            {
                temp = x1max/xabs;
                s1 = 1 + s1*SQR(temp);
                x1max = xabs;
            }
            else
            {
                temp = xabs/x1max;
                s1 += SQR(temp);
            }
            continue;
        }
// **  sum for small components.
        if (xabs > x3max)
        {
            temp = x3max/xabs;
            s3 = 1 + s3*SQR(temp);
            x3max = xabs;
        }
        else
        {
            if (xabs != 0.)
            {
                temp = xabs/x3max;
                s3 += SQR(temp);
            }
        }
    }

// *** calculation of norm.

    if (s1 != 0)
        return x1max*sqrt(s1 + (s2/x1max)/x1max);
    if (s2 != 0)
    {
        if (s2 >= x3max)
            return sqrt( s2*(1+(x3max/s2)*(x3max*s3)) );
        else
            return sqrt( x3max*((s2/x3max)+(x3max*s3)) );
    }

    return x3max*sqrt(s3);
}


// Random numbers ==========================================================

// uniform point
double randuni(double& x)
{ double m35=34359738368., m36=68719476736., m37=137438953472.;
  float a=0.,b=1.;
  if( x < 0 ) // Initialize. process
  {
    int j;
    double R;
    j=rand();
    R = ceil(24359738368.*j/RAND_MAX + 10000000000.);
    if( !fmod(R,2) )
      R=R+1.;
    x = R;
  }
  x=x*5.;
  if(x>=m37) x=x-m37;
  if(x>=m36) x=x-m36;
  if(x>=m35) x=x-m35;
 return(x/m35*(b-a)+a);
}

// normal point
double randnorm(double& x)
{ double R1=0.;
  int j;
  for(j=0;j<101;j++)
    R1+=randuni(x);
      R1=(R1-101./2.)/pow(101./12.,0.5);
           R1=1./6.*(R1-(-3.0));
           if(R1<0.) R1=0.;
           if(R1>1.) R1=1.;
           return(R1);
/*return(1./9.*(R1-(-4.5)));*/
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS  1.2e-7
#define RNMX  (1.0-EPS)

// Long period (> 2 � 1018) random number generator of L�Ecuyer with Bays-Durham shuffle
// and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
// the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
// idum between successive deviates in a sequence. RNMX should approximate the largest floating
// value that is less than 1.
float ran2(long& idum)
{
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[NTAB];
   float temp;
   if (idum <= 0)
   { // Initialize.
      if ( -idum < 1)
         idum=1; // Be sure to prevent idum = 0.
      else
         idum = -idum;
      idum2= idum;
     for (j=NTAB+7;j>=0;j--) // Load the shuffle table (after 8 warm-ups).
     {
       k = idum/IQ1;
       idum = IA1*(idum-k*IQ1)-k*IR1;
       if (idum < 0)
         idum += IM1;
       if (j < NTAB)
         iv[j] = idum;
     }
   iy=iv[0];
  }
   k = idum/IQ1;       //    Start here when not initializing.
   idum = IA1 * (idum-k*IQ1) - k*IR1;  //Compute idum=(IA1*idum) % IM1 without
                                       // overflows by Schrage�s  method.
   if ( idum < 0 )
      idum += IM1;
   k = idum2/IQ2;
   idum2 = IA2*(idum2-k*IQ2)-k*IR2;  //Compute idum2=(IA2*idum) % IM2 likewise.
   if (idum2 < 0)
      idum2 += IM2;
   j = iy/NDIV;                      // Will be in the range 0..NTAB-1.
   iy = iv[j]-idum2;                 // Here idum is shuffled, idum and idum2 are
   iv[j] = idum;                     // combined to generate output.
   if (iy < 1)
      iy += IMM1;
   if ( ( temp=AM*iy) > RNMX)
       return (float)RNMX;                //  Because users don�t expect endpoint values.
   return temp;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// According to Knuth, any large MBIG, and any smaller (but still large) MSEED
// can be substituted for the above values.
// Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to
// initialize or reinitialize the sequence.
float ran3(long& idum)
{
  static int inext,inextp;
  static long ma[56];     // The value 56 (range ma[1..55]) is special and
  static int iff=0;       // should not be modified; see  Knuth.
  long mj,mk;
  int i,ii,k;
  if ( idum < 0 || iff == 0) // Initialization.
  {
     iff=1;
     mj = labs( MSEED - labs(idum));  // Initialize ma[55] using the seed idum
     mj %= MBIG;                      // and the large number MSEED.
     ma[55] = mj;
     mk = 1;
     for (i=1;i<=54;i++)   // Now initialize the rest of the table,
     {  ii=(21*i) % 55;    // in a slightly random order,
        ma[ii]=mk;         // with numbers that are not especially random.
        mk=mj-mk;
        if (mk < MZ)
          mk += MBIG;
        mj=ma[ii];
      }
      for (k=1;k<=4;k++)    // We randomize them by �warming upthe generator.�
        for(i=1;i<=55;i++)
        {
         ma[i] -= ma[1+(i+30) % 55];
         if (ma[i] < MZ)
            ma[i] += MBIG;
         }
       inext=0;     // Prepare indices for our first generated number.
       inextp=31;  //  The constant 31 is special; see Knuth.
       idum=1;
   }
   //  Here is where we start, except on initialization.
    if (++inext == 56)
       inext=1;           // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56)
       inextp=1;
    mj = ma[inext]-ma[inextp];  // Generate a new random number subtractively.
    if (mj < MZ)
         mj += MBIG;    //  Be sure that it is in range.
    ma[inext]=mj;       // Store it,
    return mj*FAC;      // and output the derived uniform deviate.
}

//-----------------------End of num_methods.cpp--------------------------

