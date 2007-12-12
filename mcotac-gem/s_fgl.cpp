//-------------------------------------------------------------------
// $Id: s_fgl.cpp 895 2007-03-20 15:26:05Z wagner $
//
// Copyright (C) 2004-2007  S.Churakov, Th.Wagner, D.Kulik
//
// Implementation of TFGLcalc class (see s_fgl.h)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// and part of the GEMIPM2K standalone code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------

#include <math.h>
#include <stdio.h>

#include "s_fgl.h"
#ifndef IPMGEMPLUGIN
  #include "m_const.h"
#endif

//--------------------------------------------------------------------//

double TCGFcalc::DIntegral(double T, double ro,unsigned IType)
{
  static double TOld,roOld;
  static double a,b,c,d,e;
   static double data[][6]=
      {{-0.257431, 0.439229,  0.414783,  -0.457019, -0.145520,  0.299666},
      {-0.396724, 0.690721,  0.628935,  -0.652622, -0.201462, -0.23163 },
      {-0.488498, 0.863195,  0.761344,  -0.750086, -0.218562, -0.538463},
      {-0.556600, 0.995172,  0.852903,  -0.804710, -0.214736, -0.761700},
      {-0.611295, 1.103390,  0.921359,  -0.838804, -0.197999, -0.940714},
      {-0.657866, 1.196189,  0.975721,  -0.862346, -0.172526, -1.091678},
      {-0.698790, 1.278054,  1.020604,  -0.880027, -0.140749, -1.222733},
      {-0.735855, 1.351533,  1.058986,  -0.894024, -0.104174, -1.338626},
      {-0.769504, 1.418223,  1.092052,  -0.905347, -0.063730, -1.442391},
      {-0.800934, 1.479538,  1.121453,  -0.914864, -0.020150, -1.536070},
      {-0.829779, 1.535822,  1.147161,  -0.922381, 0.026157 , -1.621183},
      {-0.856655, 1.587957,  1.169885,  -0.928269, 0.074849 , -1.698853},
      {-0.881757, 1.636402,  1.190082,  -0.932668, 0.125590 , -1.769898},
      {-0.904998, 1.681421,  1.207610,  -0.935419, 0.178283 , -1.835070},
      {-0.926828, 1.723393,  1.223088,  -0.936667, 0.232649 , -1.894899},
      {-0.946773, 1.762571,  1.236007,  -0.936403, 0.288687 , -1.949858},
      {-0.965248, 1.799170,  1.246887,  -0.934650, 0.346207 , -2.000344}};
  //static double dt12[]=
//      {-2.139734,1.971553, 0.945513, -1.901492,-0.588630,-5.390941};
//      {-0.637684, 0.708107,  0.222086,  -0.481116, -0.332141, -3.492213};

   unsigned n;
   double *dtmp,rez;

  if ( (T!=TOld) || (ro!=roOld) )
  {
    TOld=T;
    roOld=ro;
    e=log(T);
    b=ro*ro;
    d=ro;
    c=ro*e;
    a=b*e;
  }

  // special case
  /*
  if ( IType==12 )
  {
    rez=(dt12[0]*T + dt12[1])*b +
    (dt12[2]*T + dt12[3])*ro + dt12[4]*T + dt12[5];
    return exp(rez);
  }
    */

  n=IType-4;
  dtmp=data[n];
  rez=dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
  return exp(rez);
}


double TCGFcalc::LIntegral(double T, double ro,unsigned IType)
{
  static double TOld,roOld;
  static double a,b,c,d,e;
  static double data[][6]=
  {{ -1.010391, 1.628552,  2.077476,  -2.30162 , -0.689931, -2.688117},
   { -1.228611, 2.060090,  2.463396,  -2.453303, -0.573894, -3.350638},
   { -1.354004, 2.402034,  2.718124,  -2.462814, -0.412252, -4.018632}};

   double *dtmp,rez;

  if ( (T!=TOld) || (ro!=roOld) )
  {
    TOld=T;
    roOld=ro;
    a=ro*ro*log(T);
    b=ro*ro;
    c=ro*log(T);
    d=ro;
    e=log(T);
  }

  switch ( IType )
  {
    case 662:
         dtmp=data[0];
         break;
    case 1262:
         dtmp=data[1];
         break;
    case 12122:
         dtmp=data[2];
         break;
    default:
         return 0;
  }

  rez=dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
  return -exp(rez);

}

double TCGFcalc::KIntegral(double T, double ro,long unsigned IType)
{
  static double TOld,roOld;
  static double a,b,c,d,e;
  static double data[][6]=
  {{ -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720},
   { -1.309550, 2.249120,  2.135877,  -2.278530, -0.773166, -3.704690},
   { -1.490116, 2.619997,  2.404319,  -2.420706, -0.829466, -3.930928},
   { -1.616385, 2.881007,  2.577600,  -2.484990, -0.828596, -4.175589},
   { -1.940503, 3.552034,  2.940925,  -2.593808, -0.724353, -4.899975}};

   double *dtmp,rez;

  if ( (T!=TOld) || (ro!=roOld) )
  {
    TOld=T;
    roOld=ro;
    a=ro*ro*log(T);
    b=ro*ro;
    c=ro*log(T);
    d=ro;
    e=log(T);
  }

  switch ( IType )
  {
    case 222333:
         dtmp=data[0];
         break;
    case 233344:
         dtmp=data[1];
         break;
    case 334445:
         dtmp=data[2];
   rez=dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
         return -exp(rez);
    case 444555:
         dtmp=data[3];
         break;
    case 666777:
         dtmp=data[4];
         break;
    default:
         return 0;

  }

   rez=dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
   return exp(rez);
}


double TCGFcalc::K23_13(double T, double ro)
{
  static double TOld,roOld,KOLD;
  static double a,b,c,d,e;
  static double dtmp[]=
  { -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720};

  if ( (T!=TOld) || (ro!=roOld) )
  {
    TOld=T;
    roOld=ro;
    a=ro*ro*log(T);
    b=ro*ro;
    c=ro*log(T);
    d=ro;
    e=log(T);
  }
  else return KOLD;

   KOLD=dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
   KOLD=exp(KOLD/3.);
   return KOLD;

  }

/////////////////////////////////////////////////////////////////////////////
// Implementation of TCGFcalc class

int TCGFcalc::CGFugacityPT( float *EoSparam, float *EoSparPT, double &Fugacity,
        double &Volume, double &DeltaH, double &DeltaS, double P, double T )
{
      int iRet=0; double ro;
      double X[1]={1.};
      double FugPure[1];

		// modification to simplify CG database structure, TW 20/03/2007
        EoSparPT[0] = EoSparam[0]+EoSparam[4]*(float)exp(T*EoSparam[5]);
        EoSparPT[1] = EoSparam[1]+EoSparam[6]*(float)exp(T*EoSparam[7]);
        EoSparPT[2] = EoSparam[2]+EoSparam[8]/((float)T+EoSparam[9]);
        EoSparPT[3] = EoSparam[3]+EoSparam[10]/((float)T+EoSparam[11]);

      /*switch (int(EoSparam[4]))
      {
       case 0:
        EoSparPT[0]=EoSparam[0];
        EoSparPT[1]=EoSparam[1];
        EoSparPT[2]=EoSparam[2];
        EoSparPT[3]=EoSparam[3];
       break;
       case 1:  // H2O type
        EoSparPT[0]=EoSparam[0]+EoSparam[5]/((float)T+EoSparam[6]);
        EoSparPT[1]=EoSparam[1]+EoSparam[7]/((float)T+EoSparam[8]);
        EoSparPT[2]=EoSparam[2]+EoSparam[9]/((float)T+EoSparam[10]);
        EoSparPT[3]=EoSparam[3]+EoSparam[11]/((float)T+EoSparam[12]);
        break;
       case 2:  // CO2 type
        EoSparPT[0]=EoSparam[0]+EoSparam[5]*(float)exp(T*EoSparam[6]);
        EoSparPT[1]=EoSparam[1]+EoSparam[7]*(float)exp(T*EoSparam[8]);
        EoSparPT[2]=EoSparam[2]+EoSparam[9]*(float)exp(T*EoSparam[10]);
        EoSparPT[3]=EoSparam[3]+EoSparam[11]*(float)exp(T*EoSparam[12]);
        break;
        default:

        return 1;// Error: Wrong type of equation
      };*/


 // returns density!
      ro = CGActivCoefPT( X, EoSparPT, FugPure, 1, P, T );
      if( ro < 0.  )
      {
          return -1;
      };
      Fugacity= FugPure[0];
      ro = DENSITY( X, EoSparPT, 1, P, T );
      if( ro < 0 )
      {  // error - density could not be calculated
         iRet = -2; ro = 1.0;
      }
      Volume=0.1/ro;  // in J/bar
//
      DeltaH=0.;  // To be completed
      DeltaS=0.;  // To be completed
//
      return iRet;
  }


   double TCGFcalc::CGActivCoefPT(double *X,float *param, double *act, unsigned NN,
        double Pbar, double T )
   {
      //double act[MAXPARAM];
      //unsigned ncmp;
       double *xtmp,*Fx;
       double P=Pbar/10.;
      //ncmp=unsigned((nn-2)/5);
//      try
//      {
        xtmp=new double [NN];
          Fx=new double [NN];


//      }
//      catch(xalloc)
//      {
//        printf("Can't allocate memory\n");
//        exit(1);
//      }

      EOSPARAM paar(X,param,NN);
      double F0,Z,F1,fideal;
      //double e[4],s3[4],m,a,xnonp;
      double ro,delta=DELTA,ax,dx /*,tmp*/;
      /* unsigned */ int i;

       norm(paar.XX0,paar.NCmp());
       copy(paar.XX0,xtmp,paar.NCmp());

        paar.ParamMix(xtmp);

        ro=ROTOTALMIX(P,T,paar);
if( ro < 0.0 ) //  Too low pressure - no corrections will be done
  return ( -1. );
        Z=P/(R*T*ro);
        F0=FTOTALMIX(T,ro,paar);


  //       fideal=log(R*T*ro/BARMPA);
         fideal=log(R*T*ro/0.1);
         ax= Z - 1.+fideal;

       for ( i=0;i<paar.NCmp();i++)
       {
         if ( xtmp[i]>0. )
         {
          copy(paar.XX0,xtmp,paar.NCmp());
          dx=xtmp[i]*delta;
          xtmp[i]+=dx;
          norm(xtmp,paar.NCmp());

          paar.ParamMix(xtmp);
          F1=FTOTALMIX(T,ro,paar)*(1.+dx);

          Fx[i]=(F1-F0)/(dx);
         }
         else Fx[i]=0.;
       };

      // GMix=0.;
       for ( i=0;i<paar.NCmp();i++)
       {
         if ( xtmp[i]>0. && Fx[i]< 100. )
         {
   //       tmp=log(paar.XX0[i]);
      //    GMix+=tmp*paar.XX0[i];
          act[i] = exp(ax+Fx[i]);
         }
        else
         {
          act[i]=0.;
         }
       };
    //   GMix+=F0 + ax;

         //MLPutRealList(stdlink,act,paar.NCmp());
        delete [] xtmp;
        delete [] Fx;

         return ro;
   };


   //void ACTDENS(double *data,long nn, double *act )
   int TCGFcalc::CGActivCoefRhoT(double *X,float *param,double *act, unsigned NN,
       double ro, double T )
   {

      //double  act[MAXPARAM];
     // unsigned ncmp;
      //ncmp=unsigned((nn-2)/5);

      //double  T = data[nn - 2];
      //double  ro = data[nn - 1];

       double *Fx,*xtmp;
//      try
//      {
        xtmp=new double [NN];
        Fx=new double [NN];

//      }
//      catch(xalloc)
//      {
//        printf("Cannot allocate memory\n");
//        exit(1);
//      }
      EOSPARAM paar(X,param,NN);

       double   F0,Z,F1,GMix,fideal;
      double delta=DELTA,ax,dx,tmp;
      int i;

       norm(paar.XX0,paar.NCmp());
       copy(paar.XX0,xtmp,paar.NCmp());

        paar.ParamMix(xtmp);
        Z=ZTOTALMIX(T,ro,paar);

        F0=FTOTALMIX(T,ro,paar);
         fideal=log(R*T*ro/0.1);
         ax= Z - 1.+fideal;

       for ( i=0;i<paar.NCmp();i++)
       {
         if ( xtmp[i]>0. )
         {
          copy(paar.XX0,xtmp,NN);
          if ( xtmp[i]>DELTAMOLLIM )
          {
            dx=xtmp[i]*delta;
          }
          else
          {
            dx=DELTAMOLLIM*delta;
          }

          xtmp[i]+=dx;
          norm(xtmp,paar.NCmp());

          paar.ParamMix(xtmp);
          F1=FTOTALMIX(T,ro,paar)*(1.+dx);

          Fx[i]=(F1-F0)/(dx);
         }
         else Fx[i]=0.;
       };

       GMix=0.;
       for ( i=0;i<paar.NCmp();i++)
       {
         if ( xtmp[i]>0. )
         {
          tmp=log(paar.XX0[i]);
          GMix+=tmp*paar.XX0[i];
          act[i] = exp(ax+Fx[i]);
         }
        else
         {
          act[i]=0.;
         }
       };

        delete [] xtmp;
        delete [] Fx;
        return 0;

    //   MLPutRealList(stdlink,act,paar.NCmp());
   };

   double TCGFcalc::DENSITY(double *X,float *param, unsigned NN ,double Pbar, double T )
   {
      double P = Pbar * 0.1;
      //unsigned ncmp;
      //ncmp=unsigned((nn-2)/5);

     // double  P = data[nn - 2];
     // double  T = data[nn - 1];

       double *xtmp;
//      try
//      {
        xtmp=new double [NN];


//      }
//      catch(xalloc)
//      {
//        printf("Can't allocate memory\n");
//        exit(1);
//      }
      EOSPARAM paar(X,param,NN);

      double ro;

       norm(paar.XX0,paar.NCmp());
       copy(paar.XX0,xtmp,paar.NCmp());

        paar.ParamMix(xtmp);
        ro=ROTOTALMIX(P,T,paar);
        delete [] xtmp;
    if( ro < 0. )
       printf(" Error - density cannot be found at this T,P" );
        return ro;
   };

   double TCGFcalc::PRESSURE(double *X,float *param,unsigned NN,double ro, double T)
   {
   //   unsigned ncmp;
   //   ncmp=unsigned((nn-2)/5);

   //   double  T = data[nn - 2];
   //   double ro = data[nn - 1];


      double *xtmp;
//      try
//      {
        xtmp=new double [NN];

//      }
//      catch(xalloc)
//      {
//        printf("Can't allocate memory\n");
//        exit(1);
//      }
      EOSPARAM paar(X,param,NN);

       norm(paar.XX0,paar.NCmp());
       copy(paar.XX0,xtmp,paar.NCmp());

        paar.ParamMix(xtmp);
        double P=PTOTALMIX(T,ro,paar);
        delete [] xtmp;
        return P*10.;
   };


//------------------------------------------------------------ private



void TCGFcalc::copy(double* sours,double *dest,unsigned num)
 {
  unsigned i;
       for ( i=0;i<num;i++)
       {
        dest[i]=sours[i];
       };
 }

void TCGFcalc::norm(double *X,unsigned mNum)
 {
  double tmp=0.;
  unsigned i;
  for ( i=0;i<mNum;i++ )
  {
    tmp+=X[i];
  }
  tmp=1./tmp;
  for ( i=0;i<mNum;i++ )
  {
    X[i]*=tmp;
  }
 }

double TCGFcalc::RPA(double beta,double nuw)
{
  double fi1,fi2;
 fi1=(1.20110+(0.064890+(-76.860+(562.686+(-2280.090+(6266.840+(-11753.40+(14053.8+(-9491.490 +2731.030*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw;
 fi2=(0.588890+(-7.455360+(40.57590+(-104.8970+(60.25470+(390.6310+(-1193.080+(1576.350+(-1045.910+283.7580*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw*nuw;
 return  (-12.*fi1 + 192.*fi2*beta)*beta*beta/PI;
}


double TCGFcalc::dHS(double beta,double ro )
{
/// service constants
   double DV112=1./12.;
   double DV712=7./12.;
/// local variables
   double T12,T112,T712,B13,
          dB,delta,d;
   double a0,a1,a6,a3,a4,a7,a9,a12;
   double p0,p2,p6,p3,p5,p8,p11;
   double dbdl,ri6ro,ri6ro2,d3,d2,dnew,F0,F1;
   unsigned i;

   T12=sqrt(beta);
   T112=exp(DV112*log(beta));
   T712=exp(DV712*log(beta));

   B13=(1+beta);
   B13=B13*B13*B13;

   dB= (P1*T112+PP2*T712+(P3+(P4+P5*beta)*beta)*beta)/B13;
   delta=(P6+P7*T12)/(1.+(P8+(P9+P10*T12)*T12)*T12);

   dbdl=dB*delta;
   ri6ro=PISIX*ro;
   ri6ro2=ri6ro*ri6ro;

   a0=dB+dbdl;
   a1=-1.;
   a3= (-1.5*dB -3.75*dbdl)*ri6ro;
   a4= (1.5*ri6ro);
   a6= (2.*dB + dbdl)*0.25*ri6ro2;
   a7= -0.5*ri6ro2;
   a9= -2.89325*ri6ro2*ri6ro*dbdl;
   a12= -0.755*ri6ro2*ri6ro2*dbdl;

   p0  = -1.;
   p2  = a3*3.;
   p3  = a4*4.;
   p5  = a6*6.;
   p6  = a7*7.;
   p8  = a9*9.;
   p11 = a12*12.;

   d=dB;
   i=0;
   while ( i++<21 )
   {
      d2=d*d;
      d3=d*d*d;

      F0= a0+(a1+(a3+(a4+(a6+(a7+(a9+a12*d3)*d2)*d)*d2)*d)*d2)*d;
      F1= p0+(p2+(p3+(p5+(p6+(p8+p11*d3)*d2)*d)*d2)*d)*d2;

      dnew=d-F0/F1;
      if ( fabs(dnew-d)<1.E-7 )
      {
        return dnew;
      }
      d=dnew;
   }

     if ( i>=20 )
     {
      return dB;
     }

   return dnew;

};


double TCGFcalc::FWCA(double T,double ro)
{
  static double TOld,roOld,F;
  double d,beta,nu,nuw;
  double nu1w1,nu1w2,nu1w3,nu1w4,nu1w5;
  double a0,a1,a2,a3;
  double I2;
  double I1_6,I1_12;
  double dW,dW12,dW6;
  double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
  double F0,F1,FA;
  double rm,rmdw1,rmdw2,rmdw3,rmdw4,rmdw5;


  if ((T==TOld) && (ro==roOld))
  {
    return F;
  }
  else
  {
   TOld=T;
   roOld=ro;
  }

      rm=TWOPOW1SIX;

  beta=1./T;

  d=dHS( beta, ro );

  tmp2=PISIX*d*d*d;
  nu=tmp2*ro;
  tmp1=(1. - nu/16.);
  nuw=nu*tmp1;
  dW=d*exp(1./3.*log(tmp1));


  nu1w1=(1.-nuw);
  nu1w2=nu1w1*nu1w1;
  nu1w3=nu1w2*nu1w1;
  nu1w4=nu1w2*nu1w2;
  nu1w5=nu1w2*nu1w3;

  tmp1=(1-nu);
  tmp1=tmp1*tmp1;
  F0= ((4.-3.*nu)*nu)/tmp1;

  a0=fa0( nuw , nu1w2);
  a1=fa1( nuw , nu1w3);
  a2=fa2( nuw , nu1w4);
  a3=fa3( nuw , nu1w5);

  I1_6 =fI1_6( nuw );
  I1_12=fI1_12( nuw );

  rmdw1=rm/dW;
  rmdw2=rmdw1*rmdw1;
  rmdw3=rmdw1*rmdw2;
  rmdw4=rmdw2*rmdw2;
  rmdw5=rmdw3*rmdw2;

  dW6=dW*dW*dW;
  dW6=1./(dW6*dW6);
  dW12=dW6*dW6;

  tmp1=(a0/4.+ a1/12. + a2/24. + a3/24.)*dW6;
  tmp2=(a0/10.+ a1/90. + a2/720. + a3/5040.)*(-dW12);
  tmp3=(a0 - a1/3. + a2/12 - a3/60)/8.;
  tmp4=(a0 - a1 + a2/2. - a3/6.)*rmdw2*(-9.)/40.;
  tmp5=(a1 - a2 + a3/2)*rmdw3*(-2.)/9.;
  tmp6=(a2 - a3)*rmdw4*(-9.)/64.;
  tmp7=a3*(-3.)/35.*rmdw5;

  I2 = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

  F1=48.*nuw*(I1_12*dW12-I1_6*dW6 + I2)*beta;
  FA=RPA(beta,nuw);


  F=F0+F1+FA;

 return F;
}

 double TCGFcalc::ZWCANum(double T,double ro)
 {
  double delta=DELTA;
  double a0,a1;
  a1=FWCA(T,ro*(1.+delta));
  a0=FWCA(T,ro);
  return 1.+(a1-a0)/delta;
 }

 double TCGFcalc::UWCANum(double T,double ro)
 {
  double delta=DELTA;
  double a0,a1,beta0,beta1;
  beta0=1./T;
  beta1=beta0*(1.+delta);
  a1=FWCA(1./beta1,ro);
  a0=FWCA(T,ro);
  return (a1-a0)/(beta1-beta0);
 }

 double TCGFcalc::FDipPair(double T,double ro,double m2)
 {
  double kappa,Z,U,beta,F;
   kappa=m2*m2/(24.*T);
   beta=1./T;
   Z=ZWCANum(T,ro);
   U=UWCANum(T,ro);
   F=kappa*(4.*beta*U-Z+1.);
   return F;
 }

 double TCGFcalc::J6LJ(double T,double ro)
 {
  double kappa,Z,U,beta,F;
   beta=1./T;
   Z=ZWCANum(T,ro);
   kappa=-16.*PI*ro*beta;
   U=UWCANum(T,ro);
   F=(4.*beta*U-Z+1.)/kappa;
   return F;
 }

 double TCGFcalc::FTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param)
  {
    double FF,A0,A2,A3,AP,A1;
    //unsigned iall,inopol;
    double emix,s3mix,rotmp,T2R;
    double Jind,Jdp;
    int /*itmp,jtmp,ktmp,*/ i,j,k;
    double s3tmp,mtmp,IK /*,atmp*/;
    double imtmp,jmtmp,iatmp,jatmp;
    double m2i,m2j,m2k;
    double s3tmpij,s3tmpik,s3tmpjk;
    double IKtmpij,IKtmpik,IKtmpjk;

    //iall=param.inonzero();
    //inopol=param.inonpolar();
    emix=param.EMIX();
    s3mix=param.S3MIX();

      rotmp=NA*ro_Real;
      T2R=T_Real*T_Real;

      A0=FWCA(T_Real/emix,s3mix*rotmp);


     // if ( inopol< iall )
      {
        /// dipole part
        A2=0.;
        for ( i=0;i<param.NCmp()-1;i++ )
        {
          //itmp=param.ind(i);
          for ( j=i+1; j<param.NCmp(); j++ )
          {
           // jtmp=param.ind(j);

            s3tmp=param.MIXS3(i,j);
            Jdp=J6LJ(T_Real*s3tmp/param.MIXES3(i,j) , s3tmp*rotmp);
            A2+=param.M2R(i)*param.M2R(j)*Jdp*
                           param.X(i)*param.X(j)/s3tmp;
          }
        }
          A2*=2.;
          for ( i=0; i<param.NCmp(); i++ )
          {
            //itmp=param.ind(i);
            mtmp=param.M2R(i);
            s3tmp=param.SIG3(i);
            Jdp=J6LJ(T_Real/param.EPS(i),s3tmp*rotmp);
            A2+=mtmp*mtmp*Jdp*param.X(i)*param.X(i)/s3tmp;
          }

         A2=-A2*TWOPI*rotmp/(3.*T2R);
         ///// A2 done ////

         if ( A2!=0. )
         {

          A3=0.;
          for ( i=0; i<param.NCmp(); i++ )
          {
            //itmp=param.ind(i);
            m2i=param.M2R(i);

            for ( j=0; j<param.NCmp(); j++  )
            {
             // jtmp=param.ind(j);
              m2j=param.M2R(j);

              s3tmpij=param.MIXS3(i,j);
              IKtmpij=K23_13(T_Real*s3tmpij/param.MIXES3(i,j),
                                                    s3tmpij*rotmp);
              for ( k=0; k<param.NCmp(); k++  )
              {
               //ktmp=param.ind(k);
               m2k=param.M2R(k);

               s3tmpik=param.MIXS3(i,k);
               s3tmpjk=param.MIXS3(j,k);

               IKtmpik=K23_13(T_Real*s3tmpik/param.MIXES3(i,k),
                                                    s3tmpik*rotmp);
               IKtmpjk=K23_13(T_Real*s3tmpjk/param.MIXES3(j,k),
                                                    s3tmpjk*rotmp);

               IK=IKtmpij*IKtmpik*IKtmpjk;
               A3+= m2i*m2j*m2k*IK*pow(s3tmpij*s3tmpik*s3tmpjk,-1./3.)*
               param.X(i)*param.X(j)*param.X(k);
              }
            }
          }
            A3=A3*32.*sqrt(14.*PI/5.)*
                  rotmp*rotmp*PI*PI*PI/(135.*T_Real*T2R);
            AP= A2/(1. - A3/A2);
         }
         else AP=0.;

        /// induced interaction
        A1=0.;
        for ( i=0;i<param.NCmp();i++ )
        {
         // itmp=param.ind(i);
          iatmp=param.A(i);
          imtmp=param.M2R(i);
          for ( j=0;j<param.NCmp();j++ )
          {
            //jtmp=param.ind(j);
            jatmp=param.A(j);
            jmtmp=param.M2R(j);

            s3tmp=param.MIXS3(i,j);
            Jind=J6LJ(T_Real*s3tmp/param.MIXES3(i,j),s3tmp*rotmp);

           A1+= (iatmp*jmtmp + jatmp*imtmp)
                  *Jind*param.X(i)*param.X(j)/s3tmp;
          }
        }
        A1=-A1*TWOPI*rotmp/T_Real;
//        A1=-A1*FOURPI*rotmp/T_Real;
      //  A1=0.;
        ///
      }/// end of polar contribution

     FF=A0 + A1 + AP;
     //printf("%g %g %g %g %g",A0,A1,A2,A3,AP);
     //exit(1) ;

    return FF;
  }


double TCGFcalc::UTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param)
{
  double T /*,ro,s3 */;
  double delta=DELTA;
  double a0,a1,beta0,beta1,eps;
  eps=param.EMIX();
  T=T_Real/eps;

  beta0=1./T;
  beta1=beta0*(1.+delta);
  a1=FTOTALMIX((1./beta1)*eps,ro_Real,param);
  a0=FTOTALMIX(T_Real,ro_Real,param);
  return (a1-a0)/(beta1-beta0);
 }


 double TCGFcalc::ZTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param)
 {
  double delta=DELTA;
  double a0,a1;
  a1=FTOTALMIX(T_Real,ro_Real*(1.+delta),param);
  a0=FTOTALMIX(T_Real,ro_Real,param);

  return 1.+(a1-a0)/delta;
 }


 double TCGFcalc::PTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param)
 {
  double Z;
    Z = ZTOTALMIX(T_Real,ro_Real,param);
    return Z*R*T_Real*ro_Real;
 }


 /*  melting density  */
 double Melt(double T)
 {

  return T*0.+.9;

 };

 double Melt2(double T)
 {
  return T*0.+3.;

 };

 #define FIRSTSEED (15)
 #define ROMIN (1.E-2)
 #define NPOINT (5)

 void  choose(double *pres, double P,unsigned &x1,unsigned &x2)
 {
  unsigned i;
  double deltam=-10000000.,tmp;
  double deltap=10000000.;

  for ( i=0; i<NPOINT;i++ )
  {
    tmp=P-pres[i];

    if ( tmp>0. )
    {
       if ( tmp<deltap )
       {
        deltap=tmp;
        x1=i;
       }
    }
    else
    {
       if ( tmp>deltam )
       {
        deltam=tmp;
        x2=i;
       }
    }
  }
     return ;

 }

/////////////////////////////////////////////////////////////////
 double TCGFcalc::ROTOTALMIX(double P,double TT,EOSPARAM& param)
 {
     unsigned i;
     double T /*,ro*/;
     double fact, fact0, romax, dro, roarr[FIRSTSEED];
     double Ptmp[FIRSTSEED], ro0, ro1, rotest, PP0, PP1 /* ,Ptest */;
     double a,b;
     double inttofloat;
     double f[4],x[4],ff,dens[5],pres[5];
     unsigned x1,x2;
//     double ptmp;

     T=TT/param.EMIX();
     fact0=1./(param.S3MIX()*NA);
     fact=R*TT*fact0;

     romax=Melt(T);
     inttofloat=FIRSTSEED-1;
     dro=(romax-ROMIN)/inttofloat;
     roarr[0]=ROMIN;
     roarr[1]=2.*ROMIN;

     for ( i=2;i<FIRSTSEED;i++)
     {
       inttofloat=i;
       roarr[i]=ROMIN+inttofloat*dro;
     }

     for ( i=0;i<FIRSTSEED;i++)
     {
      Ptmp[i]=ZTOTALMIX(TT,roarr[i]*fact0,param);
      Ptmp[i] *= roarr[i] * fact;
      if ( Ptmp[i] > P )
      {
        break;
      }
     }

     if ( i==0 )// Uses aproximation of ideal gas
     {
            return P/(R*TT);
     }


     // additional high pressure inteval
     if ( i==FIRSTSEED )
     {

     //roarr[0]=romax-0.0001;
     roarr[0]=roarr[FIRSTSEED-1];
     Ptmp[0]=Ptmp[FIRSTSEED-1];

     romax=Melt2(T);
     inttofloat=FIRSTSEED-1;
     dro=(romax-ROMIN)/inttofloat;
     for ( i=1;i<FIRSTSEED;i++)
     {
       inttofloat=i;
       roarr[i]=ROMIN+inttofloat*dro;
     }

     for ( i=1;i<FIRSTSEED;i++)
     {
      Ptmp[i]=ZTOTALMIX(TT,roarr[i]*fact0,param)*roarr[i]*fact;
      if ( Ptmp[i]>P )
      {
        break;
      }
     }

     if ( i==FIRSTSEED || i==0 )
     {
         printf("Input pressure is too high!\n");
//         exit(1);
         return (-1.0);
     }
     }

     ro0=roarr[i-1];
     ro1=roarr[i];
     PP0=Ptmp[i-1];
     PP1=Ptmp[i];
     i=0;

   while ( i++<20 )
   {
     //Start interp
     ff = ro0;
     dens[0]=ro0;
     dens[1]=ro1;
     pres[0]=PP0;
     pres[1]=PP1;
     //first order
     x[0]=P-pres[0];
     f[0]=(dens[1]-dens[0])/(pres[1]-pres[0]);
     ff+=f[0]*x[0];

     //second order
     dens[2]=ff;
     pres[2]=ZTOTALMIX(TT,ff*fact0,param)*ff*fact;

     if ( fabs(pres[2]-P)<1E-5 )
     {
       return ff*fact0;
     }

     x[1]=x[0]*(P-pres[1]);
     f[1]=(dens[2]-dens[1])/(pres[2]-pres[1]);

     f[0]=(f[1]-f[0])/(pres[2]-pres[0]);
     ff+=f[0]*x[1];

     //third order
     dens[3]=ff;
     pres[3]=ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[3]-P)<1E-6 )
     {
      return ff*fact0;
     }
     x[2]=x[1]*(P-pres[2]);
     f[2]=(dens[3]-dens[2])/(pres[3]-pres[2]);
     f[1]=(f[2]-f[1])/(pres[3]-pres[1]);
     f[0]=(f[1]-f[0])/(pres[3]-pres[0]);
     ff+=f[0]*x[2];
     dens[4]=ff;
     pres[4]=ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[4]-P)<1e-6 )
     {
      return ff*fact0;
     }

     choose(pres,P,x1,x2);

     ro0=dens[x1];
     ro1=dens[x2];
     PP0=pres[x1];
     PP1=pres[x2];

      if ( fabs((ro1-ro0))<0.001 )
      {
          a=(PP1-PP0)/(ro1-ro0);
          b=PP1-a*ro1;
          rotest=(P-b)/a;
          return rotest*(fact0);
      }
   }
        //return 10.;
         /// bad result

          a=(PP1-PP0)/(ro1-ro0);
          b=PP1-a*ro1;
          rotest=(P-b)/a;
          return rotest*(fact0);

 }


// Implementation of EOSPARAM class

EOSPARAM::~EOSPARAM()
{
  unsigned  i;
  if ( isize > 0)
  {
         for ( i=0;i<isize;i++ )
         {
           delete [] mixpar[i];
         }

   delete [] epspar;
   delete [] sig3par;
   delete [] XX;
   delete [] eps;
   delete [] eps05;
   delete [] sigpar;
   delete [] mpar;
   delete [] apar;
   delete [] aredpar;
   delete [] m2par;
   delete [] XX0;

   delete [] mixpar;

  }
}


void EOSPARAM::allocate(unsigned inew)
{
   unsigned i;
  if ( (isize > 0) && (inew > isize) )
  {
         for ( i=0;i<isize;i++ )
         {
           delete [] mixpar[i];
         }

   delete [] epspar;
   delete [] sig3par;
   delete [] XX;
   delete [] eps;
   delete [] eps05;
   delete [] sigpar;
   delete [] mpar;
   delete [] apar;
   delete [] aredpar;
   delete [] m2par;
   delete [] XX0;

   delete [] mixpar;

  }

  if ( (inew > isize) )
  {

//    try{
           mixpar=new   double* [inew];
         for ( i=0;i<inew;i++ )
         {
           mixpar[i]=new   double [inew];
         }

   epspar =new double [inew];
   sig3par=new double [inew];
   XX     =new double [inew];
   eps    =new double [inew];
   eps05  =new double [inew];
   sigpar =new double [inew];
   mpar   =new double [inew];
   apar   =new double [inew];
   aredpar=new double [inew];
   m2par  =new double [inew];
   XX0    =new double [inew];
//  }
//  catch (xalloc)
//  {
//    printf("Can't allocate mamory\n");
//    exit(1);
//  }

  isize=inew;
  }
}

unsigned EOSPARAM::ParamMix(double *Xin)
  {
    /* unsigned */ int j,i;
    double tmp,tmp1,tmp2;
    for ( i=0; i<NComp; i++ ) XX[i]=Xin[i];

    emix=0.;
    s3mix=0.;

    for ( i=0;i<NComp-1;i++ )
    {
      for ( j=i+1;j<NComp;j++ )
      {
          tmp=XX[i]*XX[j];
          tmp2=mixpar[j][i]; //eps
          tmp1=mixpar[i][j]; //signa
          s3mix+= tmp1*tmp;
          emix+= tmp2*tmp;
      }
    }

    s3mix*=2.;
    emix*=2.;

    for ( i=0;i<NComp;i++ )
    {
          tmp=XX[i]*XX[i];

          s3mix+= sig3par[i]*tmp;
          emix +=  epspar[i]*tmp;
    }
    emix=emix/s3mix;
    return NComp;
  }

  void EOSPARAM::PureParam(double* e,double* s,double* m,double* a)
  {
    /* unsigned */ int i;

    for ( i=0;i<NComp;i++ )
    {
      e[i]=eps[i];
      s[i]=sigpar[i];
      m[i]=mpar[i];
      a[i]=aredpar[i];
    }

  }

 ////////////////////////////////////////////////////////////////////
void EOSPARAM::init(double *Xinp, float * data,unsigned ncmp)
{
  int i,j;
  double tmp;

  //if ( ncmp>MAXPARAM ) NComp=MAXPARAM;
   allocate(ncmp);
   NComp=(int)ncmp;

  for ( i=0;i<NComp;i++ )
  {
       XX0[i] = Xinp[i];

    sigpar[i] = data[i*4    ];
       eps[i] = data[i*4 + 1];
      mpar[i] = data[i*4 + 2];
      apar[i] = data[i*4 + 3];

  }

    for ( i=0;i<NComp;i++ )
    {
      tmp=sigpar[i];
      tmp=tmp*tmp*tmp;
      sig3par[i]=tmp;
      eps05[i]=sqrt(eps[i]);
      epspar[i]=tmp*eps[i];
      m2par[i]=mpar[i]*mpar[i]/(1.38048E-4);
      aredpar[i]=apar[i]/tmp;
    }

  /// calculation of mixing properties ///
    for ( i=0;i<NComp-1;i++ )
    {
     for ( j=i+1;j<NComp;j++ )
     {
       tmp=(sigpar[i]+sigpar[j])*0.5;
       tmp=tmp*tmp*tmp;
       mixpar[i][j]=tmp;
       mixpar[j][i]=tmp*eps05[i]*eps05[j];
     }
    }
};

/////////////////////////////////////////////////////////////////////

void EOSPARAM::copy(double* sours,double *dest,unsigned num)
 {
  unsigned i;
       for ( i=0;i<num;i++)
       {
        dest[i]=sours[i];
       };
 }

void EOSPARAM::norm(double *X,unsigned mNum)
 {
  double tmp=0.;
  unsigned i;
  for ( i=0;i<mNum;i++ )
  {
    tmp+=X[i];
  }
  tmp=1./tmp;
  for ( i=0;i<mNum;i++ )
  {
    X[i]*=tmp;
  }
 }

// -------------------------------------------------------------------------

// Implementation
// TPRSVcalc class - private methods
int
TPRSVcalc::PureParam( double *Eos2parPT )
{ // calculates a and b arrays
	// calculates a and b parameters of pure species
   int i;
   double Tcrit, Pcrit, omg, k1, k2, k3, apure, bpure;

   for (i=0; i<NComp; i++)
   {
      Tcrit = Eosparm[i][0];
      Pcrit = Eosparm[i][1];
      omg = Eosparm[i][2];
      k1 = Eosparm[i][3];
      k2 = Eosparm[i][4];
      k3 = Eosparm[i][5];

      apure = A(Tcrit, omg, k1, k2, k3, Pcrit);
      bpure = B(Tcrit, Pcrit);
      Pureparm[i][0] = apure;
      Pureparm[i][1] = bpure;
      Eos2parPT[0] = apure;
      Eos2parPT[1] = bpure;
   }
   return 0;

}

double
TPRSVcalc::A(double Tcrit, double omg, double k1, double k2, double k3, double Pcrit )
{
 // calculates a term of cubic EoS
  double Tred, k0, k, alph, aprsv;

  Tred = Tk/Tcrit;
  k0 = 0.378893 + 1.4897153*omg - 0.17131848*pow(omg,2.) + 0.0196554*pow(omg,3.);
  if(Tk < Tcrit)
  	k = k0 + (k1 + k2*(k3-Tred)*(1.-sqrt(Tred))) * (1.+sqrt(Tred)) * (0.7-Tred);
  else
  	k = k0;
  alph = pow(1. + k*(1.-sqrt(Tred)), 2.);
  aprsv = alph*(0.457235*pow(R_CONSTANT,2.)*pow(Tcrit,2.) / Pcrit);
  return aprsv;
}

double
TPRSVcalc::B(double Tcrit, double Pcrit)
{
    double bprsv;

    bprsv = 0.077796*R_CONSTANT*Tcrit/Pcrit;
    return bprsv;
}

int
TPRSVcalc::FugacityPure( )
{ // Calculates the fugacity of pure species
// calculates fugacity and state functions of pure species
    int i;
	double Tcrit, Pcrit, Tred, aprsv, bprsv, alph, k, aa, bb, a2, a1, a0,
               z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, gdep, hdep, sdep, g, h, s, fugpure;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONSTANT*log(P);
	gig = hig - Tk*sig;

	for (i=0; i<NComp; i++)
	{
		// calculate a and b terms of cubic EoS
		Tcrit = Eosparm[i][0];
		Pcrit = Eosparm[i][1];
		Tred = Tk/Tcrit;
		aprsv = Pureparm[i][0];
		bprsv = Pureparm[i][1];
		// solve cubic equation
		aa = aprsv*P/(pow(R_CONSTANT,2.)*pow(Tk,2.));
		bb = bprsv*P/(R_CONSTANT*Tk);
		a2 = bb - 1.;
		a1 = aa - 3.*pow(bb,2.) - 2.*bb;
		a0 = pow(bb,3.) + pow(bb,2.) - aa*bb;
		Cardano(a2, a1, a0, z1, z2, z3);

		// find stable roots
		vol1 = z1*R_CONSTANT*Tk/P;
		vol2 = z2*R_CONSTANT*Tk/P;
		vol3 = z3*R_CONSTANT*Tk/P;
		if (z1 > bb)
			lnf1 = (-1.)*log(z1-bb)
    - aa/(bb*sqrt(8.))*log((z1+(1.+sqrt(2.))*bb)/(z1+(1.-sqrt(2.))*bb))+z1-1.;
		else
			lnf1 = 1000.;
		if (z2 > bb)
			lnf2 = (-1.)*log(z2-bb)
    - aa/(bb*sqrt(8.))*log((z2+(1.+sqrt(2.))*bb)/(z2+(1.-sqrt(2.))*bb))+z2-1.;
		else
			lnf2 = 1000.;
		if (z3 > bb)
			lnf3 = (-1.)*log(z3-bb)
    - aa/(bb*sqrt(8.))*log((z3+(1.+sqrt(2.))*bb)/(z3+(1.-sqrt(2.))*bb))+z3-1.;
		else
			lnf3 = 1000.;

		if (lnf2 < lnf1)
		{
			z = z2; vol = vol2; lnf = lnf2;
		}
		else
		{
			z = z1; vol = vol1; lnf = lnf1;
		}
		if (lnf3 < lnf)
		{
			z = z3; vol = vol3; lnf = lnf3;
		}
		else
		{
			z = z; vol = vol; lnf = lnf;
		}
		// calculate thermodynamic properties
		alph = aprsv/(0.457235*pow(R_CONSTANT,2.)*pow(Tcrit,2.) / Pcrit);
		k = (sqrt(alph)-1.)/(1.-sqrt(Tred));
		gdep = R_CONSTANT*Tk*(z-1.-log(z-bb)-aa/(bb*sqrt(8.))
                       *log((z+(1+sqrt(2.))*bb)/(z+(1-sqrt(2.))*bb)));
		hdep = R_CONSTANT*Tk*(z-1.-log((z+(1+sqrt(2.))*bb)/(z+(1-sqrt(2.))*bb))
                       *aa/(bb*sqrt(8.))*(1+k*sqrt(Tred)/sqrt(alph)));
		sdep = (hdep-gdep)/Tk;
		g = gig + gdep;
		h = hig + hdep;
		s = sig + sdep;
		fugpure = exp(lnf);
		Fugpure[i][0] = fugpure;
		Fugpure[i][1] = g;
		Fugpure[i][2] = h;
		Fugpure[i][3] = s;
                Fugpure[i][4] = vol;
	}
        return 0;
}

int
TPRSVcalc::Cardano(double a2, double a1, double a0, double &z1, double &z2, double &z3)
{
   // finds roots of cubic equation
   double q, rc, q3, rc2, theta, ac, bc;

   q = (pow(a2,2.) - 3.*a1)/9.;
   rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
   q3 = pow(q,3.);
   rc2 = pow(rc,2.);
   if (rc2 < q3)   // three real roots
   {
      theta = acos(rc/sqrt(q3));
       z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
       z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
       z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
   }
   else   // one real root
   {
  	ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
  	if (ac != 0.)
        	bc = q/ac;
   	else
  		bc = 0.;
    	z1 = ac+bc-a2/3.;
   	z2 = ac+bc-a2/3.;
   	z3 = ac+bc-a2/3.;
   }
   return 0;
}

int
TPRSVcalc::MixParam( double &amix, double &bmix)
{;
	// calculates a and b parameters of the mixture
	int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// K = KK0ij[i][j] + KK1ij[i][j]*Tk;
                        K = KK0ij[i][j];
			AAij[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}
	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + Wx[i]*Wx[j]*AAij[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + Wx[i]*Pureparm[i][1];
	}
	return 0;
}

int
TPRSVcalc::FugacityMix( double amix, double bmix,
    double &fugmix, double &zmix, double &vmix)
{
	// calculates fugacity of the mixture
    double aa, bb, a2, a1, a0, z1, z2, z3;
	double vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	aa = amix*P/(pow(R_CONSTANT,2.)*pow(Tk,2.));
	bb = bmix*P/(R_CONSTANT*Tk);
	a2 = bb - 1.;
	a1 = aa - 3.*pow(bb,2.) - 2.*bb;
	a0 = pow(bb,3.) + pow(bb,2.) - aa*bb;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONSTANT*Tk/P;
	vol2 = z2*R_CONSTANT*Tk/P;
	vol3 = z3*R_CONSTANT*Tk/P;
	if (z1 > bb)
		lnf1 = (-1.)*log(z1-bb)
    - aa/(bb*sqrt(8.))*log((z1+(1.+sqrt(2.))*bb)/(z1+(1.-sqrt(2.))*bb))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > bb)
		lnf2 = (-1.)*log(z2-bb)
    - aa/(bb*sqrt(8.))*log((z2+(1.+sqrt(2.))*bb)/(z2+(1.-sqrt(2.))*bb))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > bb)
		lnf3 = (-1.)*log(z3-bb)
    - aa/(bb*sqrt(8.))*log((z3+(1.+sqrt(2.))*bb)/(z3+(1.-sqrt(2.))*bb))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}
	fugmix = exp(lnf);
        PhVol = vmix;
	return 0;
}

// #define MAXPUREPARAM 7
int
TPRSVcalc::FugacitySpec( double *fugpure, float *binpar, float *params  )
{
    // calculates fugacity and activity of species
    int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double au, bu, lnfci, fci;

    // Reload params to Pureparm
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/P;
//      for( i=0; i<4; i++ )
//        Pureparm[j][i] = (double)params[j*4+i];
      for( i=0; i<4; i++ )  // Temporary
      {
           Pureparm[j][i] = (double)params[j*10+i+6];
      }
    }

/*
    for( j=0; j<NComp; j++ )
      for( i=0; i<NComp; i++ )
        KK0ij[j][i] = (double)binpar[j*NComp+i];
*/
	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix);
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		au = amix*P/(pow(R_CONSTANT, 2.)*pow(Tk, 2.));
		bu = bmix*P/(R_CONSTANT*Tk);
		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + Wx[j]*AAij[i][j];
		}
		lnfci = Pureparm[i][1]/bmix*(zmix-1.) - log(zmix-bu)
		      + au/(sqrt(8.)*bu)*(2.*sum/amix-Pureparm[i][1]/bmix)
                      * log((zmix+bu*(1.-sqrt(2.)))/(zmix+bu*(1.+sqrt(2.))));
		fci = exp(lnfci);
		Fugci[i][0] = fci;	// fugacity coefficient using engineering convention
		Fugci[i][1] = Wx[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (Wx[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/Wx[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;

	}
	return iRet;
}

int
TPRSVcalc::GetEosParam( float *EoSparam )
{
   // reads EoS parameters from database into work array
   int i;
   if( !EoSparam )
     return -1;  // Memory alloc error

   for (i=0; i<NComp; i++)
   {
      Eosparm[i][0] = (double)EoSparam[0];   // critical temperature in K
      Eosparm[i][1] = (double)EoSparam[1];   // critical pressure in bar
      Eosparm[i][2] = (double)EoSparam[2];   // Pitzer acentric factor omega
      Eosparm[i][3] = (double)EoSparam[3];   // empirical EoS parameter k1
      Eosparm[i][4] = (double)EoSparam[4];   // empirical EoS parameter k2
      Eosparm[i][5] = (double)EoSparam[5];   // empirical EoS parameter k3
    }
    return 0;
}

int
TPRSVcalc::GetMoleFract( double *X )
{
  ; // Loads mole fractions for NComp species
  for( int j = 0; j< NComp; j++ )
    Wx[j] = X[j];
  return 0;
}

double
TPRSVcalc::ObtainResults( double *ActCoef )
{
    int j;

    for(j=0; j<NComp; j++)
       ActCoef[j] = Fugci[j][3];

    return PhVol;
}


// - - - - - - -

// TPRSVcalc class - high-level methods
// Constructor

TPRSVcalc::TPRSVcalc( int NCmp, double Pp, double Tkp )
    {
       int i;

       NComp = NCmp;
       P = Pp;
       Tk = Tkp;
       R_CONSTANT = 8.31451;

       Wx = new double [NComp];
       Eosparm = new double *[NComp];
       Pureparm = new double *[NComp];
       Fugpure = new double *[NComp];
       Fugci = new double *[NComp];

       for (i=0; i<NComp; i++)
       {
         Eosparm[i] = new double [6];
       }
       for (i=0; i<NComp; i++)
       {
         Pureparm[i] = new double [4];
       }
       for (i=0; i<NComp; i++)
       {
         Fugpure[i] = new double [5];
       }
       for (i=0; i<NComp; i++)
       {
         Fugci[i] = new double [4];
       }
//
       KK0ij = new double *[NComp];
       KK1ij = new double *[NComp];
       AAij = new double *[NComp];

       for (i=0; i<NComp; i++)
       {
		KK0ij[i] = new double[NComp];
       }
       for (i=0; i<NComp; i++)
       {
        	KK1ij[i] = new double[NComp];
       }
       for (i=0; i<NComp; i++)
       {
		AAij[i] = new double[NComp];
       }
    }

TPRSVcalc::~TPRSVcalc()
    {
        int i;

        delete[]Wx;

       for (i=0; i<NComp; i++)
       {
         delete[]Eosparm[i];
       }
       for (i=0; i<NComp; i++)
       {
         delete[]Pureparm[i];
       }
       for (i=0; i<NComp; i++)
       {
         delete[]Fugpure[i];
       }
       for (i=0; i<NComp; i++)
       {
         delete[]Fugci[i];
       }

       for (i=0; i<NComp; i++)
       {
	  delete[]KK0ij[i];
       }
       for (i=0; i<NComp; i++)
       {
          delete[]KK1ij[i];
       }
       for (i=0; i<NComp; i++)
       {
	  delete[]AAij[i];
       }

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
        delete[]KK0ij;
	delete[]KK1ij;
	delete[]AAij;
    }

 // Implementation of the TPRSVcalc class

int
TPRSVcalc::PRFugacityPT( double P, double Tk, float *EoSparam, double *Eos2parPT,
        double &Fugacity, double &Volume, double &DeltaH, double &DeltaS )
 {

      int iRet = 0;

      iRet = GetEosParam( EoSparam );
      if( iRet)
        return iRet;

      iRet = PureParam( Eos2parPT ); // Calculates A and B
                                     // - attraction and repulsion terms
      if( iRet)
        return iRet;

      iRet = FugacityPure();
      if( iRet)
        return iRet;

      Fugacity = Fugpure[0][0]; // Fugacity coefficient
      DeltaH = Fugpure[0][2];   //
      DeltaS = Fugpure[0][3];   // Entropy includes ideal gas contribution
      Volume = Fugpure[0][4];   //  J/bar

      return iRet;
 }

 // Called from IPM-Gamma() where activity coefficients are computed
int
TPRSVcalc::PRActivCoefPT( int NComp, double Pbar, double Tk, double *X,
    double *fugpure, float *binpar, float *param, double *act, double &PhaseVol,
    int NPar, int NPcoef, int MaxOrd, short *aIPx  )
{

   int iRet;
   int j, i, ip;
   int index1, index2;

   if( NPcoef > 0 )
   {
      // fill internal array of interaction parameters with standard value
      for( j=0; j<NComp; j++ )
        for( i=0; i<NComp; i++ )
          KK0ij[j][i] = 0.;

      // transfer those interaction parameters that have non-standard value
      for ( ip=0; ip<NPar; ip++ )
      {
         index1 = (int)aIPx[MaxOrd*ip];
         index2 = (int)aIPx[MaxOrd*ip+1];
	 KK0ij[index1][index2] = binpar[NPcoef*ip];
	 KK0ij[index2][index1] = binpar[NPcoef*ip];	// symmetric case
      }
    }

    GetMoleFract( X );

    iRet = FugacitySpec( fugpure, binpar, param );

    PhaseVol = ObtainResults( act );

    return iRet;
}

//--------------------- End of s_fgl.cpp ---------------------------
