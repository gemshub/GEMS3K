//-------------------------------------------------------------------
// $Id: s_fgl.h 858 2007-02-16 17:11:38Z gems $
//
// Copyright (C) 2003-2007  S.Churakov, Th.Wagner, D.Kulik, S.Dmitrieva
//
// Declaration of new implementation of Fluid EoS classes
// (PRSV and CG EoS models)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//

#ifndef _s_fgl_h_
#define _s_fgl_h_

// Added 09 May 2003
// Definition of a class for CG EOS calculations for fluids
// Incorporates a C++ program written by Sergey Churakov (CSCS ETHZ)
// implementing a paper Churakov & Gottschalk 2003 GCA v.   p.


/*-----------------07.11.99 12:14-------------------
 Structure of parameers for LJ and Shokmayer potential
 sig(A)=F1+F5*T*1E-4
 epsilon(K)=F2+F3*Exp(-F4*T*1E-3)
 polaris(A)=F6
 mu(debue)=F7+F8*Exp(F9*T)
--------------------------------------------------*/
//#define MOLNUM (12)
//#define MOLPOL (6)
//#define MOLNONPOL (6)
/* Ar CO2 CH4 O2 N2 H2 CO SO2 HCl H2O H2S NH3 */

class EOSPARAM
{
  //static double Told;
  unsigned isize;  // int isize;
  double emix,s3mix;
  int NComp;
   double *epspar,*sig3par;
   double *XX;
   double *eps;
   double *eps05;
   double *sigpar;
   double *mpar;
   double *apar;
   double *aredpar;
   double *m2par;
   double **mixpar;
  void init(double*,float *, unsigned );
  void allocate(unsigned );

  void copy(double* sours,double *dest,unsigned num);
  void norm(double *X,unsigned mNum);

  public:
  double *XX0;

  EOSPARAM():isize(0),emix(0),s3mix(0),NComp(0){};
    //EOSPARAM(double*data, unsigned nn):isize(0){allocate(nn);init(data,nn);};

    EOSPARAM(double *Xtmp, float *data, unsigned nn)
    :isize(0),emix(0),s3mix(0),NComp(0)
         {allocate(nn);init(Xtmp,data,nn);};

  ~EOSPARAM();

  double EPS05(unsigned i)           {return eps05[i];};
  double X(unsigned i)               {return XX[i];};
  double EPS(unsigned i)             {return eps[i];};
  double EMIX(void){return emix;};
  double S3MIX(void){return s3mix;};
  double MIXS3(unsigned i,unsigned j)
  {
    if (i==j) return sig3par[i];
    if (i<j) return mixpar[i][j]; else return mixpar[j][i];
  };

  double MIXES3(unsigned i,unsigned j)
  {
    if ( i==j ) return epspar[i];
    if (i<j) return mixpar[j][i]; else return mixpar[i][j];
  };

  double SIG3(unsigned i){return sig3par[i];};
  double M2R(unsigned i) {return m2par[i];};
  double A(unsigned i)  {return apar[i];};

  unsigned ParamMix( double *Xin);
  void PureParam(double* ,double*,double*,double* );
  int NCmp(){return NComp;};
};


class TCGFcalc // Churakov & Gottschalk (2003) EOS calculations
{
  private:

  const double
          PI,    // pi
          TWOPI,    // 2.*pi
          PISIX,    // pi/6.
          TWOPOW1SIX,   // 2^(1/6)
          DELTA,
          DELTAMOLLIM,
          R,  NA,  P1,
          PP2, P3, P4,
          P5,  P6, P7,
          P8,  P9, P10,
          AA1, AA2, AA3,
          A4, A5, A6,
          BB1, BB2, BB3,
          B4,  B5,  B6,
          A00, A01, A10,
          A11, A12, A21,
          A22, A23, A31,
          A32, A33, A34;

  public:

        TCGFcalc():
         PI( 3.141592653589793120),    // pi
         TWOPI( 6.283185307179586230),    // 2.*pi
         PISIX( 0.523598775598298927),    // pi/6.
         TWOPOW1SIX(1.12246204830937302),   // 2^(1/6)
         DELTA (0.00001),
         DELTAMOLLIM (0.0000001),
         R(8.31439),
         NA(0.6023),
         P1(1.186892378996),
         PP2(-0.4721963005527),
         P3(3.259515855283),
         P4(3.055229342609),
         P5(1.095409321023),
         P6(1.282306659774E-2),
         P7(9.55712461425E-2),
         P8(13.67807693107),
         P9(35.75464856619),
         P10(16.04724381643),
         AA1(-0.120078459237),
         AA2(-.808712488307),
         AA3(.321543801337),
         A4(1.16965477132),
         A5(-.410564939543),
         A6(-.516834310691),
         BB1(-2.18839961483),
         BB2(1.59897428009),
         BB3(-.392578806128),
         B4(-.189396607904),
         B5(-.576898496254),
         B6(-0.0185167641359),
         A00(.9985937977069455),
         A01(.5079834224407451),
         A10(1.021887697885469),
         A11(-5.136619463333883),
         A12(-5.196188074016755),
         A21(-6.049240839050804),
         A22(18.67848155616692),
         A23(20.10652684217768),
         A31(9.896491419756988),
         A32(14.6738380473899),
         A33(-77.44825116542995),
         A34(-4.82871082941229)
      {}

protected:

   void copy(double* sours,double *dest,unsigned num);
   void norm(double *X,unsigned mNum);
   double RPA(double beta,double nuw);
   double dHS(double beta,double ro );
   inline double fI1_6(double nuw)
{
   return (1.+(A4+(A5+A6*nuw)*nuw)*nuw)/
   ((1.+(AA1+(AA2+AA3*nuw)*nuw)*nuw)*3.);
};

inline double fI1_12(double nuw)
{
   return (1.+(B4+(B5+B6*nuw)*nuw)*nuw)/
   ((1.+(BB1+(BB2+BB3*nuw)*nuw)*nuw)*9.);
};


inline double fa0(double nuw ,double nu1w2)
{
     return (A00 + A01*nuw)/nu1w2;
};

inline double fa1(double nuw ,double nu1w3)
{
     return (A10+(A11+A12*nuw)*nuw)/nu1w3;
};

inline double fa2(double nuw ,double nu1w4)
{
     return ((A21+(A22+A23*nuw)*nuw)*nuw)/nu1w4;
};

inline double fa3(double nuw ,double nu1w5)
{
     return ((A31+(A32+(A33+A34*nuw)*nuw)*nuw)*nuw)/nu1w5;
};

   double DIntegral(double T, double ro, unsigned IType);
   double LIntegral(double T, double ro, unsigned IType);
   double KIntegral(double T, double ro, long unsigned IType);
   double K23_13(double T, double ro);
   double J6LJ(double T,double ro);
   double FDipPair(double T,double ro,double m2);
   double UWCANum(double T,double ro);
   double ZWCANum(double T,double ro);

   double FWCA(double T,double ro);
   double FTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param);
   double UTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param);
   double ZTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param);
   double PTOTALMIX(double T_Real,double ro_Real,EOSPARAM& param);
   double ROTOTALMIX(double P,double TT,EOSPARAM& param);

public:
    //
   double PRESSURE(double *X, float *param, unsigned NN, double ro, double T );
   double DENSITY(double *X,float *param, unsigned NN ,double Pbar, double T );
   int CGActivCoefRhoT(double *X,float *param, double *act, unsigned NN,
     double ro, double T );
   double CGActivCoefPT(double *X,float *param,double *act, unsigned NN,
     double Pbar, double T );

    int CGcalcFug( void );  // Calc. fugacity for 1 species at X=1
    int CGFugacityPT( float *EoSparam, float *EoSparPT, double &Fugacity,
        double &Volume, double &DeltaH, double &DeltaS, double P, double T );

};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Added 19 July 2006 by Th.Wagner and D.Kulik
// Definition of a class for PRSV EOS calculations for fluids
// Incorporates a C++ program written by Thomas Wagner (Univ. Tuebingen)
//
//
class TPRSVcalc // Peng-Robinson-Styjek-Vera EOS calculations
{
  private:
     double R_CONSTANT;
     int NComp; // number of species;
//     int i;
     double P, Tk, PhVol;   // bar, T Kelvin, phase volume in cm3
     // main work arrays
     double *Wx;         // mole fractions of species
     double **Eosparm;   // EoS parameters
     double **Pureparm;  // Parameters a and b for cubic EoS
     double **Fugpure;   // Fugacity parameters of pure gas species
     double **Fugci;     // Fugacity parameters of species in the mixture

     double **KK0ij;    //  Constant term of the binary interaction parameter
     double **KK1ij;    //  T-dependent term
     double **AAij;     //  binary a terms in the mixture

    public:

    TPRSVcalc( int NCmp, double Pp, double Tkp );

    ~TPRSVcalc();

   // Called from IPM-Gamma() where activity coefficients are computed
   int PRActivCoefPT( int NComp, double Pbar, double Tk, double *X,
        double *fugpure, float *binpar, float *param, double *act, double &PhaseVol );

   int CalcFugPure( void );
   // Calc. fugacity for 1 species at X=1
   int PRFugacityPT( double P, double Tk, float *EoSparam, double *Eos2parPT,
        double &Fugacity, double &Volume, double &DeltaH, double &DeltaS );

protected:

int PureParam( double *params ); // calculates a and b arrays
double A(double Tcrit, double omg, double k1, double k2, double k3, double Pcrit);
double B(double Tcrit, double Pcrit);
int FugacityPure( void ); // Calculates the fugacity of pure species
int Cardano(double a2, double a1, double a0, double &z1, double &z2, double &z3);
int MixParam( double &amix, double &bmix);
int FugacityMix( double amix, double bmix,
     double &fugmix, double &zmix, double &vmix);
int FugacitySpec( double *fugpure, float *binpar, float *params  );

int GetEosParam( float *params ); // Loads EoS parameters for NComp species
int GetMoleFract( double *Wx ); // Loads mole fractions for NComp species
double ObtainResults( double *ActCoef ); // returns activity coeffs and phase volume

};



#endif  // _v_fgl_h
