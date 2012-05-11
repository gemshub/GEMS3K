//-------------------------------------------------------------------
// $Id: s_fgl.h 1445 2009-10-15 12:42:54Z gems $
//
// Copyright (C) 2003-2010  T.Wagner, D.Kulik, S.Dmitrieva, S.Churakov
//
// Declaration of new versions of fluid, liquid, aquous
// and solid-solution models
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// and part of the standalone GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//------------------------------------------------------------------
//

#include <string.h>
#ifndef _s_fgl_h_
#define _s_fgl_h_

enum fluid_mix_rules {  // Codes to identify specific mixing rules in EoS models
    MR_WAAL_ = 'W',		// Basic Van der Waals mixing rules in cubic EoS models
    MR_CONST_ = 'C',    // Constant one-term interaction parameter kij
    MR_TEMP_ = 'T'		// Temperature-dependent one-term interaction parameter kij (Jaubert et al. 2005)
};

enum dc_class_codes {  // codes to identify fluid types in EoS models
    DC_GAS_H2O_ = 'V',
    DC_GAS_CO2_ = 'C',
    DC_GAS_H2_ = 'H',
    DC_GAS_N2_ = 'N',
    DC_GAS_COMP_ = 'G'
};


// ------------------------------------------------------------------
// base class for subclasses of built-in mixing models
// (c) March 2007 DK/TW

#define MAXPHASENAME 16

class TSolMod
{
	protected:
		char ModCode;   // Code of the mixing model
		char MixCode;	// Code for specific EoS mixing rules
        char PhaseName[MAXPHASENAME+1];    // Phase name (for specific built-in models)

        long int NComp;    	// Number of components in the solution phase
        long int NPar;     	// Number of non-zero interaction parameters
        long int NPcoef;   	// Number of coeffs per parameter (columns in the aIPc table)
        long int MaxOrd;   	// max. parameter order (or number of columns in aIPx)
        long int NP_DC;    	// Number of coeffs per one DC in the phase (columns in aDCc)
        long int NPTP_DC;   // Number of properties per one DC at T,P of interest (columns in aDC)
        long int *aIPx;  	// Pointer to list of indexes of non-zero interaction parameters

        double R_CONST; // R constant
        double Tk;    	// Temperature, K
        double Pbar;  	// Pressure, bar

        double *aIPc;  // Table of interaction parameter coefficients
        double *aIP;   // Vector of interaction parameters corrected to T,P of interest
        double *aDCc;  // End-member properties coefficients
        double **aDC;  // Table of corrected end member properties at T,P of interest
        double *x;     // Pointer to mole fractions of end members (provided)
        double *phVOL; // phase volumes, cm3/mol (now obsolete)

        // Results
        double Gam;   	// work cell for activity coefficient of end member
        double lnGamRT;
        double lnGam;
        double Gex, Hex, Sex, CPex, Vex, Aex, Uex;   // molar excess properties
        double Gid, Hid, Sid, CPid, Vid, Aid, Uid;   // molar ideal mixing properties
        double Gdq, Hdq, Sdq, CPdq, Vdq, Adq, Udq;   // molar Darken quadratic terms
        double Grs, Hrs, Srs, CPrs, Vrs, Ars, Urs;   // molar residual functions (fluids)
        double *lnGamma;   // Pointer to ln activity coefficients of end members
                           // (memory must be provided from the calling program)

	public:

		// Generic constructor
		TSolMod( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
					long int NPperDC, long int NPTPperDC, char Mod_Code, char Mix_Code,
					long int* arIPx, double* arIPc, double* arDCc, double *arWx,
					double *arlnGam, double *aphVOL, double T_k, double P_bar );

		// Destructor
		virtual ~TSolMod();

		virtual long int PureSpecies()
		{
			return 0;
		};

		virtual long int PTparam()
		{
			return 0;
		};

		virtual long int MixMod()
		{
			return 0;
		};

		virtual long int ExcessProp( double *Zex )
		{
			return 0;
		};

		virtual long int IdealProp( double *Zid )
		{
			return 0;
		};

		// set new system state
		long int UpdatePT ( double T_k, double P_bar );

		bool testSizes( long int NSpecies, long int NParams, long int NPcoefs,
				long int MaxOrder, long int NPperDC, char Mod_Code, char Mix_Code );

		// getting phase name
		void GetPhaseName( const char *PhName );

};



// Churakov & Gottschalk (2003) EOS calculations
// declaration of EOSPARAM class (used by the TCGFcalc class)
class EOSPARAM
{
	private:
		//static double Told;
		// unsigned long int isize;  // int isize = NComp;
		long int NComp;
		double emix, s3mix;
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

		void allocate();
		void free();

	public:

		double *XX0;

		//EOSPARAM():isize(0),emix(0),s3mix(0),NComp(0){};
		//EOSPARAM(double*data, unsigned nn):isize(0){allocate(nn);init(data,nn);};

		EOSPARAM( double *Xtmp, double *data, long int nn )
			:NComp(nn), emix(0),s3mix(0)
			{ allocate(); init(Xtmp,data,nn);};

		~EOSPARAM()
			{ free(); }

		void init( double*,double *, long int );
		long int NCmp()   {return NComp;};

		double EPS05( long int i)
			{return eps05[i];};
		double X( long int i)
			{return XX[i];};
		double EPS( long int i)
			{return eps[i];};
		double EMIX(void)
			{return emix;};
		double S3MIX(void)
			{return s3mix;};

		double MIXS3( long int i, long int j)
		{
			if (i==j) return sig3par[i];
			if (i<j) return mixpar[i][j]; else return mixpar[j][i];
		};

		double MIXES3( long int i, long int j)
		{
			if ( i==j ) return epspar[i];
			if (i<j) return mixpar[j][i]; else return mixpar[i][j];
		};

		double SIG3( long int i){return sig3par[i];};
		double M2R( long int i) {return m2par[i];};
		double A( long int i)   {return apar[i];};

		long int ParamMix( double *Xin);
};



// -------------------------------------------------------------------------------------
// Churakov and Gottschalk (2003) EOS calculations
// Added 09 May 2003
// Declaration of a class for CG EOS calculations for fluids
// Incorporates a C++ program written by Sergey Churakov (CSCS ETHZ)
// implementing papers by Churakov and Gottschalk (2003a, 2003b)

class TCGFcalc: public TSolMod
{
	private:

		double
        PI_1,    // pi
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

		//  double PhVol;  // phase volume in cm3
		double *Pparc;     // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
		double *phWGT;
		double *aX;        // DC quantities at eqstate x_j, moles - primal IPM solution [L]
		double *aGEX;      // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
		double *aVol;      // DC molar volumes, cm3/mol [L]

		// main work arrays
		EOSPARAM *paar;
		EOSPARAM *paar1;
		double *FugCoefs;
		double *EoSparam;
		double *EoSparam1;

		// internal functions
		void alloc_internal();
		void free_internal();
		void set_internal();

		void choose( double *pres, double P,unsigned long int &x1,unsigned long int &x2 );
		double Melt2( double T );
		double Melt( double T );

		void copy( double* sours,double *dest,unsigned long int num );
		void norm( double *X,unsigned long int mNum );
		double RPA( double beta,double nuw );
		double dHS( double beta,double ro );

		inline double fI1_6( double nuw )
		{
			return (1.+(A4+(A5+A6*nuw)*nuw)*nuw)/
				((1.+(AA1+(AA2+AA3*nuw)*nuw)*nuw)*3.);
		};

		inline double fI1_12( double nuw )
		{
			return (1.+(B4+(B5+B6*nuw)*nuw)*nuw)/
				((1.+(BB1+(BB2+BB3*nuw)*nuw)*nuw)*9.);
		};

		inline double fa0( double nuw ,double nu1w2 )
		{
			return (A00 + A01*nuw)/nu1w2;
		};

		inline double fa1( double nuw ,double nu1w3 )
		{
			return (A10+(A11+A12*nuw)*nuw)/nu1w3;
		};

		inline double fa2( double nuw ,double nu1w4 )
		{
			return ((A21+(A22+A23*nuw)*nuw)*nuw)/nu1w4;
		};

		inline double fa3( double nuw ,double nu1w5 )
		{
			return ((A31+(A32+(A33+A34*nuw)*nuw)*nuw)*nuw)/nu1w5;
		};

		double DIntegral( double T, double ro, unsigned long int IType ); // not used
		double LIntegral( double T, double ro, unsigned long int IType ); // not used
		double KIntegral( double T, double ro, unsigned long int IType ); // not used
		double K23_13( double T, double ro );
		double J6LJ( double T,double ro );
		double FDipPair( double T,double ro,double m2 ); // not used
		double UWCANum( double T,double ro );
		double ZWCANum( double T,double ro );

		double FWCA( double T,double ro );
		double FTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double UTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param ); // not used
		double ZTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double PTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double ROTOTALMIX( double P,double TT,EOSPARAM* param );

		double PRESSURE( double *X, double *param, unsigned long int NN, double ro, double T ); // not used
		double DENSITY( double *X,double *param, unsigned long int NN ,double Pbar, double T );
		long int CGActivCoefRhoT( double *X,double *param, double *act, unsigned long int NN,
				double ro, double T ); // not used

		long int CGActivCoefPT(double *X,double *param,double *act, unsigned long int NN,
				double Pbar, double T, double &roro );

	public:

		// Constructor
		TCGFcalc( long int NCmp, double Pp, double Tkp );
		TCGFcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *aPparc, double *aphWGT, double *arX,
				double *arGEX, double *arVol, double T_k, double P_bar );

		// Destructor
		~TCGFcalc();

		// calculates of pure species properties (pure fugacities)
		long int PureSpecies( );

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		// CGofPureGases, calculates fugacity for 1 species at (X=1)
		long int CGcalcFugPure( double Tmin, float *Cemp, double *FugProps );  // called from DCthermo
		long int CGFugacityPT( double *EoSparam, double *EoSparPT, double &Fugacity,
				double &Volume, double P, double T, double &roro );

		// calculates departure functions
		long int CGResidualFunct( double *X, double *param, double *param1, unsigned long int NN,
				double ro, double T );

		double GetDELTA( void )
		{
			return DELTA;
		};
};



// -------------------------------------------------------------------------------------
// Peng-Robinson-Stryjek-Vera (PRSV) model for fluid mixtures
// References: Stryjek and Vera (1986)
// (c) TW July 2006

class TPRSVcalc: public TSolMod

{
	private:

		double PhVol;   // phase volume in cm3
		double *Pparc;  // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
		double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
		double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
		double (*Eosparm)[6];   // EoS parameters
		double (*Pureparm)[4];  // Parameters a, b, da/dT, d2a/dT2 for cubic EoS
		double (*Fugpure)[6];   // fugacity parameters of pure gas species
		double (*Fugci)[4];     // fugacity parameters of species in the mixture

		double **a;		// arrays of generic parameters
		double **b;
		double **KK;     // binary interaction parameter
		double **dKK;    // derivative of interaction parameter
		double **d2KK;   // second derivative
		double **AA;     // binary a terms in the mixture



		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double k1, double k2, double k3,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

		// Constructor
		TPRSVcalc( long int NCmp, double Pp, double Tkp );
		TPRSVcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *aPparc, double *arGEX, double *arVol,
				double T_k, double P_bar );

		// Destructor
		~TPRSVcalc();

		// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

		// Calculates T,P corrected interaction parameters
		long int PTparam();

		// Calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		// Calculates pure species properties (called from DCthermo)
		long int PRSVCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
// Soave-Redlich-Kwong (SRK) model for fluid mixtures
// References: Soave (1972); Soave (1993)
// (c) TW December 2008

class TSRKcalc: public TSolMod

{
	private:

		double PhVol;   // phase volume in cm3
		double *Pparc;  // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
		double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
		double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
		double (*Eosparm)[4];   // EoS parameters
		double (*Pureparm)[4];  // Parameters a, b, da/dT, d2a/dT2 for cubic EoS
		double (*Fugpure)[6];   // Fugacity parameters of pure gas species
		double (*Fugci)[4];     // Fugacity parameters of species in the mixture

		double **a;		// arrays of generic parameters
		double **b;
		double **KK;    // binary interaction parameter
		double **dKK;   // derivative of interaction parameter
		double **d2KK;  // second derivative
		double **AA;    // binary a terms in the mixture

		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double N,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

		// Constructor
		TSRKcalc( long int NCmp, double Pp, double Tkp );
		TSRKcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *aPparc, double *arGEX, double *arVol,
				double T_k, double P_bar );

		// Destructor
		~TSRKcalc();

		// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

		// Calculates T,P corrected interaction parameters
		long int PTparam();

		// Calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		// Calculates pure species properties (called from DCthermo)
		long int SRKCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
// Peng-Robinson (PR78) model for fluid mixtures
// References: Peng and Robinson (1976); Peng and Robinson (1978)
// (c) TW July 2009

class TPR78calc: public TSolMod

{
	private:

		double PhVol;   // phase volume in cm3
		double *Pparc;  // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
		double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
		double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
		double (*Eosparm)[4];   // EoS parameters
		double (*Pureparm)[4];  // Parameters a, b, da/dT, d2a/dT2 for cubic EoS
		double (*Fugpure)[6];   // Fugacity parameters of pure gas species
		double (*Fugci)[4];     // Fugacity parameters of species in the mixture

		double **a;		// arrays of generic parameters
		double **b;
		double **KK;    // binary interaction parameter
		double **dKK;   // derivative of interaction parameter
		double **d2KK;  // second derivative
		double **AA;    // binary a terms in the mixture

		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double N,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

		// Constructor
		TPR78calc( long int NCmp, double Pp, double Tkp );
		TPR78calc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *aPparc, double *arGEX, double *arVol,
				double T_k, double P_bar );

		// Destructor
		~TPR78calc();

		// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

		// Calculates T,P corrected interaction parameters
		long int PTparam();

		// Calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		// Calculates pure species properties (called from DCthermo)
		long int PR78CalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
// Compensated Redlich-Kwong (CORK) model for fluid mixtures
// References: Holland and Powell (1991)
// (c) TW May 2010

class TCORKcalc: public TSolMod

{
        private:

                double PhVol;   // phase volume in cm3
                double *Pparc;  // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
                double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
                double *aVol;   // DC molar volumes, cm3/mol [L]

                // main work arrays
                double (*Eosparm)[2];   // EoS parameters
                double (*Fugpure)[6];   // Fugacity parameters of pure gas species
                double (*Fugci)[4];     // Fugacity parameters of species in the mixture

                char *EosCode;    // identifier of EoS routine
                double *Phi;        // phi parameters
                double **A;         // binary interaction parameters
                double **W;         // volume scaled interaction parameters (derivatives)
                double **dW;
                double **d2W;
                double **dWp;

                double RR;    // gas constant in kbar
                double Pkb;   // pressure in kbar

                // internal functions
                void alloc_internal();
                void free_internal();
                long int FugacityPT( long int j, double *EoSparam );
                long int FugacityH2O( long int j );
                long int FugacityCO2( long int j );
                long int FugacityCorresponding( long int j );
                long int VolumeFugacity( long int phState, double pp, double p0, double a, double b, double c,
                        double d, double e, double &vol, double &fc );
                long int Cardano( double cb, double cc, double cd, double &v1, double &v2, double &v3 );
                long int FugacityMix();
                long int ResidualFunct();


        public:

                // Constructor
                TCORKcalc( long int NCmp, double Pp, double Tkp, char *Eos_Code );
                TCORKcalc( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
                                long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
                                double *arIPc, double *arDCc, double *arWx, double *arlnGam,
                                double *aphVOL, double *aPparc, double *arGEX, double *arVol,
                                double T_k, double P_bar, char *Eos_Code );

                // Destructor
                ~TCORKcalc();

                // Calculates pure species properties (pure fugacities)
                long int PureSpecies();

                // Calculates T,P corrected interaction parameters
                long int PTparam();

                // Calculates activity coefficients
                long int MixMod();

                // calculates excess properties
                long int ExcessProp( double *Zex );

                // calculates ideal mixing properties
                long int IdealProp( double *Zid );

                // Calculates pure species properties (called from DCthermo)
                long int CORKCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
// Van Laar model for solid solutions
// References:  Holland and Powell (2003)
// (c) TW March 2007

class TVanLaar: public TSolMod
{
	private:
		double *Wu;
		double *Ws;
		double *Wv;
		double *Wpt;   // Interaction coeffs at P-T
		double *Phi;   // Mixing terms
		double *PsVol; // End member volume parameters

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TVanLaar( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TVanLaar();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates of activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Regular model for multicomponent solid solutions
// References: Holland and Powell (1993)
// (c) TW March 2007

class TRegular: public TSolMod
{
	private:
		double *Wu;
		double *Ws;
		double *Wv;
		double *Wpt;   // Interaction coeffs at P-T

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TRegular( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TRegular();

		// calculates T,P corrected interaction parameters
		long int PTparam( );

		// calculates of activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Redlich-Kister model for multicomponent solid solutions
// References: Hillert (1998)
// (c) TW March 2007

class TRedlichKister: public TSolMod
{
	private:
		double (*Lu)[4];
		double (*Ls)[4];
		double (*Lcp)[4];
		double (*Lv)[4];
		double (*Lpt)[4];

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TRedlichKister( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TRedlichKister();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Non-random two liquid (NRTL) model for liquid solutions
// References: Renon and Prausnitz (1968), Prausnitz et al. (1997)
// (c) TW June 2008

class TNRTL: public TSolMod
{
	private:
		double **Tau;
		double **dTau;
		double **d2Tau;
		double **Alp;
		double **dAlp;
		double **d2Alp;
		double **G;
		double **dG;
		double **d2G;

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TNRTL( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TNRTL();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Wilson model for liquid solutions
// References: Prausnitz et al. (1997)
// (c) TW June 2008

class TWilson: public TSolMod
{
	private:
		double **Lam;
		double **dLam;
		double **d2Lam;

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TWilson( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TWilson();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// SIT model reimplementation for aqueous electrolyte solutions
// (c) DK/TW June 2009

class TSIT: public TSolMod
{
	private:

		// data objects copied from MULTI
		double *z;    // vector of species charges (for aqueous models)
		double *m;    // vector of species molalities (for aqueous models)
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties

		// internal work objects
		double I;	// ionic strength
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double **E0;  // interaction parameter
		double **E1;
		double **dE0;
		double **dE1;
		double **d2E0;
		double **d2E1;

		// internal functions
		double IonicStrength();
		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TSIT( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TSIT();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		// Calculation of internal tables (at each GEM iteration)
		long int PTparam();

};

// -------------------------------------------------------------------------------------
// Pitzer model, Harvie-Moller-Weare (HMW) version, with explicit temperature dependence
// References:
// (c) SD/FH February 2009

class TPitzer: public TSolMod
{

private:
	long int Nc;	 // Number of cations
	long int Na;     // Number of anions
	long int Nn;     // Number of neutral species
	long int Ns;     // Total number of aqueous species (without H2O); index of H2O in aq phase
					 // Conversion of species indexes between aq phase and Pitzer parameter tables
	long int *xcx;   // list of indexes of Nc cations in aqueous phase
	long int *xax;   // list of indexes of Na anions in aq phase
	long int *xnx;   // list of indexes of Nn neutral species in aq phase
	double *aZ;    // Vector of species charges (for aqueous models)
	double *zc;
	double *za;
	double *aM;    // Vector of species molality (for aqueous models)
	double *mc;
	double *ma;
	double *mn;
	double *RhoW;  // water density properties
	double *EpsW;  // water dielectrical properties

    double Aphi, dAphidT, d2AphidT2, dAphidP;  // Computing A-Factor
	double I;  // Ionic Strength
	double Is;  // Ionic Strength square root
	double Ffac; // F-Factor
	double Zfac; // Z-Term

	// Input parameter arrays
			//for Gex and activity coefficient calculation
	double **Bet0;     // Beta0 table for cation-anion interactions [Nc][Na]
	double **Bet1;	   // Beta1 table for cation-anion interactions [Nc][Na]
	double **Bet2;	   // Beta2 table for cation-anion interactions [Nc][Na]
	double **Cphi;     // Cphi  table for cation-anion interactions [Nc][Na]
	double **Lam;      // Lam table for neutral-cation interactions [Nn][Nc]
	double **Lam1;     // Lam1 table for neutral-anion interactions [Nn][Na]
	double **Theta;    // Theta table for cation-cation interactions [Nc][Nc]
	double **Theta1;   // Theta1 table for anion-anion interactions [Na][Na]
	double ***Psi;     // Psi array for cation-cation-anion interactions [Nc][Nc][Na]
	double ***Psi1;    // Psi1 array for anion-anion-cation interactions [Na][Na][Nc]
	double ***Zeta;    // Zeta array for neutral-cation-anion interactions [Nn][Nc][Na]


        // Work parameter arrays
		// double *B1;      // B' table for cation-anion interactions corrected for IS [Nc][Na]
		// double *B2;      // B table for cation-anion interactions corrected for IS [Nc][Na]
		// double *B3;      // B_phi table for cation-anion interactions corrected for IS [Nc][Na]
		// double *Phi1;    // Phi' table for anion-anion interactions corrected for IS [Na][Na]
		// double *Phi2;    // Phi table for cation-cation interactions corrected for IS [Nc][Nc]
		// double *Phi3;    // PhiPhi table for anion-anion interactions corrected for IS [Na][Na]
		// double *C;       // C table for cation-anion interactions corrected for charge [Nc][Na]
		// double *Etheta;  // Etheta table for cation-cation interactions [Nc][Nc]
		// double *Ethetap; // Etheta' table for anion-anion interactions [Na][Na]
		// double bk[21];   // work space
		// double dk[21];   // work space

	// McInnes parameter array and gamma values
	double *McI_PT_array;
	double *GammaMcI;

	enum eTableType
	{
		bet0_ = -10, bet1_ = -11, bet2_ = -12, Cphi_ = -20, Lam_ = -30, Lam1_ = -31,
		Theta_ = -40,  Theta1_ = -41, Psi_ = -50, Psi1_ = -51, Zeta_ = -60
	};

	// internal setup
	void calcSizes();
	void alloc_internal();
	void free_internal();

	// build conversion of species indexes between aq phase and Pitzer parameter tables
	void setIndexes();
	double setvalue(long int ii, int Gex_or_Sex);

	// internal calculations
	// Calculation of Etheta and Ethetap values
	void Ecalc( double z, double z1, double I, double DH_term,
					double& Etheta, double& Ethetap );
	inline long int getN() const
	{
		return Nc+Na+Nn;
	}

	double Z_Term( );
	double IonicStr( double& I );
	void getAlp( long int c, long int a, double& alp, double& alp1 );
	double get_g( double x_alp );
	double get_gp( double x_alp );
	double G_ex_par5( long int ii );
	double G_ex_par8( long int ii );
	double S_ex_par5( long int ii );
	double S_ex_par8( long int ii );
	double CP_ex_par5( long int ii );
	double CP_ex_par8( long int ii );
	double F_Factor( double DH_term );
	double lnGammaN( long int N );
	double lnGammaM( long int M, double DH_term );
	double lnGammaX( long int X, double DH_term );
	double lnGammaH2O( double DH_term );

	// calc vector of interaction parameters corrected to T,P of interest
	void PTcalc( int Gex_or_Sex );

	// calculation KCl activity coefficients for McInnes scaling
	double McInnes_KCl();

	inline long int getIc( long int jj )
    {
		for( long int ic=0; ic<Nc; ic++ )
			if( xcx[ic] == jj )
				return ic;
		return -1;
    }

	inline long int getIa( long int jj )
    {
		for( long int ia=0; ia<Na; ia++ )
			if( xax[ia] == jj )
				return ia;
		return -1;
    }

	inline long int getIn( long int jj )
    {
		for( long int in=0; in<Nn; in++ )
			if( xnx[in] == jj )
				return in;
		return -1;
    }

	inline double sum( double* arr, long int *xx, long int Narr )
    {
		double sum_ =0.;
		for( long int i=0; i<Narr; i++ )
          sum_ += arr[xx[i]];
		return sum_;
    }

public:

    // Constructor
	TPitzer( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
	         long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
	         double *arIPc, double *arDCc, double *arWx, double *arlnGam,
	         double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
	         double *dW, double *eW );

	// Destructor
	~TPitzer();

	// Calculation of T,P corrected interaction parameters
	long int PTparam();


    long int MixMod();

	// calculates activity coefficients
	long int Pitzer_calc_Gamma();
	long int Pitzer_McInnes_KCl();

    // calculates excess properties
    long int ExcessProp( double *Zex );

	// calculates ideal mixing properties
	long int IdealProp( double *Zid );

	void Pitzer_test_out( const char *path, double Y );

};



// -------------------------------------------------------------------------------------
// Extended universal quasi-chemical (EUNIQUAC) model for aqueous electrolyte solutions
// References: Nicolaisen et al. (1993), Thomsen et al. (1996), Thomsen (2005)
// (c) TW/FH May 2009

class TEUNIQUAC: public TSolMod
{
	private:

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties

		// internal work objects
		double *R;   // volume parameter
		double *Q;   // surface parameter
		double *Phi;
		double *Theta;
		double **U;   // interaction energies
		double **dU;   // first derivative
		double **d2U;   // second derivative
		double **Psi;
		double **dPsi;
		double **d2Psi;
		double IS;  // ionic strength
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)

		// objects needed for debugging output
		double gammaDH[200];
		double gammaC[200];
		double gammaR[200];

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

		// Constructor
		TEUNIQUAC( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TEUNIQUAC();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

		void Euniquac_test_out( const char *path );

};



// -------------------------------------------------------------------------------------
// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Helgesons variant
// References: Helgeson et al. (1981); Oelkers and Helgeson (1990); Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW July 2009

class THelgeson: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // flag for water
		long int flagNeut;  // flag for neutral species
		long int flagElect;  // flag for selection of background electrolyte model

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties
		double *an;  // individual ion size-parameters
		double *bg;  // individual extended-term parameters
		double ac;  // common ion size parameters
		double bc;  // common extended-term parameter

		// internal work objects
		double ao, daodT, d2aodT2, daodP;  // ion-size parameter (TP corrected)
		double bgam, dbgdT, d2bgdT2, dbgdP;  // extended-term parameter (TP corrected)
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double molZ;  // total molality of charged species
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)
		double B, dBdT, d2BdT2, dBdP;  // B term of DH equation (and derivatives)
		double Gf, dGfdT, d2GfdT2, dGfdP;  // g function (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

		// Constructor
		THelgeson( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~THelgeson();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Davies variant
// References: Langmuir (1997)
// (c) TW July 2009

class TDavies: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // flag for water
		long int flagNeut;  // flag for neutral species
		long int flagMol;  // flag for molality correction

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties

		// internal work objects
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

		// Constructor
		TDavies( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TDavies();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Debye-Hueckel (DH) limiting law for aqueous electrolyte solutions
// References: Langmuir (1997)
// (c) TW July 2009

class TLimitingLaw: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // flag for water
		long int flagNeut;  // flag for neutral species

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties

		// internal work objects
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

		// Constructor
		TLimitingLaw( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TLimitingLaw();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Two-term Debye-Hueckel (DH) model for aqueous electrolyte solutions
// References: Helgeson et al. (1981)
// uses individual ion-size parameters, optionally individual salting-out coefficients
// (c) TW July 2009

class TDebyeHueckel: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // flag for water
		long int flagNeut;  // flag for neutral species

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties
		double *an;  // individual ion size-parameters
		double *bg;  // individual extended-term parameters
		double ac;  // common ion size parameters
		double bc;  // common extended-term parameter

		// internal work objects
		double ao;  // average ion-size parameter
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)
		double B, dBdT, d2BdT2, dBdP;  // B term of DH equation (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

		// Constructor
		TDebyeHueckel( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TDebyeHueckel();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Karpovs variant
// References: Karpov et al. (1997); Helgeson et al. (1981); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW July 2009

class TKarpov: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // flag for water
		long int flagNeut;  // flag for neutral species
		long int flagElect;  // flag for selection of background electrolyte model

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties
		double *an;  // individual ion size-parameters at T,P
		double *bg;  // individual extended-term parameters
		double ac;  // common ion size parameters
		double bc;  // common extended-term parameter

		// internal work objects
		double ao;  // average ion-size parameter
		double bgam, dbgdT, d2bgdT2, dbgdP;  // extended-term parameter (TP corrected)
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double molZ;  // total molality of charged species
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)
		double B, dBdT, d2BdT2, dBdP;  // B term of DH equation (and derivatives)
		double Gf, dGfdT, d2GfdT2, dGfdP;  // g function (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

		// Constructor
		TKarpov( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TKarpov();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Shvarov variant
// References: Shvarov (2007); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW July 2009

class TShvarov: public TSolMod
{
	private:

		// status flags copied from MULTI
		long int flagH2O;  // new flag for water
		long int flagNeut;  // new flag for neutral species
		long int flagElect;  // flag for selection of background electrolyte model

		// data objects copied from MULTI
		double *z;   // species charges
		double *m;   // species molalities
		double *RhoW;  // water density properties
		double *EpsW;  // water dielectrical properties
		double *bj;  // individual ion parameters
		double ac;  // common ion size parameters
		double bc;  // common extended-term parameter

		// internal work objects
		double ao, daodT, d2aodT2, daodP;  // ion-size parameter (TP corrected)
		double bgam, dbgdT, d2bgdT2, dbgdP;  // extended-term parameter (TP corrected)
		double *LnG;  // activity coefficient
		double *dLnGdT;  // derivatives
		double *d2LnGdT2;
		double *dLnGdP;
		double IS;  // ionic strength
		double molT;  // total molality of aqueous species (except water solvent)
		double A, dAdT, d2AdT2, dAdP;  // A term of DH equation (and derivatives)
		double B, dBdT, d2BdT2, dBdP;  // B term of DH equation (and derivatives)
		double Gf, dGfdT, d2GfdT2, dGfdP;  // g function (and derivatives)

		// internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

		// Constructor
		TShvarov( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double *arM, double *arZ, double T_k, double P_bar,
				double *dW, double *eW );

		// Destructor
		~TShvarov();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Class for hardcoded models for solid solutions
// (c) TW January 2009

class TModOther: public TSolMod
{
	private:
		double PhVol;   // phase volume in cm3
		double *Pparc;  // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
		double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
		double *aVol;   // DC molar volumes, cm3/mol [L]
		double *Gdqf;	// DQF correction terms
		double *Hdqf;
		double *Sdqf;
		double *CPdqf;
		double *Vdqf;

		void alloc_internal();
		void free_internal();

	public:

		// Constructor
		TModOther( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar, double *dW, double *eW );

		// Destructor
		~TModOther();

		// calculates pure species properties (pure fugacities, DQF corrections)
		long int PureSpecies();

		// calculates T,P corrected interaction parameters
		long int PTparam();

		// calculates activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

                // functions for individual models (under construction)
                long int Amphibole1();
                long int Biotite1();
                long int Chlorite1();
                long int Clinopyroxene1();
                long int Feldspar1();
                long int Feldspar2();
                long int Garnet1();
                long int Muscovite1();
                long int Orthopyroxene1();
                long int Staurolite1();
                long int Talc();

};



// -------------------------------------------------------------------------------------
// Ternary Margules (regular) model for solid solutions
// References: Anderson and Crerar (1993); Anderson (2006)
// (c) TW/DK June 2009

class TMargules: public TSolMod
{
	private:

		double WU12, WS12, WV12, WG12;
		double WU13, WS13, WV13, WG13;
		double WU23, WS23, WV23, WG23;
		double WU123, WS123, WV123, WG123;

	public:

		// Constructor
		TMargules( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TMargules();

		// calculates T,P corrected interaction parameters
		long int PTparam( );

		// calculates of activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Binary Margules (subregular) model for solid solutions
// References: Anderson and Crerar (1993); Anderson (2006)
// (c) TW/DK June 2009

class TSubregular: public TSolMod
{
	private:

		double WU12, WS12, WV12, WG12;
		double WU21, WS21, WV21, WG21;

	public:

		// Constructor
		TSubregular( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TSubregular();

		// calculates T,P corrected interaction parameters
		long int PTparam( );

		// calculates of activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
// Binary Guggenheim (Redlich-Kister) model for solid solutions
// References: Anderson and Crerar (1993); Anderson (2006)
// uses normalized (by RT) interaction parameters
// (c) TW/DK June 2009

class TGuggenheim: public TSolMod
{
	private:

		double a0, a1, a2;

	public:

		// Constructor
		TGuggenheim( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
				long int NPperDC, char Mod_Code, char Mix_Code, long int *arIPx,
				double *arIPc, double *arDCc, double *arWx, double *arlnGam,
				double *aphVOL, double T_k, double P_bar );

		// Destructor
		~TGuggenheim();

		// calculates T,P corrected interaction parameters
		long int PTparam( );

		// calculates of activity coefficients
		long int MixMod();

		// calculates excess properties
		long int ExcessProp( double *Zex );

		// calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



#endif

// _s_fgl_h
