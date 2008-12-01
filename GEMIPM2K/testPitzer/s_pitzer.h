
class TSolMod
{

protected:
        char ModCode;   // Code of the mixing model

        long int NComp;    	// Number of components in the solution phase
        long int NPar;     	// Number of non-zero interaction parameters
        long int NPcoef;   	// Number of coeffs per parameter (columns in the aIPc table)
        long int MaxOrd;   	// max. parameter order (or number of columns in aIPx)
        long int NP_DC;    	// Number of coeffs per one DC in the phase (columns in aDCc)
        long int *aIPx;  	// Pointer to list of indexes of non-zero interaction parameters

        double R_CONST; // R constant
        double RhoW;	// Density of liquid water, added 04.06.2008 (TW)
        double EpsW;	// Dielectric constant of liquid water
        double IonStr;	// Ionic strength
        double Tk;    	// Temperature, K
        double Pbar;  	// Pressure, bar

double *aIPc;  	// Table of interaction parameter coefficients
double *aIP;    // Vector of interaction parameters corrected to T,P of interest
        double *aDCc;  	// End-member parameter coefficients
        double *x;    	// Pointer to mole fractions of end members (provided)
double *aZ;    // Vector of species charges (for aqueous models)
double *aM;    // Vector of species molality (for aqueous models)

// Results
        double Gam;   	// work cell for activity coefficient of end member
        double lnGamRT;
        double lnGam;
        double Gex;   	// Molar excess Gibbs energy
        double Vex;   	// Excess molar volume
        double Hex;   	// Excess molar enthalpy
        double Sex;   	// Excess molar entropy
        double CPex;  	// Excess heat capacity
        double *lnGamma;   // Pointer to ln activity coefficients of end members
                           // (memory must be provided from the calling program)
        
public:
// Generic constructor
    TSolMod( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
         long int NPperDC, double T_k, double P_bar, char Mod_Code,
         long int* arIPx, double* arIPc, double* arDCc,
         double *arWx, double *arlnGam, double *arM, double *arZ, 
         double dW, double eW, double iS );
    ~TSolMod();


};



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

//	double I;        // Ionic strength
//	double *z;    // Vector of species charges (for aqueous models)
//	double *m;    // Vector of species molality (for aqueous models)
//	double *aIP;  // Vector of interaction parameters corrected to T,P of interest

	 double Aphi; //------------ Computing A- Factor
	 double I;  //----------- Ionic Strength
	 double Is;  //----------- Ionic Strength
	 double Ffac; //----------- ------- F-Factor
	 double	 Zfac; //----------- Z- Term

				// Input parameter arrays
	double *abet0;    // Beta0 table for cation-anion interactions [Nc][Na]
	double *abet1;    // Beta1 table for cation-anion interactions [Nc][Na]
	double *abet2;    // Beta2 table for cation-anion interactions [Nc][Na]
	double *aCphi;    // Cphi  table for cation-anion interactions [Nc][Na]
	double *aLam;     // Lam table for neutral-cation interactions [Nn][Nc]
	double *aLam1;    // Lam1 table for neutral-anion interactions [Nn][Na]
	double *aTheta;   // Theta table for cation-cation interactions [Nc][Nc]
	double *aTheta1;  // Theta1 table for anion-anion interactions [Na][Na]
	double *aPsi;     // Psi array for cation-cation-anion interactions [Nc][Nc][Na]
	double *aPsi1;    // Psi1 array for anion-anion-cation interactions [Na][Na][Nc]
	double *aZeta;    // Zeta array for neutral-cation-anion interactions [Nn][Nc][Na]
	
	// internal values
	
                  // Work parameter arrays
/*	double *B1;      // B' table for cation-anion interactions corrected for IS [Nc][Na]
	double *B2;      // B table for cation-anion interactions corrected for IS [Nc][Na]
	double *B3;      // B_phi table for cation-anion interactions corrected for IS [Nc][Na]
	double *Phi1;    // Phi' table for anion-anion interactions corrected for IS [Na][Na]
	double *Phi2;    // Phi table for cation-cation interactions corrected for IS [Nc][Nc]
	double *Phi3;    // PhiPhi table for anion-anion interactions corrected for IS [Na][Na]
	double *C;       // C table for cation-anion interactions corrected for charge [Nc][Na]
	double *Etheta;  // Etheta table for cation-cation interactions [Nc][Nc]
	double *Ethetap; // Etheta' table for anion-anion interactions [Na][Na]
	double bk[21];   // work space
	double dk[21];   // work space
*/
	enum eTableType
	{ bet0_ = -10, bet1_ = -11, bet2_ = -12, Cphi_ = -20, Lam_ = -30, Lam1_ = -31, 
	  Theta_ = -40,  Theta1_ = -41, Psi_ = -50, Psi1_ = -51, Zeta_ = -60
	};
	
 // internal calculations
	// Calculation of Etheta and Ethetap values
	void Ecalc( double z, double z1, double I, double Aphi, 
			double& Etheta, double& Ethetap);
	inline long int getN() const
	 { return Nc+Na+Nn; }
	double Z_Term();
	double A_Factor( double T );
	double IonicStr( double& I );
	void getAlp( long int c, long int a, double& alp, double& alp1 );
	double F_Factor( double Aphi, double I, double Is );

	double lnGammaN(  long int N );
	double lnGammaM(  long int M );
	double lnGammaX(  long int X );
	double lnGammaH2O();

  // internal setup	
	void calcSizes();
	void alloc_internal();
	void free_internal();
  // calc vector of interaction parameters corrected to T,P of interest
	void P_T( double T ); 
  // build conversion of species indexes between aq phase and Pitzer parameter tables
	void setupIndexes();
	void setupValues();

    
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
          sum_ += arr[ xx[i]];
		return sum_;
    }

public:
    // Constructor

	TPitzer( long int NSpecies, long int NParams, long int NPcoefs, long int MaxOrder,
	         long int NPperDC, double T_k, double P_bar, char Mod_Code,
	         long int* arIPx, double* arIPc, double* arDCc,
	         double *arWx, double *arlnGam, double *arM, double *arZ, 
	         double dW, double eW, double iS );
	
	// Destructor
	~TPitzer();
	
	// Calculation of internal tables (at each GEM iteration)
	long int Pitzer_Init();
	
	// Calculation of activity coefficients
	long int Pitzer_Gamma( long int index );
	long int Pitzer_calc_Gamma( );

	void Pitzer_test_out( const char *path );

};

#define IPc( ii, jj )  ( aIPc[ (ii) * NPcoef + jj ])
#define IPx( ii, jj )  ( aIPx[ (ii) * MaxOrd + jj ])

#define mc( ii ) (aM[ xcx[ii] ])
#define ma( ii ) (aM[ xax[ii] ])
#define mn( ii ) (aM[ xnx[ii] ])
#define zc( ii ) (aZ[ xcx[ii] ])
#define za( ii ) (aZ[ xax[ii] ])

#define bet0( c,a ) ( abet0[ ((c)*Na+a) ])
#define bet1( c,a ) ( abet1[ ((c)*Na+a) ])
#define bet2( c,a ) ( abet2[ ((c)*Na+a) ])
#define Cphi( c,a ) ( aCphi[ ((c)*Na+a) ])
#define Lam( n,c )  ( aLam[ ((n)*Nc+c) ])
#define Lam1( n,a )  ( aLam[ ((n)*Na+a) ])
#define Theta( c,c1 )  ( aTheta[ ((c)*Nc+c1) ])
#define Theta1( a,a1 )  ( aTheta1[ ((a)*Na+a1) ])

#define Psi( c,c1,a )  ( aPsi[ (( (c) * Nc + c1  ) * Na + a) ])
#define Psi1( a,a1,c )  ( aPsi1[ (( (a) * Na + a1  ) * Nc + c) ])
#define Zeta( n,c,a )  ( aZeta[ (( (n) * Nc + c  ) * Na + a) ])

 
//--------------------- End of s_pitzer.h ---------------------------
