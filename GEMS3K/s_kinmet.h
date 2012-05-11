//-------------------------------------------------------------------
// $Id: s_kinmet.h 1445 2012-05-09 12:42:54Z gems $
//
// Copyright (C) 2012  D.Kulik, B.Thien, G.Kosakowski
//
// Declaration of kinetics and metastability models (TKinMet class)
//
// This file is part of a GEM-Selektor (GEMS) v.3.x program
// environment for thermodynamic modeling in geochemistry
// and part of the standalone GEMS3K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//------------------------------------------------------------------
//
#ifndef S_KINMET_H
#define S_KINMET_H

//#include "s_fgl.h"

const int   MAXDCNAME =      16, MAXPHNAME = 16;

// Class describing a set of rate constants for a mineral precipitation/dissolution (Lasaga's form)
class TRateCon
{
protected:

char  SitT;   // Site type code, see SITETYPECODES
char  SATT;   // SACT equation code, see SACTCODES

long int NSpec;   // Number of species on this site (at least 1; 0 if site to be ignored0

// double MASDJ[DFCN]; Parameters of surface species in surface complexation models

double  qC,          //  site capacity parameter in mol/kg (CEC in eq/kg)
        GamC,        //  site density parameter in mol/m2 (CEC in eq/m2)
        AlphaF,      //  Frumkin interaction parameter;
        BETp,        //  BET isotherm parameter
        BETq;        //  BET isotherm parameter

 double dent[];      //  Species denticity or coordination number [NSpec]

// enum {
//	//[0] - max site density in mkmol/(g sorbent); [1] - species charge allocated to 0 plane;
//	//[2] - surface species charge allocated to beta -or third plane; [3] - Frumkin interaction parameter;
//	//[4] species denticity or coordination number; [5]  - reserved parameter (e.g. species charge on 3rd EIL plane)]
//   XL_ST = 0, XL_EM, XL_SI, XL_SP
// };

   double nxs[];       // moles of surface species (picked up) [NSpec]
   double spDUL[];    // temporary upper constraint on species amount for SAT calculations [NSpec]

// results
//   double (*D)[MST];  // Reserved; new work array for calc. surface act.coeff.
   double XSsT;         // current total moles of species on this site type
// double lnSAC[4];   // former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]
   double ISAT[];     // ISAT for each species (in SACT calculations) [NSpec]
   double eF[];       // Frumkin exponent acting on species [NSpec]
   double SACT[];     // Surface activity coefficient terms (confid.entropy terms) [NSpec]
   double lnSACT[];   // ln of SACT

public:

    // Generic constructor
    TRateCon(  );

    // Destructor
    virtual ~TRateCon();

};



//  Data for kinetics and metastability in MULTI [Kulik,2006])


class TKinMet  // Base class for kinetics and metastability models
{
    protected:
        char ModCode;   // Code of the sorption phase- see SORPPHASECODES
        char PhaseName[MAXPHASENAME+1];    // Phase name (for specific built-in models)

        long int NComp;       // Number of components in the sorption phase
        long int NSC;         // Number of components in the sorption part of the phase; if 0 then no sorption for this phase
        long int NemS;        // Number of components in the sorbent or carrier ( < NComp, can be 0 (then site balances must be used)
        long int NsurT;       // Number of surface patch types (0 to 6), 0 means a bulk sorption capacity model

        double Tk;    	  // Temperature, K
        double Pbar;  	  // Pressure, bar
        double R_CONST;   // R constant
        double F_CONST;   // F (Faraday's) constant
        double N_C_st;    // Rsp1 Standard surface number density, 1/m2"
        double Gam_C_st;  // Standard surface density, mol/m2
        double q_C_st;    // Standard sorption capacity, mol/kg(sorbent)

        double Xcond, 	// conductivity of phase carrier, sm/m2, reserved
               Xeps,  	// rel. diel. permeability of phase carrier (solvent), reserved
               Aalp,  	// Specific surface area of phase carrier (sorbent) (m2/g) !!! m2/kg
               Sigw,  	// Specific surface free energy for phase-water interface (J/m2)
               Sigg,  	// Specific surface free energy for phase-gas interface (J/m2) (not yet used)
               Xr0,   // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved) Xr0h0[2]
               Xh0;   // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)

        char  (*SM3)[MAXDCNAME];  // pointer to the list of DC names in the sorption phase [NComp] read-only
        char  DCC[];   // pointer to the classifier of DCs involved in sorption phase [NComp] read-only
        char  (*SCM)[MST]; // Classifier of built-in EIL models applied to surface types on sorption phase [MST=6] read-only
        char  *SATT[];  // Classifier of applied SACT equations (isotherm corrections) [NComp] read-only

        TSurfPatchMod* patchMod[]; // Pointer to array of TSurfPatch instances - size NsurT

        double *nx;     // Pointer to mole amounts of sorption phase components (provided) [NComp] read-only

        // Results
        double Gex, Hex, Sex, CPex, Vex, Aex, Uex;   // molar excess properties for surface species
        double Gid, Hid, Sid, CPid, Vid, Aid, Uid;   // molar ideal mixing properties for surface species

        double *lnSACT;    // Pointer to ln SACT for surface species [NSC]
        double *CoulTerms; // Pointer to Coulombic correction terms (electrostatic activity coefficients) [NSC]
        double *lnGamma;   // Pointer to ln activity coefficients of sorption phase components
                           // (memory under pointers must be provided from the calling program)

        public:
                // Generic constructor
                TSorpMod( long int NSpecies, long int NSurSpecies, long int NSorbentEMs, long int NSurfTypes,
                          char Mod_Code, char **SM3_, char *DCC_, char **SCM_, char **SATT_,
                          double Aalp_, double *arnx_, double *arlnGam, double *arlnSACT, double *arCoulT,
                          double T_k, double P_bar );

                // Destructor
                virtual ~TSorpMod();

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

                virtual long int EILMod()
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

                bool testSizes( long int NSpecies, long int NSurSpecies, long int NSorbentEMs,
                                long int NSurfTypes, char Mod_Code, char EIL_Code );

                // getting phase name
                void GetPhaseName( const char *PhName );

};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - derived classes - - - - - - - - - - - - - -

class TPalandri: public TKinMet  // Generic dissolution/precipitation models following Lasaga 1998
{
    protected:

        double Xcond, 	// conductivity of phase carrier, sm/m2, reserved
               Xeps,  	// rel. diel. permeability of phase carrier (solvent), reserved
               Aalp,  	// Specific surface area of phase carrier (sorbent) (m2/g) !!! m2/kg
               Sigw,  	// Specific surface free energy for phase-water interface (J/m2)
               Sigg,  	// Specific surface free energy for phase-gas interface (J/m2) (not yet used)
               Xr0,   // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved) Xr0h0[2]
               Xh0;   // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)

//        TSurfPatchMod* patchMod[]; // Pointer to array of TSurfPatch instances - size NsurT

        // internal functions
        void alloc_internal();
        void free_internal();


public:

        // Constructor
        TLasaga( long int NCmp, double Pp, double Tkp );
        TLasaga( long int NSpecies, long int NSurSpecies, long int NSorbentEMs, long int NSurfTypes,
                  char Mod_Code, char *Phase_Name, char **SM3_, char *DCC_, char **SCM_, char **SATT_,
                  double Aalp_, double *arnx_, double *arlnGam, double *arlnSACT, double T_k, double P_bar  );
        // Destructor
        ~TLasaga();

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


#endif // S_KINMET_H
