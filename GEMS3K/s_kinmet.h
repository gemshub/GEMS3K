//-------------------------------------------------------------------
// $Id$
//
// Declaration of TKinMet class for kinetics/metastability models
//
// Copyright (C) 2012-2013  D.Kulik, B.Thien, G.Kosakowski
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------
//
#ifndef S_KINMET_H
#define S_KINMET_H

//#include "s_fgl.h"

const int   MAXDCNAME_ =      16, MAXPHNAME_ = 16;   // see also v_mod.h

enum kinmet_controls {  /// Re-declared codes to control kinetic rate models (see also m_phase.h)

    KM_UNDEF_ = 'N',     /// not defined, no account
    KM_RATE_SURF_ = 'A', /// surface-scaled rate model (k in mol/m2/s)
    KM_RATE_PV_ = 'V',   /// pore-volume-scaled model (k in mol/m3/s)

    KM_DIR_DISSOL_ = 'D', ///	dissolution
    KM_DIR_GROWTH_ = 'G', /// growth on seed particles
    KM_DIR_NUCPREC_ ='P', ///	nucleation and precipitation
    KM_DIR_BOTH_ = 'B',   /// bidirectional (D or G depending on stability index)

    KM_UPT_ENTRAP_ = 'E',  ///	unified entrapment model
    KM_UPT_ENTRDIF_ = 'M', ///	unified entrapment model with in/out diffusion term

    KM_LNK_SURF_ = 'S',  /// link to (fraction of) solid substrate surface
    KM_LNK_PVOL_ = 'P',  ///	link to (fraction of) solid substrate (pore) volume
    KM_LNK_MASS_ = 'M',  ///	link to (fraction of) solid substrate mass

    KM_SIZED_UNI_ = 'U', /// 	uniform
    KM_SIZED_BIN_ = 'B', /// 	binodal
    KM_SIZED_FUN_ = 'F'  ///    distribution function

};

enum affin_term_op_codes {
    ATOP_CLASSIC_ = 0,       /// classic TST affinity term (see .../Doc/Descriptions/KinetParams.pdf)
    ATOP_CLASSIC_REV_ = 1,   /// classic TST affinity term, reversed
    ATOP_SCHOTT_ = 2,        /// Schott et al. 2012 fig. 1e
    ATOP_HELLMANN_ = 3,      /// Hellmann Tisserand 2006 eq 9
    ATOP_TENG1_ = 4,         /// Teng et al. 2000, eq 13
    ATOP_TENG2_ = 5,         /// Teng et al. 2000, Fig 6
    ATOP_FRITZ_ = 6          /// Fritz et al. 2009, eq 6 nucleation and growth

};

// This data structure should be passed in order to create an instance of TKinMet derived class for a phase
struct KinMetData {

    char  KinRateCod_;  /// Type code of the kinetic/metastability model, see enum kinmet_controls
    char  KinDirCod_;   /// Code of direction of the kinetic process, see enum kinmet_controls
    char  KinUptCod_;   /// Type code of Tr uptake model (solution/sorption phases only), see enum kinmet_controlsS
    char  KinLnkCod_;   /// Type code of metastability links of this phase to other phases,see enum kinmet_controls
    char  KinSizedCod_; /// Type of particle/pore size distribution and specific surface area correction, see enum kinmet_controls
    char  KinRezdCod_;  /// Reserved { N }
    char  PhasNam_[MAXPHNAME_];      /// Phase name (for specific built-in models)

    long int NComp_;   /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    long int nlPh_;    /// Number of linked phases (cf. lPh), default 0
    long int nlPc_;    /// TKinMet, TSorpMod: number of parameters per linked phase, default 0.

    long int nPRk_;  /// number of «parallel reactions» that affect amount constraints for k-th phase (1, 2, 3, ...), 1 by default
    long int nSkr_;  /// number of (aqueous or gaseous or surface) species from other reacting phases involved, 0 by default
    long int nrpC_;  /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12 + 3res.)
    long int naptC_; /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    long int nAscC_; /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    long int numpC_; /// number of uptake model parameter coefficients (per end member)
    long int iRes4_;  // reserved

    double T_k_;     /// Temperature, K (initial)
    double P_bar_;   /// Pressure, bar (initial)
    double kTau_;    /// current time, s (initial)
    double kdT_;     /// current time step (initial)
  //
    double IS_;      /// Effective molal ionic strength of aqueous electrolyte
    double pH_;      /// pH of aqueous solution
    double pe_;      /// pe of aqueous solution
    double Eh_;      /// Eh of aqueous solution, V
  //
    double nPh_;     /// current amount of this phase, mol (read-only)
    double mPh_;     /// current mass of this phase, g (read-only)
    double vPh_;     /// current volume of this phase, cm3 (read-only)
    double sAPh_;    /// current surface of this phase, m2
    double LaPh_;    /// phase stability index (log scale)
    double OmPh_;    /// phase stability index (activity scale) 10^LaPh_
  //
    double sSA_;    /// Specific surface area of the phase, m2/g, default: 0.
    double sgw_;    /// Standard mean surface energy of solid-aqueous interface, J/m2
    double sgg_;    /// Standard mean surface energy of gas-aqueous interface, J/m2
    double rX0_;    /// Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
    double hX0_;    /// Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
    double sVp_;    /// Specific pore volume of phase, m3/g (default: 0)
    double sGP_;    /// surface free energy of the phase, J (YOF*PhM)
    double nPul_;   /// upper restriction to this phase amount, mol (calculated here)
    double nPll_;   /// lower restriction to this phase amount, mol (calculated here)

    double *arlPhc_;   /// TsolMod, TKinMet, TSorpMod: pointer to input array of phase link parameters [nlPh*nlPc]
    double *arfeSAr_;  /// Pointer to input fractions of surface area of the solid related to different parallel reactions [nPRk] read-only
    double *arrpCon_;  /// Pointer to input array of kinetic rate constants for faces and 'parallel reactions' [nPRk*nrpC] read-only
    double *arapCon_;  /// Pointer to array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC] read-only
    double *arAscp_;   /// Pointer to array of parameter coefficients of equation for correction of specific surface area [nAscC] read-only
    double *arUmpCon_; /// Pointer to input array of uptake model coefficients [nComp*numpC] read-only
    // new:new: array of nucleation model parameters (A.Testino?)

    char  (*SM_)[MAXDCNAME_];  /// pointer to the classifier of DCs involved in sorption phase [NComp] read-only
    char  *arDCC_;       /// pointer to the classifier of DCs involved in the phase [NComp] read-only
    char  *arlPhC_;      /// TSolMod, TKinMet: Phase linkage type codes [nlPh] { TBA  }

    long int *arocPRk_; /// pointer to operation codes for kinetic parallel reaction affinity terms [nPRk] read-only
    long int *arxSKr_;  /// pointer to input array of DC indexes used in activity products [nSKr_] read-only
    long int *arxlPh_;  /// pointer to input array of linked phase indexes  [nlPh_] read-only

    double *arym_;    /// Pointer to molalities of all species in MULTI (provided), read-only
    double *arla_;    /// Pointer to lg activities of all species in MULTI (provided), read-only

    double *arnx_;    /// Pointer to mole amounts of phase components (provided) [NComp] read-only

    double *arnxul_;  /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    double *arnxll_;  /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    double *arWx_;    /// Species (end member) mole fractions ->NComp
    double *arVol_;   /// molar volumes of end-members (species) cm3/mol ->NSpecies

};

// Class describing a 'parallel reaction' region kinetic rate law data (Schott ea 2012 general form)
class TKinReact
{
protected:

     // Input data
     long int xPR;   /// index of this parallel reaction
     long int iRes; // reserved

     long int ocPRk; /// operation code for this kinetic parallel reaction affinity term
     long int *xSKr;  /// pointer to input array of DC indexes used in activity products [nSKr] (copy)

     double feSAr;   /// input fraction of surface area of the solid related to this parallel reaction
     double *rpCon;  /// input array of kinetic rate constants for this 'parallel reaction' [nrpC]
     double **apCon; /// input array of parameters per species involved in 'activity product' terms [nSkr*naptC]

     // work data: unpacked rpCon[nrpC]
     double ko,  /// rate constant at standard temperature (mol/m2/s)
            Ap,  /// Arrhenius parameter
            Ea,  /// activation energy at st.temperature J/mol
            bI,
            bpH,
            bpe,
            bEh,
            pPR,
            qPR,
            mPR,
            OmEff,

            Omg; /// Input stability index non-log (d-less)

     // Results of rate term calculation

     double
        arf,  // Arrhenius factor (temperature correction on kappa)
        cat,  // catalytic product term (f(prod(a))
        aft,  // affinity term (f(Omega))

        kPR,   // rate constant (involving all corrections) in mol/m2/s
        rPR,   // rate for this region (output) in mol/s
        rmol,   // rate for the whole face (output) in mol/s
//        velo,   // velocity of face growth (positive) or dissolution (negative) nm/s
        ;

public:
    // Generic constructor
    TKinReact( long int xPR_p, long int ocPRk_p, long int *xSKr_p, double feSAr_p, double Omg_p,
               double *rpCon_p, double **apCon_p );

    // Destructor
    virtual ~TKinReact();

    double PRrateCon( long int xPR_p );  // calculates the rate constant for xPR parallel reaction

};


class TKinMet  // Base class for MWR kinetics and metastability models
{
    protected:
    char  KinRateCode;  /// Type code of the kinetic/metastability model, see enum kinmet_controls
    char  KinDirCode;   /// Code of direction of the kinetic process, see enum kinmet_controls
    char  KinUptCode;   /// Type code of Tr uptake model (solution/sorption phases only), see enum kinmet_controlsS
    char  KinLnkCode;   /// Type code of metastability links of this phase to other phases,see enum kinmet_controls
    char  KinSizedCode; /// Type of particle/pore size distribution and specific surface area correction, see enum kinmet_controls
    char  KinRezdCode;  /// Reserved { N }
    char  PhasName[MAXPHNAME_];      /// Phase name (for specific built-in models)

    long int NComp;   /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    long int nlPh;    /// Number of linked phases (cf. lPh), default 0
    long int nlPc;    /// TKinMet, TSorpMod: number of parameters per linked phase, default 0.

  long int nPRk;      /// number of «parallel reactions» that affect amount constraints for k-th phase (1, 2, 3, ...), 1 by default
  long int nSkr;      /// number of (aqueous or gaseous) species from other reacting phases involved, 0 by default
  long int nrpC;      /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12 + 3res.)
  long int naptC;     /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    long int nAscC;   /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    long int numpC;   /// number of uptake model parameter coefficients (per end member)
    long int iRes4;   // reserved

    double R_CONST;     /// Gas constant, 8.31451 J/K/mol
    double T_k;         /// Temperature, K
    double P_bar;       /// Pressure, bar
    double kTau;        /// current time, s
    double kdT;         /// current time step

    // These values will be corrected after GEM converged at each time step
    double IS;          /// Effective molal ionic strength of aqueous electrolyte
    double pH;          /// pH of aqueous solution
    double pe;          /// pe of aqueous solution
    double Eh;          /// Eh of aqueous solution, V
    double nPh;     /// current amount of this phase, mol
    double mPh;     /// current mass of this phase, g
    double vPh;     /// current volume of this phase, cm3
    double sAPh;    /// current surface of this phase, m2
    double LaPh;    /// phase stability index (log scale)
    double OmPh;    /// phase stability index (activity scale) 10^LaPh

    // These values may be corrected inside of TKinMet class instance over time steps
    double sSA;    /// Specific surface area of the phase, m2/g, default: 0.
    double sgw;    /// Standard mean surface energy of solid-aqueous interface, J/m2
    double sgg;    /// Standard mean surface energy of gas-aqueous interface, J/m2
    double rX0;    /// Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
    double hX0;    /// Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
    double sVp;    /// Specific pore volume of phase, m3/g (default: 0)
    double sGP;    /// surface free energy of the phase, J (YOF*PhM)
    double nPul;   /// upper restriction to this phase amount, mol (calculated here)
    double nPll;   /// lower restriction to this phase amount, mol (calculated here)

    double **arlPhc; /// TsolMod, TKinMet, TSorpMod: pointer to input array of phase link parameters [nlPh*nlPc]
  double *arfeSAr;   /// input fractions of surface area of the solid related to different parallel reactions [nPRk]
  double **arrpCon;  /// input array of kinetic rate constants for faces and 'parallel reactions' [nPRk*nrpC]
  double ***arapCon; /// input array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC]
    double *arAscp;  /// input array of parameter coefficients of equation for correction of specific surface area [nAscC]
    double **arUmpCon; /// input array of uptake model coefficients [NComp*numpC] read-only
    // new:new: array of nucleation model parameters (A.Testino?)

    char  (*SM)[MAXDCNAME_];  /// pointer to the list of DC names in the phase [NComp] read-only
    char  *arDCC;       /// pointer to the classifier of DCs involved in the phase [NComp] read-only
    char  *arlPhC;      /// TSolMod, TKinMet: Phase linkage type codes [nlPh] { TBA  }

  long int *arocPRk; /// input operation codes for kinetic parallel reaction affinity terms [nPRk]
  long int *arxSKr;  /// pointer to input array of DC indexes used in activity products [nSKr]
    long int *arxlPh;  /// pointer to input array of linked phase indexes  [nlPh]

    double *arym;    /// Pointer to molalities of all species in MULTI (provided), read-only
    double *arla;    /// Pointer to lg activities of all species in MULTI (provided), read-only

    double *arnx;    /// Pointer to mole amounts of phase components (provided) [NComp] read-only

    double *arnxul;  /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    double *arnxll;  /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    double *arWx;    /// Species (end member) mole fractions ->NComp
    double *arVol;   /// molar volumes of end-members (species) cm3/mol ->NComp

// Work data and kinetic law calculation results

    TKinReact arPRt[]; /// work array of parameters and results for 'parallel reaction' terms [nPRk]

    double spcfu[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]
    double spcfl[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]

    double kTot;   /// Total rate constant (per m2 phase surface area)
    double rTot;   /// Current total MWR rate (mol/s)
    double vTot;   /// Total surface propagation velocity (nm/s)

    double sSAcor; /// Corrected specific surface area (m2/g)
    double sAph_c; /// Corrected surface area of the phase (m2/g)

    // Uptake model output

    // SS dissolution

    // SS precipitation

    // functions for allocation and initialization of kinetic rate tables
    void alloc_kinrtabs();
    long int init_kinrtabs( double *p_arlPhc, double *p_arrpCon,  double *p_arapCon,  double *p_arUmpCon );
    void free_kinrtabs();

  public:
    // Generic constructor
    TKinMet( const KinMetData *kmd );

    // Destructor
    virtual ~TKinMet();

    virtual long int SetTime()
    {
        return 0;
    };

    virtual long int PTparam()
    {
        return 0;
    };

    virtual long int RateMod()
    {
        return 0;
    };

    virtual long int SplitMod()
    {
        return 0;
    };

    virtual long int UptakeMod()
    {
        return 0;
    };

    // sets the specific surface area of the phase and 'parallel reactions' area fractions
    long int UpdateFSA( const double *fSAf_p, const double As );

    // returns modified specific surface area of the phase and 'parallel reactions' area fractions
    double ModifiedFSA ( double *fSAf_p );

    // sets new system TP state
    long int UpdatePT ( const double T_k, const double P_bar );

    // sets new time and time step
    bool UpdateTime( const double Tau, const double dTau );

    // Checks dimensions in order to re-allocate class instance, if necessary
    bool testSizes( const KinMetData *kmd );

};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - derived classes - - - - - - - - - - - - - -

class TKinTST: public TKinMet // Generic dissolution/precipitation models following Lasaga 1998, Shott ea 2012
{

    private:

    public:

    // Constructor
    TKinTST( KinMetData *kmd /*, specific params */ );

    // Destructor
     ~TKinTST();

    // Sets initial time
    long int SetTime();

    // Calculates T,P corrected kinetic parameters
    long int PTparam();

    // Calculates kinetic rates
    long int RateMod();

    long int SplitMod();

    // Calculates uptake rates
    long int UptakeMod();

};

class TPalandri: public TKinMet  // Dissolution/precipitation models following Palandri 2004
{
    private:

 // specific stuff for Palandri form

        // internal functions
//        void alloc_internal();
//        void free_internal();

public:
        // Constructor
        TPalandri( KinMetData *kmd /*, specific params */ );

        // Destructor
        ~TPalandri();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        long int SplitMod();

        // Calculates uptake rates
        long int UptakeMod();

};


class TWolthers: public TKinMet  // Calcite growth and uptake model
{
    private:

 // specific stuff for Wolthers form

        // internal functions
//        void alloc_internal();
//        void free_internal();

public:
        // Constructor
        TWolthers( KinMetData *kmd /*, specific params */ );

        // Destructor
        ~TWolthers();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        long int SplitMod();

        // Calculates uptake rates
        long int UptakeMod();

};

class TLaidler: public TKinMet  // Laidler d/p and uptake model
{
    private:

 // specific stuff for Laidler form

        // internal functions
//        void alloc_internal();
//        void free_internal();

public:

        // Constructor
        TLaidler( KinMetData *kmd /*, specific params */ );

        // Destructor
        ~TLaidler();

        // Sets initial time
        long int SetTime();

        // Calculates T,P corrected kinetic parameters
        long int PTparam();

        // Calculates kinetic rates
        long int RateMod();

        long int SplitMod();

        // Calculates uptake rates
        long int UptakeMod();


};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif // S_KINMET_H
