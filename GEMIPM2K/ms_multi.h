//-------------------------------------------------------------------
// $Id: ms_multi.h 1388 2009-08-07 14:44:31Z gems $
//
// Declaration of TMulti class, configuration, and related functions
// based on the IPM work data structure MULTI that represents chemical
// thermodynamic multisystem work data for GEM IPM-2 algorithm
//
// Rewritten from C to C++ by S.Dmytriyeva
// Copyright (C) 1995,2008 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of the
// standalone GEMIPM2K code (define IPMGEMPLUGIN).
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _ms_multi_h_
#define _ms_multi_h_

#ifdef Use_qd_real
// QD_real is enabled only if the above compiler key is used (experimental)
#include <qd/qd_real.h>
#include <qd/fpu.h>
#endif

#ifndef IPMGEMPLUGIN

#include "m_param.h"
#include "v_ipnc.h"
// Internal subroutine for ET_translate() to process Phase scripts
typedef int (tget_ndx)( int nI, int nO, int Xplace );

#else
#include <time.h>
#include "m_const.h"

#endif

#include "s_fgl.h"

typedef struct
{  // MULTI is base structure to Project (local values)
  char
    stkey[EQ_RKLEN+1],   // Record key identifying IPM minimization problem
    // NV_[MAXNV], nulch, nulch1, // Variant Nr for fixed b,P,T,V; index in a megasystem
    PunE,         // Units of energy  { j;  J c C N reserved }
    PunV,         // Units of volume  { j;  c L a reserved }
    PunP,         // Units of pressure  { b;  B p P A reserved }
    PunT;         // Units of temperature  { C; K F reserved }
  long int
    N,        	// N - number of IC in IPM problem
    NR,       	// NR - dimensions of R matrix
    L,        	// L -   number of DC in IPM problem
    Ls,       	// Ls -   total number of DC in multi-component phases
    LO,       	// LO -   index of water-solvent in IPM DC list
    PG,       	// PG -   number of DC in gas phase
    PSOL,     	// PSOL - number of DC in liquid hydrocarbon phase
    Lads,     	// Total number of DC in sorption phases included into this system.
    FI,       	// FI -   number of phases in IPM problem
    FIs,      	// FIs -   number of multicomponent phases
    FIa,      	// FIa -   number of sorption phases
    FI1,     // FI1 -   number of phases present in eqstate
    FI1s,    // FI1s -   number of multicomponent phases present in eqstate
    FI1a,    // FI1a -   number of sorption phases present in eqstate
    IT,      // It - number of completed IPM iterations
    E,       // PE - flag of electroneutrality constraint { 0 1 }
    PD,      // PD - mode of calling GammaCalc() { 0 1 2 3 4 }
    PV,      // Flag for the volume balance constraint (on Vol IC) - for indifferent equilibria at P_Sat { 0 1 }
    PLIM,    // PU - flag of activation of DC/phase restrictions { 0 1 }
    Ec,    // GammaCalc() return code: 0 (OK) or 1 (error)
    K2,    // Number of IPM loops performed ( >1 up to 3 because of PhaseSelection() )
    PZ,    // Indicator of PhaseSelection() status (since r1594): 0 untouched, 1 phase(s) inserted
              // 2 insertion done after 3 major IPM loops
    pNP, //Mode of FIA selection: 0-auto-SIMPLEX,1-old eqstate,-1-user's choice
    pESU,  // Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pIPN,  // State of IPN-arrays:  0-create; 1-available; -1 remake
    pBAL,  // State of reloading CSD:  1- BAL only; 0-whole CSD
tMin,  // Type of thermodynamic potential to minimize
    pTPD,  // State of reloading thermod data: 0- all  1 - G0 only  2 - no
    pULR,  // Start recalc kinetic constraints (0-do not, 1-do )internal
    ITaia,  // Number of IPM iterations completed in AIA mode (renamed from pRR1)
    FIat,   // max. number of surface site types
    MK,     // PM return code: 0 - continue;  1 - converged
    W1,     // Indicator of CleanupSpeciation() status (since r1594) 0 untouched, -1 phase(s) removed, 1 some DCs inserted
    is,     // is - index of IC for IPN equations ( GammaCalc() )
    js,     // js - index of DC for IPN equations ( GammaCalc() )
    next,
    sitNcat,    // SIT: number of cations
    sitNan      // SIT: number of anions
    ;
  double
    TC,TCc, 	// Temperature T, min.-max. (0,2000 C)
    T,Tc,   	// T, min.-max. K
    P,Pc,   	// Pressure P, min.-max.(0,10000 bar)
    VX_,VXc,    // V(X) - volume of the system, min.-max., cm3
    GX_,GXc,    // Gibbs potential of the system G(X), min.-max. (J)
    AX_,AXc,    // Helmholtz potential of the system F(X), reserved
    UX_,UXc,  	// Internal energy of the system U(X), reserved
    HX_,HXc, 	// Total enthalpy of the system H(X), reserved
    SX_,SXc, 	// Total entropy of the system S(X), reserved
    CpX_,CpXc,  // reserved
    CvX_,CvXc,  // reserved
    TMols,      // Input total moles in b vector before rescaling
    SMols,      // Standart total moles (upscaled) {1000}
    MBX,        // Total mass of the system, kg
    FX,    	// Current Gibbs potential of the system in IPM, moles
    IC,         // Effective molal ionic strength of aqueous electrolyte
    pH,         // pH of aqueous solution
    pe,         // pe of aqueous solution
    Eh,         // Eh of aqueous solution, V
    DHBM,       // balance (relative) precision criterion
    DSM,        // min amount of phase DS
    GWAT,       // used in ipm_gamma()
    YMET,       // reserved
    PCI,        // Current value of Dikin criterion of IPM convergence DK>=DX
    DXM,        // IPM convergence criterion threshold DX (1e-5)
    lnP,        // log Ptotal
    RT,         // RT: 8.31451*T (J/mole/K)
    FRT,        // F/RT, F - Faraday constant = 96485.309 C/mol
    Yw,         // Current number of moles of solvent in aqueous phase
    ln5551,     // ln(55.50837344)
    aqsTail,    // v_j asymmetry correction factor for aqueous species
    lowPosNum,  // Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
    logXw,      // work variable
    logYFk,     // work variable
    YFk,        // Current number of moles in a multicomponent phase
    FitVar[5];  // internal; FitVar[0] is total mass (g) of solids in the system (sum over the BFC array)
                //      FitVar[1], [2] reserved
                //       FitVar[4] is the AG smoothing parameter;
                //       FitVar[3] is the actual smoothing coefficient
double
  denW[5],   // Density of water, first T, second T, first P, second P derivative for Tc,Pc
  denWg[5],  // Density of steam for Tc,Pc
  epsW[5],   // Diel. constant of H2O(l)for Tc,Pc
  epsWg[5];  // Diel. constant of steam for Tc,Pc

  long int
    *L1,    // l_a vector - number of DCs included into each phase [Fi]
    *LsMod, // Number of interaction parameters, max. parameter order (cols in IPx),
        // and number of coefficients per parameter in PMc table [3*FIs]
    *LsMdc, // Number of parameters per component of the phase for the non-ideal mixing models [FIs]
    *IPx,   // Collected indexation table for interaction parameters of non-ideal solutions
            // ->LsMod[k,0] x LsMod[k,1]   over FIs
    *mui,   // IC indices in RMULTS IC list [N]
    *muk,   // Phase indices in RMULTS phase list [FI]
    *muj;   // DC indices in RMULTS DC list [L]
  long int  (*SATX)[4]; // Setup of surface sites and species (will be applied separately within each sorption phase) [Lads]
             // link indexes to surface type [XL_ST]; sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
  double
    *PMc,    // Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
    *DMc,    // Non-ideality coefficients f(TPX) for DC -> LsMdc[k]
    *A,      // DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,    // IC atomic (molar) mass, g/mole [0:N-1]

 // Reconsider usage
     *Wb,     //Relative Born factors (HKF, reserved) [0:Ls-1]
     *Wabs,   // Absolute Born factors (HKF, reserved) [0:Ls-1]
     *Rion,   // Ionic or solvation radii, A (reserved) [0:Ls-1]
     *HYM,    // reserved
     *ENT,    // reserved no object

 // Convert H0, A0, U0, S0, Cp0 to double
     *H0,     // DC pmolar enthalpies, reserved [L]
     *A0,     // DC molar Helmholtz energies, reserved [L]
     *U0,     // DC molar internal energies, reserved [L]
     *S0,     // DC molar entropies, reserved [L]
     *Cp0,    // DC molar heat capacity, reserved [L]
     *Cv0,    // DC molar Cv, reserved [L]

    *VL,        // ln mole fraction of end members in phases-solutions
    *Xcond, 	// conductivity of phase carrier, sm/m2   [0:FI-1], reserved
    *Xeps,  	// diel.permeability of phase carrier (solvent) [0:FI-1], reserved
    *Aalp,  	// Full vector of specific surface areas of phases (m2/g) [0:FI-1]
    *Sigw,  	// Specific surface free energy for phase-water interface (J/m2)   [0:FI-1]
    *Sigg  	// Specific surface free energy for phase-gas interface (J/m2) (not yet used)  [0:FI-1], reserved
    ;


//  Data for surface comlexation and sorption models (new variant [Kulik,2006])
  double  (*Xr0h0)[2];   // mean r & h of particles (- pores), nm  [0:FI-1][2], reserved
  double  (*Nfsp)[MST];  // Fractions of the sorbent specific surface area allocated to surface types  [FIs][FIat]
  double  (*MASDT)[MST]; // Total maximum site  density per surface type (mkmol/g)  [FIs][FIat]
  double  (*XcapF)[MST]; // Capacitance density of Ba EDL layer F/m2 [FIs][FIat]
  double  (*XcapA)[MST]; // Capacitance density of 0 EDL layer, F/m2 [FIs][FIat]
  double  (*XcapB)[MST]; // Capacitance density of B EDL layer, F/m2 [FIs][FIat]
  double  (*XcapD)[MST]; // Eff. cap. density of diffuse layer, F/m2 [FIs][FIat]
  double  (*XdlA)[MST];  // Effective thickness of A EDL layer, nm [FIs][FIat], reserved
  double  (*XdlB)[MST];  // Effective thickness of B EDL layer, nm [FIs][FIat], reserved
  double  (*XdlD)[MST];  // Effective thickness of diffuse layer, nm [FIs][FIat], reserved
  double  (*XlamA)[MST]; // Factor of EDL discretness  A < 1 [FIs][FIat], reserved
  double  (*Xetaf)[MST]; // Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [FIs][FIat]
  double  (*MASDJ)[DFCN];  // Parameters of surface species in surface complexation models [Lads][DFCN]
                          // Contents defined in the enum below this structure
// Other data

  double
    *XFs,    // Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,  // Current Karpov criteria of phase stability  F_a [0:FI-1]
    *Fug,    // Demo partial fugacities of gases [0:PG-1]
    *Fug_l,  // Demo log partial fugacities of gases [0:PG-1]
    *Ppg_l,     // Demo log partial pressures of gases [0:PG-1]

    *DUL,     // VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,     // NG Vector of lower kinetic restrictions to x_j, moles [L]
    *GEX,     // Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
    *PUL,  // Vector of upper restrictions to phases amounts X_a (reserved)[FIs]
    *PLL,  // Vector of lower restrictions to phases amounts X_a (reserved)[FIs]
    *YOF,     // Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]
    *Vol,     // DC molar volumes, cm3/mol [L]
    *MM,      // DC molar masses, g/mol [L]
    *Pparc,   // Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
    *Y_m,     // Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    // log activity of DC in multi-component phases (mju-mji0) [0:Ls-1]
    *Y_w,     // Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
    *Gamma,   // DC activity coefficients in molal or other phase-specific scale [0:L-1]
    *lnGmf,   // ln of initial DC activity coefficients for correcting G0 [0:L-1]
    *lnGmM,   // ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]
    *EZ,      // Formula charge of DC in multi-component phases [0:Ls-1]
    *FVOL,    // phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *FWGT,    // phase (carrier) masses, g                [0:FI-1]
//
    *G,    // Normalized DC energy function c(j), mole/mole [0:L-1]
    *G0,   // Input normalized g0_j(T,P) for DC at unified standard scale[L]
    *lnGam, // ln of DC activity coefficients in unified (mole-fraction) scale [0:L-1]
    *lnGmo; // Copy of lnGam from previous IPM iteration (reserved)
  double  (*lnSAC)[4]; // former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]
  double  *B,  // Input bulk chem. compos. of the system - b vector, moles of IC[N]
    *U,  // IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *U_r,  // IC chemical potentials u_i (J/mole) [0:N-1]
    *C,    // Calculated IC mass-balance deviations (moles) [0:N-1]
    *IC_m, // Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    *IC_lm,	// log total IC molalities in aqueous phase [0:N-1]
    *IC_wm,	// Total dissolved IC concentrations in g/kg_soln [0:N-1]
    *BF,    //Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC,   //Total output bulk compositions of all solid phases[1][N]
    *XF,    // Output total number of moles of phases Xa[0:FI-1]
    *YF,    // Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   // Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   // Approximation of XFA in the next IPM iteration [0:FIs-1]
    *Falp;  // Karpov phase stability criteria F_a [0:FI-1] or phase stability index (PC==2)

  double (*VPh)[MIXPHPROPS],     // Volume properties for mixed phases [FIs]
         (*GPh)[MIXPHPROPS],     // Gibbs energy properties for mixed phases [FIs]
  		 (*HPh)[MIXPHPROPS],     // Enthalpy properties for mixed phases [FIs]
         (*SPh)[MIXPHPROPS],     // Entropy properties for mixed phases [FIs]
         (*CPh)[MIXPHPROPS],     // Heat capacity Cp properties for mixed phases [FIs]
         (*APh)[MIXPHPROPS],     // Helmholtz energy properties for mixed phases [FIs]
         (*UPh)[MIXPHPROPS];     // Internal energy properties for mixed phases [FIs]

// EDL models (data for electrostatic activity coefficients)
   double (*XetaA)[MST]; // Total EDL charge on A (0) EDL plane, moles [FIs][FIat]
   double (*XetaB)[MST]; // Total charge of surface species on B (1) EDL plane, moles[FIs][FIat]
   double (*XetaD)[MST]; // Total charge of surface species on D (2) EDL plane, moles[FIs][FIat]
   double (*XpsiA)[MST]; // Relative potential at A (0) EDL plane,V [FIs][FIat]
   double (*XpsiB)[MST]; // Relative potential at B (1) EDL plane,V [FIs][FIat]
   double (*XpsiD)[MST]; // Relative potential at D (2) plane,V [FIs][FIat]
   double (*XFTS)[MST];  // Total number of moles of surface DC at surface type [FIs][FIat]
//
   double *X,  // DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,   // Copy of x_j from previous IPM iteration [0:L-1]
    *XY,  // Copy of x_j from previous loop of Selekt2() [0:L-1]
    *Qp,  // Work variables related to non-ideal phases FIs*(QPSIZE=60)
    *Qd,  // Work variables related to DC in non-ideal phases FIs*(QDSIZE=60)
    *MU,  // mu_j values of differences between dual DC chem.potentials [L]
    *EMU, // Exponents of DC increment to F_a criterion for phase [L]
    *NMU, // DC increments to F_a criterion for phase [L]
    *W,   // Weight multipliers for DC (incl restrictions) in IPM [L]
    *Fx,  // Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,  // Mole fractions Wx of DC in multi-component phases [L]
    *F,   //Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0;  // Excess Gibbs energies for (metastable) DC, mole/mole [L]
   double (*D)[MST];  // Reserved; new work array for calc. surface act.coeff.
// Name lists
  char  (*sMod)[6];   // Codes for built-in mixing models of multicomponent phases [FIs]
  char  (*SB)[MAXICNAME+MAXSYMB]; // List of IC names in the system [N]
  char  (*SB1)[MAXICNAME]; // List of IC names in the system [N]
  char  (*SM)[MAXDCNAME];  // List of DC names in the system [L]
  char  (*SF)[MAXPHNAME+MAXSYMB];  // List of phase names in the system [FI]
  char  (*SM2)[MAXDCNAME];  // List of multicomp. phase DC names in the system [Ls]
  char  (*SM3)[MAXDCNAME];  // List of adsorption DC names in the system [Lads]
  char  *DCC3;   // Classifier of DCs involved in sorption phases [Lads]
  char  (*SF2)[MAXPHNAME+MAXSYMB]; // List of multicomp. phase names in the syst [FIs]
  char  (*SFs)[MAXPHNAME+MAXSYMB];
    // List of phases currently present in non-zero quantities [FI]
  char  *pbuf, 	// Text buffer for table printouts
// Class codes
    *RLC,   // Code of metastability constraints for DCs [L] enum DC_LIMITS
    *RSC,   // Units of metastability/kinetic constraints for DCs  [L]
    *RFLC,  // Classifier of restriction types for XF_a [FIs]
    *RFSC,  // Classifier of restriction scales for XF_a [FIs]
    *ICC,   // Classifier of IC { e o h a z v i <int> } [N]
    *DCC,   // Classifier of DC { TESWGVCHNIJMDRAB0123XYZPQO } [L]
    *PHC;   // Classifier of phases { a g f p m l x d h } [FI]
  char  (*SCM)[MST]; // Classifier of built-in electrostatic models applied to surface types in sorption phases [FIs][FIat]
  char  *SATT,  // Classifier of applied SACT equations (isotherm corrections) [Lads]
    *DCCW;  // internal DC class codes [L]
//  long int
//     *sitXcat, // SIT: indices of cations (may be changed soon)
//     *sitXan;  // SIT: indices of anions
//  double
//     *sitE;    // pointer to SIT coeff. table (may be changed soon)
  long int ITF,       // Number of completed IA EFD iterations
         ITG;        // Number of completed GEM IPM iterations
  clock_t t_start, t_end;
  double t_elap_sec;  // work variables for determining IPM calculation time
#ifdef IPMGEMPLUGIN
  double *Guns;  //  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
  double *Vuns;  //  mu.L work vector of uncertainty space increments to tp->Vm
  double *tpp_G; // Partial molar(molal) Gibbs energy g(TP) (always), J/mole
  double *tpp_S;    // Partial molar(molal) entropy s(TP), J/mole/K
  double *tpp_Vm;   // Partial molar(molal) volume Vm(TP) (always), J/bar
#endif

  // additional arrays for internal calculation in ipm_main
  double *XU; //dual-thermo calculation of DC amount X(j) from A matrix and u vector [L]
  double *Uc; // Internal copy of IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
  double *Uefd; // Internal copy of IC chemical potentials u_i (mole/mole) - EFD function [N]
  char errorCode[100]; //  code of error in IPM      (Ec number of error)
  char errorBuf[1024]; // description of error in IPM
  double logCDvalues[5]; // Collection of lg Dikin crit. values for the new smoothing equation
//  qd_real qdFX;    	// Current Gibbs potential of the system in IPM, moles

  double // Iterators for MTP interpolation (do not load/unload for IPM)
Pai[4],  // Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
Tai[4];  // Temperature T, C: start, end, increment for MTP array in DataCH , Ttol

  // Experimental: modified cutoff and insertion values (DK 28.04.2010)
  double
// cutoffs (rescaled to system size)
  XwMinM,// Cutoff mole amount for elimination of water-solvent { 1e-13 }
  ScMinM,// Cutoff mole amount for elimination of solid sorbent { 1e-13 }
  DcMinM,// Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
  PhMinM,// Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
// insertion values (re-scaled to system size)
  DFYwM, // Insertion mole amount for water-solvent { 1e-6 }
  DFYaqM,// Insertion mole amount for aqueous and surface species { 1e-6 }
  DFYidM,// Insertion mole amount for ideal solution components { 1e-6 }
  DFYrM, // Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
  DFYhM, // Insertion mole amount for minor solution components { 1e-6 }
  DFYcM, // Insertion mole amount for single-component phase { 1e-6 }
  DFYsM, // Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
  SizeFactor; // factor for re-scaling the cutoffs/insertions to the system size
}
MULTI;

enum {
	//[0] - max site density in mkmol/(g sorbent); [1] - species charge allocated to 0 plane;
	//[2] - surface species charge allocated to beta -or third plane; [3] - Frumkin interaction parameter;
	//[4] species denticity or coordination number; [5]  - reserved parameter (e.g. species charge on 3rd EIL plane)]
   XL_ST = 0, XL_EM, XL_SI, XL_SP
};

// Data of MULTI
class TMulti
#ifndef IPMGEMPLUGIN
  : public TSubModule
#endif
{
    MULTI pm;
    MULTI *pmp;

// Internal arrays for the performance optimization  (since version 2.0.0)
   long int sizeN; /*, sizeL, sizeAN;*/
   double *AA;
   double *BB;
#ifdef Use_qd_real
   qd_real *qdAA;
   qd_real *qdBB;
#endif
   long int *arrL;
   long int *arrAN;

   void Alloc_A_B( long int newN );
   void Free_A_B();
   void Build_compressed_xAN();
   void Free_compressed_xAN();
   void Free_internal();

   long int sizeFIs;     // current size of phSolMod
   TSolMod* (*phSolMod); // size current FIs -   number of multicomponent phases

   void Alloc_TSolMod( long int newFIs );
   void Free_TSolMod();


#ifndef IPMGEMPLUGIN
// These pointers and methods are only used in GEMS-PSI
    SYSTEM *syp;
    MTPARM *tpp;
    RMULTS *mup;

    void MultiSystemInit();
    void multi_sys_dc();
    void multi_sys_ph();
    void ph_sur_param( int k, int kk );
    void ph_surtype_assign( int k, int kk, int jb, int je,
                            int car_l[], int car_c, int Cjs );
    void sm_text_analyze( int nph, int Type, int JB, int JE, int jb, int je );
    void SolModLoad();
    bool CompressPhaseIpxt( int kPH );
    gstring PressSolMod( int nP );
    char *ExtractEG( char *Etext, int jp, int *EGlen, int Nes );
    int find_icnum( char *name, int LNmode );
    int find_dcnum( char *name, int jb, int je, int LNmode, char *stmt  );
    int find_phnum( char *name, int LNmode );
    int find_acnum( char *name, int LNmode );

#else

   char PAalp_; // Flag for using (+) or ignoring (-) specific surface areas of phases
   char PSigm_; // Flag for using (+) or ignoring (-) specific surface free energies

#endif

// Internal functions for SCMs
   void getLsModsum( long int& LsModSum, long int& LsIPxSum );
   void getLsMdcsum( long int& LsMdcSum );

// ipm_chemical.cpp
    void XmaxSAT_IPM2();
    void XmaxSAT_IPM2_reset();
    double DualChemPot( double U[], double AL[], long int N, long int j );
    void Set_DC_limits( long int Mode );
    void TotalPhases( double X[], double XF[], double XFA[] );
    double Ej_init_calc( double, long int j, long int k);
    double  PrimalDC_ChemPot(  double G,  double logY,  double logYF,
                           double asTail,  double logYw,  char DCCW );
    void PrimalChemicalPotentials( double F[], double Y[],
                                  double YF[], double YFA[] );
    double KarpovCriterionDC( double *dNuG, double logYF, double asTail,
                 double logYw, double Wx,  char DCCW );
    void f_alpha();
    void  StabilityIndexes( );   // added 01.05.2010 DK
    double FreeEnergyIncr(   double G,  double x,  double logXF,
                             double logXw,  char DCCW );
    double GX( double LM  );
    long int  Mol_u( double Y[], double X[], double XF[], double XFA[] );
    void ConvertDCC();
    long int  getXvolume();

// ipm_chemical2.cpp
    void GasParcP();
    void phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] );
    void phase_bfc( long int k, long int jj );
    double bfc_mass( void );
    void ConCalcDC( double X[], double XF[], double XFA[],
                    double Factor, double MMC, double Dsur, long int jb, long int je, long int k );
    void ConCalc( double X[], double XF[], double XFA[]);
    long int GouyChapman(  long int jb, long int je, long int k );

//  Surface activity coefficient terms
    long int SurfaceActivityCoeff( long int jb, long int je, long int jpb, long int jdb, long int k );
//    void SurfaceActivityTerm( long int jb, long int je, long int k );  // Obsolete / deleted
    double PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag = 0L ); // Added 26.06.08

// ipm_chemical3.cpp
    void IS_EtaCalc();
    void pm_GC_ods_link( long int k, long int jb, long int jpb, long int jdb, long int ipb );
    double SmoothingFactor( );
    void SetSmoothingFactor( long int mode ); // new smoothing function (3 variants)
// Main call for calculation of activity coefficients on IPM iterations
    long int GammaCalc( long int LinkMode );
// Built-in activity coefficient models
// Generic solution model calls
    void SolModCreate( long int jb, long int je, long int jpb, long int jdb, long int k, long int ipb,
    		char ModCode, char MixCode );
    void SolModParPT( long int k, char ModCode );
    void SolModActCoeff( long int k, char ModCode );
    void SolModExcessProp( long int k, char ModCode );
    void SolModIdealProp ( long int jb, long int k, char ModCode );
    void SolModStandProp ( long int jb, long int k, char ModCode );
    void SolModDarkenProp ( long int jb, long int k, char ModCode );
// Specific phase property calculation functions
    void IdealGas( long int jb, long int k, double *Zid );
    void IdealOneSite( long int jb, long int k, double *Zid );
    void IdealMultiSite( long int jb, long int k, double *Zid );

// ipm_main.cpp - numerical part of GEM IPM-2
    void MultiCalcMain( long int rLoop );
    long int EnterFeasibleDomain( long int WhereCalledFrom );
    long int InteriorPointsMethod( long int &status, long int rLoop );
    void SimplexInitialApproximation( );

// ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
   void MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                               double *B, double *C );
//   long int CheckMassBalanceResiduals(double *Y );
   double LMD( double LM );
   void ZeroDCsOff( long int jStart, long int jEnd, long int k=-1L );
   void RaiseZeroedOffDCs( long int jStart, long int jEnd, /*double sfactor,*/ long int k=-1L );
   double RaiseDC_Value( const long int j );
   //   void LagrangeMultiplier();
   long int MetastabilityLagrangeMultiplier();
   void WeightMultipliers( bool square );
   long int SolverLinearEquations( long int N, bool initAppr );
   double calcDikin(  long int N, bool initAppr );
   double calcLM(  bool initAppr );
   void Restoring_Y_YF();
//   double calcSfactor();
   double RescaleToSize( bool standard_size ); // replaced calcSfactor() 30.08.2009 DK
   long int CleanupSpeciation( double AmountThreshold, double ChemPotDiffCutoff ); // added 25.03.10 DK
   long int PhaseSelection( long int &k_miss, long int &k_unst, long int rLoop );  // added 01.05.10 DK
   long int PhaseSelect( long int &k_miss, long int &k_unst, long int rLoop );
   bool AutoInitialApprox();

   // IPM_SIMPLEX.CPP Simplex method with two-sided constralong ints (Karpov ea 1997)
    void Simplex(long int M, long int N, long int T, double GZ, double EPS,
                 double *UND, double *UP, double *B, double *U,
                 double *AA, long int *STR, long int *NMB );
    void SPOS( double *P, long int STR[],long int NMB[],long int J,long int M,double AA[]);
    void START( long int T,long int *ITER,long int M,long int N,long int NMB[],
                double GZ,double EPS,long int STR[],long int *BASE,
                double B[],double UND[],double UP[],double AA[],double *A,
                double *Q );
    void NEW(long int *OPT,long int N,long int M,double EPS,double *LEVEL,long int *J0,
             long int *Z,long int STR[], long int NMB[], double UP[],
             double AA[], double *A);
    void WORK(double GZ,double EPS,long int *I0, long int *J0,long int *Z,long int *ITER,
              long int M, long int STR[],long int NMB[],double AA[],
              long int BASE[],long int *UNO,double UP[],double *A,double Q[]);
    void FIN(double EPS,long int M,long int N,long int STR[],long int NMB[],
             long int BASE[],double UND[],double UP[],double U[],
             double AA[],double *A,double Q[],long int *ITER);
    void GibbsMinimization();
    double calcTotalMoles( );
    void ScaleMulti(  double ScFact );
    void RescaleMulti(  double ScFact );
    void MultiConstInit(); // from MultiRemake
    void MultiCalcInit();


#ifdef Use_qd_real
// QD_real
    // ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
       void qdMassBalanceResiduals( long int N, long int L, double *A, double *Y,
                                   double *B, double *C );
       qd_real qdGX( double LM  );
       long int qdSolverLinearEquations( long int N, bool initAppr );
       double qdLMD( double LM );
       double qdcalcDikin(  long int N, bool initAppr );
#endif

public:

    void set_def( long int i=0);

#ifndef IPMGEMPLUGIN
// This is used only in GEMS-PSI
    TIArray<IPNCalc> qEp;
    TIArray<IPNCalc> qEd;

    TMulti( int nrt, SYSTEM* sy_, MTPARM *tp_, RMULTS *mu_ );
    ~TMulti()
    {  Free_internal(); };


    void ods_link( int i=0);
    void dyn_set( int i=0);
    void dyn_kill( int i=0);
    void dyn_new( int i=0);
//    void set_def( int i=0);
//    void sit_dyn_new();

    // ms_muleq.cpp
    void packData();
    void packData( TCIntArray PHon, TCIntArray DCon );
    void setSizes();
    void loadData( bool newRec );
    void unpackData();

    void MultiKeyInit( const char*key );
    void EqstatExpand( const char *key );
    void ET_translate( int nOet, int nOpex, int JB, int JE, int jb, int je,
     tget_ndx *get_ndx = 0 );
    void getNamesList( int nO, TCStringArray& lst );

   class UserCancelException {};
#else
// this allocation is used only in standalone GEMIPM2K
   TMulti()
   {
	 pmp = &pm;
     sizeN = 0;
     AA = 0;
     BB = 0;
     arrL = 0;
     arrAN = 0;

     sizeFIs = 0;
     phSolMod = 0;

     pmp->Guns = 0;
     pmp->Vuns = 0;
     pmp->tpp_G = 0;
     pmp->tpp_S = 0;
     pmp->tpp_Vm = 0;
   }

    ~TMulti()
    {  multi_free(); }

    void multi_realloc( char PAalp, char PSigm );
    void multi_free();

#endif

    MULTI* GetPM()
    { return &pm; }

    const char* GetName() const
    {  return "Multi";  }

   //connection to mass transport
    void to_file( GemDataStream& ff );
    void to_text_file( const char *path, bool append=false  );
    void from_file( GemDataStream& ff );
    void to_text_file_gemipm( const char *path, bool addMui,
    		bool with_comments = true, bool brief_mode = false );
    void from_text_file_gemipm( const char *path );

    // EXTERNAL FUNCTIONS
    // MultiCalc
    void Alloc_internal();
    double calcEqustat( long int typeMin, long int& NumIterFIA, long int& NumIterIPM );
    void MultiInit();
    void CompG0Load();
    void setErrorMessage( long int num, const char *code, const char * msg);
    void addErrorMessage( const char * msg);

   long int CheckMassBalanceResiduals(double *Y );
   double Cj_init_calc( double g0, long int j, long int k );

// connection to UnSpace
    double pb_GX( double *Gxx  );
};

// ???? syp->PGmax
typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP    =  'G',   // Gibbs energy minimization G(T,P)
    A_TV    =  'A',   // Helmholts energy minimization A(T,V)
    U_SV    =  'U',   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS    =  'H',   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH   =  '1',   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV   =  '2'    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} THERM_POTENTIALS;

typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP_    =  0,   // Gibbs energy minimization G(T,P)
    A_TV_    =  1,   // Helmholts energy minimization A(T,V)
    U_SV_    =  2,   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS_    =  3,   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH_   =  4,   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV_   =  5    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} NUM_POTENTIALS;


#endif   //_ms_multi_h

