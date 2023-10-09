//-------------------------------------------------------------------
/// \file tsolmod_multi.h
/// Declaration of subset of TMulti class, configuration, and related functions
/// based on the IPM work data structure MULTI
//
// Copyright (c) 2023 S.Dmytriyeva, D.Kulik, T.Wagner
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
//-------------------------------------------------------------------
//
#ifndef MS_MULTI_BASE_H
#define MS_MULTI_BASE_H

#include "m_const_base.h"
#include "gems3k_impex.h"
#include "solmodcalc.h"
#include "gems3k_impex.h"
#include "datach.h"
#include "databr.h"

#ifdef USE_THERMOFUN
#include "ThermoFun/ThermoFun.h"
#endif

class GemDataStream;

const int  QPSIZE = 180, // earlier 20, 40 SD oct 2005
           QDSIZE = 60;

// Physical constants - see m_param.cpp or ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;

struct BASE_PARAM /// Flags and thresholds for numeric modules
{
   short
           PC,   ///< Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
           PD,   ///< abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }.
                 ///< Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 ///< if PD < 0 then use test qd_real accuracy mode
           PRD,  ///< Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
           PSM,  ///< Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   ///< Maximum allowed number of iterations in the MassBalanceRefinement() procedure {  30 }
           DW,   ///< Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   ///< Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, ///< IPM tolerance for detecting divergence in dual solution { 10; range 1 to 1000; 0 disables the detection }
           PE,   ///< Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
           IIM   ///< Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 1000 }
           ;
         double DG,   ///< Standart total moles { 1e5 }
           DHB,  ///< Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-10 }
           DS,   ///< Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
           DK,   ///< IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-4) { 1e-5 }
           DF,   ///< Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
           DFM,  ///< Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.1 }
           DFYw, ///< Insertion mole amount for water-solvent { 1e-6 }
           DFYaq,///< Insertion mole amount for aqueous species { 1e-6 }
           DFYid,///< Insertion mole amount for ideal solution components { 1e-6 }
           DFYr, ///< Insertion mole amount for major solution components { 1e-6 }
           DFYh, ///< Insertion mole amount for minor solution components { 1e-6 }
           DFYc, ///< Insertion mole amount for single-component phase { 1e-6 }
           DFYs, ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
           DB,   ///< Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
           AG,   ///< Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { -1 }
           DGC,  ///< Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { -0.99 }
           GAR,  ///< Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
           GAH,  ///< Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
           GAS,  ///< Since r1583/r409: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }.
                 ///< before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
           DNS,  ///< Standard surface density (nm-2) for calculating activity of surface species (12.05)
           XwMin,///< Cutoff mole amount for elimination of water-solvent { 1e-9 }
           ScMin,///< Cutoff mole amount for elimination of solid sorbent {1e-7}
           DcMin,///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
           PhMin,///< Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
           ICmin,///< Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
           EPS,  ///< Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
           IEPS, ///< Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
           DKIN; ///< Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
    char *tprn;       ///< internal
};


typedef struct
{  // MULTI is base structure to Project (local values)
    char
    stkey[EQ_RKLEN+5],   ///< Record key identifying IPM minimization problem
    // NV_[MAXNV], nulch, nulch1, ///< Variant Nr for fixed b,P,T,V; index in a megasystem
    PunE,         ///< Units of energy  { j;  J c C N reserved }
    PunV,         ///< Units of volume  { j;  c L a reserved }
    PunP,         ///< Units of pressure  { b;  B p P A reserved }
    PunT;         ///< Units of temperature  { C; K F reserved }
    long int
    N,        	///< N - number of IC in IPM problem
    NR,       	///< NR - dimensions of R matrix
    L,        	///< L -   number of DC in IPM problem
    Ls,       	///< Ls -   total number of DC in multi-component phases
    LO,       	///< LO -   index of water-solvent in IPM DC list
    PG,       	///< PG -   number of DC in gas phase
    PSOL,     	///< PSOL - number of DC in liquid hydrocarbon phase
    Lads,     	///< Total number of DC in sorption phases included into this system.
    FI,       	///< FI -   number of phases in IPM problem
    FIs,      	///< FIs -   number of multicomponent phases
    FIa,      	///< FIa -   number of sorption phases
    FI1,     ///< FI1 -   number of phases present in eqstate
    FI1s,    ///< FI1s -   number of multicomponent phases present in eqstate
    FI1a,    ///< FI1a -   number of sorption phases present in eqstate
    IT,      ///< It - number of completed IPM iterations
    E,       ///< PE - flag of electroneutrality constraint { 0 1 }
    PD,      ///< PD - mode of calling CalculateActivityCoefficients() { 0 1 2 3 4 }
    PV,      ///< Flag for the volume balance constraint (on Vol IC) - for indifferent equilibria at P_Sat { 0 1 }
    PLIM,    ///< PU - flag of activation of DC/phase restrictions { 0 1 }
    Ec,     ///< CalculateActivityCoefficients() return code: 0 (OK) or 1 (error)
    K2,     ///< Number of IPM loops performed ( >1 up to 6 because of PSSC() )
    PZ,     ///< Indicator of PSSC() status (since r1594): 0 untouched, 1 phase(s) inserted
    ///< 2 insertion done after 5 major IPM loops
    pNP,    ///< Mode of FIA selection: 0-automatic-LP AIA, 1-smart SIA, -1-user's choice
    pESU,   ///< Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pIPN,   ///< State of IPN-arrays:  0-create; 1-available; -1 remake
    pBAL,   ///< State of reloading CSD:  1- BAL only; 0-whole CSD
    tMin,   ///< Type of thermodynamic potential to minimize
    pTPD,   ///< State of reloading thermod data: 0-all  -1-full from database   1-new system 2-no
    pULR,   ///< Start recalc kinetic constraints (0-do not, 1-do )internal
    pKMM, ///< new: State of KinMet arrays: 0-create; 1-available; -1 remake
    ITaia,  ///< Number of IPM iterations completed in AIA mode (renamed from pRR1)
    FIat,   ///< max. number of surface site types
    MK,     ///< IPM return code: 0 - continue;  1 - converged
    W1,     ///< Indicator ofSpeciationCleanup() status (since r1594) 0 untouched, -1 phase(s) removed, 1 some DCs inserted
    is,     ///< is - index of IC for IPN equations ( CalculateActivityCoefficients() )
    js,     ///< js - index of DC for IPN equations ( CalculateActivityCoefficients() )
    next,   ///< for IPN equations (is it really necessary? TW please check!
    sitNcat,    //< Can be re-used
    sitNan      // Can be re-used
    ;
    double
    TC,  	///< Temperature T, min. (0,2000 C)
    TCc, 	///< Temperature T, max. (0,2000 C)
    T,   	///< T, min. K
    Tc,   	///< T, max. K
    P,      ///< Pressure P, min(0,10000 bar)
    Pc,   	///< Pressure P, max.(0,10000 bar)
    VX_,    ///< V(X) - volume of the system, min., cm3
    VXc,    ///< V(X) - volume of the system, max., cm3
    GX_,    ///< Gibbs potential of the system G(X), min. (J)
    GXc,    ///< Gibbs potential of the system G(X), max. (J)
    AX_,    ///< Helmholtz potential of the system F(X)
    AXc,    ///<  reserved
    UX_,  	///< Internal energy of the system U(X)
    UXc,  	///<  reserved
    HX_,    ///< Total enthalpy of the system H(X)
    HXc, 	///<  reserved
    SX_,    ///< Total entropy of the system S(X)
    SXc,   ///<	 reserved
    CpX_,  ///< reserved
    CpXc,  ///< 20 reserved
    CvX_,  ///< reserved
    CvXc,  ///< reserved
    // TKinMet stuff
    kTau,  ///< current time, s (kinetics)
    kdT,   ///< current time step, s (kinetics)

    TMols,      ///< Input total moles in b vector before rescaling
    SMols,      ///< Standart total moles (upscaled) {1000}
    MBX,        ///< Total mass of the system, kg
    FX,    	    ///< Current Gibbs potential of the system in IPM, moles
    IC,         ///< Effective molal ionic strength of aqueous electrolyte
    pH,         ///< pH of aqueous solution
    pe,         ///< pe of aqueous solution
    Eh,         ///< Eh of aqueous solution, V
    DHBM,       ///< balance (relative) precision criterion
    DSM,        ///< min amount of phase DS
    GWAT,       ///< used in ipm_gamma()
    YMET,       ///< reserved
    PCI,        ///< Current value of Dikin criterion of IPM convergence DK>=DX
    DXM,        ///< IPM convergence criterion threshold DX (1e-5)
    lnP,        ///< log Ptotal
    RT,         ///< RT: 8.31451*T (J/mole/K)
    FRT,        ///< F/RT, F - Faraday constant = 96485.309 C/mol
    Yw,         ///< Current number of moles of solvent in aqueous phase
    ln5551,     ///< ln(55.50837344)
    aqsTail,    ///< v_j asymmetry correction factor for aqueous species
    lowPosNum,  ///< Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
    logXw,      ///< work variable
    logYFk,     ///< work variable
    YFk,        ///< Current number of moles in a multicomponent phase
    FitVar[5];  ///< Internal. FitVar[0] is total mass (g) of solids in the system (sum over the BFC array)
    ///<      FitVar[1], [2] reserved
    ///<       FitVar[4] is the AG smoothing parameter;
    ///<       FitVar[3] is the actual smoothing coefficient
    double
    denW[5],   ///< Density of water, first T, second T, first P, second P derivative for Tc,Pc
    denWg[5],  ///< Density of steam for Tc,Pc
    epsW[5],   ///< Diel. constant of H2O(l)for Tc,Pc
    epsWg[5];  ///< Diel. constant of steam for Tc,Pc

    long int
    *L1,    ///< l_a vector - number of DCs included into each phase [Fi]
    // TSolMod stuff
    *LsMod, ///< Number of interaction parameters. Max parameter order (cols in IPx),
    ///< and number of coefficients per parameter in PMc table [3*FIs]
    *LsMdc, ///<  for multi-site models: [3*FIs] - number of nonid. params per component;
    /// number of sublattices nS; number of moieties nM
    *LsMdc2, ///<  new: [3*FIs] - number of DQF coeffs; reciprocal coeffs per end member;
    /// reserved
    *IPx,   ///< Collected indexation table for interaction parameters of non-ideal solutions
    ///< ->LsMod[k,0] x LsMod[k,1]   over FIs
    *LsPhl,  ///< new: Number of phase links; number of link parameters; [Fi][2]
    (*PhLin)[2];  ///< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over Fi)

    double
    // TSolMod stuff
    *PMc,    ///< Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
    *DMc,    ///< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
    *MoiSN,  ///< End member moiety- site multiplicity number tables ->  L1[k] x LsMdc[k][1] x LsMdc[k][2]
    *SitFr,  ///< Tables of sublattice site fractions for moieties -> LsMdc[k][1] x LsMdc[k][2]

    // Stoichiometry basis
    *A,   ///< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,    ///< IC atomic (molar) mass, g/mole [0:N-1]

    // Reconsider usage
    *Wb,     ///< Relative Born factors (HKF, reserved) [0:Ls-1]
    *Wabs,   ///< Absolute Born factors (HKF, reserved) [0:Ls-1]
    *Rion,   ///< Ionic or solvation radii, A (reserved) [0:Ls-1]

    *H0,     ///< DC pmolar enthalpies, reserved [L]
    *A0,     ///< DC molar Helmholtz energies, reserved [L]
    *U0,     ///< DC molar internal energies, reserved [L]
    *S0,     ///< DC molar entropies, reserved [L]
    *Cp0,    ///< DC molar heat capacity, reserved [L]

    *VL,      ///< ln mole fraction of end members in phases-solutions
    // TSolMod stuff
    *lPhc,  ///< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
    *DQFc  ///< new: Collected array of DQF parameters for DCs in phases -> L1[k] x LsMdc2[k][0]
    ;
    // until here move to --> datach.h

    // Contents defined in the enum below this structure
    // Other data
    double
    *XFs,    ///< Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,  ///< Current Karpov criteria of phase stability  F_a [0:FI-1]
    *DUL,     ///< VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,     ///< NG Vector of lower kinetic restrictions to x_j, moles [L]
    *fDQF,    ///< Increments to molar G0 values of DCs from pure gas fugacities or DQF terms, normalized [L]
    *YOF,     ///< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]
    *Vol,     ///< DC molar volumes, cm3/mol [L]
    *MM,      ///< DC molar masses, g/mol [L]
    *Pparc,   ///< Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
    *Y_m,     ///< Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    ///< log activity of DC in multi-component phases (mju-mji0) [0:L-1]
    *Y_w,     ///< Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
    *Gamma,   ///< DC activity coefficients in molal or other phase-specific scale [0:L-1]
    *lnGmf,   ///< ln of initial DC activity coefficients for correcting G0 [0:L-1]
    *lnGmM,   ///< ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]
    *EZ,      ///< Formula charge of DC in multi-component phases [0:Ls-1]
    *FVOL,    ///< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *FWGT,    ///< phase (carrier) masses, g                [0:FI-1]
    //
    *G,       ///< Normalized DC energy function c(j), mole/mole [0:L-1]            --> activities.h
    *G0,      ///< Input normalized g0_j(T,P) for DC at unified standard scale[L]   --> activities.h
    *lnGam,   ///< ln of DC activity coefficients in unified (mole-fraction) scale [0:L-1] --> activities.h
    *lnGmo;   ///< Copy of lnGam from previous IPM iteration (reserved)

    // TSolMod stuff (detailed output on partial energies of mixing)   --> activities.h
    double *lnDQFt; ///< new: DQF terms adding to overall activity coefficients [Ls_]
    double *lnRcpt; ///< new: reciprocal terms adding to overall activity coefficients [Ls_]
    double *lnExet; ///< new: excess energy terms adding to overall activity coefficients [Ls_]
    double *lnCnft; ///< new: configurational terms adding to overall activity [Ls_]
    double *CTerms;   ///< new: Coulombic terms (electrostatic activity coefficients) [Ls_]


    double  *B,  ///< Input bulk chem. compos. of the system - b vector, moles of IC[N]
    *U,  ///< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *U_r,  ///< IC chemical potentials u_i (J/mole) [0:N-1]
    *C,    ///< Calculated IC mass-balance deviations (moles) [0:N-1]
    *IC_m, ///< Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    *IC_lm,	///< log total IC molalities in aqueous phase [0:N-1]
    *IC_wm,	///< Total dissolved IC concentrations in g/kg_soln [0:N-1]
    *BF,    ///< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC,   ///< Total output bulk composition of all solid phases [1][N]
    *XF,    ///< Output total number of moles of phases Xa[0:FI-1]
    *YF,    ///< Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   ///< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   ///< Approximation of XFA in the next IPM iteration [0:FIs-1]
    *Falp;  ///< Karpov phase stability criteria F_a [0:FI-1] or phase stability index (PC==2)

    double (*VPh)[MIXPHPROPS],     ///< Volume properties for mixed phases [FIs]
    (*GPh)[MIXPHPROPS],     ///< Gibbs energy properties for mixed phases [FIs]
    (*HPh)[MIXPHPROPS],     ///< Enthalpy properties for mixed phases [FIs]
    (*SPh)[MIXPHPROPS],     ///< Entropy properties for mixed phases [FIs]
    (*CPh)[MIXPHPROPS],     ///< Heat capacity Cp properties for mixed phases [FIs]
    (*APh)[MIXPHPROPS],     ///< Helmholtz energy properties for mixed phases [FIs]
    (*UPh)[MIXPHPROPS];     ///< Internal energy properties for mixed phases [FIs]

    double *X,  ///< DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,   ///< Copy of x_j from previous IPM iteration [0:L-1]
    *XY,  ///< Copy of x_j from previous loop of Selekt2() [0:L-1]
    *Qp,  ///< Work variables related to non-ideal phases FIs*(QPSIZE=180)
    *Qd,  ///< Work variables related to DC in non-ideal phases FIs*(QDSIZE=60)
    *MU,  ///< mu_j values of differences between dual and primal DC chem.potentials [L]
    *EMU, ///< Exponents of DC increment to F_a criterion for phase [L]
    *NMU, ///< DC increments to F_a criterion for phase [L]
    *W,   ///< Weight multipliers for DC (incl restrictions) in IPM [L]
    *Fx,  ///< Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,  ///< Mole fractions Wx of DC in multi-component phases [L]
    *F,   /// <Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0;  ///< Excess Gibbs energies for (metastable) DC, mole/mole [L]
    // Name lists
    char (*sMod)[8];   ///< new: Codes for built-in mixing models of multicomponent phases [FIs]
    char  (*dcMod)[6];   ///< Codes for PT corrections for dependent component data [L]
    char  (*SB)[MAXICNAME+MAXSYMB]; ///< List of IC names in the system [N]
    char  (*SM)[MAXDCNAME];  ///< List of DC names in the system [L]
    char  (*SF)[MAXPHNAME+MAXSYMB];  ///< List of phase names in the system [FI]
    // Class codes
    char *RLC,   ///< Code of metastability constraints for DCs [L] enum DC_LIMITS
    *RSC,   ///< Units of metastability/kinetic constraints for DCs  [L]
    *RFLC,  ///< Classifier of restriction types for XF_a [FIs]
    *RFSC,  ///< Classifier of restriction scales for XF_a [FIs]
    *ICC,   ///< Classifier of IC { e o h a z v i <int> } [N]
    *DCC,   ///< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
    *PHC;   ///< Classifier of phases { a g f p m l x d h } [FI]
    char *DCCW;  ///< internal DC class codes [L]

    long int ITF,       ///< Number of completed IA EFD iterations
    ITG,         ///< Number of completed GEM IPM iterations
    ITau,    /// new: Time iteration for TKinMet class calculations
    IRes1;
    clock_t t_start, t_end;
    double t_elap_sec;  ///< work variables for determining IPM calculation time
    double *Guns;     ///<  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
    double *Vuns;     ///<  mu.L work vector of uncertainty space increments to tp->Vm
    double *tpp_G;    ///< Partial molar(molal) Gibbs energy g(TP) (always), J/mole
    double *tpp_S;    ///< Partial molar(molal) entropy s(TP), J/mole/K
    double *tpp_Vm;   ///< Partial molar(molal) volume Vm(TP) (always), J/bar

    // additional arrays for internal calculation in ipm_main
    double *XU;      ///< dual-thermo calculation of DC amount X(j) from A matrix and u vector [L]
    double (*Uc)[2]; ///< Internal copy of IC chemical potentials u_i (mole/mole) at r-1 and r-2 [N][2]
    double *Uefd;    ///< Internal copy of IC chemical potentials u_i (mole/mole) - EFD function [N]
    char errorCode[100]; ///<  code of error in IPM      (Ec number of error)
    char errorBuf[1024]; ///< description of error in IPM
    double logCDvalues[5]; ///< Collection of lg Dikin crit. values for the new smoothing equation
    double *GamFs;   ///< Copy of activity coefficients Gamma before the first enter in PhaseSelection() [L] new

    double // Iterators for MTP interpolation (do not load/unload for IPM)
    Pai[4],    ///< Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
    Tai[4],    ///< Temperature T, C: start, end, increment for MTP array in DataCH , Ttol
    Fdev1[2],  ///< Function1 and target deviations for  minimization of thermodynamic potentials
    Fdev2[2];  ///< Function2 and target deviations for  minimization of thermodynamic potentials

    // Experimental: modified cutoff and insertion values (DK 28.04.2010)
    double
    // cutoffs (rescaled to system size)
    XwMinM, ///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
    ScMinM, ///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
    DcMinM, ///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
    PhMinM, ///< Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
    ///< insertion values (re-scaled to system size)
    DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
    DFYaqM, ///< Insertion mole amount for aqueous and surface species { 1e-6 }
    DFYidM, ///< Insertion mole amount for ideal solution components { 1e-6 }
    DFYrM,  ///< Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
    DFYhM,  ///< Insertion mole amount for minor solution components { 1e-6 }
    DFYcM,  ///< Insertion mole amount for single-component phase { 1e-6 }
    DFYsM,  ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
    SizeFactor; ///< factor for re-scaling the cutoffs/insertions to the system size
}
TSOLMOD_MULTI;

extern const BASE_PARAM pa_p_;

// Data of MULTI
class TSolModMulti
{
    /// Internal TSolMod decorator
    std::vector<SolModCalc> phase_models;
    /// Names of multi-component phases
    std::vector<std::string> phase_names;
    void InitalizeTSolMod();

    // Stuff from TNode
    DATACH* CSD;  ///< Pointer to chemical system data structure CSD (DATACH)
    DATABR* CNode;  ///< Pointer to a work node data bridge structure (node)
    bool load_thermodynamic_data = false;
    std::string current_output_set_name;
    std::string current_input_set_name;
    std::string ipmlog_error;

#ifdef USE_THERMOFUN
    std::unique_ptr<ThermoFun::ThermoEngine> thermo_engine;
#endif
    std::string thermo_json_string;

    /// Clear thermodynamic data from ThermoEngine
    void clear_ThermoEngine();
    /// Read ThermoEngine
    bool load_ThermoEngine(const std::string& thermo_file_or_string);
    void allocMemory();
    void freeMemory();
    void clearipmLogError() {
        ipmlog_error.clear();
    }

public:

    explicit TSolModMulti();
    virtual ~TSolModMulti();

    /// Initialization of GEM IPM2 data structures in coupled RMT-GEM programs
    ///  that use GEMS3K module. Also reads in the IPM, DCH and DBR text input files
    ///  in key-value, json or binary format. Parameters:
    ///  ipmfiles_lst_name - name of a text file that contains:
    ///    " -j | -t |-b <DCH_DAT file name> <IPM_DAT file name> [<>] <dataBR file name>
    ///    or " -f <DCH_DAT file name> <IPM_DAT file name> <ThermoFun JSON format file> <dataBR file name>
    ///  dbfiles_lst_name - name of a text file that contains:
    ///    <dataBR  file name1>, ... , <dataBR file nameN> "
    ///    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
    ///    exist in the same directory where the ipmfiles_lst_name file is located.
    ///    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
    ///    and must contain valid initial chemical systems (of the same structure
    ///    as described in the DCH_DAT file) to set up the initial state of the FMT
    ///    node array.
    ///  If -t flag or no flag is specified then all data files must be in key-value text
    ///    (ASCII) format (and file names must have .dat extension);
    ///  If -j and -f flag is specified then all data files must be in JSON format (and file names
    ///    must have .json extension);
    ///  If -f flag is specified then the use of ThermoFun along with GEMS3K in place of the interpolation
    ///  of lookup arrays for standard thermodynamic data for substances;
    ///  if -b flag is specified then all data files are assumed to be binary (little-endian)
    ///    files.
    ///  @returns:
    ///    0: OK;
    ///    1: GEM IPM read file or input file format error;
    ///   -1: System error (e.g. memory allocation).
    long int GEM_init(const std::string& ipmfiles_lst_name);

    /// Initialization of GEM IPM3 data structures in coupled programs that use GEMS3K module.
    /// Also reads the input data from the IPM, DCH and one DBR JSON input strings
    /// (e.g. exported from GEM-Selektor or retrieved from JSON database).
    /// Parameters:
    ///  @param dch_json -  DATACH - the Data for CHemistry data structure as a json/key-value string
    ///  @param ipm_json -  Parameters and settings for GEMS3K IPM-3 algorithm as a json/key-value string
    ///  @param dbr_json -  DATABR - the data bridge structure as a json/key-value string
    ///  @param fun_json -  ThermoFun data structure as a json string
    ///  @returns:
    ///    0: OK;
    ///    1: GEM IPM read file or input file format error;
    ///   -1: System error (e.g. memory allocation).
    long GEM_init(const std::string& dch_json, const std::string& ipm_json,
                  const std::string& dbr_json, const std::string& fun_json);

    /// Update thermodynamic data according new TK and P
    void UpdateThermodynamic(double TK, double PPa);

    /// Access to idxs phase model
    SolModCalc &get_phase(std::size_t idx);
    /// Access to phase model by phase name
    SolModCalc &get_phase(const std::string& name);

    /// Trace output full internal structure
    void to_text_file(const char *path, bool append=false);

protected:

    // Internal MULTI subset (could be remove unused in future)
    char PAalp_; ///< Flag for using (+) or ignoring (-) specific surface areas of phases
    char PSigm_; ///< Flag for using (+) or ignoring (-) specific surface free energies
    std::shared_ptr<BASE_PARAM> pa_standalone;
    TSOLMOD_MULTI pm;
    TSOLMOD_MULTI *pmp;

    BASE_PARAM* base_param() const
    {
        return pa_standalone.get();
    }

    // reading multi structure
    template<typename TIO>
    void to_text_file_gemipm( TIO& out_format, bool addMui,
                              bool with_comments = true, bool brief_mode = false );
    template<typename TIO>
    void from_text_file_gemipm( TIO& in_format,  DATACH  *dCH );

    /// Writes Multi to a json/key-value string
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    std::string gemipm_to_string( bool addMui, const std::string& test_set_name, bool with_comments = true, bool brief_mode = false );
    /// Reads Multi structure from a json/key-value string
    bool gemipm_from_string( const std::string& data,  DATACH  *dCH, const std::string& test_set_name );

    ///  Reads the contents of the work instance of the DATABR structure from a stream.
    ///   \param stream    string or file stream.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    void  read_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f, DATACH  *dCH, const std::string& test_set_name );

    /// Writes the contents of the work instance of the DATABR structure into a stream.
    ///   \param stream    string or file stream.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    ///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
    ///                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
    ///                 If   false (0), comments will not be written;
    ///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
    ///                                   with_comments=off as "condensed JSON"
    ///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
    void  write_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f,
                                   bool addMui, bool with_comments, bool brief_mode, const std::string& test_set_name );

    void multi_realloc(char PAalp, char PSigm);
    void multi_kill();
    void set_def(int i=0);

    // New functions for TSolMod parameter arrays
    void getLsModsum( long int& LsModSum, long int& LsIPxSum );
    void getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum );
    /// Get dimensions from LsPhl array
    void getLsPhlsum( long int& PhLinSum,long int& lPhcSum );
    /// Get dimensions from LsMdc2 array
    void getLsMdc2sum( long int& DQFcSum,long int& rcpcSum );

    void get_PAalp_PSigm(char &PAalp, char &PSigm);
    void alloc_IPx( long int LsIPxSum );
    void alloc_PMc( long int LsModSum );
    void alloc_DMc( long int LsMdcSum );
    void alloc_MoiSN( long int LsMsnSum );
    void alloc_SitFr( long int LsSitSum );
    void alloc_DQFc( long int DQFcSum );
    void alloc_PhLin( long int PhLinSum );
    void alloc_lPhc( long int lPhcSum );

    // Fill multi arrays
    void InitalizeGEM_IPM_Data();
    void MultiConstInit();
    long getXvolume();
    void load_all_thermodynamic_from_grid(double TK, double PPa);
    bool load_all_thermodynamic_from_thermo(double TK, double PPa);
    void LoadThermodynamicData(double TK, double PPa);
    void ConvertDCC();
    double ConvertGj_toUniformStandardState(double g0, long j, long k);
};

// syp->PGmax
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

typedef enum {  // Field index into outField structure
    f_pa_PE = 0,  f_PV,  f_PSOL,  f_PAalp,  f_PSigm,
    f_Lads,  f_FIa,  f_FIat
} MULTI_STATIC_FIELDS;

typedef enum {  // Field index into outField structure
    f_sMod = 0,  f_LsMod,  f_LsMdc,  f_B,  f_DCCW,
    f_Pparc,  f_fDQF,  f_lnGmf,  f_RLC,  f_RSC,
    f_DLL,  f_DUL,  f_Aalp,  f_Sigw,  f_Sigg,
    f_YOF,  f_Nfsp,  f_MASDT,  f_C1,  f_C2,
    f_C3,  f_pCh,  f_SATX,  f_MASDJ,  f_SCM,
    f_SACT,  f_DCads,
    // static
    f_pa_DB,  f_pa_DHB,  f_pa_EPS,  f_pa_DK,  f_pa_DF,
    f_pa_DP,  f_pa_IIM,  f_pa_PD,  f_pa_PRD,  f_pa_AG,
    f_pa_DGC,  f_pa_PSM,  f_pa_GAR,  f_pa_GAH,  f_pa_DS,
    f_pa_XwMin,  f_pa_ScMin,  f_pa_DcMin,  f_pa_PhMin,  f_pa_ICmin,
    f_pa_PC,  f_pa_DFM,  f_pa_DFYw,  f_pa_DFYaq,  f_pa_DFYid,
    f_pa_DFYr,  f_pa_DFYh,  f_pa_DFYc,  f_pa_DFYs,  f_pa_DW,
    f_pa_DT,  f_pa_GAS,  f_pa_DG,  f_pa_DNS,  f_pa_IEPS,
    f_pKin,  f_pa_DKIN,  f_mui,  f_muk,  f_muj,
    f_pa_PLLG,  f_tMin,  f_dcMod,
    //new
    f_kMod, f_LsKin, f_LsUpt, f_xICuC, f_PfFact,
    f_LsESmo, f_LsISmo, f_SorMc, f_LsMdc2, f_LsPhl

} MULTI_DYNAMIC_FIELDS;


enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

#endif   //tsolmod_multi_h

