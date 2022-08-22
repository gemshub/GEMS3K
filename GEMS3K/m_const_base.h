//-------------------------------------------------------------------
// $Id$
//
/// \file m_const_base.h
/// Declarations of enums and constants from GEM-Selektor code that
/// are used in GEMS3K standalone code (this file is not used otherwise).
//
// Copyright (c) 1995-2012 S.Dmytriyeva, D.Kulik
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

#ifndef M_CONST_BASE_H
#define M_CONST_BASE_H

const char S_OFF = '-',
           S_ON = '+',
           S_REM = '*',
           A_NUL ='?';

const long int MST =  6,
               DFCN = 6; // number of columns in MASDJ table

const unsigned long int
            MAXICNAME =      6,
            MAXSYMB =        4,
            TDBVERSION =     64;

const long int
           MAXDCNAME =      16,
           MAXPHNAME =      16,
           EQ_RKLEN = 58;

const int 	MPP_TOT = 0,       // index of column with total mixed phase property
            MPP_STD = 1,       // index of column with standard property sum for mixed phases
            MPP_RCP = 2,       // index of column with reciprocal (Darken) property sum for mixed phases
            MPP_IDL = 3,       // index of column with ideal mixing property for the phases
            MPP_EXS = 4,       // index of column with excess mixing property for the phases
            MIXPHPROPS = 5;    // Number of columns in the property table for mixed phases


typedef enum {  // Limits on DC and phases
    // type of lmits
    NO_LIM = 'O', LOWER_LIM ='L', UPPER_LIM = 'U', BOTH_LIM ='B',
    // mode recalc of limits Set_DC_Limits()
    DC_LIM_INIT = 0, DC_LIM_CURRENT = 1
} DC_LIMITS;

// Work DC classifier codes  pm->DCCW
enum SolDCodes {

    DC_SINGLE = 'U',        // This DC is a single-component phase
    DC_SYMMETRIC = 'I',     // This DC is in symmetric solution phase
    DC_ASYM_SPECIES = 'S',  // This is DC-solute(sorbate) in asymmetric phase
    DC_ASYM_CARRIER = 'W'   // This is carrier(solvent) DC in asymmetric phase
};

enum solmod_switches { // indexes of keys of phase (solution, sorption, kinetic) models
    SPHAS_TYP,
    DCOMP_DEP,
    SPHAS_DEP,
    SGM_MODE,
    DCE_LINK,
    MIX_TYP,

    // Link state of CalculateActivityCoefficients()
    LINK_UX_MODE,
    LINK_TP_MODE,
    LINK_PP_MODE,  // LINK_PHP_MODE,
    LINK_IN_MODE, // Initialization mode for kinetics and other time-dependent processes
    SORP_MOD,  // new, see also enum sorption_control
    KINR_MOD,  /// see also enum kinmet_controls
    // Posible modes of calculation of activity coefficients (private, public)
    SM_UNDEF = 'N',
    SM_TPDEP = 'T',
    SM_UXDEP = 'X',
    SM_PRIVATE_ = 'P',
    SM_PUBLIC = 'U',

    // Posible modes of calculation of activity coefficients (built-in or scripted models)
    SM_STNGAM = 'S',        // Built-in function for activity coefficients
    SM_NOSTGAM = 'N',

    // Codes to identify the mixing models used (during IPM iterations)
    SM_IDEAL =  'I',	// ideal solution or single-component phase
    SM_BERMAN = 'B',    // built-in multicomponent multisite (a)symmetric solid-solution model
    SM_CEF = '$',    //     built-in multicomponent multisite solid-solution model (CALPHAD)
    SM_MBW = '#',    //     built-in Modified Bragg-Williams model, [Vinograd et al. 2018]
    SM_REDKIS = 'G', 	// built-in binary Guggenheim (Redlich-Kister) solid-solution model
    SM_MARGB = 'M',	// built-in binary Margules solid-solutions (subregular)
    SM_MARGT = 'T',	// built-in ternary Margules solid-solution (regular)
    SM_VANLAAR = 'V',	// built-in multi-component Van Laar solid-solution model
    SM_GUGGENM = 'K',	// built-in multi-component Guggenheim solid-solution model
    SM_REGULAR = 'R',	// built-in multi-component Regular solid-solution model
    SM_NRTLLIQ = 'L',	// built-in multi-component NRTL model for liquid solutions
    SM_WILSLIQ = 'W',	// built-in multi-component Wilson model for liquid solutions
    SM_CGFLUID = 'F',	// built-in multi-component Churakov-Gottschalk (CG) fluid EoS model
    SM_PRFLUID = 'P',	// built-in Peng-Robinson-Stryjek-Vera (PRSV) fluid EoS model
    SM_PCFLUID = '5',   // built-in perturbed-chain statistical-association (PCSAFT) fluid EoS model (reserved)
    SM_STFLUID = '6',   // built-in Sterner-Pitzer (STP) fluid EoS model
    SM_PR78FL = '7',	// built-in Peng-Robinson (PR78) fluid EoS model
    SM_CORKFL = '8',    // built-in compensated Redlich-Kwong (CORK) fluid EoS model
    SM_REFLUID = '9',   // built-in reference EoS fluid model (reserved)
    SM_SRFLUID = 'E',	// built-in Soave-Redlich-Kwong (SRK) fluid EoS model
    SM_AQDAV = 'D',	// built-in Davies model (with 0.3) for aqueous electrolytes
    SM_AQDH1 = '1',	// built-in Debye-Hueckel limiting law for aqueous electrolytes
    SM_AQDH2 = '2',	// built-in 2-term Debye-Hueckel model for aqueous electrolytes
    SM_AQDH3 = '3',	// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Karpov version)
    SM_AQDHH = 'H',	// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Helgeson version)
    SM_AQDHS = 'Y',	// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Shvarov version)
    SM_AQSIT = 'S',	// built-in SIT model for aqueous electrolytes
    SM_AQEXUQ = 'Q',    // built-in extended UNIQUAC model for aqueous electrolytes
    SM_AQPITZ = 'Z',    // built-in Pitzer HMW model for aqueous electrolytes
    SM_AQMIX = 'C',     // built-in mixed-solvent aqueous Debye-Hueckel model (reserved)
    SM_AQELVIS = 'J',   // built-in modified extended UNIQUAC model (ELVIS) for aqueous electrolyte
    SM_DONNAN = 'X',    // ion exchange (Donnan volume model) (reserved)
    SM_SURCOM = 'A',	// models of surface complexation at solid-aqueous interface
    SM_USERDEF = 'U',	// user-defined mixing model (scripts in Phase record)
    SM_OTHER = 'O',	// other built-in phase-specific models of non-ideal solutions (selected by phase name)

    // Codes to identify specific mixing rules and temperature functions in EoS and activity models
    MR_UNDEF = 'N', // Default mixing rule or form of interaction parameter coefficients; NEM for adsorption 'A'
    MR_WAAL = 'W',	// Basic Van der Waals mixing rules in cubic EoS models
    MR_CONST = 'C',	// Constant one-term interaction parameter kij; CCM for sorption 'A'
    MR_TEMP = 'T',	// Temperature-dependent one-term interaction parameter kij (Jaubert et al. 2005); TLM for adsorption 'A'
    MR_LJ = 'J',        // Lemmon-Jacobsen mixing rule (Lemmon and Jacobsen, 1999)
    MR_KW1 = 'K',       // Kunz-Wagner mixing rule (Kunz and Wagner, 2007)
    MR_PITZ5 = '5',     // 5-term Pitzer model temperature dependence (TOUGHREACT variant)
    MR_PITZ6 = '6',     // 6-term Pitzer model temperature dependence (FREZCHEM variant)
    MR_PITZ8 = '8',     // 8-term Pitzer model temperature dependence
    MR_B_RCPT = 'R',    // Use CEF reciprocal non-ideality terms in Berman multi-site ss model
    MR_A_DLM  = 'D',    // Diffuse-layer electrostatic model (DLM) for adsorption 'A'
    MR_A_BSM  = 'B',    // Basic Stern electrostatic model (BSM) for adsorption 'A'
    MR_A_CDLM  = 'M',    // CD-MUSIC (3-layer) electrostatic model (DLM) for adsorption 'A'
    MR_A_ETLM  = 'E'     // Extended TLM electrostatic model (ETLM) for adsorption 'A'
};


typedef enum {  // classes of independent components IC, used in ccIC code list
    IC_ELEMENT  =  'e',  // chemical element (except oxygen and hydrogen)
    IC_OXYGEN   =  'o',  // oxygen
    IC_HYDROGEN =  'h',  // hydrogen (natural mixture of isotopes) H
    IC_PROTIUM   = 'p',  // protium (reserved) Hp
    IC_DEYTERIUM = 'd',  // deuterium (reserved) D
    IC_TRITIUM  =  't',  // tritium (reserved) T
    IC_FORMULA  =  'f',  // formula unit (eg. for Sio - a symbol of SiO2)
    IC_METALION =  'm',  // metal ion (cation), reserved
    IC_LIGAND   =  'l',  // ligand (anion), reserved
    IC_ADDIT    =  'a',  // IC with unknown stoichiometry (eg: Hum, humic acid)
    IC_ISOTOPE  =  'i',  // isotope of chemical element from 1 to 250
    IC_OXYGEN16 =  'q',  // q: oxygen 16O (reserved)
    IC_OXYGEN18 =  'r',  // r: oxygen 18O (reserved)
    IC_CHARGE   =  'z',  // z: electrical charge
    IC_VOLUME   =  'v',  // volume
    IC_SITE     =  's'   // sorption site for site balance constraint (reserved)
} IC_CLASSES;


typedef enum {  // Classifications of DC
    // Type of input data for
    SRC_DCOMP = 'd',        // the key points to existing PDB record in DCOMP chain
    SRC_REACDC = 'r',       // the key points to existing PDB record in REACDC chain
    SRC_NEWDC = 'n',        // the key new reaction-defined component
    SRC_NEWISO = 'i',       // the same as n, but this component is an isotopic form
    SRC_FICT = 'f',         // fictive species
    // Aqueous electrolyte phase
    DC_AQ_PROTON = 'T',     // hydrogen ion H+
    DC_AQ_ELECTRON = 'E',   // electron (as a DC)
    DC_AQ_SPECIES  = 'S',   // other aqueous species (ions, complexes and ion pairs)
    DC_AQ_SURCOMP = 'K',    // Surface complex represented as aqueous species
    DC_AQ_SOLVENT  = 'W',   // water H2O (major solvent)
    DC_AQ_SOLVCOM  = 'L',   // other components of a solvent (eg. alcohol)
    // Gas phase ( G code can be used for all gases; V,C,H,N codes are reserved)
    DC_GAS_COMP  = 'G',     // other gases
    DC_GAS_H2O = 'V',       // H2O steam
    DC_GAS_CO2 = 'C',       // CO2 (carbon dioxide)
    DC_GAS_H2 = 'H',        // H2 hydrogen
    DC_GAS_N2 = 'N',        // N2 nitrogen
    // Solid/liquid non-electrolyte multicomponent phases
    DC_SOL_IDEAL = 'I',     // end-member of ideal solution
    DC_SOL_MINOR = 'J',     // junior independent end member (for initial approximation)
    DC_SOL_MAJOR = 'M',     // major independent end member (for initial approximation)
    DC_SOL_MINDEP = 'F',    // junior dependent end member (for initial approximation)
    DC_SOL_MAJDEP = 'D',    // major dependent end member (for initial approximation)
    // Sorption phases and poly(oligo)electrolytes
    DC_SUR_CARRIER = 'Q',   // Principal end-member of solid carrier
    DC_SUR_MINAL = 'P',     // Minor end-member of solid carrier
    DC_PEL_CARRIER = 'R',   // Carrier of poly(oligo)electrolyte
    DC_SSC_A0 = '0',        // Strong surface complex on site type 0 (A plane)
    DC_SSC_A1 = '2',        // Strong surface complex on site type 2 (A plane)
    DC_SSC_A2 = '4',        // Strong surface complex on site type 4 (A plane)
    DC_SSC_A3 = '6',        // Strong surface complex on site type 6 (A plane)
    DC_SSC_A4 = '8',        // Strong surface complex on site type 8 (A plane)
    DC_WSC_A0 = '1',        // Weak surface complex on site type 1 (B plane)
    DC_WSC_A1 = '3',        // Weak surface complex on site type 3 (B plane)
    DC_WSC_A2 = '5',        // Weak surface complex on site type 5 (B plane)
    DC_WSC_A3 = '7',        // Weak surface complex on site type 7 (B plane)
    DC_WSC_A4 = '9',        // Weak surface complex on site type 9 (B plane)
    DC_IESC_A  = 'A',       // Strong exchange ion const-charge plane
    DC_IEWC_B  = 'B',       // Weak exchange ion const-charge plane

    // Aliaces for 1-site model
    DC_SUR_GROUP = 'X',      // Surface site A plane -> '0'
    DC_SUR_COMPLEX = 'Y',    // Strong sur. complex A plane -> '0'
    DC_SUR_IPAIR = 'Z',      // Weak sur complex B plane -> '1'

    // Single-component phases
    DC_SCP_CONDEN = 'O',     // DC forming a single-component phase

    // New surface complexation models (added 16.11.2017 by DK)
    DC_SCM_SPECIES = 'U'

} DC_CLASSES;

#define SORPTION_DC "0246813579ABXYZPQ"

// This code defines standard state and reference scale of concentrations
// for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES {  // Possible values
    PH_AQUEL = 'a',  	// aqueous electrolyte
    PH_GASMIX = 'g',  	// mixture of gases
    PH_FLUID = 'f',  	// fluid phase
    PH_PLASMA = 'p',  	// plasma
    PH_LIQUID = 'l',  	// non-electrolyte liquid (melt)
    PH_SIMELT = 'm',  	// silicate (magmatic) melt or non-aqueous electrolyte
    PH_SORPTION = 'x',  // Sorption phase (sorbent+sorbates) in aqueous electrolyte
    PH_POLYEL = 'y',  	// colloidal poly- (oligo)electrolyte e.g. Donnan volume phase
    PH_IONEX = 'i',     // ion exchange on permanent charge ligand e.g. B&B Clay
    PH_ADSORPT = 'z',   // surface complexation (adsorption) on hydrated amphoteric surface
    PH_SINCOND = 's',  	// condenced solid phase, also multicomponent (solid solution)
    PH_SINDIS = 'd',  	// dispersed solid phase, also multicomponent
    PH_HCARBL = 'h'   	// mixture of condensed hydrocarbons
};

enum sorption_control {
    // EDL interface models - separate for site types in v. 3.1
    SC_DDLM = 'D',  SC_CCM = 'C',     SC_TLM = 'T',   SC_MTL = 'M',
    SC_MXC = 'E',   SC_NNE = 'X',     SC_IEV  = 'I',  SC_BSM = 'S',
    SC_3LM = '3', SC_NOT_USED = 'N',

    // Methods of Surface Activity Terms calculation
    SAT_COMP = 'C', SAT_NCOMP = 'N', SAT_SOLV = 'S', SAT_INDEF = 'I',

    // New methods for surface activity coefficient terms (2004)
    SAT_L_COMP = 'L', SAT_QCA_NCOMP = 'Q', SAT_QCA1_NCOMP = '1',
    SAT_QCA2_NCOMP = '2', SAT_QCA3_NCOMP = '3', SAT_FREU_NCOMP = 'f',
    SAT_QCA4_NCOMP = '4', SAT_BET_NCOMP = 'B', SAT_VIR_NCOMP = 'W',
    SAT_FRUM_NCOMP = 'F', SAT_FRUM_COMP = 'R', SAT_PIVO_NCOMP = 'P',

    // Assignment of surtype to carrier (end-member)
    CCA_VOL = 'V',
    CCA_0 = '0', CCA_1, CCA_2, CCA_3, CCA_4, CCA_5,
    CCA_6, CCA_7, CCA_8, CCA_9, SPL_0='0', SPL_1, SPL_2, SPL_3,
    SPL_B = 'b', SPL_D = 'd', SPL_C = 'c',
    SDU_N = 'n', SDU_m = 'm', SDU_M = 'M', SDU_g = 'g',
    CST_0 = '0', CST_1, CST_2, CST_3, CST_4, CST_5, // surface type index
    CSI_0 = '0', CSI_1, CSI_2, CSI_3, CSI_4, CSI_5, // surface site index

    // Number of parameters per surface species in the MaSdj array
    // MCAS = 6 = DFCN ; position index    added by KD 25.10.2004
    // Column indices of surface species allocation table MCAS
    SA_MCA=0, SA_EMX, SA_STX, SA_PLAX, SA_SITX, SA_UNIT,

    // Column indices of MaSdj table of adsorption parameters
    PI_DEN=0, PI_CD0, PI_CDB, PI_P1, PI_P2, PI_P3
};

enum ph_kinmet_controls { /// TKinMet: codes to control kinetic rate models

    KM_UNDEF = 'N',       /// not defined, no account for
    KinProCode = 2,
    KM_PRO_MWR = 'M',     /// Kinetics of generic dissolution/precipitation (no uptake, ionex, adsorption)
    KM_PRO_UPT = 'U',     /// Kinetics of uptake/entrapment (of minor/trace element) into solid solution
    KM_PRO_IEX = 'X',     /// Kinetics of ion exchange (clays, C-S-H, zeolites, ...)
    KM_PRO_ADS = 'A',     /// Kinetics of adsorption (on MWI), redox
    KM_PRO_NUPR = 'P',    /// Kinetics of nucleation and precipitation
    KinModCode = 3,
    KM_MOD_TST = 'T',     /// Generic TST dissolution/precipitation model following Shott ea 2012
    KM_MOD_PAL = 'P',     /// Dissolution/precipitation model of the form (Palandri 2004)
    KM_MOD_WOL = 'W',     /// Carbonate growth model following (Wolthers 2012)
    KM_MOD_NUGR = 'U',    /// Mineral nucleation and growth model with nuclei/particle size distr. (TBD)
    KinSorpCode = 4,
    KM_UPT_ENTRAP = 'E',  ///	Unified entrapment model (Thien,Kulik,Curti 2012)
    KM_UPT_UPDP = 'M',    ///	DePaolo (2011) uptake kinetics model
    KM_UPT_SEMO = 'G',    ///  Growth (surface) entrapment model (Watson 2004)
    KM_IEX_FAST = 'F',    ///  Fast ion exchange kinetics (e.g. montmorillonite, CSH)
    KM_IEX_SLOW = 'L',    ///  Slow ion exchange kinetics (e.g. illite, zeolites)
    KM_ADS_INHIB = 'I',   ///  Adsorption inhibition
    KM_NUCL_SSMP  = 'P',  ///  Solid solution nucleation model (Prieto 2013)
    KinLinkCode = 5,
    KM_LNK_SURF = 'S',    ///   Link to (fraction of) solid substrate surface
    KM_LNK_PVOL = 'P',    ///    Link to (fraction of) solid substrate (pore) volume
    KM_LNK_MASS = 'M',    ///	Link to (fraction of) solid substrate mass
    KinSizedCode = 6,
    KM_SIZED_ETM = 'T',   ///  Empirical f(time) cubic polynomial f = a + bt +ct^2 + dt^3 (default)
    KM_SIZED_ESI = 'S',   ///  Empirical f(lgSI) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_ESA = 'A',   ///  Empirical f(sarea-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_EVOL = 'V',  ///  Empirical f(volume-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_MASS = 'M',  ///  Empirical f(mass-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_MOL = 'X',   ///  Empirical f(amount-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_UNI = 'U',   /// 	Uniform particle/pore size distribution
    KM_SIZED_BIN = 'B',   /// 	Binodal particle/pore size distribution
    KM_SIZED_FUN = 'F',   ///    Empirical distribution function
    KinResCode = 7,
    KM_RES_SURF_N = 'A',   /// surface-scaled rate constant (k in mol/m2/s), default
    KM_RES_SURF_M = 'M',   /// surface-scaled rate constant (k in kg/m2/s)
    KM_RES_PVS_N  = 'V',   /// pore-volume-scaled rate constant (k in mol/m3/s)
    KM_RES_PVS_M  = 'W',   /// pore-volume-scaled rate constant (k in kg/m3/s)
    KM_RES_ABS_N  = 'F',   /// absolute (unscaled) rate constant (k in mol/s)
    KM_RES_ABS_M  = 'G',   /// absolute (unscaled) rate constant (k in kg/s)
    KM_LIN_RATE   = 'L'   /// linear growth/dissolution rate constant (v in m/s)

};

typedef enum { // Units of measurement of quantities and concentrations
    // number of components and phases
    QUAN_MKMOL = 'Y',  QUAN_MMOL = 'h',  QUAN_MOL = 'M',  // NUMBER OF MOLES
    QUAN_MGRAM = 'y',  QUAN_GRAM = 'g',  QUAN_KILO = 'G', // MASS
    // concentrations of components and phases
    CON_MOLFR = 'n', CON_MOLPROC = 'N', CON_pMOLFR = 'f', // MOLE FRACTION
    CON_VOLFR = 'v', CON_VOLPROC = 'V', CON_pVOLFR = 'u', // VOLUME FRACTION
    CON_MOLAL = 'm', CON_MMOLAL =  'i', CON_pMOLAL = 'p', // MOLALITY
    CON_MOLAR = 'L', CON_MMOLAR =  'j', CON_pMOLAR = 'q', // MOLARITY
    CON_WTFR  = 'w', CON_WTPROC =  '%', CON_PPM =    'P', // MASS FRACTION
    CON_AQWFR = 'C', CON_AQWPROC = 'c', CON_AQPPM =  'a', // CONCENTRATION
    // aqueous species
    CON_AQGPL = 'd', CON_AQMGPL = 'e', CON_AQMKGPL = 'b', //VOLUME CONCENTRATION

    //Units of measurement of pressure Pr, P  { b B p P A }'
    PVT_BAR =  'b', // bar - default, 1 bar = 1e5 Pa
    PVT_KBAR = 'B', // kbar, 1 kbar = 1000 bar
    PVT_PASC = 'p', // Pascal (Pa)
    PVT_KPASC = 'P',// MPa, 1 MPa = 10 bar = 1e6 Pa
    PVT_ATM =  'A', // atm, 1 atm = 1.013 bar
    //Attention: Only b code can be used in this version!

    //Units of measurement of molar volume  { c j a L m }'
    PVT_CM3 =  'c',  // cm3, cm3/mole
    PVT_LITR =  'L', // liters (L) - volume of the system only, 1 L = 1000 cm3
    PVT_JBAR =  'j', // J/bar, 10 cm3/mole = 1 J/bar
    PVT_CBAR = 'a',  // (cal/bar), 41.84 cm3/mole = 1 cal/bar
    // m  - reserved.
    //Attention: only j code can be used in this version!

    //Units of measurement of reference temperature Tr { C K F }'
    PVT_CELS = 'C',   // degrees Celsius (C)
    PVT_KELVIN = 'K', // Kelvins (K), 0 C = 273.15 K
    PVT_FAREN = 'F',  // degrees Fahrenheit (F)
    //Attention: Only C code can be used in this version of GEM IPM algorithm.

    // Units of measurement of energy values { j c J C n N }
    TDAT_JOUL = 'j',  // Joules (J/mole)
    TDAT_KJOUL = 'J', // kilojoules (kJ/mole)
    TDAT_CAL = 'c',   // calories (cal/mole); 1 cal = 4.184 J;
    TDAT_KCAL = 'C',  // kilocalories (kcal/mole)
    TDAT_NORM = 'N'   // normalized (mole/mole, J/mole/RT, cal/mole/RT)
    // Attention: Only j code can be used in this version!
} SPPA_UNITS;

#endif  // M_CONST_BASE_H


