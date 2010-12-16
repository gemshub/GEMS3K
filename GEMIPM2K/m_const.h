//-------------------------------------------------------------------
// $Id: m_const.h 1431 2009-08-28 16:28:04Z gems $
//
// Copyright (C) 2006,2009  S.Dmitrieva, D.Kulik
//
// Codes and parameters used in GEM IPM work structure (standalone version)
//
// This file is part of the standalone GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//
#ifndef _m_const_h
#define _m_const_h

#include <ctype.h>
#include <fstream>

using namespace std;
#include "verror.h"
#include "v_user.h"

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
	MPP_RES = 2,       // index of column with residual property sum for mixed phases
	MPP_ID = 3,        // index of column with ideal mixing property for the phases
	MPP_EX = 4,        // index of column with excess mixing property for the phases
	MIXPHPROPS = 5;    // Number of columns in the property table for mixed phases

enum solmod_switches { // indexes of keys of model solution
    SPHAS_TYP,
    DCOMP_DEP,
    SPHAS_DEP,
    SGM_MODE,
    DCE_LINK,
    MIX_TYP,
    // link state
    LINK_UX_MODE,
    LINK_TP_MODE,
    LINK_FIA_MODE,
    LINK_PHP_MODE,
    // Posible values of �of keys of model solution - DCOMP_DEP, SPHAS_DEP
    SM_UNDEF = 'N',
    SM_TPDEP = 'T',
    SM_UXDEP = 'X',
    SM_PRIVATE_ = 'P',
    SM_PUBLIC = 'U',
    // Posible modes calculating of activity coefficients SGM_MODE
    SM_STNGAM = 'S',
    SM_NOSTGAM = 'N',
    //  Possible values: (SPHAS_TYP)
    // Code to identify the mixing models used (during IPM iterations)
    SM_IDEAL =  'I',	// ideal solution or single-component phase;
    SM_REDKIS = 'G', 	// built-in binary Guggenheim (Redlich-Kister) solid-solution model
    SM_MARGB = 'M',		// built-in binary Margules solid-solutions (subregular)
    SM_MARGT = 'T',		// built-in ternary Margules solid-solution (regular)
    SM_VANLAAR = 'V',	// built-in multi-component Van Laar solid-solution model
    SM_GUGGENM = 'K',	// built-in multi-component Guggenheim solid-solution model
    SM_REGULAR = 'R',	// built-in multi-component Regular solid-solution model
    SM_NRTLLIQ = 'L',	// built-in multi-component NRTL model for liquid solutions
    SM_WILSLIQ = 'W',	// built-in multi-component Wilson model for liquid solutions
    SM_CGFLUID = 'F',	// built-in multi-component Churakov-Gottschalk (CG) fluid EOS model
    SM_PRFLUID = 'P',	// built-in Peng-Robinson-Stryjek-Vera (PRSV) fluid EOS model
    SM_SRFLUID = 'E',	// built-in Soave-Redlich-Kwong (SRK) fluid EOS model
    SM_PR78FL = '7',	// built-in Peng-Robinson (PR78) fluid EoS model (under construction)
    SM_CORKFL = '8',    // built-in compensated Redlich-Kwong (CORK) fluid EoS model (under construction)
    SM_AQDAV = 'D',		// built-in Davies model (with 0.3) for aqueous electrolytes
    SM_AQDH1 = '1',		// built-in Debye-Hueckel limiting law for aqueous electrolytes
    SM_AQDH2 = '2',		// built-in 2-term Debye-Hueckel model for aqueous electrolytes
    SM_AQDH3 = '3',		// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Karpov version)
    SM_AQDHH = 'H',		// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Helgeson version)
    SM_AQDHS = 'Y',		// built-in 3-term Debye-Hueckel model for aqueous electrolytes (Shvarov version)
    SM_AQSIT = 'S',		// built-in SIT model for aqueous electrolytes
    SM_AQEXUQ = 'Q',    // built-in EUNIQUAC model for aqueous electrolytes
    SM_AQPITZ = 'Z',    // built-in Pitzer HMW model for aqueous electrolytes
		// SM_IONEX = 'X',		// ion exchange (Donnan, Nikolskii) (reserved)
    SM_SURCOM = 'A',	// models of surface complexation at solid-aqueous interface
    SM_USERDEF = 'U',	// user-defined mixing model (scripts in Phase record)
    SM_OTHER = 'O'		// other built-in phase-specific models of non-ideal solutions
    	                //    (selected through phase name)
};

#ifndef _chbr_classes_h_
#define _chbr_classes_h_

typedef enum {  // classes of independent components IC, used in ccIC code list
    IC_ELEMENT  =  'e',  // chemical element (except oxygen and hydrogen)
    IC_OXYGEN   =  'o',  // oxygen
    IC_HYDROGEN =  'h',  // hydrogen (natural mixture of isotopes) H
    IC_PROTIUM   = 'p',  // protium Hp (reserved)
    IC_DEYTERIUM = 'd',  // deuterium D (reserved)
    IC_TRITIUM  =  't',  // tritium T (reserved)
    IC_FORMULA  =  'f',  // formula unit (eg. for Sio - a symbol of SiO2)
    IC_METALION =  'm',  // metal ion (cation), reserved
    IC_LIGAND   =  'l',  // ligand (anion), reserved
    IC_ADDIT    =  'a',  // IC with unknown stoichiometry (eg: Hum - humate ligand)
    IC_ISOTOPE  =  'i',  // isotope of chemical element (mass from 1 to 250)
    IC_OXYGEN16 =  'q',  // q  - oxygen 16O (reserved)
    IC_OXYGEN18 =  'r',  // r  - oxygen 18O (reserved)
    IC_CHARGE   =  'z',  // z  - electrical charge
    IC_VOLUME   =  'v',  // volume (for the volume balance constraint)
        IC_SITE     =  's'   // sorption site for site balance constraint (reserved)
} ICL_CLASSES;

typedef enum {  // Classifications of DC
    // Type of input data for
    SRC_DCOMP =  'd',  // the key points to existing PDB record in DCOMP chain
    SRC_REACDC = 'r',  // the key points to existing PDB record in REACDC chain
    SRC_NEWDC =  'n',  // the key new reaction-defined component
    SRC_NEWISO = 'i',  // the same as n, but this component is an isotopic form
    SRC_FICT =   'f',  // fictive species
    // Aqueous electrolyte phase:
    DC_AQ_PROTON   = 'T',      // hydrogen ion H+
    DC_AQ_ELECTRON = 'E',      // electron (as a DC)
    DC_AQ_SPECIES  = 'S',      // other aqueous species (ions, complexes and ion pairs)
DC_AQ_SURCOMP = 'K',     // Surface complex represented as aqueous species
    DC_AQ_SOLVENT  = 'W',      // water H2O (major solvent)
    DC_AQ_SOLVCOM  = 'L',      // other components of a solvent (eg. alcohol)
    // Gas phase ( G code can be used for all gases; V,C,H,N codes are reserved
    // for future use of the built-in equations of state in FGL module):
    DC_GAS_COMP    = 'G',   // other gases
    DC_GAS_H2O     = 'V',   // H2O steam
    DC_GAS_CO2     = 'C',   // CO2 (carbon dioxide)
    DC_GAS_H2      = 'H',   // H2 hydrogen
    DC_GAS_N2      = 'N',   // N2 nitrogen
    // Solid/liquid non-electrolyte multicomponent phases:
    DC_SOL_IDEAL   = 'I',   // end-member component with ideal behaviour
    DC_SOL_MINOR   = 'J',   // junior component (Henry's Law)
    DC_SOL_MAJOR   = 'M',   // major component (Raoult's Law)
    // Sorption phases and poly(oligo)electrolytes
    DC_SUR_CARRIER = 'Q',   // Principal end-member of solid carrier
    DC_SUR_MINAL   = 'P',   // Minor end-member of solid carrier
    DC_PEL_CARRIER = 'R',   // Carrier of poly(oligo)electrolyte

    DC_SSC_A0 = '0', DC_SSC_A1 = '2', DC_SSC_A2 = '4', DC_SSC_A3 = '6',
    DC_SSC_A4 = '8', // Strong surface complex on site type 0,1,2,3,4 - A plane
    DC_WSC_A0 = '1', DC_WSC_A1 = '3', DC_WSC_A2 = '5', DC_WSC_A3 = '7',
    DC_WSC_A4 = '9', // Weak surface complex on site type 0,1,2,3,4 - B plane
    DC_IESC_A  = 'A', // Strong exchange ion const-charge plane
    DC_IEWC_B  = 'B', // Weak exchange ion const-charge plane

    // Aliaces for 1-site model
    DC_SUR_GROUP    = 'X',  // Surface group on A plane -> '0'
    DC_SUR_COMPLEX = 'Y',   // Strong sur. complex A plane -> '0'
    DC_SUR_IPAIR   = 'Z',   // Weak sur complex B plane -> '1'

    // Single-component phases:
    DC_SCP_CONDEN  = 'O'   // DC forming a single-component phase

} DC_CLASSES;


//    This code defines standard state and reference scale of concentra-
// tions for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES{  // Possible values
    PH_AQUEL    = 'a',  // aqueous electrolyte
    PH_GASMIX   = 'g',  // mixture of gases
    PH_FLUID    = 'f',  // fluid phase
    PH_PLASMA   = 'p',  // plasma
    PH_LIQUID   = 'l',  // non-electrolyte liquid (melt)
    PH_SIMELT   = 'm',  // silicate (magmatic) melt or non-aqueous electrolyte
    PH_SORPTION = 'x',  // dilspersed solid with adsorption (ion exchange) in aqueous
    PH_POLYEL = 'y',    // colloidal poly- (oligo)electrolyte
    PH_SINCOND  = 's',  // condenced solid phase, also multicomponent
    PH_SINDIS   = 'd',  // dispersed solid phase, also multicomponent
    PH_HCARBL   = 'h'   // mixture of condensed hydrocarbons
};

#else

// This code defines standard state and reference scale of concentra-
// tions for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES2{  // Possible values
    PH_PLASMA   = 'p',  // plasma
    PH_SIMELT   = 'm',  // silicate (magmatic) melt or non-aqueous electrolyte
    PH_HCARBL   = 'h'   // mixture of condensed hydrocarbons
};

#endif

typedef enum {  // Limits on DC and phases
    // type of lmits
    NO_LIM = 'O', LOWER_LIM ='L', UPPER_LIM = 'U', BOTH_LIM ='B',
    // mode recalc of limits Set_DC_Limits()
    DC_LIM_INIT = 0, DC_LIM_CURRENT = 1
} DC_LIMITS;


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
    CCA_VOL = 'V', CCA_0 = '0', CCA_1, CCA_2, CCA_3, CCA_4, CCA_5,
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

typedef enum { // Units of measurement of quantities and concentrations
    // amounts of components and phases
    QUAN_MKMOL = 'Y',  QUAN_MMOL = 'h',  QUAN_MOL = 'M', // MOLES
    QUAN_MGRAM = 'y',  QUAN_GRAM = 'g',  QUAN_KILO = 'G',// MASSES
    // concentrations of components and phases
    CON_MOLFR = 'n', CON_MOLPROC = 'N', CON_pMOLFR = 'f', // MOLE FRACTION
    CON_VOLFR = 'v', CON_VOLPROC = 'V', CON_pVOLFR = 'u', // VOLUME FRACTION
    CON_MOLAL = 'm', CON_MMOLAL =  'i', CON_pMOLAL = 'p', // MOLALITY
    CON_MOLAR = 'L', CON_MMOLAR =  'j', CON_pMOLAR = 'q', // MOLARITY
    CON_WTFR  = 'w', CON_WTPROC =  '%', CON_PPM =    'P', // MASS FRACTION
    CON_AQWFR = 'C', CON_AQWPROC = 'c', CON_AQPPM =  'a', // CONCENTRATION
    // aqueous species
    CON_AQGPL = 'd', CON_AQMGPL = 'e', CON_AQMKGPL = 'b',//VOLUME CONCENTRATION

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

const char S_OFF = '-',
                   S_ON = '+',
                          S_REM = '*',
                                  A_NUL ='?';

#endif
// m_const.h in GEMIPM2K
