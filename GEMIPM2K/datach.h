//-------------------------------------------------------------------
// $Id: datach.h 1430 2009-08-27 17:02:11Z gems $
//
// DataCHemistry - contains chemical system definitions common to all
// nodes for the exchange between the GEM IPM and the FMT code parts.
// Contains dimensions and index lists for ICs, DCs, Phases in DATABR structure.
// Also contains thermodynamic data as grid arrays for interpolation over T,P.
// Used in TNode and TNodeArray classes
//
//      CH: chemical structure in GEM IPM
//      FMT: fluid mass transport

// Copyright (C) 2003,2009 by D.Kulik, S.Dmytriyeva, F.Enzmann, W.Pfingsten
// This file is part of GEMIPM2K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization
// developed in the Laboratory for Waste Management, Paul Scherrer Institute

// This file may be distributed together with GEMIPM2K source code
// under the licence terms defined in GEMIPM2K.QAL
//
// See also http://gems.web.psi.ch/
// E-mail: gems2.support@psi.ch
//------------------------------------------------------------------------------
//
#ifndef _DataCh_H_
#define _DataCh_H_

const long int
    MaxICN =      6,      // IC name length
    MaxDCN =      16,     // DC name length
    MaxPHN =      16;     // PH name length

typedef struct   // Structure DataCH
{
  long int     // Dimensionalities chemical system definition
//  These dimensionalities should be the same as in the GEMIPM work structure (MULTI)
    nIC,    // Number of Independent Components (stoichiometry units, usually chemical elements and charge)
    nDC,    // Total number of Dependent Components (chemical species made of Independent Components)
    nPH,    // Number of phases (into which Dependent Components are grouped)
    nPS,    // Number of phases-solutions (multicomponent phases) in the chemical system definition, nPS <= nPH
    nDCs,   // Number of Dependent Components in phases-solutions (multicomponent phases)
    nTp,    // Number of temperature grid points in interpolation lookup arrays, 1 or more
    nPp,    // Number of pressure grid points in interpolation lookup arrays, 1 or more
    iGrd,   // Flag for selection Diffusition coefficients array provided in the DCH file. (0 or 1)
    nAalp,  // Flag for keeping specific surface areas of phases in DATABR structure (1) or ignoring them (0)

  // These dimensionalities define sizes of packed arrays in DATABR structures
  // describing nodes. They are needed to save on the storage demand for nodes.
  // Connection between any node and DATACH occurs through the xIC, xPH and xDC
  // index lists (see below)
    nICb,   // Number of Independent Components kept in the DBR file and DATABR memory structure (<= nIC)
    nDCb,  	// Number of Dependent Components kept in the DBR file and DATABR memory structure (<=nDC)
    nPHb,   // Number of Phases to be kept in the DBR file and DATABR structure (<= nPH)
    nPSb,   // Number of Phases-solutions (multicomponent phases) to be kept in the DBR file and DATABR memory structure (<= nPS)
    uRes1,  // reserved

// Lists, vectors and matrices
    *nDCinPH,  // This vector tells how many Dependent Components is included in each phase [nPH]

  // Indices connecting the lists used in nodes (DATABR structure), see
  //    databr.h, with the lists in this (DATACH) structure
    *xic,   // DATACH access index list for IC kept in the DATABR structure and in DBR files [nICb]
    *xdc,   // DATACH access index list of DC kept in the DATABR  structure and in DBR files [nDCb]
    *xph;   // DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]

  double
    Ttol,    // Tolerance for the temperature interpolation (K)
    Ptol,    // Tolerance for the pressure interpolation (Pa)
    dRes1,   // reserved
    dRes2,   // reserved

// Data vectors - must be loaded before calling GEMIPM2K
    *TKval,  // Temperature values for the interpolation grid (Kelvin) for the lookup arrays of thermodynamic data [nTp]
    *Pval,   // Pressure values for the interpolation grid (Pa) for the lookup arrays of thermodynamic data [nPp]
    *A,      // Stoichiometry matrix A for Dependent Components. [nIC][nDC] elements

   // Values for IC (independent components)
    *ICmm,   // Atomic (molar) masses of Independent Components  (kg/mol) [nIC]

    // DC - related values
    *DCmm,   // Molar masses of Dependent Components (kg/mol) [nDC]
    *DD,     // Lookup array for diffusion coefficients of DCs (reserved) [nDC][nPp][nTp]  for now constant

    // Look-up grid arrays of thermodynamic data require a Lagrange interpolation subroutine to extract data
    // for a given P,T point (new interpolation is done when P or T differs
    // from the previous P,T by more than Ptol, Ttol)
    *denW,  // Lookup array for the density of water-solvent (kg/m3) [5][nPp][nTp]
    *denWg, // Optional lookup array for the density of water vapour (kg/m3) [5][nPp][nTp]
//  *visW, // Optional lookup array for the viscosity of liquid water (units?) [5][nPp][nTp] reserved
    *epsW,  // Lookup array for the dielectric constant of water-solvent (dimensionless) [5][nPp][nTp]
    *epsWg, // Optional lookup array for the dielectric constant of water vapour [5][nPp][nTp]
    *G0,    // Obligatory lookup array for DC molar Gibbs energy function g(T,P) (J/mol) [nDC][nPp][nTp]
    *V0,    // Obligatory lookup array for (standard) molar volumes of DC V(T,P) (J/Pa) [nDC][nPp][nTp]
    *S0,    // Optional lookup array for the DC absolute entropy function S(T,P) (J/K/mol) [nDC][nPp][nTp]
    *H0,    // Optional lookup array for DC molar enthalpy h(T,P) (J/mol) [nDC][nPp][nTp]
    *Cp0,   // Optional lookup array for DC heat capacity function Cp(T,P) (J/K/mol) [nDC][nPp][nTp]
    *A0,    // Optional lookup array for Helmholtz energy of DC (J/mol) reserved, [nDC][nPp][nTp]
    *U0;    // Optional lookup array for Internal energy of DC (J/K/mol) [nDC][nPp][nTp]

// Name lists
  char (*ICNL)[MaxICN]; // List of IC names in the system, [nIC]  of MaxICN length
  char (*DCNL)[MaxDCN]; // List of DC names in the system, [nDC] of MaxDCN length
  char (*PHNL)[MaxPHN]; // List of Phase names  [nPH]  of MaxPHN length

// Class code lists
   char *ccIC,   // Class codes of IC, see  enum ICL_CLASSES  [nIC]
        *ccDC,   // Type codes of DC, see  enum DCL_CLASSES  [nDC]
        *ccPH;   // Class codes of phases, see enum PHL_CLASSES [nPH]
}
DATACH;

#ifdef IPMGEMPLUGIN

#ifndef _chbr_classes_h_
#define _chbr_classes_h_

typedef enum {  // classes of independent components IC used in ccIC code list
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

typedef enum {  // Classes of dependent components DC used in ccDC code list

	// Single-component (pure) condensed phases:
    DC_SCP_CONDEN  = 'O',       // DC forming a single-component phase

	// Solid/liquid non-electrolyte multicomponent phases:
    DC_SOL_IDEAL   = 'I',   // ideal end-member component (Raoult)
    DC_SOL_MINOR   = 'J',   // junior (minor) component (Henry)
    DC_SOL_MAJOR   = 'M',   // major component (Raoult)

	// Aqueous electrolyte phase:
    DC_AQ_PROTON   = 'T',      // hydrogen ion H+
    DC_AQ_ELECTRON = 'E',      // electron (as a DC)
    DC_AQ_SPECIES  = 'S',      // other aqueous species (ions, complexes and ion pairs)
    DC_AQ_SURCOMP = 'K',     // Surface complex represented as aqueous species
    DC_AQ_SOLVENT  = 'W',      // water H2O (major solvent)
    DC_AQ_SOLVCOM  = 'L',      // other components of a solvent (eg. alcohol)

	// Gas phase ( G code can be used for all gases; V,C,H,N codes are reserved
    // for future use in the built-in equations of state):
    DC_GAS_COMP    = 'G',      // other gases
    DC_GAS_H2O     = 'V',      // H2O steam
    DC_GAS_CO2     = 'C',      // CO2 (carbon dioxide)
    DC_GAS_H2      = 'H',      // H2 hydrogen
    DC_GAS_N2      = 'N',      // N2 nitrogen

	// Sorption phases and poly(oligo)electrolytes
    DC_SUR_CARRIER = 'Q',   // Principal end-member of solid carrier (sorbent)
    DC_SUR_MINAL   = 'P',   // Minor end-member of solid carrier (sorbent)
    DC_PEL_CARRIER = 'R',   // Carrier of poly(oligo)electrolyte (for future use)

    // GEM CD-MUSIC and NE surface complexation models
    DC_SUR_GROUP   = 'X',   // Surface group (surface solvent), also fictive
    DC_SUR_COMPLEX = 'Y',   // Inner-sphere (strong) surface complex, the same as '0' code
    DC_SUR_IPAIR   = 'Z',   // Outer-sphere (weak) surface complex, surface ion pair,
	                          // exchange ion (the same as '1')

    // Obsolete codes for old GEM SCMs - usage in newly created models is not recommended
    DC_SSC_A0 = '0', DC_SSC_A1 = '2', DC_SSC_A2 = '4', DC_SSC_A3 = '6',
    DC_SSC_A4 = '8',        // Strong surface complex on site type 0,1,2,3,4 - A plane
    DC_WSC_A0 = '1', DC_WSC_A1 = '3', DC_WSC_A2 = '5', DC_WSC_A3 = '7',
    DC_WSC_A4 = '9',        // Weak surface complex on site type 0,1,2,3,4 - B plane
    DC_IESC_A  = 'A',       // Strong exchange ion const-charge plane
    DC_IEWC_B  = 'B',       // Weak exchange ion const-charge plane

    // Special class codes for diffusing species etc. (reserved)
    DCaquoCATION   = 'c',
    DCaquoANION    = 'n',
    DCaquoLIGAND   = 'l',
    DCaquoCOMPLEX  = 'x',
    DCaquoIONPAIR  = 'p',
    DCaquoGAS      = 'g',

} DCL_CLASSES;

typedef enum {  // Classes of Phases used in ccPH code list
    PH_AQUEL    = 'a',  // aqueous electrolyte (also with HKF EoS)
    PH_GASMIX   = 'g',  // mixture of gases (also corresponding states theory)
    PH_FLUID    = 'f',  // supercritical fluid phase with special EoS
    PH_LIQUID   = 'l',  // non-electrolyte liquid (melt)
    PH_SORPTION = 'x',  // dispersed solid with adsorption (ion exchange) in aqueous system
    PH_POLYEL   = 'y',  // colloidal poly- (oligo)electrolyte (reserved)
    PH_SINCOND  = 's',  // condenced solid phase, also multicomponent (solid solution)
    PH_SINDIS   = 'd',  // dispersed solid phase, also multicomponent
} PHL_CLASSES;

#endif
#endif
/*
// Codes allowed in Generic DC code list ccDCW (structure MULTI)
enum SolDCLcodes {
    DCl_SINGLE = 'U',        // This DC is a single-component (pure) phase
    DCl_SYMMETRIC = 'I',     // This DC is symmetric component (end member)
	                         // of a solution phase or a gas mixture
    DCl_ASYM_SPECIES = 'S',  // This DC is asymmetric component
                             // (solute, sorbate species) in a solution or sorption phase
    DCl_ASYM_CARRIER = 'W'   // This is a solvent or carrier in a solution
                             // (sorption) phase
};
*/
#endif
// -----------------------------------------------------------------------------
// End of datach.h


