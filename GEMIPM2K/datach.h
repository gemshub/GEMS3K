//-------------------------------------------------------------------
// $Id: datach.h 1066 2008-05-16 14:16:59Z gems $
//
// DataCHemistry - contains chemical system definitions common to all
// nodes for the exchange between the GEM IPM and the FMT code parts.
// Contains dimensions and index lists for ICs, DCs, Phases in DATABR structure.
// Also contains thermodynamic data as grid arrays for interpolation over T,P.
// Used in TNode and TNodeArray classes
//
//      CH: chemical structure in GEM IPM
//      FMT: fluid mass transport

// Copyright (C) 2003,2008 by D.Kulik, S.Dmytriyeva, W.Pfingsten, F.Enzmann
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

const unsigned int
    MaxICN =      6,      // IC name length
    MaxDCN =      16,     // DC name length
    MaxPHN =      16;     // PH name length

typedef struct
{  // Structure DataCH
// Dimensionalities
  short
//  These dimensionalities should be the same as in the GEMIPM work structure (MULTI)
    nIC,    // Total number of IC (independent components) in the reactive part
    nDC,    // Total number of DC (chemical species) in the reactive part
    nPH,    // Total number of phases included into the GEM IPM problem
    nPS,    // Number of multicomponent phases, nPS <= nPH
    nDCs,   // Total Number of DC in phases-solutions
    nTp,    // Number of temperature points in grid arrays
	    //       for the interpolation of thermodynamic data
    nPp,    // Number of pressure points in grid arrays
            //       for the interpolation of thermodynamic data
    iGrd,   // flag for grid array setup: 0 - only V0 and G0; 1 - plus H0;
            //    2 - plus S0; 3 - plus Cp0; 4 - plus A0 (Helmholtz);
            //    -1 - V0, G0 and DD (diffusion coefficients in aq, reserved)
    nAalp,  // Flag for considering surface areas of phases

  // These dimensionalities define sizes of packed arrays in DATABR structures
  // describing nodes. They are needed to save on the storage demand for nodes.
  // Connection between any node and DATACH occurs through the xIC, xPH and xDC
  // index lists (see below)
    nICb,       // number of IC (stoichiometry units) (<= nIC) used in nodes
    nDCb,      	// number of DC (chemical species, <= nDC) used in nodes
    nPHb,     	// number of Phases (<= nPH) used in nodes
    nPSb,       // number of Phases-solutions (<= nPS) used in nodes
    uRes1,      // reserved

// Lists, vectors and matrices
    *nDCinPH,  // number of DC included into each phase, [nPH] elements

// Indices connecting the lists used in nodes (DATABR structure), see
//    databr.h, with the lists in this (DATACH) structure
    *xIC,   // IC name indices in DATABR IC vectors, [nICb] elements
    *xDC,   // DC name indices in DATABR DC vectors, [nDCb] elements
    *xPH;   // PH name indices in DATABR Phase vectors, [nPHb] elements
            // see below definitions of the ICNL, DCNL and PHNL lists

  float
    *TCval,   // discrete values of Temperature (C), [nTp] elements,
 // that correspond to grid arrays for the interpolation of thermodynamic data
    *Pval,   // discrete values of Pressure (bar), [nPp] elements,
 // that correspond to grid arrays for the interpolation of thermodynamic data
    *A;      // Stoichiometry matrix A containing elemental stoichiometries
             // of Dependent Components, [nIC][nDC] elements
 double
    Ttol,    // Temperature tolerance (K) for interpolation of thermodynamic data
    Ptol,    // Pressure tolerance (bar) for interpolation of thermodynamic data
    dRes1,   // reserved
    dRes2,   // reserved

// Data vectors - must be loaded before calling GEMIPM2K

// Values for IC (independent components)
    *ICmm,   // IC atomic (molar) mass, g/mol, [nIC] elements

// DC - related values
    *DCmm,   // DC molar mass, g/mol, [nDC] elements
    *DD,     // Diffusition coefficients, [nDC][nPp][nTp] elements, for now constant

// Look-up grid arrays of thermodynamic data
// Require a Lagrange interpolation subroutine to extract data
// for a given P,T point (new interpolation is done when P or T differs
// from the previous P,T by more than Ptol, Ttol)
    *roW,   // density of water-solvent, g/cm3, [ nPp][nTp] elements
    *epsW,  // dielectric  constant of water-solvent, [nPp][nTp] elements
    *G0,    // G0 standard molar Gibbs energy of DC, J/mol, [nDC][nPp][nTp] elements
    *V0,    // V0 standard molar volume of DC, J/bar, [nDC][nPp][nTp] elements
    *S0,    // S0 standard molar entropy of DC, J/K/mol, [nDC][nPp][nTp] elements
    *H0,    // H0 standard molar enthalpy of DC, J/mol, reserved, [nDC][nPp][nTp] elements
    *Cp0;   // Cp0 molar heat capacity of DC, J/K/mol, [nDC][nPp][nTp] elements

// Name lists
   // List of IC names in the system, [nIC] elements of MaxICN length
 char (*ICNL)[MaxICN];
   // List of DC names in the system, [nDC] elements of MaxDCN length
 char (*DCNL)[MaxDCN];
   // List of phase names in the system, [nPH] elements of MaxPHN length
 char (*PHNL)[MaxPHN];

// Class code lists
 char   *ccIC,   // Class codes of IC, see  enum ICL_CLASSES  ([nIC] elements)
        *ccDC,   // Class codes of DC, see  enum DCL_CLASSES  ([nDC] elements)
        *ccPH;   // Class codes of phases, see enum PHL_CLASSES ([nPH] elements)
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

#endif
// -----------------------------------------------------------------------------
// End of datach.h


