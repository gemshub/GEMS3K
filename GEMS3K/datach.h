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


#endif

// -----------------------------------------------------------------------------
// End of datach.h


