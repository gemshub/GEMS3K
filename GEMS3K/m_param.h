//-------------------------------------------------------------------
// $Id: m_param.h 1406 2009-08-21 08:37:56Z gems $
//
// Declaration of TProfil class, config and calculation functions
//
// Rewritten from C to C++ by S.Dmytriyeva
// Copyright (C) 1995-2007 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#ifndef _m_param_h_
#define _m_param_h_

#include <math.h>

// Physical constants - see m_param.cpp or ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;


#include "gdatastream.h"
// #include "ms_multi.h"
#include "verror.h"


typedef struct
{ // Flags and thresholds for numeric modules
   long int
           PC,   // Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
           PD,   // abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }
                 // Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 // if PD < 0 then use test qd_real accuracy mode
           PRD,  // Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
           PSM,  // Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   // Maximum allowed number of iterations in the MassBalanceRefinement() procedure {  30 }
           DW,   // Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   // Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, // IPM tolerance for detecting divergence in dual solution { 10; range 1 to 1000; 0 disables the detection }
           PE,   // Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
           IIM   // Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 1000 }
           ;
         double DG,   // Standart total moles { 1e5 }
           DHB,  // Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-10 }
           DS,   // Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
           DK,   // IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-4) { 1e-5 }
           DF,   // Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
           DFM,  // Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.1 }
           DFYw, // Insertion mole amount for water-solvent { 1e-6 }
           DFYaq,// Insertion mole amount for aqueous species { 1e-6 }
           DFYid,// Insertion mole amount for ideal solution components { 1e-6 }
           DFYr, // Insertion mole amount for major solution components { 1e-6 }
           DFYh, // Insertion mole amount for minor solution components { 1e-6 }
           DFYc, // Insertion mole amount for single-component phase { 1e-6 }
           DFYs, // Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
           DB,   // Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
           AG,   // Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { -1 }
           DGC,  // Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { -0.99 }
           GAR,  // Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
           GAH,  // Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
           GAS,  // Since r1583/r409: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }
                 // before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
           DNS,  // Standard surface density (nm-2) for calculating activity of surface species (12.05)
           XwMin,// Cutoff mole amount for elimination of water-solvent { 1e-9 }
           ScMin,// Cutoff mole amount for elimination of solid sorbent {1e-7}
           DcMin,// Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
           PhMin,// Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
           ICmin,// Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
           EPS,  // Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
           IEPS, // Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
           DKIN; // Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
 //   char *tprn;       // internal
  void write(fstream& oss);
} BASE_PARAM ;

typedef struct   
{   // Base Parametres of SP
    char *ver; // Version & Copyright 64
    BASE_PARAM p; //
   void write(fstream& oss);
}SPP_SETTING;


// Module TParam ( +MULTY )
//class TProfil //: public TCModule
//{

//public:

 //   static TProfil* pm;

 //   TMulti* multi;
 //   MULTI *pmp;
 //   SPP_SETTING pa;

 //   TProfil( TMulti* amulti );



//};

enum QpQdSizes {   // see m_phase.h
   QPSIZE = 60,    // This enum is for GEMIPM2K only!
   QDSIZE = 60
};



// Work DC classifier codes  pm->DCCW
enum SolDCodes {

    DC_SINGLE = 'U',        // This DC is a single-component phase
    DC_SYMMETRIC = 'I',     // This DC is in symmetric solution phase
    DC_ASYM_SPECIES = 'S',  // This is DC-solute(sorbate) in asymmetric phase
    DC_ASYM_CARRIER = 'W'   // This is carrier(solvent) DC in asymmetric phase
};


#endif  // _m_param_h


