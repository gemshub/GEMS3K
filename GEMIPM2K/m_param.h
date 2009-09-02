//-------------------------------------------------------------------
// $Id: m_param.h 1387 2009-08-07 12:31:14Z gems $
//
// Copyright (C) 2006,2007  S.Dmitrieva, D.Kulik
//
// Basic codes and numerical settings used in GEM IPM-2 operation
// (standalone version)
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
#ifndef _m_param_h_
#define _m_param_h_


#ifndef IPMGEMPLUGIN
 #define IPMGEMPLUGIN
#endif

#include <math.h>
#include "gdatastream.h"
#include "ms_multi.h"
#include "verror.h"


// Physical constants - see ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;
//

struct BASE_PARAM
{ // Flags and thresholds for numeric modules
   long int
    PC,   // Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
    PD,   // Mode of execution of the built-in DebyeHueckel()/Davies() functions { 3 }
          // Mode of DHH():0-invoke,1-at FIA only,2-last IPM it. 3-every IPM it.
    PRD,  //  Negative number (from -1 to -50): the number |PRD| of additional full IPM-2 loops
          //  to improve the GEM final solution, else  no additional loops, { 3 }
    PSM,  // Level of diagnostic messages: 0- disabled (no ipmlog file); 1- normal; 2-including warnings { 1 }
    DP,   // Maximum allowed number of iterations in the EnterFeasibleDomain() procedure {  150 }
    DW,   // Maximum number of allowed IPM-2 mass balance accuracy improvement loops { 0- disable >=1  - enable ; default 15}
    DT,   // IPM-2 cutoff exponent to recover small elements of the primal x vector using the dual solution vector u { -3 }
    PLLG, // TIPM-2 tolerance for checking change in dual solution among GEM refinement loops { 200 }
          //      { 0 to 1000 mol/mol, default 0 or 32000 means no check }
    PE,   // Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
    IIM   // Maximum allowed number of iterations in the MainIPM_Descent() procedure { 500 }
    ;
  double DG,   // Threshold for minimum descent step size Lambda in EntryFeasibleDomain() { 1e-5 }
    DHB,  // Maximum allowed mass balance residual (moles) for major Independent Components { 1e-8 }
    DS,   // Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
    DK,   // IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-3) { 1e-4 }
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
    AG,   // Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations {0.7}
    DGC,  // Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { 0.07 }
    GAR,  // Initial activity coefficient value for major (M) species in a solution phase before Simplex() approximation { 1 }
    GAH,  // Initial activity coefficient value for minor (J) species in a solution phase before Simplex() approximation { 1000 }
    GAS,  // IPM-2 balance accuracy control ratio DHBM[i]/b[i], determines the maximum allowed mass balance residual for minor IC { 1e-3 }
    DNS,  // Standard surface density (nm-2) for calculating activity of surface species (12.05)
    XwMin,// Cutoff mole amount for elimination of water-solvent { 1e-9 }
    ScMin,// Cutoff mole amount for elimination of solid sorbent {1e-7}
    DcMin,// Cutoff mole amount for elimination of solution- or surface species { 1e-20 }
    PhMin,// Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
    ICmin,// Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
    EPS,  // Precision criterion of the simplex() procedure to obtain the automatic initial approximation { 1e-6 to 1e-14, default 1e-7 }
    IEPS, // Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
    DKIN; // Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
    char *tprn;       // internal

    void write(fstream& oss);
};

struct SPP_SETTING
{   // Base Parametres of SP
    char ver[TDBVERSION]; // Version & Copyright 64
    BASE_PARAM p; //
    void write(fstream& oss);
};

extern SPP_SETTING pa_;

// Module TParam ( +MULTY )
class TProfil //: public TCModule
{

public:

    static TProfil* pm;

    TMulti* multi;
    MULTI *pmp;
    SPP_SETTING pa;

    TProfil( TMulti* amulti );

    const char* GetName() const
    {
        return "Project";
    }

   void outMulti( GemDataStream& ff, gstring& path  );
   void outMultiTxt( const char *path, bool append=false  );
   void readMulti( GemDataStream& ff );
   void readMulti( const char* path );

   double calcMulti( long int& PrecLoops_, long int& NumIterFIA_, long int& NumIterIPM_ );
   long int testMulti( );
   void test_G0_V0_H0_Cp0_DD_arrays( long int nT, long int nP );
};

// Work DC classifier codes  pm->DCCW
enum SolDCodes {

    DC_SINGLE = 'U',        // This DC is a single-component phase
    DC_SYMMETRIC = 'I',     // This DC is in symmetric solution phase
    DC_ASYM_SPECIES = 'S',  // This is DC-solute(sorbate) in asymmetric phase
    DC_ASYM_CARRIER = 'W'   // This is carrier(solvent) DC in asymmetric phase
};

enum QpQdSizes {   // see m_phase.h
   QPSIZE = 60,    // This enum is for GEMIPM2K only!
   QDSIZE = 60
};

#endif
// _m_param_h
