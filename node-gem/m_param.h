//-------------------------------------------------------------------
// $Id: m_param.h 921 2007-09-10 12:04:49Z gems $
//
// Copyright (C) 2006-2007  S.Dmitrieva, D.Kulik
//
// Base parameters and codes used in GEM IPM work structure
// (standalone version)
//
// This file is part of the standalone GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
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


// Physical constants - see m_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg;
//

struct BASE_PARAM
{ // Flags and thresholds for numeric modules
    short PC,    // "Mode of Selekt2() subroutine operation"
    PD,  // Mode of DHH():0-invoke,1-at FIA only,2-last IPM it. 3-every IPM it.
    PRD, //Positive: mode GammaCalc(): 0-inactive, 1-FIA only, 2-SELEKT2 only 3-every IPM iteration
         // Negative: number of additional EFD-IPM improvement loops (-1 to -30 )
    PSM, //  Mode of setting FIA Gamma in phases: 0-off; 1 - set activ.coeff. GAN
    DP,  //  Max number of EnterFeasibleDomain() iterations { default 144 }
    DW,  // IPM-2 precision mode: 0- disable >=1  - enable (number of loops, default 15)
    DT,  // Exponent for restoring x_j from dual solution {-5 to +2, default -3}
    PLLG, // IPM view debug time on Progress Dialog screen from (sec) or 0
    PE,  // Include electroneutrality condition? { 0 1 }
    IIM  // Maximum number of iterations in the main IPM algorithm { 500 }
    ;
    double DG,   // Precision of LEM IPM solver (Jordan) { 1e-15 }
    DHB,  // Precision of mass balance deviations at EFD { 1e-8 }
    DS,   // Cutoff number of moles of a phase { 1e-10 }
    DK,   // Threshold of Dikin criterion of IPM convergence { 1e-5 }
    DF,   // Threshold of Karpov' criteria for phase stability { 0.01 }
    DFM,  // Threshold for Karpov' criteria for insertion of phases { -0.1 }
    DFYw, // N of moles of H2O-solvent for phase insertion { 2e-5 }
    DFYaq,// N of moles of aqueous DC for phase insertion { 1e-6 }
    DFYid,// N of moles of DC in ideal solution for phase insertion { 1e-6 }
    DFYr, // Cutoff moles of major DC in non-ideal solution { 1e-6 }
    DFYh, // Cutoff moles of junior DC in non-ideal solution { 1e-6 }
    DFYc, // N of moles for insertion of single-component phase { 1e-6 }
    DFYs, // Moles of DC to boost phases-solutions (SELEKT2) { 1e-6 }
    DB,   // Cutoff number of moles of IC in the system 1e-9
    AG,   //Smoothing parameter for free energy increments on iterations {0.7}
    DGC,  // Power to raise Pa_AG on steps of smoothing { 0.07 }
    GAR,  // Initial activity coefficient of major DC in a phase (FIA) { 1 }
    GAH,  // Initial activity coefficient of minor DC in phase (FIA) { 1000 }
    GAS,  // IPM-2 precision factor (DHB[i] / b[i]) { 1e-3 }
    DNS,  // Default standard surface site density 1/nm2 (12.05)
    XwMin, // Cutoff N of moles of water in aqueous phase { 1e-5 }
    ScMin, // Cutoff N of moles of main DC in solid carrier for sorption {1e-7}
    DcMin, //  Cutoff N of moles of solute (sorbate) { 1e-19 }
    PhMin, // Cutoff N of moles of a phase { 1e-12 }
    ICmin, // Cutoff value of ionic strength to run DH equation { 1e-5 }
    EPS,   // Precision of FIA by SIMPLEX-method { 1e-15 }
    IEPS,  // Accuracy of cutoff limits 1+eps for surface activity terms {1e-8}
    DKIN;  // Range of lower to upper kinetic restrictions for DC { 0.001 }
    char *tprn;       // internal

    void write(ostream& oss);
};

struct SPP_SETTING
{   // Base Parametres of SP
    char ver[TDBVERSION]; // Version & Copyright 64
    BASE_PARAM p; //
    void write(ostream& oss);
};


// Module TParam (RMULTS+MTPARM+SYSTEM+MULTY see SubModules)
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
   void outMultiTxt( const char *path  );
   void readMulti( GemDataStream& ff );
   void readMulti( const char* path );
   double calcMulti();

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
