//-------------------------------------------------------------------
// $Id$
//
/// \file node.cpp
/// Implementation of TNode class functionality including initialization
/// and execution of the GEM IPM 3 kernel
/// Works with DATACH and DATABR structures
//
// Copyright (c) 2005-2012 S.Dmytriyeva, D.Kulik, G.Kosakowski, F.Hingerl
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include "node.h"
#include "num_methods.h"
#include "kinetics.h"
#include "v_service.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>


// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> TNode::node_logger = spdlog::stdout_color_mt("tnode");
// Thread-safe logger to file
std::shared_ptr<spdlog::logger> TNode::ipmlog_file = spdlog::rotating_logger_mt("ipmlog", "ipmlog.txt", 1048576, 3);

// Conversion factors
const double bar_to_Pa = 1e5,
m3_to_cm3 = 1e6,
kg_to_g = 1e3;

void TNode::write_ThermoFun_format_stream(std::iostream &stream, bool compact) {
    stream << thermo_json_string;
}

// Constructor of the class instance in memory for standalone GEMS3K or coupled program
TNode::TNode()
{
    CSD = NULL;
    CNode = NULL;
    allocMemory();
    dbr_file_name = "dbr_file_name";
    load_thermodynamic_data = false;
}

TNode::TNode( const TNode& otherNode )
{
    CSD = 0;
    CNode = 0;

    allocMemory();
    dbr_file_name = otherNode.dbr_file_name;
    // copy data from otherNode
    datach_copy( otherNode.CSD );
    databr_copy( otherNode.CNode );

    multi_ptr()->copyMULTI( *otherNode.multi_base );

    // copy intervals for minimizatiom
    pmm->Pai[0] = CSD->Pval[0]/bar_to_Pa;
    pmm->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
    pmm->Pai[2] = getStep( pmm->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
    pmm->Pai[3] = CSD->Ptol/bar_to_Pa;

    pmm->Tai[0] = CSD->TKval[0]-C_to_K;
    pmm->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
    pmm->Tai[2] = getStep( pmm->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
    pmm->Tai[3] = CSD->Ttol;

    pmm->Fdev1[0] = 0.;
    pmm->Fdev1[1] = 1e-6;   // 24/05/2010 must be copy from GEMS3 structure
    pmm->Fdev2[0] = 0.;
    pmm->Fdev2[1] = 1e-6;

    node_logger->info("copy constructor...");
}

TNode::~TNode()
{
    freeMemory();
}

void TNode::allocMemory()
{
    // memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;
    // mem_set( CSD, 0, sizeof(DATACH) );
    datach_reset();
    // mem_set( CNode, 0, sizeof(DATABR) );
    databr_reset( CNode, 2 );

    // allocation internal structures
    internal_multi.reset(new TMultiBase(this));
    multi_base = internal_multi.get();
    multi_base->set_def();
    pmm = multi_base->GetPM();
    atp.reset( new TActivity( CSD, CNode, this ) );
    //    atp->set_def();
    kip.reset( new TKinetics( CSD, CNode, this ) );
    kip->set_def();
}

void TNode::freeMemory()
{
    datach_free();
    // CSD = 0;
    delete CSD;
    CNode = databr_free( CNode );
}

std::string TNode::system_id() const
{
    return char_array_to_string(multi_ptr()->GetPM()->stkey, EQ_RKLEN);
}

// Checks if given temperature TK and pressure P fit within the interpolation
// intervals of the DATACH lookup arrays (returns true) or not (returns false)
bool  TNode::check_TP( double TK, double P ) const
{
    bool okT = true, okP = true;
    double T_=TK, P_=P;

    if( CSD->mLook == 1 )
    {
        for(long int  jj=0; jj<CSD->nPp; jj++)
            if( (fabs( P - CSD->Pval[jj] ) < CSD->Ptol ) && ( fabs( TK - CSD->TKval[jj] ) < CSD->Ttol ) )
            {
                return true;
            }
        Error( "check_TP: ", std::string("Temperature ")+std::to_string(TK)+
               " and pressure "+std::to_string(P)+" out of range");
        //return false;
    }
    else
    {
        if( TK <= CSD->TKval[0] - CSD->Ttol )
        { 				// Lower boundary of T interpolation interval
            okT = false;
            T_ = CSD->TKval[0] - CSD->Ttol;
        }
        if( TK >= CSD->TKval[CSD->nTp-1] + CSD->Ttol )
        {
            okT = false;
            T_ = CSD->TKval[CSD->nTp-1] + CSD->Ttol;
        }
        if( okT == false ) {
            ipmlog_file->info("In node {},  Given TK={} is beyond the interpolation "
                              "range for thermodynamic data near boundary T_= {}", CNode->NodeHandle, TK, T_);
        }

        if( P <= CSD->Pval[0] - CSD->Ptol )
        {
            okP = false;
            P_ = CSD->Pval[0] - CSD->Ptol;
        }
        if( P >= CSD->Pval[CSD->nPp-1] + CSD->Ptol )
        {
            okP = false;
            P_ = CSD->Pval[CSD->nPp-1] + CSD->Ptol;
        }
        if( !okP ) {
            ipmlog_file->info("In node {},  Given P={} is beyond the interpolation "
                              "range for thermodynamic data near boundary P_= {}", CNode->NodeHandle, P, P_);
        }
        return ( okT && okP );
    }
    return ( okT && okP );
}

//-------------------------------------------------------------------------------------------------------------------------------
// (2) Main call for GEM IPM calculations using the input bulk composition, temperature, pressure
//   and metastability constraints provided in the work instance of DATABR structure.
//   Actual calculation will be performed only when dBR->NodeStatusCH == NEED_GEM_SIA (5) or dBR->NodeStatusCH = NEED_GEM_AIA (1).
//   By other values of NodeStatusCH, no calculation will be performed and the status will remain unchanged.
//  In "smart initial approximation" (SIA) mode, the program can automatically switch into the "automatic initial
//  approximation" (AIA) mode and return  OK_GEM_AIA instead of OK_GEM_SIA.
//  Parameter:
//   uPrimalSol  flag to define the mode of GEM smart initial approximation
//               (only if dBR->NodeStatusCH = NEED_GEM_SIA has been set before GEM_run() call).
//               false  (0) -  use speciation and activity coefficients from previous GEM_run() calculation
//               true  (1)  -  use speciation provided in the DATABR memory structure (e.g. after reading the DBR file)
//  Return values:    NodeStatusCH  (the same as set in dBR->NodeStatusCH). Possible values (see "databr.h" file for the full list)
long int TNode::GEM_run( bool uPrimalSol )
{
    CalcTime = 0.0;
    PrecLoops = 0; NumIterFIA = 0; NumIterIPM = 0;
    clearipmLogError();

    try
    {
        ipmlog_file->debug(" GEM_run() begin Mode= {}", CNode->NodeStatusCH);
        // Checking T and P  for interpolation intervals
        check_TP( CNode->TK, CNode->P);
        // Unpacking work DATABR structure into MULTI (GEM IPM structure): uses DATACH
        // setting up up PIA or AIA mode
        if( CNode->NodeStatusCH == NEED_GEM_SIA )
        {
            pmm->pNP = 1;
            unpackDataBr( uPrimalSol );
        }
        else if( CNode->NodeStatusCH == NEED_GEM_AIA )
        {
            pmm->pNP = 0; // As default setting AIA mode
            if (CNode->dt > 0.)
                uPrimalSol = true;
            unpackDataBr( uPrimalSol );
        }
        else
            return CNode->NodeStatusCH;

        // added 18.12.14 DK : setting chemical kinetics time counter and variables
        node_logger->debug("GEM_run dTime:{}, TimeStep: {} Time: {}", CNode->dt, CNode->NodeStatusFMT, CNode->Tm);
        if( CNode->dt <= 0. )
        {  // no kinetics to consider
            pmm->kTau = 0.;
            pmm->kdT = 0.;
            pmm->ITau = -1;
            pmm->pKMM = 2;  // no need to allocate TKinMet instances
        }
        else {   // considering kinetics
            pmm->kTau = CNode->Tm;
            pmm->kdT = CNode->dt;
            if( pmm->ITau < 0 && CNode->Tm/CNode->dt < 1e-9 )
            {   // we need to initialize TKinMet
                pmm->pKMM = -1;
                pmm->ITau = -1;
            }
            else  // TKinMet exists, simulation continues
                pmm->pKMM = 1; // pmm->ITau = CNode->Tm/CNode->dt;
        }

        // GEM IPM calculation of equilibrium state
        CalcTime = multi_ptr()->CalculateEquilibriumState( /*RefineLoops,*/ NumIterFIA, NumIterIPM  );

        // Extracting and packing GEM IPM results into work DATABR structure
        packDataBr();
        CNode->IterDone = NumIterFIA+NumIterIPM;
        //**************************************************************
        // only for testing output results for files
        //    GEM_write_dbr( "calculated_dbr.dat",  2 /*json*/ );
        //    GEM_print_ipm( "calc_multi.ipm" );
        // *********************************************************

        // test error result GEM IPM calculation of equilibrium state in MULTI
        long int erCode =  multi_ptr()->testMulti();

        if( erCode )
        {
            if( CNode->NodeStatusCH  == NEED_GEM_AIA )
                CNode->NodeStatusCH = BAD_GEM_AIA;
            else
                CNode->NodeStatusCH = BAD_GEM_SIA;

            // internal multi error
            ipmlog_error = pmm->errorCode+ std::string(": ") +  pmm->errorBuf;
        }
        else
        {
            if( CNode->NodeStatusCH  == NEED_GEM_AIA )
                CNode->NodeStatusCH = OK_GEM_AIA;
            else
                CNode->NodeStatusCH = OK_GEM_SIA;
        }

        return CNode->NodeStatusCH;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
        if( CNode->NodeStatusCH  == NEED_GEM_AIA )
            CNode->NodeStatusCH = ERR_GEM_AIA;
        else
            CNode->NodeStatusCH = ERR_GEM_SIA;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
        CNode->NodeStatusCH = T_ERROR_GEM;
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        CNode->NodeStatusCH = T_ERROR_GEM;
    }
    ipmlog_file->error("Error Node:{}  time:{}  dt:{}", CNode->NodeHandle, CNode->Tm, CNode->dt);
    if( multi_ptr()->base_param()->PSM >= 2  )
        ipmlog_file->error("{}", ipmlog_error);

    return CNode->NodeStatusCH;
}

// Returns GEMIPM2 calculation time in seconds elapsed during the last call of GEM_run() -
// can be used for monitoring the performance of calculations.
// Return value:  double number, may contain 0.0 if the calculation time is less than the
//                internal time resolution of C/C++ function
double TNode::GEM_CalcTime( long int& NumK2, long int& NumIterFIA1, long int& NumIterIPM1)
{
    NumK2 = pmm->K2;
    NumIterFIA1 = pmm->ITF;
    NumIterIPM1 = pmm->ITG;
    return CalcTime;
}

// To obtain the number of GEM IPM2 iterations performed during the last call of GEM_run() e.g. for monitoring the
// performance of GEMS3K in AIA or SIA modes, or for problem diagnostics.
// Parameters:  long int variables per reference (must be allocated before calling GEM_Iterations(), previous values will be lost. See Return values.
// Return values:
//   Function         Total number of EFD + IPM iterations from the last call to GEM_run()
//   PrecLoops        Number of performed IPM-2 precision refinement loops
//   NumIterFIA       Total number of performed MBR() iterations to obtain a feasible initial approximation for the IPM algorithm.
//   NumIterIPM       Total number of performed IPM main descent algorithm iterations.
long int TNode::GEM_Iterations( long int& PrecLoops_, long int& NumIterFIA_, long int& NumIterIPM_ )
{
    PrecLoops_ = PrecLoops;
    NumIterFIA_ = NumIterFIA;
    NumIterIPM_ = NumIterIPM;
    return NumIterFIA+NumIterIPM;
}

// Extracting and packing GEM IPM results into work DATABR structure
void TNode::packDataBr()
{
    long int ii;

    // set default data to DataBr
    //   CNode->NodeStatusCH = NEED_GEM_AIA;
    if( pmm->pNP == 0 )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    else
        CNode->NodeStatusCH = NEED_GEM_SIA;

    CNode->TK = pmm->TCc+C_to_K; //25
    CNode->P = pmm->Pc*bar_to_Pa; //1
    //   CNode->IterDone = pmm->IT;
    CNode->IterDone = /*pmm->ITF+*/pmm->IT;   // Now complete number of FIA and IPM iterations
    // values
    CNode->Vs = pmm->VXc*1.e-6; // from cm3 to m3
    CNode->Gs = pmm->FX;
    CNode->Hs = pmm->HXc;
    CNode->IC = pmm->IC;
    CNode->pH = pmm->pH;
    CNode->pe = pmm->pe;
    //  CNode->Eh = pmm->FitVar[3];  Bugfix 19.12.2006  KD
    CNode->Eh = pmm->Eh;
    CNode->Ms = pmm->MBX;

    // arrays
    for( ii=0; ii<CSD->nPHb; ii++ )
    {  CNode->xPH[ii] = pmm->XF[ CSD->xph[ii] ];
        if( CSD->nAalp >0 )
            CNode->aPH[ii] = pmm->Aalp[ CSD->xph[ii] ]*kg_to_g;
        CNode->omPH[ii] = pmm->Falp[ CSD->xph[ii] ];
    }
    for( ii=0; ii<CSD->nPSb; ii++ )
    {   CNode->vPS[ii] = pmm->FVOL[ CSD->xph[ii] ]/m3_to_cm3;
        CNode->mPS[ii] = pmm->FWGT[ CSD->xph[ii] ]/kg_to_g;
        CNode->xPA[ii] = pmm->XFA[ CSD->xph[ii] ];
        CNode->amru[ii] = pmm->PUL[ CSD->xph[ii] ];
        CNode->amrl[ii] = pmm->PLL[ CSD->xph[ii] ];
    }
    for( ii=0; ii<CSD->nPSb; ii++ )
        for(long int jj=0; jj<CSD->nICb; jj++ )
        { long int   new_ndx= (ii*CSD->nICb)+jj,
                    mul_ndx = ( CSD->xph[ii]*CSD->nIC )+ CSD->xic[jj];
            CNode->bPS[new_ndx] = pmm->BF[ mul_ndx ];
        }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->xDC[ii] = pmm->X[ CSD->xdc[ii] ];
        CNode->gam[ii] = pmm->Gamma[ CSD->xdc[ii] ];
        CNode->dul[ii] = pmm->DUL[ CSD->xdc[ii] ];// always for GEM2MT init
        CNode->dll[ii] = pmm->DLL[ CSD->xdc[ii] ];// always for GEM2MT init
    }
    for( ii=0; ii<CSD->nICb; ii++ )
    {  CNode->bIC[ii] = pmm->B[ CSD->xic[ii] ]; // always for GEM2MT  init
        CNode->rMB[ii] = pmm->C[ CSD->xic[ii] ];
        CNode->uIC[ii] = pmm->U[ CSD->xic[ii] ];
        CNode->bSP[ii] = pmm->BFC[ CSD->xic[ii] ];
    }
}

// Unpacking work DATABR structure into MULTI
//(GEM IPM work structure): uses DATACH
//  if uPrimalSol is true then the primal solution (vectors x, gamma, IC etc.)
//  will be unpacked - as an option for PIA mode with previous GEM solution from
//  the same node.
//  If uPrimalSol = false then the primal solution data will not be unpacked
//  into the MULTI structure (AIA mode or SIA mode with primal solution retained
//    in the MULTI structure from previous IPM calculation)
void TNode::unpackDataBr( bool uPrimalSol )
{
    long int ii;

    CheckMtparam(); // T or P change detection - moved to here from InitalizeGEM_IPM_Data() 11.10.2012
    pmm->kTau = CNode->Tm;  // added 18.12.14 DK
    pmm->kdT = CNode->dt;   // added 18.12.14 DK

    pmm->TCc = CNode->TK-C_to_K;
    pmm->Tc = CNode->TK;
    pmm->Pc  = CNode->P/bar_to_Pa;
    pmm->VXc = CNode->Vs/1.e-6; // from cm3 to m3
    // Obligatory arrays - always unpacked!
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        pmm->DUL[ CSD->xdc[ii] ] = CNode->dul[ii];
        pmm->DLL[ CSD->xdc[ii] ] = CNode->dll[ii];
        if( pmm->DUL[ CSD->xdc[ii] ] < pmm->DLL[ CSD->xdc[ii] ] )
        {
            Error("unpackDataBr", std::string("Upper kinetic restriction less than the lower one for DC&RC")
                  +char_array_to_string( pmm->SM[CSD->xdc[ii]], MAXDCNAME));
        }
    }
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        pmm->B[ CSD->xic[ii] ] = CNode->bIC[ii];
        if( ii < CSD->nICb-1 && pmm->B[ CSD->xic[ii] ] < multi_ptr()->base_param()->DB )
        {
            Error("unpackDataBr", std::string("Bulk mole amount of IC ")+
                  char_array_to_string(pmm->SB[CSD->xic[ii]], 6)+" is "+
                    std::to_string(pmm->B[ CSD->xic[ii] ])+" - out of range" );
        }
    }
    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        if( CSD->nAalp > 0 )
            pmm->Aalp[ CSD->xph[ii] ] = CNode->aPH[ii]/kg_to_g;
        pmm->Falp[ CSD->xph[ii] ] = CNode->omPH[ii];
    }

    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        pmm->XF[ CSD->xph[ii] ] = CNode->xPH[ii];
        //pmm->YF[ CSD->xph[ii] ] = CNode->xPH[ii];
    }

    if( !uPrimalSol )
    {
        //  Using primal solution retained in the MULTI structure instead -
        ; // the primal solution data from the DATABR structure are not unpacked

        //   pmm->IT = 0;
    }
    else {   // Unpacking primal solution provided in the node DATABR structure
        pmm->IT = CNode->IterDone; // ?  pmm->ITF+pmm->IT;
        //pmm->IT = 0;
        pmm->MBX = CNode->Ms;
        pmm->IC = CNode->IC;
        pmm->Eh = CNode->Eh;
        for( ii=0; ii<CSD->nDCb; ii++ )
            /*    pmm->X[ CSD->xdc[ii] ] = */
            pmm->Y[ CSD->xdc[ii] ] = CNode->xDC[ii];

        for( ii=0; ii<CSD->nPSb; ii++ )
        {
            //      pmm->FVOL[ CSD->xph[ii] ] = CNode->vPS[ii]*m3_to_cm3;
            //      pmm->FWGT[ CSD->xph[ii] ] = CNode->mPS[ii]*kg_to_g;
            pmm->PUL[ CSD->xph[ii] ] = CNode->amru[ii];
            pmm->PLL[ CSD->xph[ii] ] = CNode->amrl[ii];
        }

        for( ii=0; ii<CSD->nPHb; ii++ )
        {
            pmm->XF[ CSD->xph[ii] ] =
                    pmm->YF[ CSD->xph[ii] ] = this->Ph_Mole(ii);
            pmm->FVOL[ CSD->xph[ii] ] = this->Ph_Volume(ii)*m3_to_cm3;
            pmm->FWGT[ CSD->xph[ii] ] = this->Ph_Mass(ii)*kg_to_g;
        }

        for( long int k=0; k<CSD->nPSb; k++ )
            for(long int i=0; i<CSD->nICb; i++ )
            { long int dbr_ndx= (k*CSD->nICb)+i,
                        mul_ndx = ( CSD->xph[k]*CSD->nIC )+ CSD->xic[i];
                pmm->BF[ mul_ndx ] = CNode->bPS[dbr_ndx];
            }

        for( ii=0; ii<CSD->nPSb; ii++ )
            pmm->XFA[ CSD->xph[ii] ] = pmm->YFA[ CSD->xph[ii] ] = CNode->xPA[ii];

        for( ii=0; ii<CSD->nICb; ii++ )
            pmm->C[ CSD->xic[ii] ] = CNode->rMB[ii];
        for( ii=0; ii<CSD->nICb; ii++ )
            pmm->U[ CSD->xic[ii] ] = CNode->uIC[ii];

        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            pmm->Gamma[ CSD->xdc[ii] ] = CNode->gam[ii];
        }

        long int jb, je = 0;
        for( long int k=0; k<pmm->FIs; k++ )
        { // loop on solution phases
            jb = je;
            je += pmm->L1[ k ];
            // Load activity coeffs for phases-solutions
            for( ii=jb; ii<je; ii++ )
            {
                pmm->lnGam[ii] =  multi_ptr()->PhaseSpecificGamma( ii, jb, je, k, 1L );
            } // ii
        }
    }
    //  End
}


// (9) Optional, for passing the current mass transport time and time step into the work
// DATABR structure (for using it in TKinMet, or tracing/debugging, or in writing DBR files for nodes)
// This call should be used instead of obsolete GEM_set_MT() (provided below for compatibility with older codes)
//
void TNode::GEM_from_MT_time(
        double p_Tm,      // actual total simulation time, s          +       -      -
        double p_dt       // actual time step, s                      +       -      -
        )
{
    CNode->Tm = p_Tm;
    CNode->dt = p_dt;
}

void TNode::GEM_set_MT( // misleading name of the method - use GEM_from_MT_time() instead, see above
                        double p_Tm,      // actual total simulation time, s          +       -      -
                        double p_dt       // actual time step, s                      +       -      -
                        )
{
    CNode->Tm = p_Tm;
    CNode->dt = p_dt;
}

// (7a) Optional, to check if the time step in the work DATABR structure was o.k. for TKinMet calculations,
//  compared with the time step p_dt given before the GEM calculation. Checks the criteria for the validity
//  of time step. If time step was acceptable by a TKinMet model used, returns the actual time step after
//  copying (changed) AMRs into p_dul and p_dll vectors, as well as (changed) specific surface areas of
//  some (kinetically controlled) phases. Otherwise, returns a (smaller) suggested time step, while the
//  p_dul, p_pll, and p_asPH vectors remain unchanged.
//  Returns 0 or a negative number (unchanged p_dul and p_dll), if TKinMet calculations failed.
//
double GEM_to_MT_time(
        double p_dt,       ///< Actual time step, s                                     -       -     (+)   (+)
        double *p_dul,    ///< Upper AMR restrictions to amounts of DC [nDCb]          -       -      +     -
        double *p_dll     ///< Lower AMR restrictions to amounts of DC [nDCb]          -       -      +     -
        );

// (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
//  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
//   chemical data for all nodes used in the mass-transport model.
void TNode::GEM_restore_MT(
        long int  &p_NodeHandle,   // Node identification handle
        long int  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
        //                                    GEM input output  FMT control
        double &p_TK,      // Temperature T, Kelvin                       +       -      -
        double &p_P,      // Pressure P,  Pa                              +       -      -
        double &p_Vs,     // Volume V of reactive subsystem,  m3         (+)      -      +
        double &p_Ms,     // Mass of reactive subsystem, kg               -       -      +
        double *p_bIC,    // Bulk mole amounts of IC  [nICb]              +       -      -
        double *p_dul,    // Upper restrictions to amounts of DC [nDCb]   +       -      -
        double *p_dll,    // Lower restrictions to amounts of DC [nDCb]   +       -      -
        double *p_asPH    // Specific surface areas of phases,m2/kg[nPHb] +       -      -
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_TK = CNode->TK;
    p_P = CNode->P;
    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
        p_bIC[ii] = CNode->bIC[ii];
    for( ii=0; ii<CSD->nDCb; ii++ )
    {  p_dul[ii] = CNode->dul[ii];
        p_dll[ii] = CNode->dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            p_asPH[ii] = CNode->aPH[ii];
}

// (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
//  This call is useful after the GEM_init() (1) and GEM_run() (2) calls to initialize the arrays which keep the
//   chemical data for all nodes used in the mass-transport model.
void TNode::GEM_restore_MT(
        long int  &p_NodeHandle,   // Node identification handle
        long int  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
        //                                    GEM input output  FMT control
        double &p_TK,      // Temperature T, Kelvin                       +       -      -
        double &p_P,      // Pressure P,  Pa                              +       -      -
        double &p_Vs,     // Volume V of reactive subsystem,  m3         (+)      -      +
        double &p_Ms,     // Mass of reactive subsystem, kg               -       -      +
        double *p_bIC,    // Bulk mole amounts of IC  [nICb]              +       -      -
        double *p_dul,    // Upper restrictions to amounts of DC [nDCb]   +       -      -
        double *p_dll,    // Lower restrictions to amounts of DC [nDCb]   +       -      -
        double *p_asPH,   // Specific surface areas of phases,m2/kg[nPHb] +       -      -
        double *p_omPH,  // Stability indices of phases,log10 scale [nPHb] (+)     +      -
        double *p_amru,   // Upper AMRs to masses of sol. phases [nPSb]   +       -      -
        double *p_amrl    // Lower AMRs to masses of sol. phases [nPSb]   +       -      -
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_TK = CNode->TK;
    p_P = CNode->P;
    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
        p_bIC[ii] = CNode->bIC[ii];
    for( ii=0; ii<CSD->nDCb; ii++ )
    {  p_dul[ii] = CNode->dul[ii];
        p_dll[ii] = CNode->dll[ii];
    }
    if( CSD->nAalp >0 )
    {  for( ii=0; ii<CSD->nPHb; ii++ )
            p_asPH[ii] = CNode->aPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        p_omPH[ii] = CNode->omPH[ii];
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        p_amru[ii] = CNode->amru[ii];
        p_amrl[ii] = CNode->amrl[ii];
    }
}

// (7)  Retrieves the GEMIPM2 chemical speciation calculation results from the work DATABR structure instance
//   into memory provided by the mass transport part. Dimensions and order of elements in the arrays must correspond
//   to those in currently existing DATACH memory structure.
void TNode::GEM_to_MT(
        long int &p_NodeHandle,    // Node identification handle
        long int &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
        long int &p_IterDone,      // Number of iterations performed in the last GEM IPM calculation
        //                                                  GEM input output  FMT control
        // Chemical scalar variables
        double &p_Vs,    // Total volume V of reactive subsystem at given P,T, m3    -      -      +     +
        double &p_Ms,    // Total mass of the reactive subsystem, kg                 -      -      +     +
        double &p_Gs,    // Total Gibbs energy of the reactive subsystem, J          -      -      +     +
        double &p_Hs,    // Total enthalpy of reactive subsystem, J (reserved)       -      -      +     +
        double &p_IC,    // Effective aqueous ionic strength, molal                  -      -      +     +
        double &p_pH,    // pH of aqueous solution                                   -      -      +     +
        double &p_pe,    // pe of aqueous solution                                   -      -      +     +
        double &p_Eh,    // Eh of aqueous solution, V                                -      -      +     +
        // Dynamic data - dimensions see in DATACH.H structure
        double  *p_rMB,  // Mole balance residuals for Independent Components [nICb] -      -       +     +
        double  *p_uIC,  // Dual solution: IC chemical potentials, mol/mol [nICb]    -      -       +     +
        double  *p_xDC,  // Primal solution: DC mole amounts  [nDCb]                 -      -       +     +
        double  *p_gam,  // External activity coefficients of DC [nDCb]              -      -       +     +
        double  *p_xPH,  // Total mole amounts of all phases [nPHb]                  -      -       +     +
        double  *p_vPS,  // Total volumes of multicomponent phases, m3   [nPSb]      -      -       +     +
        double  *p_mPS,  // Total mass of multicomponent phase (carrier),kg [nPSb]   -      -       +     +
        double  *p_bPS,  // Bulk compositions of phases  [nPSb][nICb]                -      -       +     +
        double  *p_xPA,  //Amount of carrier in a multicomponent asymmetric phase[nPSb]-    -       +     +
        double  *p_aPH,  //surface area for phases, m2                           [nPHb]-    -       +     +
        double  *p_bSP   //Bulk composition of all solids, moles [nICb]                -    -       +     +
        )
{
    long int ii;
    p_NodeHandle = CNode->NodeHandle;
    p_NodeStatusCH = CNode->NodeStatusCH;
    p_IterDone = CNode->IterDone;

    p_Vs = CNode->Vs;
    p_Ms = CNode->Ms;
    p_Gs = CNode->Gs;
    p_Hs = CNode->Hs;
    p_IC = CNode->IC;
    p_pH = CNode->pH;
    p_pe = CNode->pe;
    p_Eh = CNode->Eh;

    for( ii=0; ii<CSD->nICb; ii++ )
    {
        p_rMB[ii] = CNode->rMB[ii];
        p_uIC[ii] = CNode->uIC[ii];
        p_bSP[ii] = CNode->bSP[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        p_xDC[ii] = CNode->xDC[ii];
        p_gam[ii] = CNode->gam[ii];
    }
    for( ii=0; ii<CSD->nPHb; ii++ )
    {
        p_xPH[ii] = CNode->xPH[ii];
        p_aPH[ii] = CNode->aPH[ii]*Ph_Mass( ii );  // correction 9.10.2013 DK
    }
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        p_vPS[ii] = CNode->vPS[ii];
        p_mPS[ii] = CNode->mPS[ii];
        p_xPA[ii] = CNode->xPA[ii];
    }
    for( ii=0; ii<CSD->nPSb*CSD->nICb; ii++ )
        p_bPS[ii] = CNode->bPS[ii];
}

// (8) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH   // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
}

//(8a) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//This overloaded variant uses the xDC speciation vector for setting the
// new bulk chemical composition to be used in the next GEM_run() calculation.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,  // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
        double *p_xDC    // Mole amounts of DCs [nDCb] - will be convoluted
        // and added to the bIC GEM input vector (if full speciation
        // and not just increments then p_bIC vector must be zeroed off -
        // it will be calculated from p_xDC and stoichiometry matrix A
        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-simplex IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is ordered, leave if simplex is ordered (KD)

    // Optional part - convolution of xDC vector into bIC vector
    if( p_xDC )
    {  long int jj;
        // Correction of bIC vector by convoluting the amounts of DCs
        for( jj=0; jj<CSD->nDCb; jj++ )
            if( p_xDC[jj] > 0.0 )
                for( ii=0; ii<CSD->nICb; ii++ )
                    CNode->bIC[ii] += p_xDC[jj] * DCaJI( jj, ii );
    }
}

//(8b) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
//  structure contains a full list of Dependent Components used in GEM IPM2 calculations.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,   // Specific surface areas of phases, m2/kg [nPHb]  +       -      -
        double *p_xDC,  // Mole amounts of DCs [nDCb] - old primal soln.     +       -      -
        double *p_gam   // DC activity coefficients [nDCb] - old primal s.   +       -      -
        )
{
    long int ii;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
        for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii];

    // Optional part - copying old primal solution from p_xDC and p_gam vectors
    if( p_xDC && p_gam )
    {
        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            CNode->xDC[ii] = p_xDC[ii];
            CNode->gam[ii] = p_gam[ii];
        }
    }
    else if( CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal solution provided!

    //  Discuss the policy!
    //   if( p_xDC )
    //   {  long int jj;
    //      // Correction of bIC vector by convoluting the amounts of DCs
    //      for( jj=0; jj<CSD->nDCb; jj++ )
    //        if( p_xDC[jj] )
    //          for( ii=0; ii<CSD->nICb; ii++ )
    //            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
    //   }

}


// (8c) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(long int  p_NodeHandle,   // Node identification handle
                        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                        //                                              GEM input output  FMT control
                        double p_TK,     // Temperature T, Kelvin                            +       -      -
                        double p_P,      // Pressure P, Pa                                   +       -      -
                        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
                        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_dll    // Lower restrictions to amounts of DC [nDCb]       +       -      -
                        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
}

// (8d) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//     This call is usually preceeding the GEM_run() call
void TNode::GEM_from_MT(long int  p_NodeHandle,   // Node identification handle
                        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
                        //                                              GEM input output  FMT control
                        double p_TK,     // Temperature T, Kelvin                            +       -      -
                        double p_P,      // Pressure P, Pa                                   +       -      -
                        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
                        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
                        double *p_asPH,  // Specific surface areas of phases, m2/kg [nPHb]   +       -      -
                        double *p_omPH,   // Stability indices of phases,log10 scale [nPHb]  (+)      +      -
                        double *p_amru,   // Upper AMR to masses of sol. phases [nPSb]        +       -      -
                        double *p_amrl    // Lower AMR to masses of sol. phases [nPSb]        +       -      -
                        )
{
    long int ii;
    bool useSimplex = false;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {  //  SD 11/02/05 for test
        //if( fabs(CNode->bIC[ii] - p_bIC[ii] ) > CNode->bIC[ii]*1e-4 ) // bugfix KD 21.11.04
        //     useSimplex = true;
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( useSimplex && CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;
    // Switch only if SIA is selected, leave if LPP AIA is prescribed (KD)
    if( CSD->nAalp >0 )
    { for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        CNode->omPH[ii] = p_omPH[ii];
    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        CNode->amru[ii] = p_amru[ii];
        CNode->amrl[ii] = p_amrl[ii];
    }
}

//(8e) Loads the GEMS3K input data for a given mass-transport node into the work instance of DATABR structure.
//In addition, provides access to speciation vector p_xDC and DC activity coefficients p_gam that will be used in
// GEM "smart initial approximation" SIA mode if dBR->NodeStatusCH == NEED_GEM_SIA (5) and
// uPrimalSol = true are set for the GEM_run() call (see Section 2) . This works only when the DATACH
//  structure contains a full list of Dependent Components used in GEM IPM2 calculations.
void TNode::GEM_from_MT(
        long int  p_NodeHandle,   // Node identification handle
        long int  p_NodeStatusCH, // Node status code (NEED_GEM_SIA or NEED_GEM_AIA)
        //                                              GEM input output  FMT control
        double p_TK,     // Temperature T, Kelvin                            +       -      -
        double p_P,      // Pressure P, Pa                                   +       -      -
        double p_Vs,     // Volume V of reactive subsystem, m3               -       -      +
        double p_Ms,     // Mass of reactive subsystem, kg                   -       -      +
        double *p_bIC,   // Bulk mole amounts of IC [nICb]                   +       -      -
        double *p_dul,   // Upper restrictions to amounts of DC [nDCb]       +       -      -
        double *p_dll,   // Lower restrictions to amounts of DC [nDCb]       +       -      -
        double *p_asPH,   // Specific surface areas of phases, m2/kg [nPHb]  +       -      -
        double *p_omPH,   // Stability indices of phases,log10 scale [nPHb]  (+)      +      -
        double *p_amru,   //< Upper AMR to masses of sol. phases [nPSb]      +       -      -
        double *p_amrl,   //< Lower AMR to masses of sol. phases [nPSb]      +       -      -
        double *p_xDC,  // Mole amounts of DCs [nDCb] - old primal soln.     +       -      -
        double *p_gam   // DC activity coefficients [nDCb] - old primal s.   +       -      -
        )
{
    long int ii;

    CNode->NodeHandle = p_NodeHandle;
    CNode->NodeStatusCH = p_NodeStatusCH;
    CNode->TK = p_TK;
    CNode->P = p_P;
    CNode->Vs = p_Vs;
    CNode->Ms = p_Ms;
    // Checking if no-LPP IA is Ok
    for( ii=0; ii<CSD->nICb; ii++ )
    {
        CNode->bIC[ii] = p_bIC[ii];
    }
    for( ii=0; ii<CSD->nDCb; ii++ )
    {
        CNode->dul[ii] = p_dul[ii];
        CNode->dll[ii] = p_dll[ii];
    }
    if( CSD->nAalp >0 )
    { for( ii=0; ii<CSD->nPHb; ii++ )
            CNode->aPH[ii] = p_asPH[ii]; }
    for( ii=0; ii<CSD->nPHb; ii++ )
        CNode->omPH[ii] = p_omPH[ii];
    // Optional part - copying old primal solution from p_xDC and p_gam vectors
    if( p_xDC && p_gam )
    {
        for( ii=0; ii<CSD->nDCb; ii++ )
        {
            CNode->xDC[ii] = p_xDC[ii];
            CNode->gam[ii] = p_gam[ii];
        }
    }
    else if( CNode->NodeStatusCH == NEED_GEM_SIA )
        CNode->NodeStatusCH = NEED_GEM_AIA;   // no complete old primal solution provided!

    for( ii=0; ii<CSD->nPSb; ii++ )
    {
        CNode->amru[ii] = p_amru[ii];
        CNode->amrl[ii] = p_amrl[ii];
    }
    //  Discuss the policy!
    //   if( p_xDC )
    //   {  long int jj;
    //      // Correction of bIC vector by convoluting the amounts of DCs
    //      for( jj=0; jj<CSD->nDCb; jj++ )
    //        if( p_xDC[jj] )
    //          for( ii=0; ii<CSD->nICb; ii++ )
    //            CNode->bIC[ii] += p_xDC[jj] * nodeCH_A( jj, ii );
    //   }

}

//-----------------------End of node.cpp--------------------------



