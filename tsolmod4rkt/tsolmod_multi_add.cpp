//-------------------------------------------------------------------
// $Id$
//
/// \file tsolmod_multi_add.cpp
/// Addition functions from TMultiBase class
//
// Copyright (c) 2023 S.Dmytriyeva,D.Kulik
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

#include "v_detail.h"
#include "io_template.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "io_keyvalue.h"
#include "tsolmod_multi.h"
#include "datach_api.h"
#include "s_solmod.h"
#include "num_methods.h"

// Conversion factors
const double bar_to_Pa = 1e5,
m3_to_cm3 = 1e6,
kg_to_g = 1e3;


const BASE_PARAM pa_p_ =
        {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     30000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
}; // BASE_PARAM

// New --------------------------------------------------------------------

// Constructor of the class instance in memory for standalone GEMS3K or coupled program
TSolModMulti::TSolModMulti()
{
    pa_standalone.reset( new BASE_PARAM() );
    *pa_standalone = pa_p_;

    pmp = &pm;
    sizeN = 0;
    AA = nullptr;
    BB = nullptr;
    arrL = nullptr;
    arrAN = nullptr;

    U_mean = nullptr;
    U_M2 = nullptr;
    U_CVo = nullptr;
    U_CV = nullptr;
    ICNud = nullptr;
    pm.errorCode[0] ='\0';
    pm.errorBuf[0] ='\0';

    sizeFIs = 0;
    phSolMod = nullptr;

    pmp->Guns = nullptr;
    pmp->Vuns = nullptr;
    pmp->tpp_G = nullptr;
    pmp->tpp_S = nullptr;
    pmp->tpp_Vm = nullptr;
    set_def();

    // from node
    CSD = NULL;
    CNode = NULL;
    allocMemory();
    load_thermodynamic_data = false;
}

TSolModMulti::~TSolModMulti()
{
    clear_ThermoEngine();
    freeMemory();
}

void TSolModMulti::allocMemory()
{
    // memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;
    // mem_set( CSD, 0, sizeof(DATACH) );
    dbr_dch_api::datach_reset(CSD);
    // mem_set( CNode, 0, sizeof(DATABR) );
    dbr_dch_api::databr_reset(CNode, 2);
}

void TSolModMulti::freeMemory()
{
    dbr_dch_api::databr_free_internal(CNode);
    delete CNode;
    CNode = nullptr;
    dbr_dch_api::datach_free(CSD);
    delete CSD;
    CSD = nullptr;
}

// (1) Initialization of GEM IPM2 data structures in coupled RMT-GEM programs
//  that use GEMS3K module. Also reads in the IPM, DCH and DBR text input files
//  in key-value, json or binary format. Parameters:
//  ipmfiles_lst_name - name of a text file that contains:
//    " -f | -j | -t |-b <DCH_DAT file name> <IPM_DAT file name> <dataBR file name>
//  dbfiles_lst_name - name of a text file that contains:
//    <dataBR  file name1>, ... , <dataBR file nameN> "
//    These files (one DCH_DAT, one IPM_DAT, and at least one dataBR file) must
//    exist in the same directory where the ipmfiles_lst_name file is located.
//    the DBR_DAT files in the above list are indexed as 1, 2, ... N (node handles)
//    and must contain valid initial chemical systems (of the same structure
//    as described in the DCH_DAT file) to set up the initial state of the FMT
//    node array.
//  If -t flag or no flag is specified then all data files must be in key-value text
//    (ASCII) format (and file names must have .dat extension);
//  If -j and -f flag is specified then all data files must be in JSON format (and file names
//    must have .json extension);
//  if -b flag is specified then all data files are assumed to be binary (little-endian)
//    files.
//-------------------------------------------------------------------


long int  TSolModMulti::GEM_init( const char* ipmfiles_lst_name )
{
    clearipmLogError();
    clear_ThermoEngine();

    try
    {
        //  Syntax: -t/-b  "<DCH_DAT file name>"  "<IPM_DAT file name>"
        //       "<DBR_DAT file1 name>" [ ...  "<DBR_DAT fileN name>"]
        GEMS3KGenerator generator( ipmfiles_lst_name );
        current_output_set_name = current_input_set_name = generator.get_name();

        switch( generator.files_mode() )
        {
        case GEMS3KGenerator::f_binary:
        {
            Error( generator.get_dch_path(), "binary input files not implemented");
        }
            break;
        default:
        {
            std::fstream f_ch( generator.get_dch_path(), std::ios::in );
            ErrorIf( !f_ch.good() , generator.get_dch_path(), "DCH_DAT fileopen error");
            dbr_dch_api::read_dch_format_stream(current_input_set_name, CSD, f_ch, generator.files_mode() );

            if( generator.files_mode()>=GEMS3KGenerator::f_thermofun )
            {
                load_ThermoEngine(generator.get_thermofun_path());
            }

            std::fstream ff( generator.get_ipm_path(), std::ios::in );
            ErrorIf( !ff.good() , generator.get_ipm_path(), "Fileopen error");
            read_ipm_format_stream( ff,generator.files_mode(), CSD, current_input_set_name);

        }
            break;
        }

        // copy intervals for minimization
        pmp->Pai[0] = CSD->Pval[0]/bar_to_Pa;
        pmp->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
        pmp->Pai[2] = getStep( pmp->Pai, CSD->nPp )/bar_to_Pa;//(pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
        pmp->Pai[3] = CSD->Ptol/bar_to_Pa;

        pmp->Tai[0] = CSD->TKval[0]-C_to_K;
        pmp->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
        pmp->Tai[2] = getStep( pmp->Tai, CSD->nTp );//(pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
        pmp->Tai[3] = CSD->Ttol;

        pmp->Fdev1[0] = 0.;
        pmp->Fdev1[1] = 1e-6;   // 24/05/2010 must be copied from GEMS3 structure
        pmp->Fdev2[0] = 0.;
        pmp->Fdev2[1] = 1e-6;

        // Reading DBR_DAT file into work DATABR structure from ipmfiles_lst_name
        std::string dbr_file = generator.get_dbr_path( 0 );
        ErrorIf( dbr_file.empty() , ipmfiles_lst_name, " Undefined DBR_DAT file name");
        std::fstream in_br( dbr_file, std::ios::in );
        ErrorIf( !in_br.good() , dbr_file, "DBR_DAT fileopen error");
        dbr_dch_api::read_dbr_format_stream(current_input_set_name, CSD, CNode, in_br, generator.files_mode());

        // Creating and initializing the TActivity class instance for this TNode instance
        /// init_into_gems3k();
        TSolMod::solmod_logger->info("Initialization of system {}", char_array_to_string(pmp->stkey,EQ_RKLEN));
        return 0;
    }
    catch(TError& err)
    {
        ipmlog_error = err.title + std::string(": ") + err.mess;
    }
    catch(std::exception& e)
    {
        ipmlog_error = std::string("std::exception: ") + e.what();
    }
    catch(...)
    {
        ipmlog_error = "unknown exception";
        return -1;
    }

    if( ipmfiles_lst_name ) {
        TSolMod::solmod_logger->error("GEMS3K input : file {}", ipmfiles_lst_name);
    }
    if( !ipmlog_error.empty() ) {
        TSolMod::solmod_logger->error("GEM_init error: {}", ipmlog_error);
    }
    return 1;
}


// Changed ----------------------------------------------------------------

void TSolModMulti::clear_ThermoEngine()
{
#ifdef USE_THERMOFUN
    // clear previous
    thermo_engine.reset();
    thermo_json_string="";
#endif
}

bool TSolModMulti::load_ThermoEngine(const std::string &thermo_file_or_string)
{
    if(thermo_file_or_string.find("\"elements\"")!=std::string::npos) {
        // input string
        thermo_json_string = thermo_file_or_string;
    }
    else {
        // input file
        std::ifstream f_fun(thermo_file_or_string);
        ErrorIf(!f_fun.good(), thermo_file_or_string, "ThermoFun JSON format file fileopen error");
        std::stringstream buffer;
        buffer << f_fun.rdbuf();
        thermo_json_string = buffer.str();
    }

#ifdef USE_THERMOFUN
    thermo_engine.reset(new ThermoFun::ThermoEngine(thermo_file_or_string));
    TSolMod::solmod_logger->trace("Read ThermoEngine: {}", thermo_file_or_string);
    return true;
#else
    TSolMod::solmod_logger->warn("Try read ThermoEngine not in USE_THERMOFUN mode {}", thermo_file_or_string);
    return false;
#endif
}


// Copy -------------------------------------------------------------------

/// Converting DC class codes into generic internal codes of IPM
void TSolModMulti::ConvertDCC()
{
    long int i, j, k, iRet=0;
    char DCCW;

    j=0;
    for( k=0; k< pm.FI; k++ )
    { // phase loop
        i=j+pm.L1[k];
        if( pm.L1[k] == 1 )
        {
            pm.DCCW[j] = DC_SINGLE;
            goto NEXT_PHASE;
        }
        for( ; j<i; j++ )
        { // DC loop
            switch( pm.DCC[j] ) // select v_j expression
            {
            case DC_SCP_CONDEN:
                DCCW = DC_SINGLE;
                break;
            case DC_GAS_COMP:
            case DC_GAS_H2O:
            case DC_GAS_CO2:
            case DC_GAS_H2:
            case DC_GAS_N2:
            case DC_SOL_IDEAL:
            case DC_SOL_MINOR:
            case DC_SOL_MAJOR:
            case DC_SOL_MINDEP:
            case DC_SOL_MAJDEP:
            case DC_SCM_SPECIES:
                DCCW = DC_SYMMETRIC;
                break;
            case DC_AQ_PROTON:
            case DC_AQ_ELECTRON:
            case DC_AQ_SPECIES:
            case DC_AQ_SURCOMP:
                DCCW = DC_ASYM_SPECIES;
                break;
            case DC_AQ_SOLVCOM:
            case DC_AQ_SOLVENT:
                DCCW = DC_ASYM_CARRIER;
                break;
            case DC_IESC_A:
            case DC_IEWC_B:
                DCCW = DC_ASYM_SPECIES;
                break;
                // Remapping
            case DC_SUR_GROUP:
            case DC_SUR_COMPLEX:
                DCCW = DC_ASYM_SPECIES;
                pm.DCC[j] = DC_SSC_A0;
                break;
            case DC_SUR_IPAIR:
                DCCW = DC_ASYM_SPECIES;
                pm.DCC[j] = DC_WSC_A0;
                break;
            case DC_SUR_MINAL:
            case DC_SUR_CARRIER:
            case DC_PEL_CARRIER:
                 DCCW = DC_ASYM_CARRIER;
                break;
            default:
                if( isdigit( pm.DCC[j] ))
                {
                    if( pm.PHC[k] == PH_SORPTION ||
                            pm.PHC[k] == PH_POLYEL )
                    {
                        DCCW = DC_ASYM_SPECIES;
                        break;
                    }
                }
                DCCW = DC_SINGLE;
                iRet++;  // error the class code
            }
            pm.DCCW[j] = DCCW;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    ErrorIf( iRet>0, "E19IPM: ConvertDCC()", "Invalid DC class code. Memory corruption?");
}




// tsolmod_multi_add.cpp
