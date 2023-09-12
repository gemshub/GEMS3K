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
#include "jsonconfig.h"

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
        /// ?? Used only T and P
        std::string dbr_file = generator.get_dbr_path( 0 );
        ErrorIf( dbr_file.empty() , ipmfiles_lst_name, " Undefined DBR_DAT file name");
        std::fstream in_br( dbr_file, std::ios::in );
        ErrorIf( !in_br.good() , dbr_file, "DBR_DAT fileopen error");
        dbr_dch_api::read_dbr_format_stream(current_input_set_name, CSD, CNode, in_br, generator.files_mode());

        // Creating and initializing the TActivity class instance for this TNode instance
        init_into_gems3k();
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


bool TSolModMulti::load_all_thermodynamic_from_thermo(double TK, double PPa)
{
#ifdef USE_THERMOFUN
    if( !thermo_engine.get() )
        return false;
    try{
        TSolMod::solmod_logger->info("Calc ThermoEngine T: {}  P: {}", TK, PPa);
        long int j, jj, k, jb, je=0;
        double G0, P = PPa/bar_to_Pa;

        pmp->T = pmp->Tc = TK;
        pmp->TC = pmp->TCc = TK-C_to_K;
        // new API
        double funT = TK, funP=P*bar_to_Pa;   // T in K, P in Pa

        DATACH  *dCH = CSD;
        auto water_props = thermo_engine->propertiesSolvent(funT,funP, "H2O@");
        auto water_electro = thermo_engine->electroPropertiesSolvent(funT,funP, "H2O@");

        auto water_vapor = thermo_engine->database().getSubstance("H2O@");
        water_vapor.setMethod_P( ThermoFun::MethodCorrP_Thrift::type::CPM_GAS);

        ThermoFun::Database db2;
        db2.addSubstance(water_vapor);

        ThermoFun::ThermoEngine te2(db2);
        auto water_gas_props = te2.propertiesSolvent(funT,funP, "H2O@");
        auto water_gas_electro = te2.electroPropertiesSolvent(funT,funP, "H2O@");

        pmp->Pc = P;
        if( P < 1e-5 )
        { // Pressure at saturated H2O vapour at given temperature
            P = funP/bar_to_Pa;
            //        long int xT =check_grid_T(TK);
            //        if(xT>= 0)
            //            P = dCH->Psat[xT]/bar_to_Pa;
            //        else
            //            P =  LagranInterp( &PPa, dCH->TKval, dCH->Psat, PPa, TK, dCH->nTp, 1,6 )/bar_to_Pa;
        }

        pmp->P = P;
        pmp->RT = R_CONSTANT * pmp->Tc;
        pmp->FRT = F_CONSTANT/pmp->RT;
        pmp->lnP = log( P );

        pmp->denW[0] = water_props.density.val/1e3;
        pmp->epsW[0] = water_electro.epsilon.val;
        pmp->denW[1] = water_props.densityT.val/1e3;
        pmp->epsW[1] = water_electro.epsilonT.val;
        pmp->denW[2] = water_props.densityTT.val/1e3;
        pmp->epsW[2] = water_electro.epsilonTT.val;
        pmp->denW[3] = water_props.densityP.val*1e2; // /1e3;
        pmp->epsW[3] = water_electro.epsilonP.val;
        pmp->denW[4] = water_props.densityPP.val/1e3;
        pmp->epsW[4] = water_electro.epsilonPP.val;

        pmp->denWg[0] = water_gas_props.density.val/1e3;
        pmp->epsWg[0] = water_gas_electro.epsilon.val;
        pmp->denWg[1] = water_gas_props.densityT.val/1e3;
        pmp->epsWg[1] = water_gas_electro.epsilonT.val;
        pmp->denWg[2] = water_gas_props.densityTT.val/1e3;
        pmp->epsWg[2] = water_gas_electro.epsilonTT.val;
        pmp->denWg[3] = water_gas_props.densityP.val*1e2; // /1e3;
        pmp->epsWg[3] = water_gas_electro.epsilonP.val;
        pmp->denWg[4] = water_gas_props.densityPP.val/1e3;
        pmp->epsWg[4] = water_gas_electro.epsilonPP.val;

#ifdef  USE_THERMO_LOG
        std::fstream f_log;
        if(GemsSettings::log_thermodynamic) {
            f_log.open(GemsSettings::with_directory("thermodynamic-log.csv"), std::ios::out/*|std::ios::app*/ );
            f_log << "\nCalc ThermoEngine;T;" << TK << ";P;" << PPa << "\n";
            f_log << "denW";
            for( jj=0; jj<5; jj++)
                f_log << ";" << floating_point_to_string(pmp->denW[jj]);
            f_log << "\nepsW";
            for( jj=0; jj<5; jj++)
                f_log << ";" << floating_point_to_string(pmp->epsW[jj]);
            f_log << "\ndenWg";
            for( jj=0; jj<5; jj++)
                f_log << ";" << floating_point_to_string(pmp->denWg[jj]);
            f_log << "\nepsWg";
            for( jj=0; jj<5; jj++)
                f_log << ";" << floating_point_to_string(pmp->epsWg[jj]);
        }
#endif

        for( k=0; k<pmp->FI; k++ )
        {
            jb = je;
            je += pmp->L1[k];
            // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
            // depending on the presence of these arrays in DATACH and Multi structures
            for( j=jb; j<je; j++ )
            {
                std::string symbol = std::string(CSD->DCNL[j], 0, MaxDCN);
                auto propAl    = thermo_engine->thermoPropertiesSubstance(funT,funP, symbol);

                G0 = propAl.gibbs_energy.val;
                pmp->Vol[j] = propAl.volume.val*10;
                if( dCH->S0 ) pmp->S0[j] = propAl.entropy.val;
                if( dCH->H0 ) pmp->H0[j] = propAl.enthalpy.val;
                if( dCH->Cp0 ) pmp->Cp0[j] = propAl.heat_capacity_cp.val;
                if( dCH->A0 ) pmp->A0[j] = propAl.helmholtz_energy.val;
                if( dCH->U0 ) pmp->U0[j] = propAl.internal_energy.val;

                pmp->G0[j] = ConvertGj_toUniformStandardState(G0, j, k);
#ifdef  USE_THERMO_LOG
                if(GemsSettings::log_thermodynamic) {
                    f_log << "\n" << symbol << ";" << floating_point_to_string(G0)
                          << ";" << floating_point_to_string(pmp->G0[j])
                          << ";" << floating_point_to_string(pmp->Vol[j]);
                    if( dCH->S0 ) f_log << ";" << floating_point_to_string(pmp->S0[j]);
                    if( dCH->H0 ) f_log << ";" << floating_point_to_string(pmp->H0[j]);
                    if( dCH->Cp0 ) f_log << ";" << floating_point_to_string(pmp->Cp0[j]);
                    if( dCH->A0 ) f_log << ";" << floating_point_to_string(pmp->A0[j]);
                    if( dCH->U0 ) f_log << ";" << floating_point_to_string(pmp->U0[j]);
                }
#endif
            }  // j
        } // k
    }
    catch (const std::runtime_error& exception)
    {
        TSolMod::solmod_logger->error("ThermoEngine error: {}", exception.what());
        Error( "ThermoEngine error:", exception.what());
    }
    return true;
#else
    return false;
#endif

}

void TSolModMulti::init_into_gems3k()
{
    //InitReadActivities( mult_in.c_str(),CSD ); // from DCH file in future?
    InitalizeGEM_IPM_Data();              // In future, initialize data in TActivity also
    //this->InitCopyActivities( CSD, pmm, CNode );
}

// Before Calculations
/// Calculation by IPM (preparing for calculation, unpacking data) In IPM
void TSolModMulti::InitalizeGEM_IPM_Data() // Reset internal data formerly MultiInit()
{

    MultiConstInit();
    // initalizeGEM_IPM_Data_GUI(); empty

    if( pm.pTPD < 2 )
    {
        //  pm.pTPD == 0 => changed T or P;   pm.pTPD == 1 => changed system;
        LoadThermodynamicData(CNode->TK, CNode->P); // Loading thermodynamic data into MULTI structure
    }
    Alloc_internal();
    Alloc_uDD( pm.N );      // Added 06.05.2011 DK

    // calculate mass of the system
    pm.MBX = 0.0;
    for(int i=0; i<pm.N; i++ )
        pm.MBX += pm.B[i] * pm.Awt[i];
    pm.MBX /= 1000.;
    /**
    RescaleToSize( true );  // Added to set default cutoffs/inserts 30.08.2009 DK

    if(  pm.pNP )
    {  //  Smart Initial Approximation mode
        long int j,k;

        loadData( false );  // unpack SysEq record into MULTI

        bool AllPhasesPure = true;   // Added by DK on 09.03.2010
        // checking if all phases are pure
        for( k=0; k < pm.FI; k++ )
            if( pm.L1[k] > 1 )
                AllPhasesPure = false;

        for( j=0; j< pm.L; j++ )
            pm.X[j] = pm.Y[j];
        //       pm.IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateConcentrations( pm.X, pm.XF, pm.XFA);  // 13.03.2008  DK
        // test multicomponent phases and load data for mixing models
        if( pm.FIs && AllPhasesPure == false )
        {
            // Load activity coeffs for phases-solutions
            int jb, je=0;
            for( int kk=0; kk<pm.FIs; kk++ )
            { // loop on solution phases
                jb = je;
                je += pm.L1[kk];
                // Load activity coeffs for phases-solutions
                for( j=jb; j< je; j++ )
                {
                    pm.lnGmo[j] = pm.lnGam[j];
                    if( fabs( pm.lnGam[j] ) <= 84. )
                        //                pm.Gamma[j] = exp( pm.lnGam[j] );
                        pm.Gamma[j] = PhaseSpecificGamma( j, jb, je, kk, 0 );
                    else pm.Gamma[j] = 1.0;
                } // j
            }  // kk
        }
    }*/
}

/// Setup/copy flags and thresholds for numeric modules to TMulti structure.
/// Do it before calculations
void TSolModMulti::MultiConstInit() // from MultiRemake
{
    pm.FI1 = 0;
    pm.FI1s = 0;
    pm.FI1a = 0;
    pm.ITF = 0; pm.ITG = 0;
    pm.PD = base_param()->PD;
    pm.Ec = pm.K2 = pm.MK = 0;
    pm.W1 = 0;
    pm.is = 0;
    pm.js = 0;
    pm.next  = 0;
    pm.ln5551 = log( H2O_mol_to_kg );             // constant corrected 30.08.2008
    pm.lowPosNum = Min_phys_amount;               // = 1.66e-24 mol
    pm.logXw = -16.;
    pm.logYFk = -9.;
    pm.DXM = base_param()->DK;

    //  ???????
    pm.FX = 7777777.;
    if( pm.pH < -15. || pm.pH > 16.  )   // Check for trash in pH - bugfix 19.06.2013
        pm.pH = pm.Eh = pm.pe = 0.0;
    pm.YMET = 0;                      // always 0.0 ????
    pm.PCI = 1.0;
    pm.FitVar[4] = 1.0;

    // from multiConstInit_PN();
    pm.PZ = base_param()->DW;  // in IPM
    //  pm.FitVar[0] = 0.0640000030398369;
}

/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void TSolModMulti::LoadThermodynamicData(double TK, double PPa)
{
    // try generate thermodynamic data from ThermoEngine
    if( !load_all_thermodynamic_from_thermo( TK, PPa ))
    {
        load_all_thermodynamic_from_grid(TK, PPa);
    }
    pm.pTPD = 2;
}


/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void TSolModMulti::load_all_thermodynamic_from_grid(double TK, double PPa )
{
    std::string error_msg;
    long int j, jj, k, xTP, jb, je=0;
    double Go, Gg=0., Ge=0., Vv, h0=0., S0 = 0., Cp0= 0., a0 = 0., u0 = 0.;
    double P = PPa/bar_to_Pa;
    DATACH  *dCH = CSD;

    //    ipm_logger->info("Calc Lookup T: {}  P: {}", TK, PPa);
    error_msg = dbr_dch_api::check_TP(dCH, TK, PPa);
    if( dCH->nTp <1 || dCH->nPp <1 || !error_msg.empty() )
    {
        Error("load_all_thermodynamic_from_grid: ",
              ( error_msg.empty() ? "empty lookup arrays" : error_msg) );
        return;
    }

    pm.T = pm.Tc = TK;
    pm.TC = pm.TCc = TK-C_to_K;
    pm.Pc = P;
    if( P < 1e-5 )
    { // Pressure at saturated H2O vapour at given temperature
        long int xT = dbr_dch_api::check_grid_T(dCH, TK);
        if(xT>= 0)
            P = dCH->Psat[xT]/bar_to_Pa;
        else
            P =  LagranInterp( &PPa, dCH->TKval, dCH->Psat, PPa, TK, dCH->nTp, 1,6 )/bar_to_Pa;
    }
    pm.P = P;
    pm.RT = R_CONSTANT * pm.Tc;
    pm.FRT = F_CONSTANT/pm.RT;
    pm.lnP = log( P );

    xTP = dbr_dch_api::check_grid_TP(dCH, TK, PPa);

    for( k=0; k<5; k++ )
    {
        jj =  k * dbr_dch_api::gridTP(dCH);
        if( xTP >= 0 )
        {
            pm.denW[k] = dCH->denW[jj+xTP]/1e3;
            pm.epsW[k] = dCH->epsW[jj+xTP];
            pm.denWg[k] = dCH->denWg[jj+xTP]/1e3;
            pm.epsWg[k] = dCH->epsWg[jj+xTP];
        }
        else
        {
            pm.denW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,6 )/1e3;// from test denW enough
            pm.epsW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,5 );// from test epsW enough
            pm.denWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 )/1e3;
            pm.epsWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 );
        }
    }

#ifdef  USE_THERMO_LOG
    std::fstream f_log;
    if(GemsSettings::log_thermodynamic) {
        f_log.open(GemsSettings::with_directory("thermodynamic-log-lookup.csv"), std::ios::out/*|std::ios::app*/ );
        f_log << "\nCalc ThermoEngine;T;" << TK << ";P;" << PPa << "\n";
        f_log << "denW";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pm.denW[jj]);
        f_log << "\nepsW";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pm.epsW[jj]);
        f_log << "\ndenWg";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pm.denWg[jj]);
        f_log << "\nepsWg";
        for( jj=0; jj<5; jj++)
            f_log << ";" << floating_point_to_string(pm.epsWg[jj]);
    }
#endif
    long int xVol =  getXvolume();

    for( k=0; k<pm.FI; k++ )
    {
        jb = je;
        je += pm.L1[k];
        // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
        // depending on the presence of these arrays in DATACH and Multi structures
        for( j=jb; j<je; j++ )
        {
            jj =  j * dbr_dch_api::gridTP(dCH);
            if( xTP >= 0 )
            {
                Go = dCH->G0[ jj+xTP];
                Vv = dCH->V0[ jj+xTP]*1e5;
                if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
                if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
                if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
                if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
                if( dCH->U0 ) h0 = dCH->U0[ jj+xTP];
            }
            else
            {
                Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
                Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 5 )*1e5;
                if( dCH->S0 ) S0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->S0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp, 4 ); // from test S0[Ca+2] enough
                if( dCH->H0 ) h0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->H0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->Cp0 ) Cp0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->Cp0+jj,
                                                    PPa, TK, dCH->nTp, dCH->nPp, 3 ); // from test Cp0[Ca+2] not more
                if( dCH->A0 ) a0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->A0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->U0 ) u0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->U0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
            }
            if( pm.tpp_G )
                pm.tpp_G[j] = Go;
            if( pm.Guns )
                Gg = pm.Guns[j];
            else
                Gg = 0.;

            Ge = 0.;
            pm.G0[j] = ConvertGj_toUniformStandardState( Go+Gg+Ge, j, k ); // formerly Cj_init_calc()
            // Inside this function, pm.YOF[k] can be added!

                switch( pm.PV )
                { // put molar volumes of components into A matrix or into the vector of molar volumes
                // to be checked!
                case VOL_CONSTR:
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
                    if( xVol >= 0. )
                        pm.A[j*pm.N+xVol] = Vv;
                    [[fallthrough]];
                case VOL_CALC:
                case VOL_UNDEF:
                    if( pm.tpp_Vm )
                        pm.tpp_Vm[j] = Vv;
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
                    pm.Vol[j] = Vv  * 10.;
                    break;
                }
            if( pm.S0 ) pm.S0[j] = S0;
            if( pm.H0 ) pm.H0[j] = h0;
            if( pm.Cp0 ) pm.Cp0[j] = Cp0;
            if( pm.A0 ) pm.A0[j] = a0;
            if( pm.U0 ) pm.U0[j] = u0;

#ifdef  USE_THERMO_LOG
            if(GemsSettings::log_thermodynamic) {
                f_log << "\n" << std::string(dCH->DCNL[j], 0, MaxDCN) << ";" << floating_point_to_string(Go)
                      << ";" << floating_point_to_string(pm.G0[j])
                      << ";" << floating_point_to_string(pm.Vol[j]);
                if( dCH->S0 ) f_log << ";" << floating_point_to_string(pm.S0[j]);
                if( dCH->H0 ) f_log << ";" << floating_point_to_string(pm.H0[j]);
                if( dCH->Cp0 ) f_log << ";" << floating_point_to_string(pm.Cp0[j]);
                if( dCH->A0 ) f_log << ";" << floating_point_to_string(pm.A0[j]);
                if( dCH->U0 ) f_log << ";" << floating_point_to_string(pm.U0[j]);
            }
#endif
        }  // j
    } // k
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

/// Conversion of g(T,P) value for DCs into the uniform cj scale.
/// \param k - index of phase, \param j - index DC in phase
/// \return if error code, returns 777777777.
double TSolModMulti::ConvertGj_toUniformStandardState(double g0, long int j, long int k)
{
    double G, YOF=0;

    G = g0/pm.RT;
    if( pm.YOF )
        YOF = pm.YOF[k];     // should be already normalized (J/mol/RT)
    // Calculation of standard concentration scaling terms
    switch( pm.DCC[j] )
    { // Aqueous electrolyte
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
    case DC_AQ_SURCOMP:
        G += pm.ln5551;
        // calculate molar mass of solvent
         [[fallthrough]];
    case DC_AQ_SOLVCOM:
    case DC_AQ_SOLVENT:
        //if( testTSyst() ) // empty function in standalone
          if( noZero( YOF ) )
            G += YOF;  // In GEMS3K, YOF[k] is the only way to influence G directly
        break;
    case DC_GAS_COMP: // gases except H2O and CO2
    case DC_GAS_H2O: // index to switch off?
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR: // changed by DK on 4.12.2006
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
case DC_SCM_SPECIES:
        if( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_FLUID
            || pm.PHC[k] == PH_PLASMA )
        {
//        if( pm.Pparc[j] != 1.0 && pm.Pparc[j] > 1e-30 )
//           G += log( pm.Pparc[j] ); // log partial pressure/fugacity
//        else
               G += log( pm.P ); // log general pressure (changed 04.12.2006)
        }
        // non-electrolyte condensed mixtures
         [[fallthrough]];
    case DC_SCP_CONDEN: // single-component phase
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:
    case DC_PEL_CARRIER:
        // if( testTSyst() )  // empty function in standalone
            if( noZero( YOF ) )
                G += YOF;
        break;
        // Sorption phases
    case DC_SSC_A0:
    case DC_SSC_A1:
    case DC_SSC_A2:
    case DC_SSC_A3:
    case DC_SSC_A4:
    case DC_WSC_A0:
    case DC_WSC_A1:
    case DC_WSC_A2:
    case DC_WSC_A3:
    case DC_WSC_A4:
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
        G += pm.ln5551;
        break;
    default: // error - returning 7777777
        return 7777777.;
    }
    return G;
}

/// Get the index of volume IC ("Vv") for the volume balance constraint
long int TSolModMulti::getXvolume()
{
    long int ii, ret = -1;
    for( ii = pm.N-1; ii>=0; ii--)
    {
        if( pm.ICC[ii] == IC_VOLUME )
        { ret = ii; break; }
    }
    return ret;
}

// tsolmod_multi_add.cpp
