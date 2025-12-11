//-------------------------------------------------------------------
/// \file solmodfactory.cpp
///
/// Implementation of subset of TMulti class,
//
// Copyright (c) 2023-2024 S.Dmytriyeva,D.Kulik
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

#include <fstream>
#include "solmodfactory.h"
#include "jsonconfig.h"
#include "datach_api.h"
#include "num_methods.h"
#include "v_service.h"

// Conversion factors
const double bar_to_Pa = 1e5,
m3_to_cm3 = 1e6,
kg_to_g = 1e3;

// New --------------------------------------------------------------------

// Constructor of the class instance in memory for standalone GEMS3K or coupled program
SolModFactory::SolModFactory(const std::string &ipmfiles_lst_name)
{
    alloc_main();
    // Initialization of GEMS3K internal data by reading files
    if( GEM_init(ipmfiles_lst_name) )
    {
        Error(ipmfiles_lst_name, "error occured during reading the files");
    }
}

SolModFactory::SolModFactory(const std::string &dch_json, const std::string &ipm_json,
                             const std::string &dbr_json, const std::string &fun_json)
{
    alloc_main();
    // Initialization of GEMS3K internal data by reading files
    if( GEM_init(dch_json, ipm_json, dbr_json, fun_json) )
    {
        Error("Init", "error occured during reading the strings");
    }
}

void SolModFactory::alloc_main()
{
    pmp = &pm;
    pm.errorCode[0] ='\0';
    pm.errorBuf[0] ='\0';
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

SolModFactory::~SolModFactory()
{
    clear_ThermoEngine();
    freeMemory();
}

void SolModFactory::allocMemory()
{
    // memory allocation for data bridge structures
    CSD = new DATACH;
    CNode = new DATABR;
    dbr_dch_api::datach_reset(CSD);
    dbr_dch_api::databr_reset(CNode, 2);
}

void SolModFactory::freeMemory()
{
    dbr_dch_api::databr_free_internal(CNode);
    delete CNode;
    CNode = nullptr;
    dbr_dch_api::datach_free(CSD);
    delete CSD;
    CSD = nullptr;
}

//-------------------------------------------------------------------
long int  SolModFactory::GEM_init( const std::string& ipmfiles_lst_name )
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
        pmp->Pai[2] = getStep( pmp->Pai, CSD->nPp )/bar_to_Pa;
        pmp->Pai[3] = CSD->Ptol/bar_to_Pa;

        pmp->Tai[0] = CSD->TKval[0]-C_to_K;
        pmp->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
        pmp->Tai[2] = getStep( pmp->Tai, CSD->nTp );
        pmp->Tai[3] = CSD->Ttol;

        pmp->Fdev1[0] = 0.;
        pmp->Fdev1[1] = 1e-6;
        pmp->Fdev2[0] = 0.;
        pmp->Fdev2[1] = 1e-6;

        // Reading DBR_DAT file into work DATABR structure from ipmfiles_lst_name
        /// ?? Used only T and P
        std::string dbr_file = generator.get_dbr_path( 0 );
        ErrorIf( dbr_file.empty() , ipmfiles_lst_name, " Undefined DBR_DAT file name");
        std::fstream in_br( dbr_file, std::ios::in );
        ErrorIf( !in_br.good() , dbr_file, "DBR_DAT fileopen error");
        dbr_dch_api::read_dbr_format_stream(current_input_set_name, CSD, CNode, in_br, generator.files_mode());

        // Loading thermodynamic, set default values and other (????)
        InitalizeGEM_IPM_Data();
        // Creating and initializing the TSolMod array
        InitalizeTSolMod();
        unpackDataBr(true);
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
        Error(ipmfiles_lst_name, "unknown exception");
    }

    if( !ipmfiles_lst_name.empty() ) {
        TSolMod::solmod_logger->error("GEMS3K input : file {}", ipmfiles_lst_name);
    }
    if( !ipmlog_error.empty() ) {
        TSolMod::solmod_logger->error("GEM_init error: {}", ipmlog_error);
    }
    return 1;
}

long int  SolModFactory::GEM_init( const std::string& dch_json, const std::string& ipm_json,
                                   const std::string& dbr_json, const std::string& fun_json)
{
    load_thermodynamic_data = false; // need load thermo
    clearipmLogError();
    clear_ThermoEngine();

    try
    {
        if(fun_json.empty()) {
            GEMS3KGenerator::default_type_f  = GEMS3KGenerator::f_json;
        }
        else {
            GEMS3KGenerator::default_type_f  = GEMS3KGenerator::f_thermofun;
        }

        // This check of data consistency temporarily disabled for perfomance testing
        //if( GEMS3KGenerator::default_type_f  == GEMS3KGenerator::f_json )
        //{
        current_output_set_name = current_input_set_name = extract_string_json( "set", dch_json );
        //     auto ipm_set =  extract_string_json( "set", ipm_json );
        //     auto dbr_set =  extract_string_json( "set", dbr_json );
        //     ErrorIf(  ipm_set!=current_input_set_name,  "GEM_init error", "Multi structure as a json has different set name:  "+ipm_set );
        //     ErrorIf(  dbr_set!=current_input_set_name,  "GEM_init error", "The data bridge structure as a json has different set name:  "+dbr_set );
        //}

        // Reading DCH_DAT data
        dbr_dch_api::datach_from_string(current_input_set_name, CSD, dch_json);

        if( !fun_json.empty() )
        {
            load_ThermoEngine(fun_json);
        }
        // Reading IPM_DAT file into structure MULTI (GEM IPM work structure)
        gemipm_from_string( ipm_json, CSD, current_input_set_name );

        // copy intervals for minimization
        pmp->Pai[0] = CSD->Pval[0]/bar_to_Pa;
        pmp->Pai[1] = CSD->Pval[CSD->nPp-1]/bar_to_Pa;
        pmp->Pai[2] = getStep( pmp->Pai, CSD->nPp )/bar_to_Pa;

        pmp->Tai[0] = CSD->TKval[0]-C_to_K;
        pmp->Tai[1] = CSD->TKval[CSD->nTp-1]-C_to_K;
        pmp->Tai[2] = getStep( pmp->Tai, CSD->nTp );

        pmp->Fdev1[0] = 0.;
        pmp->Fdev1[1] = 1e-6;
        pmp->Fdev2[1] = 1e-6;

        // Reading DBR_DAT file into work DATABR structure from ipmfiles_lst_name
        dbr_dch_api::databr_from_string(current_input_set_name, CSD, CNode, dbr_json);

        // Loading thermodynamic, set default values and other (????)
        InitalizeGEM_IPM_Data();
        // Creating and initializing the TSolMod array
        InitalizeTSolMod();
        unpackDataBr(true);
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
        Error("GEM_init", "unknown exception");
    }
    if( !ipmlog_error.empty() ) {
        TSolMod::solmod_logger->error("GEM_init error: {}", ipmlog_error);
    }
    return 1;
}


void SolModFactory::UpdateThermoData(double TK, double PPa)
{
    LoadThermodynamicData(TK, PPa);
    for(auto& phase_model: phase_models) {
        phase_model.UpdatePT(TK, PPa);
        phase_model.SolModPTParams();
    }
}

SolModEngine &SolModFactory::Sol_Phase(std::size_t idx)
{
    ErrorIf( idx>=phase_models.size(), "SolModFactory", "array index " + std::to_string(idx) + " is out of range" );
    return phase_models[idx];
}

SolModEngine &SolModFactory::SolPhase(const std::string &name)
{
    for(size_t kk=0; kk<phase_names.size(); ++kk ) {
        if(phase_names[kk] == name) {
            return  Sol_Phase(kk);
        }
    }
    Error( "SolModFactory", "Phase '" + name + "' not found" );
}

std::vector<std::vector<double>> SolModFactory::Get_StoichiometryMatrix()
{
    std::vector<std::vector<double>> A_matr;
    for(long int ii=0; ii<CSD->nDC; ++ii) {
        A_matr.push_back(std::vector<double>(CSD->A+ii*CSD->nIC, CSD->A+ii*CSD->nIC+CSD->nIC));
    }
    return A_matr;
}

std::vector<std::string> SolModFactory::Get_AllElementNames()
{
    std::vector<std::string> names;
    for(long int ii=0; ii<CSD->nIC; ++ii) {
        names.push_back(char_array_to_string(CSD->ICNL[ii], MAXICNAME));
    }
    return names;
}

std::vector<std::string> SolModFactory::Get_AllSpeciesNames()
{
    std::vector<std::string> names;
    for(long int ii=0; ii<CSD->nDC; ++ii) {
        names.push_back(char_array_to_string(CSD->DCNL[ii], MAXDCNAME));
    }
    return names;
}

std::vector<std::string> SolModFactory::Get_AllPhasesNames()
{
    std::vector<std::string> names;
    for(long int ii=0; ii<CSD->nPH; ++ii) {
        names.push_back(char_array_to_string(CSD->PHNL[ii], MAXPHNAME));
    }
    return names;
}

// Calculation by IPM (preparing for calculation, unpacking data) In IPM
void SolModFactory::InitalizeGEM_IPM_Data() // Reset internal data formerly MultiInit()
{
    MultiConstInit();
    LoadThermodynamicData(CNode->TK, CNode->P); // Loading thermodynamic data into MULTI structure
}

void SolModFactory::InitalizeTSolMod()
{
    phase_models.clear();
    phase_names.clear();

    long int k, j, jb, je=0, jpb, jdb;
    long int ipb, jpe=0, jde=0, ipe=0;
    long int  jdqfc=0;
    long int  jmb, jme=0, jsb, jse=0;
    char *sMod;

    for( k=0; k<pm.FI; k++ )
    {
        // loop on solution phases
        auto phase_name = char_array_to_string(pm.SF[k]+MAXSYMB, MAXPHNAME);
        strip(phase_name);
        phase_names.push_back(phase_name);

        jb = je;
        je += pm.L1[k];

        SolutionData sd;
        AddSolutionData addsd;

        sd.phaseName = phase_name;
        sd.Mod_Code = 'N';
        sd.Mix_Code = 'N';
        sd.NSpecies = pm.L1[k];          // Number of components (end members) in the phase

        // properties generic to all models
        sd.arWx = pm.Wx+jb;       // End member mole fractions
        sd.arlnGam = pm.lnGam+jb; // End member ln activity coeffs

        if( je <= pm.Ls) {
            sd.arlnDQFt = pm.lnDQFt+jb; // End member ln activity coeffs
            sd.arlnRcpt = pm.lnRcpt+jb; // End member ln activity coeffs
            sd.arlnExet = pm.lnExet+jb; // End member ln activity coeffs
            sd.arlnCnft = pm.lnCnft+jb; // End member ln activity coeffs
            sd.arCTermt = pm.CTerms+jb; // End member coulombic terms
        }
        sd.aphVOL = pm.FVOL+k;  // CalculateConcentrations
        sd.DC_Codes = pm.DCC+jb;  // pointer to Dcomp class codes (added 02.05.2010 TW)
        sd.arGEX = pm.fDQF+jb;      // DQF parameters or pure-gas fugacities
        sd.arPparc = pm.Pparc+jb;
        sd.TP_Code = &pm.dcMod[jb];
        sd.T_k = pm.Tc;
        sd.P_bar = pm.P;
        sd.arVol = pm.Vol+jb;
        sd.arSM = pm.SM+jb;

        addsd.arZ = pm.EZ+jb;
        addsd.arM = pm.Y_m+jb;
        addsd.ardenW = pm.denW;
        addsd.arepsW = pm.epsW;
        addsd.arG0 = pm.G0+jb;
        addsd.arFWGT = pm.FWGT+k;
        addsd.arX = pm.X+jb;

        if( k >= pm.FIs || ( pm.L1[k] == 1 && !( pm.PHC[k] == PH_GASMIX ||
                                pm.PHC[k] == PH_PLASMA ||
                                pm.PHC[k] == PH_FLUID ) ) )
        {
            // empty model
            phase_models.push_back(SolModEngine(k, jb, sd, addsd));
            continue;
        }
        // Indexes for extracting data from IPx, PMc and DMc arrays
        ipb = ipe;
        ipe += pm.LsMod[k*3]*pm.LsMod[k*3+1];
        jpb = jpe;
        jpe += pm.LsMod[k*3]*pm.LsMod[k*3+2];
        jdb = jde;
        jde += pm.LsMdc[k*3]*pm.L1[k];

        jmb = jme;
        jme += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2]*pm.L1[k];
        jsb = jse;
        jse += pm.LsMdc[k*3+1]*pm.LsMdc[k*3+2];
        sMod = pm.sMod[k];

        double nxk = 1./pm.L1[k];
        for( j= jb; j<je; j++ )   {
            // top if(pm.XF[k] < std::min( pm.DSM, pm.PhMinM ) ) // pm.lowPosNum )
            if(pm.Wx[j] < 1e-10 ) // ???
                pm.Wx[j] = nxk;  // need this eventually to avoid problems with zero mole fractions
            pm.fDQF[j] =0.0;  // cleaning fDQF in TP mode!
            pm.lnGmo[j] = pm.lnGam[j]; // saving activity coefficients in TP mode
        }

        // the following section probably needs to be re-written to allow more flexible
        // combinations of fluid models for pure gases with gE mixing models,
        // scheme should probably be the same as in LINK_UX_MODE, 03.06.2008 (TW)
        switch( pm.PHC[k] )
        {
        case PH_AQUEL: case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
        case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID: case PH_ADSORPT:
        case PH_IONEX:
        {
            // init model data
            sd.Mod_Code = sMod[SPHAS_TYP];
            sd.Mix_Code = sMod[MIX_TYP];
            sd.NParams = pm.LsMod[k*3];      // Number of interaction parameters
            sd.NPcoefs = pm.LsMod[k*3+2];    // and number of coefs per parameter in PMc table
            sd.MaxOrder =  pm.LsMod[k*3+1];  // max. parameter order (cols in IPx)
            sd.NPperDC = pm.LsMdc[k*3];      // Number of non-ideality coeffs per one DC in multicomponent phase
            sd.NSublat = pm.LsMdc[k*3+1];    // Number of site types (sublattices) for multi-site SS model
            sd.NMoiet = pm.LsMdc[k*3+2];     // Number of moieties for multi-site SS model
            sd.NDQFpDC = pm.LsMdc2[k*3];

            // properties generic to all models
            sd.arIPx = pm.IPx+ipb;   // Pointer to list of indexes for non-ideal solutions -> NPar x MaxOrd
            sd.arIPc = pm.PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
            sd.arDCc = pm.DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC

            sd.arMoiSN = pm.MoiSN+jmb;  // Pointer to sublattice-moiety multiplicity array
            sd.arSitFr = pm.SitFr+jsb;  // Pointer to sublattice-moiety multiplicity array
            sd.arDQFc = pm.DQFc+ jdqfc;

            phase_models.push_back(SolModEngine(k, jb, sd, addsd));
            // new solution models (TW, DK 2007)
            phase_models.back().SolModPTParams();
            break;
        }
        default:
            // empty model
            phase_models.push_back(SolModEngine(k, jb, sd, addsd));
            break;
        }
        // move pointers
        jdqfc += (pm.LsMdc2[k*3]*pm.L1[k]);

    } // k
}


// Changed ----------------------------------------------------------------

void SolModFactory::clear_ThermoEngine()
{
#ifndef NO_USE_THERMOFUN
    // clear previous
    thermo_engine.reset();
    thermo_json_string="";
#endif
}

bool SolModFactory::load_ThermoEngine(const std::string &thermo_file_or_string)
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

#ifndef NO_USE_THERMOFUN
    thermo_engine.reset(new ThermoFun::ThermoEngine(thermo_file_or_string));
    TSolMod::solmod_logger->trace("Read ThermoEngine: {}", thermo_file_or_string);
    return true;
#else
    TSolMod::solmod_logger->warn("Try read ThermoEngine not in USE_THERMOFUN mode {}", thermo_file_or_string);
    return false;
#endif
}

bool SolModFactory::load_all_thermodynamic_from_thermo(double TK, double PPa)
{
#ifndef NO_USE_THERMOFUN
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

        if( CSD->ccPH[0] == PH_AQUEL )
        {
            std::string h2o_key = dCH->DCNL[pmp->LO];
            TSolMod::solmod_logger->info("water-solvent: {}", h2o_key);
            auto water_props = thermo_engine->propertiesSolvent(funT, funP, h2o_key);
            auto water_electro = thermo_engine->electroPropertiesSolvent(funT, funP, h2o_key);

            auto water_vapor = thermo_engine->database().getSubstance(h2o_key);
            water_vapor.setMethod_P( ThermoFun::MethodCorrP_Thrift::type::CPM_GAS);

            ThermoFun::Database db2;
            db2.addSubstance(water_vapor);

            ThermoFun::ThermoEngine te2(db2);
            auto water_gas_props = te2.propertiesSolvent(funT, funP, h2o_key);
            auto water_gas_electro = te2.electroPropertiesSolvent(funT, funP, h2o_key);
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
        }
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
                auto propAl    = thermo_engine->thermoPropertiesSubstance(funT, funP, symbol);

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

/// Setup/copy flags and thresholds for numeric modules to TMulti structure.
/// Do it before calculations
void SolModFactory::MultiConstInit() // from MultiRemake
{
    pm.FI1 = 0;
    pm.FI1s = 0;
    pm.FI1a = 0;
    pm.ITF = 0; pm.ITG = 0;
    //pm.PD__ = base_param()->PD;
    pm.Ec = pm.K2 = pm.MK = 0;
    pm.W1 = 0;
    pm.is = 0;
    pm.js = 0;
    pm.next  = 0;
    pm.ln5551 = log( H2O_mol_to_kg );             // constant corrected 30.08.2008
    pm.lowPosNum = Min_phys_amount;               // = 1.66e-24 mol
    pm.logXw = -16.;
    pm.logYFk = -9.;
    //pm.DXM__ = base_param()->DK;

    //  ???????
    pm.FX = 7777777.;
    if( pm.pH < -15. || pm.pH > 16.  )   // Check for trash in pH - bugfix 19.06.2013
        pm.pH = pm.Eh = pm.pe = 0.0;
    pm.YMET = 0;                      // always 0.0 ????
    pm.PCI = 1.0;
    pm.FitVar[4] = 1.0;

    // from multiConstInit_PN();
    //pm.PZ__ = base_param()->DW;  // in IPM
    //  pm.FitVar[0] = 0.0640000030398369;
}

/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void SolModFactory::LoadThermodynamicData(double TK, double PPa)
{
    // try generate thermodynamic data from ThermoEngine
    if( !load_all_thermodynamic_from_thermo( TK, PPa ))
    {
        load_all_thermodynamic_from_grid(TK, PPa);
    }
    pm.pTPD = 2;
}


/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void SolModFactory::load_all_thermodynamic_from_grid(double TK, double PPa )
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
                Vv = dCH->V0[ jj+xTP];
                if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
                if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
                if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
                if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
                if( dCH->U0 ) u0 = dCH->U0[ jj+xTP];
            }
            else
            {
                Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
                Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 5 );
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
            if( pm.V0 ) pm.V0[j] = Vv;
            Vv *=bar_to_Pa;

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
void SolModFactory::ConvertDCC()
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
double SolModFactory::ConvertGj_toUniformStandardState(double g0, long int j, long int k)
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
long int SolModFactory::getXvolume()
{
    long int ii, ret = -1;
    for( ii = pm.N-1; ii>=0; ii--)
    {
        if( pm.ICC[ii] == IC_VOLUME )
        { ret = ii; break; }
    }
    return ret;
}


//--------------------- end of solmodfactory.cpp ---------------------------

