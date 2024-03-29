//-------------------------------------------------------------------
/// \file solmodengine.cpp
///
/// Implementation of the SolModEngine class - c++ API for phase models
/// Decorator for TSolMod and derived classes implementing built-in models
/// of mixing in fluid, liquid, aqueous, and solid-solution phases
//
// Copyright (c) 2023-2024 S.Dmytriyeva, D.Kulik
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
//

#include "solmodengine.h"
#include "v_service.h"
#include "verror.h"
#include "io_template.h"
#include "io_nlohmann.h"
#include "io_simdjson.h"
#include "jsonconfig.h"

SolModEngine::SolModEngine(long k, long jb, SolutionData &sd, const AddSolutionData &addsd):
    model_code(sd.Mod_Code), phase_name(sd.phaseName),
    phase_ndx(k), dc_ndx(jb), dc_num(sd.NSpecies)
{
    arWx = sd.arWx;
    arM = addsd.arM;
    arFWGT = addsd.arFWGT;
    arX = addsd.arX;

    arlnGam = sd.arlnGam;
    arlnDQFt = sd.arlnDQFt;
    arlnRcpt = sd.arlnRcpt;
    arlnExet = sd.arlnExet;
    arlnCnft = sd.arlnCnft;
    arGEX = sd.arGEX;
    arVol = sd.arVol;
    aphVOL = sd.aphVOL;
    arPparc = sd.arPparc;

    for(long int ii=0; ii<dc_num; ++ii) {
        dc_names.push_back(char_array_to_string(sd.arSM[ii], MAXDCNAME));
    }
    SolMod_create(sd, addsd);
    //to_json_file(std::string("solmod_")+std::to_string(phase_ndx)+".json");
}

SolModEngine::SolModEngine(long k, long jb, const std::string &aphase):
    model_code(' '), phase_name(aphase),
    phase_ndx(k), dc_ndx(jb), dc_num(1)
{
    arWx = nullptr;
    arM = nullptr;
    arFWGT = nullptr;
    arX = nullptr;

    arlnGam = nullptr;
    arlnDQFt = nullptr;
    arlnRcpt = nullptr;
    arlnExet = nullptr;
    arlnCnft = nullptr;
    arGEX = nullptr;
    arVol = nullptr;
    aphVOL = nullptr;
    arPparc = nullptr;

    dc_names.push_back(phase_name);
    solmod_task.reset();
    model_name = "undefined";
}

void SolModEngine::SolModPTParams()
{
    if(check_mode(model_code) && solmod_task) {
        solmod_task->PTparam();
    }
}

void SolModEngine::SolModActivityCoeffs()
{
    /// ???? If no clean illegal result
    if(arlnGam) {
        map2property({}, arlnGam, 0.);
    }
    if(check_mode(model_code) && solmod_task) {
        solmod_task->MixMod();
        //solmod_task->to_text_file(std::string("solmod_act_coef_")+std::to_string(phase_ndx)+".txt", true);
    }
}

std::map<std::string, double> SolModEngine::SolModExcessProps()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(model_code != SM_IDEAL && check_mode(model_code) && solmod_task) {
        solmod_task->ExcessProp(zex);
        ex_map["Gex"] = zex[0];
        ex_map["Hex"] = zex[1];
        ex_map["Sex"] = zex[2];
        ex_map["CPex"] = zex[3];
        ex_map["Vex"] = zex[4];
        ex_map["Aex"] = zex[5];
        ex_map["Uex"] = zex[6];
        TSolMod::solmod_logger->debug("Assignment of calculated excess properties of mixing [2] Gex={} Hex={} Sex={} CPex={}",
                                      ex_map["Gex"], ex_map["Hex"], ex_map["Sex"], ex_map["CPex"]);
    }
    return ex_map;
}

std::map<std::string, double> SolModEngine::SolModIdealProps()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(model_code) && solmod_task) {
        // check what solution phase (ideal gas?)
        solmod_task->IdealProp(zex);
        ex_map["Gid"] = zex[0];
        ex_map["Hid"] = zex[1];
        ex_map["Sid"] = zex[2];
        ex_map["CPid"] = zex[3];
        ex_map["Vid"] = zex[4];
        ex_map["Aid"] = zex[5];
        ex_map["Uid"] = zex[6];
        TSolMod::solmod_logger->debug("Assignment of calculated excess properties of mixing [2] Gid={} Hid={} Sid={} CPid={}",
                                      ex_map["Gid"], ex_map["Hid"], ex_map["Sid"], ex_map["CPid"]);
    }
    return ex_map;
}

std::map<std::string, double> SolModEngine::SolModDarkenProps()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(model_code) && solmod_task) {
        //solmod_task->(zex);
        ex_map["Gdq"] = zex[0];
        ex_map["Hdq"] = zex[1];
        ex_map["Sdq"] = zex[2];
        ex_map["CPdq"] = zex[3];
        ex_map["Vdq"] = zex[4];
        ex_map["Adq"] = zex[5];
        ex_map["Udq"] = zex[6];
    }
    return ex_map;
}

std::map<std::string, double> SolModEngine::SolModStandProps()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(model_code) && solmod_task) {
        // check what solution phase (ideal gas?)
        solmod_task->StandardProp(zex);
        ex_map["Gst"] = zex[0];
        ex_map["Hst"] = zex[1];
        ex_map["Sst"] = zex[2];
        ex_map["CPst"] = zex[3];
        ex_map["Vst"] = zex[4];
        ex_map["Ast"] = zex[5];
        ex_map["Ust"] = zex[6];
        TSolMod::solmod_logger->debug("Assignment of calculated excess properties of mixing [2] Gst={} Hst={} Sst={} CPst={}",
                                      ex_map["Gst"], ex_map["Hst"], ex_map["Sst"], ex_map["CPst"]);
    }
    return ex_map;
}

void SolModEngine::Get_MoleFractions(double *molfr)
{
    if(arWx) {
        for(int jj=0; jj<dc_num; ++jj) {
            molfr[jj] = arWx[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetMoleFractions()
{
    return property2map(arWx);
}

std::vector<double> SolModEngine::Get_MoleFractions()
{
    return property2vector(arWx);
}

void SolModEngine::Get_Molalities(double *molal)
{
    if(arM) {
        for(int jj=0; jj<dc_num; ++jj) {
            molal[jj] = arM[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetMolalities()
{
    return property2map(arM);
}

std::vector<double> SolModEngine::Get_Molalities()
{
    return property2vector(arM);
}


void SolModEngine::Get_lnActivityCoeffs(double *lngamma)
{
    if(solmod_task) {
        solmod_task->Get_lnGamma(lngamma);
    }
}

std::map<std::string, double> SolModEngine::GetlnActivityCoeffs()
{
    return property2map(arlnGam);
}

std::vector<double> SolModEngine::Get_lnActivityCoeffs()
{
    return property2vector(arlnGam);
}

void SolModEngine::Get_lnActivities(double* lnactiv)
{ 
    if(lnactiv && arM && arlnGam && arWx) {
        for(int jj=0; jj<dc_num; ++jj)
        {
            lnactiv[jj] = 0.0;
            if( check_molal_scale(model_code) && arM[jj] > 1e-23 ) {
                lnactiv[jj] = log(arM[jj]) + arlnGam[jj];
            }
            else {
                if( arWx[jj] > 1e-23 ) {
                    lnactiv[jj] = log(arWx[jj]) + arlnGam[jj];
                }
            }
        }
    }
}

std::map<std::string, double> SolModEngine::GetlnActivities()
{
    auto vec = Get_lnActivities();
    return property2map(vec.data());
}

std::vector<double> SolModEngine::Get_lnActivities()
{
    std::vector<double> dsc_name_vec;
    if(!arlnGam || !arM || !arWx) { // nullptr
        return dsc_name_vec;
    }
    double lna;
    for(int jj=0; jj<dc_num; ++jj) {
        lna = 0.0;
        if( check_molal_scale(model_code) && arM[jj] > 1e-23 ) {
            lna = log(arM[jj]) + arlnGam[jj];
        }
        else {
            if( arWx[jj] > 1e-23 ) {
                lna = log(arWx[jj]) + arlnGam[jj];
            }
        }
        dsc_name_vec.push_back(lna);
    }
    return dsc_name_vec;
}

void SolModEngine::Get_lnConfTerms(double *lnGamConf)
{
    if(arlnCnft) {
        for(int jj=0; jj<dc_num; ++jj) {
            lnGamConf[jj] = arlnCnft[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetlnConfTerms()
{
    return property2map(arlnCnft);
}

std::vector<double> SolModEngine::Get_lnConfTerms()
{
   return property2vector(arlnCnft);
}

void SolModEngine::Get_lnRecipTerms(double *lnGamRecip)
{
    if(arlnRcpt) {
        for(int jj=0; jj<dc_num; ++jj) {
            lnGamRecip[jj] = arlnRcpt[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetlnRecipTerms()
{
    return property2map(arlnRcpt);
}

std::vector<double> SolModEngine::Get_lnRecipTerms()
{
    return property2vector(arlnRcpt);
}

void SolModEngine::Get_lnExcessTerms(double *lnGamEx)
{
    if(arlnExet) {
        for(int jj=0; jj<dc_num; ++jj) {
            lnGamEx[jj] = arlnExet[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetlnExcessTerms()
{
    return property2map(arlnExet);
}

std::vector<double> SolModEngine::Get_lnExcessTerms()
{
    return property2vector(arlnExet);
}

void SolModEngine::Get_lnDQFTerms(double *lnGamDQF)
{
    if(arlnDQFt) {
        for(int jj=0; jj<dc_num; ++jj) {
            lnGamDQF[jj] = arlnDQFt[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetlnDQFTerms()
{
    return property2map(arlnDQFt);
}

std::vector<double> SolModEngine::Get_lnDQFTerms()
{
    return property2vector(arlnDQFt);
}

void SolModEngine::Get_G0Increments(double *aGEX)
{
    if(arGEX) {
        for(int jj=0; jj<dc_num; ++jj) {
            aGEX[jj] = arGEX[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetG0Increments()
{
    return property2map(arGEX);
}

std::vector<double> SolModEngine::Get_G0Increments()
{
    return property2vector(arGEX);
}

void SolModEngine::Get_MolarVolumes(double *aVol)
{
    if(arVol) {
        for(int jj=0; jj<dc_num; ++jj) {
            aVol[jj] = arVol[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetMolarVolumes()
{
    return property2map(arVol);
}

std::vector<double> SolModEngine::Get_MolarVolumes()
{
    return property2vector(arVol);
}

double SolModEngine::GetPhaseVolume()
{
    if(aphVOL) {
        return *aphVOL;
    }
    else {
        return 0.;
    }
}

void SolModEngine::Get_PartialPressures(double *aPparc)
{
    if(arPparc) {
        for(int jj=0; jj<dc_num; ++jj) {
            aPparc[jj] = arPparc[jj];
        }
    }
}

std::map<std::string, double> SolModEngine::GetPartialPressures()
{
    return property2map(arPparc);
}

std::vector<double> SolModEngine::Get_PartialPressures()
{
    return property2vector(arPparc);
}


void SolModEngine::Set_MoleFractions(double *aWx)
{
    if(arWx) {
        for(int jj=0; jj<dc_num; ++jj) {
            arWx[jj] = aWx[jj];
        }
    }
}

void SolModEngine::SetMoleFractions(const std::map<std::string, double> &awx_map, double defwx)
{
    if(arWx) {
        map2property(awx_map, arWx, defwx);
    }
}

void SolModEngine::Set_SpeciesMolality(double *aM)
{
    if(arM) {
        for(int jj=0; jj<dc_num; ++jj) {
            arM[jj] = aM[jj];
        }
    }
}

void SolModEngine::SetSpeciesMolality(const std::map<std::string, double> &val_map, double def_val)
{
    if(arM) {
        map2property(val_map, arM, def_val);
    }
}

void SolModEngine::Set_SpeciesAmounts(double *aX)
{
    if(arX) {
        for(int jj=0; jj<dc_num; ++jj) {
            arX[jj] = aX[jj];
        }
    }
}

void SolModEngine::SetSpeciesAmounts(const std::map<std::string, double> &val_map, double def_val)
{
    if(arX) {
        map2property(val_map, arX, def_val);
    }
}

void SolModEngine::Set_PhaseMass(double aFWGT)
{
    if(arFWGT) {
        *arFWGT = aFWGT;
    }
}

//--------------------------------------------------------------------------------

// Wrapper calls for creating multi-component mixing models for phases
// using  TSolMod class. Now including multi-site ideal and scripted models
void SolModEngine::SolMod_create(SolutionData& sd, const AddSolutionData& addsd)
{
    model_name.clear();
    TSolMod* mySM = nullptr;

    // creating instances of subclasses of TSolMod base class
    switch(sd.Mod_Code) {
    case SM_OTHER:  // Hard-coded solid solution models (selected by phase name)
    {
        model_name = "TModOther";
        TModOther* myPT = new TModOther(&sd, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_VANLAAR:  // Van Laar solid solution model (multicomponent)
    {
        model_name = "TVanLaar";
        TVanLaar* myPT = new TVanLaar(&sd);
        mySM = myPT;
        break;
    }
    case SM_REGULAR:  // Regular solid solution model (multicomponent)
    {
        model_name = "TRegular";
        TRegular* myPT = new TRegular(&sd);
        mySM = myPT;
        break;
    }
    case SM_GUGGENM:  // Redlich-Kister solid solution model (multicomponent)
    {
        model_name = "TRedlichKister";
        TRedlichKister* myPT = new TRedlichKister(&sd);
        mySM = myPT;
        break;
    }
    case SM_NRTLLIQ:  // NRTL liquid solution model (multicomponent)
    {
        model_name = "TNRTL";
        TNRTL* myPT = new TNRTL(&sd);
        mySM = myPT;
        break;
    }
    case SM_WILSLIQ:  // Wilson liquid solution model (multicomponent)
    {
        model_name = "TWilson";
        TWilson* myPT = new TWilson(&sd);
        mySM = myPT;
        break;
    }
    case SM_MARGT:  // Margules ternary (regular) solid solution model
    {
        model_name = "TMargules";
        TMargules* myPT = new TMargules(&sd);
        mySM = myPT;
        break;
    }
    case SM_MARGB:  // Margules binary (subregular) solid solution model
    {
        model_name = "TSubregular";
        TSubregular* myPT = new TSubregular(&sd);
        mySM = myPT;
        break;
    }
    case SM_REDKIS:  // Gugenheim binary (REdlich-Kister) solid solution
    {
        model_name = "TGuggenheim";
        TGuggenheim* myPT = new TGuggenheim(&sd);
        mySM = myPT;
        break;
    }
    case SM_AQPITZ:  // Pitzer aqueous electrolyte model (multicomponent)
    {
        model_name = "TPitzer";
        TPitzer* myPT = new TPitzer(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQSIT:  // SIT aqueous electrolyte model (multicomponent)
    {
        model_name = "TSIT";
        TSIT* myPT = new TSIT(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQEXUQ:  // EUNIQUAC aqueous electrolyte model (multicomponent)
    {
        model_name = "TEUNIQUAC";
        TEUNIQUAC* myPT = new TEUNIQUAC(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
        /*      case SM_AQELVIS:  // ELVIS aqueous electrolyte model (multicomponent)
        {
                model_name = "TELVIS";
                TELVIS* myPT = new TELVIS(&sd, aM, aZ, pm.denW, pm.epsW);
                mySM = myPT;
                break;
        }
*/
    case SM_AQDH3:  // extended Debye-Hueckel aqueous electrolyte model (Karpov version)
    {
        model_name = "TKarpov";
        TKarpov* myPT = new TKarpov(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQDH2:   // Debye-Hueckel aqueous electrolyte model
    {
        model_name = "TDebyeHueckel";
        TDebyeHueckel* myPT = new TDebyeHueckel(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQDH1:   // Debye-Hueckel limiting law aqueous electrolyte model
    {
        model_name = "TLimitingLaw";
        TLimitingLaw* myPT = new TLimitingLaw(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQDHS:  // extended Debye-Hueckel aqueous electrolyte model (Shvarov version)
    {
        model_name = "TShvarov";
        TShvarov* myPT = new TShvarov(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQDHH:  // extended Debye-Hueckel aqueous electrolyte model (Helgeson version)
    {
        model_name = "THelgeson";
        THelgeson* myPT = new THelgeson(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_AQDAV:  // Davies aqueous electrolyte model (in NEA TDB version)
    {
        model_name = "TDavies";
        TDavies* myPT = new TDavies(&sd, addsd.arM, addsd.arZ, addsd.ardenW, addsd.arepsW);
        mySM = myPT;
        break;
    }
    case SM_PRFLUID:  // PRSV fluid mixture (multicomponent)
    {
        model_name = "TPRSVcalc";
        TPRSVcalc* myPT = new TPRSVcalc(&sd);
        mySM = myPT;
        break;
    }
    case SM_CGFLUID:  // CG fluid mixture (multicomponent)
    {
        model_name = "TCGFcalc";
        TCGFcalc* myPT = new TCGFcalc(&sd, addsd.arFWGT, addsd.arX);
        mySM = myPT;
        break;
    }
    case SM_SRFLUID:  // SRK fluid mixture (multicomponent)
    {
        model_name = "TSRKcalc";
        TSRKcalc* myPT = new TSRKcalc(&sd);
        mySM = myPT;
        break;
    }
    case SM_PR78FL:  // PR78 fluid mixture (multicomponent)
    {
        model_name = "TPR78calc";
        TPR78calc* myPT = new TPR78calc(&sd);
        mySM = myPT;
        break;
    }
    case SM_CORKFL:  // CORK fluid mixture (multicomponent)
    {
        model_name = "TCORKcalc";
        TCORKcalc* myPT = new TCORKcalc(&sd);
        mySM = myPT;
        break;
    }
    case SM_STFLUID:  // STP fluid mixture (H2O-CO2)
    {
        model_name = "TSTPcalc";
        TSTPcalc* myPT = new TSTPcalc(&sd);
        mySM = myPT;
        break;
    }
    case SM_BERMAN:  // Non-ideal (multi-site) model
    {
        model_name = "TBerman";
        TBerman* myPT = new TBerman(&sd, addsd.arG0);
        mySM = myPT;
        break;
    }
    case SM_CEF:  // Non-ideal (multi-site) model (CALPHAD)
    {
        model_name = "TCEFmod";
        TCEFmod* myPT = new TCEFmod(&sd, addsd.arG0);
        mySM = myPT;
        break;
    }
    case SM_MBW:
    {
        model_name = "TMBWmod";
        TMBWmod* myPT = new TMBWmod(&sd, addsd.arG0);
        mySM = myPT;
        break;
    }
    case SM_IDEAL:
    {
        model_name = "TIdeal";
        TIdeal* myPT = new TIdeal(&sd);
        mySM = myPT;
        break;
    }
        // case SM_USERDEF:
    case SM_SURCOM:
    {
        model_name = "TSCM_NEM";
        TSCM_NEM* myPT = new TSCM_NEM(&sd);
        mySM = myPT;
        break;
    }
    default:
        model_name = "Undefined";
        break;
    }
    solmod_task.reset(mySM);
}

bool SolModEngine::check_mode(char ModCode)
{
    // Extended constructor to connect to params, coeffs, and mole fractions
    switch( ModCode )
    {   // solid and liquid solutions
    // case SM_USERDEF:
    case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
    case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
    case SM_BERMAN:  case SM_CEF:  /* new ss models */   case SM_MBW: /* new mbw model */ case SM_SURCOM:
        // aqueous DH models
    case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
    case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
    case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
    case SM_STFLUID:
    {
        if(!solmod_task.get()) {
            TSolMod::solmod_logger->info("SolMod the undefined instance of phase  {}", phase_name);
            return false;
        }
        return true;
    }
    default:
        break;
    }
    return false;
}

bool SolModEngine::check_molal_scale(char ModCode)
{
    // Returns true if concentrations are expressed in molality scale
    switch( ModCode )
    {    
    // aqueous DH models
    case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
    // aqueous SIT models
    case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        return true;
    default:
        return false; 
    }
    return false;      
}

std::map<std::string, double> SolModEngine::property2map(double *dcs_size_array)
{
    std::map<std::string, double> dsc_name_map;
    if(!dcs_size_array) { // nullptr
        return dsc_name_map;
    }
    for(int jj=0; jj<dc_num; ++jj) {
        dsc_name_map[dc_names[jj]] = dcs_size_array[jj];
    }
    return dsc_name_map;
}

std::vector<double> SolModEngine::property2vector(double *dcs_size_array)
{
    if(dcs_size_array) {
        return std::vector<double>(dcs_size_array, dcs_size_array+dc_num);
    }
    return {};
}

void SolModEngine::map2property(const std::map<std::string, double> &dsc_name_map, double *dcs_size_array, double def_value)
{
    if(!dcs_size_array) { // nullptr
        return;
    }
    for(int jj=0; jj<dc_num; ++jj) {
        auto it =dsc_name_map.find(dc_names[jj]);
        if(it!=dsc_name_map.end()) {
            dcs_size_array[jj] = it->second;
        }
        else {
            dcs_size_array[jj] = def_value;
        }
    }
}

void SolModEngine::to_json_stream_short(std::ostream& ff) const
{
#ifdef USE_NLOHMANNJSON
    io_formats::NlohmannJsonWrite out_format( ff, "set_name" );
#else
    io_formats::SimdJsonWrite out_format( ff, "set_name", true );
#endif
    out_format.put_head( phase_name, "tsolmod");
    io_formats::TPrintArrays  prar( 0, {}, out_format );

    prar.addField("PhaseName", phase_name);
    prar.addField("ModCode", model_code);
    prar.addField("NComp", dc_num);
    prar.writeArray( "phVOL",  aphVOL, 1L );
    prar.writeArray( "aFWGT",  arFWGT, 1L );

    prar.writeArray( "DCNL", dc_names, 10 );
    prar.writeArray( "aGEX",  arGEX, dc_num );
    prar.writeArray( "aPparc",  arPparc, dc_num );
    prar.writeArray( "aVol",  arVol, dc_num );
    prar.writeArray( "aWx",  arWx, dc_num );
    prar.writeArray( "lnGamma",  arlnGam, dc_num );
    prar.writeArray( "aM",  arM, dc_num );
    prar.writeArray( "aX",  arX, dc_num );

    out_format.dump( true );
}

void SolModEngine::to_json_file(const std::string &path) const
{
    std::fstream ff(GemsSettings::with_directory(path), std::ios::out);
    ErrorIf(!ff.good(), path, "Fileopen error");
    to_json(ff);
}

std::ostream& operator<<(std::ostream& out, const SolModEngine& obj)
{
    obj.to_json(out);
    return out;
}

//--------------------- end of solmodengine.cpp ---------------------------

