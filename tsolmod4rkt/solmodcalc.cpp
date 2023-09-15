#include "solmodcalc.h"
#include "v_service.h"
#include "verror.h"


SolModCalc::SolModCalc(long k, long jb, SolutionData &sd, const AddSolutionData &addsd):
    mod_code(sd.Mod_Code), phase_name(sd.phaseName),
    phase_ndx(k), dc_ndx(jb), dc_num(sd.NSpecies)
{
    arWx = sd.arWx;
    arlnGam = sd.arlnGam;
    for(long int ii=0; ii<dc_num; ++ii) {
        dc_names.push_back(char_array_to_string(sd.arSM[ii], MAXDCNAME));
    }
    SolMod_create(sd, addsd);
    to_json_file(std::string("solmod_")+std::to_string(phase_ndx)+".json");
}

SolModCalc::SolModCalc(long k, long jb, const std::string &aphase):
    mod_code(' '), phase_name(aphase),
    phase_ndx(k), dc_ndx(jb), dc_num(1)
{
   arWx = nullptr;
   arlnGam = nullptr;
   dc_names.push_back(phase_name);
   solmod_task.reset();
   model_name = "undefined";
}

void SolModCalc::SolModParPT()
{
    if(check_mode(mod_code) && solmod_task) {
        solmod_task->PTparam();
    }
}

void SolModCalc::SolModActCoeff()
{
    /// ???? If no clean illegal result
    if(arlnGam) {
        map2property({}, arlnGam, 0.);
    }
    if(check_mode(mod_code) && solmod_task) {
        solmod_task->MixMod();
        //solmod_task->to_text_file(std::string("solmod_act_coef_")+std::to_string(phase_ndx)+".txt", true);
    }
}

std::map<std::string, double> SolModCalc::SolModExcessProp()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(mod_code != SM_IDEAL && check_mode(mod_code) && solmod_task) {
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

std::map<std::string, double> SolModCalc::SolModIdealProp()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(mod_code) && solmod_task) {
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

std::map<std::string, double> SolModCalc::SolModDarkenProp()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(mod_code) && solmod_task) {
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

std::map<std::string, double> SolModCalc::SolModStandProp()
{
    std::map<std::string, double> ex_map;
    // order of phase properties: G, H, S, CP, V, A, U
    double zex[7];

    for(long int j=0; j<7; j++) {
        zex[j] = 0.0;
    }
    if(check_mode(mod_code) && solmod_task) {
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

void SolModCalc::Get_lnGamma(double *lngamma)
{
    if(solmod_task) {
        solmod_task->Get_lnGamma(lngamma);
    }
}

std::map<std::string, double> SolModCalc::GetlnGamma()
{
    return property2map(arlnGam);
}

void SolModCalc::Set_MoleFractionsWx(double *aWx)
{
    if(arWx) {
        for(int jj=0; jj<dc_num; ++jj) {
            arWx[jj] = aWx[jj];
        }
    }
}

void SolModCalc::SetMoleFractionsWx(const std::map<std::string, double> &awx_map, double defwx)
{
    if(arWx) {
        map2property(awx_map, arWx, defwx);
    }
}

//--------------------------------------------------------------------------------

// Wrapper calls for creating multi-component mixing models for phases
// using  TSolMod class. Now including multi-site ideal and scripted models
void SolModCalc::SolMod_create(SolutionData& sd, const AddSolutionData& addsd)
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
        break;
    }
    solmod_task.reset(mySM);
}

bool SolModCalc::check_mode(char ModCode)
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

std::map<std::string, double> SolModCalc::property2map(double *dcs_size_array)
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

void SolModCalc::map2property(const std::map<std::string, double> &dsc_name_map, double *dcs_size_array, double def_value)
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
