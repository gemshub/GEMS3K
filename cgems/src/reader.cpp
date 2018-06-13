#import "reader.h"

Reader::Reader() {}

std::vector<std::string> Reader::split(const std::string& str){
    std::vector<std::string> result;
    std::istringstream iss(str);
    for(std::string s; iss >> s; )
        result.push_back(s);
    return result;
}


bool Reader::fileExists(const std::string& filename) {
    std::ifstream infile(filename.c_str());
    return infile.good();
}

void Reader::clearFiles() {
    std::stringstream stream;
    std::fstream file;

    if (fileExists(outFileDbr)) {
        file.open(outFileDbr.c_str(), std::fstream::out);
        file.close();
    }

    if (fileExists(outFileDhc)) {
        file.open(outFileDhc.c_str(), std::fstream::out);
        file.close();
    }

    if (fileExists(outFileIpm)) {
        file.open(outFileIpm.c_str(), std::fstream::out);
        file.close();
    }
}

void Reader::openInputFile(std::string& file, SYSDATA& data) {
    data.Clear();

    std::string  line;
    std::fstream inputFile;

    static const std::string commb = "/*";
    static const std::string comme = "*/";

    std::vector<std::string> res;

    inputFile.open(file.c_str(), std::ios::in);

    if (!inputFile.is_open()) {
        std::cout << "Unfortunatelly no input file found ... " << std::endl;
        std::cin.ignore();
        return;
    }

    clearFiles();
    while (!inputFile.eof()) {
        std::getline(inputFile, line);
        res = split(line);

        std::string::const_iterator b = line.begin();
        std::string::const_iterator e = line.end();

        if (std::search(b, e, commb.begin(), commb.end()) != e) {
            while (std::search(b, e, comme.begin(), comme.end()) == e && !inputFile.eof()){
                std::getline(inputFile, line);
                b = line.begin();
                e = line.end();
                if (inputFile.eof()){
                    std::cout <<  "Looks like you forgot to put the closing comment" << std::endl;
                    std::cin.ignore();
                }
            }
            std::getline(inputFile, line);
        }

        if (line[0]!= '#' && line[0]!= '/')
            readPairs(res.begin(), res.end(), data);
    }

    inputFile.close();
    for (unsigned int i = 0; i < data.Compounds.size(); i++) {
        data.Compounds[i].atoms = parser.readFormula(data.Compounds[i].formula);
    }

    data.totElements.clear();
    for (auto& m : data.Compounds) {
        for (auto& r : m.atoms) {
            data.totElements[r.first] += r.second;
        }
    }
    std::cout << "Read the input file" << std::endl;
}

void Reader::readPairs(std::vector<std::string>::iterator begin,
                       std::vector<std::string>::iterator end,
                       SYSDATA& data) {
    std::string unit;
    std::string value;
    ValueUnit ValUn;
    std::stringstream stream;

    while (begin != end) {
        if (*begin == "Project" || *begin == "PROJECT" || *begin == "project") {
            stream.str("");
            stream.clear();
            ++begin;
            while (begin < end) {
                stream << *begin << "  ";
                ++begin;
            }
            data.project = stream.str();
            continue;
        }

        if (*begin == "P" || *begin == "Pressure" || *begin == "pressure") {
            ++begin;
            ValUn = getValue(*begin);
            data.P = ValUn.value;
        }
        if (*begin == "T" || *begin == "Temperature" || *begin == "temperature") {
            ++begin;
            ValUn = getValue(*begin);
            data.T = ValUn.value;
        }

        if (*begin == "Compound" || *begin == "Comp" || *begin=="DComp" || *begin == "compound") {
            CompoundStr compound;
            ++begin;
            compound.name = *begin;
            ++begin;
            compound.formula = *begin;
            data.Compounds.push_back(compound);
        }

        if (*begin == "Phase" || *begin == "phase") {
            PhaseStr phase;
            ++begin;
            phase.code = (*begin);
            ++begin;
            phase.name = (*begin);
            ++begin;
            while (begin < end){
                phase.compounds.push_back(*begin);
                ++begin;
            }
            data.Phases.push_back(phase);
            continue;
        }

        if (*begin == "Concentration" || *begin == "Conc" || *begin=="conc" || *begin == "concentration") {
            ++begin;
            int i = 0;
            while (begin < end){
                //data.Concentration.push_back(ToDouble(*begin));
                data.Compounds[i].initX = ToDouble(*begin);
                ++i;
                ++begin;
            }
            continue;
        }

        if (*begin == "Ammount" || *begin == "Amm" || *begin=="amm" || *begin == "ammount") {
            data.Ammount.clear();
            ++begin;
            while (begin < end){
                data.Ammount.push_back(ToDouble(*begin));
                ++begin;
            }
            continue;
        }

        ++begin;
    }
}

//void Reader::readFormula(std::string& formula, SYSDATA& data) {
//    data.Elements.push_back(parser.readFormula(formula));
//}

void Reader::generateDHCfile(SYSDATA & data, std::string& fileOut) {
    std::fstream file;
    std::string text = "";
    std::stringstream stream;

    file.open(fileOut.c_str(), std::fstream::out);

    stream.str("");
    stream.clear();

    stream << "#  GEMS3K v.3.3 r.1036 (rc) " << std::endl;
    stream << "# File: /home/kulik/DevGEMS/CalcColumn-dch.dat" << std::endl;
    stream << "# Comments can be marked with # $ ; as the first character in the line" << std::endl;
    stream << "# DCH text input file (should be read before IPM and DBR files)" << std::endl;
    stream << std::endl;

    stream << "## (1) Dimensions for memory allocation" << std::endl;
    stream << "# nIC: Number of Independent Components (usually chemical elements and charge)" << std::endl;
    stream << "<nIC>  " << data.totElements.size() << std::endl;
    stream << "# nDC: Number of Dependent Components (chemical species made of Independent Components)" << std::endl;
    stream << "<nDC>  " << data.Compounds.size() << std::endl;
    stream << "# nPH: Number of phases (into which Dependent Components are grouped)" << std::endl;
    stream << "<nPH>  " << data.Phases.size() << std::endl;
    stream << "# nPS: Number of phases-solutions (multicomponent phases) <= nPH" << std::endl;
    stream << "<nPS>  " << data.getSolutionsCount() << std::endl;
    stream << "# nDCs: Number of Dependent Components in phases-solutions <= nDC" << std::endl;
    stream << "<nDCs>  " << data.getCoumpoundsInSolutions() << std::endl;
    stream << std::endl;

    stream << "## (2) Dimensions for DBR node recipe (memory allocation)" << std::endl;
    stream << "# nICb: Number of ICs kept in the DBR file and DATABR memory structure (<= nIC)" << std::endl;
    stream << "<nICb>  " << data.totElements.size() << std::endl;
    stream << "# nDCb: Number of DCs kept in the DBR file and DATABR memory structure (<=nDC)" << std::endl;
    stream << "<nDCb>  " << data.Compounds.size() << std::endl;
    stream << "# nPHb: Number of phases kept in the DBR file and DATABR structure (<=nPH)" << std::endl;
    stream << "<nPHb>  " << data.Phases.size() << std::endl;
    stream << "# nPSb: Number of phases-solutions kept in the DBR file and DATABR structure (<=nPS)" << std::endl;
    stream << "<nPSb>  " << data.getSolutionsCount() << std::endl;
    stream << std::endl;

    stream << "## (3) Dimensions for thermodynamic data arrays" << std::endl;
    stream << "# nTp: Number of temperature grid points in lookup arrays for data interpolation, >=1" << std::endl;
    stream << "<nTp>  " << "1" << std::endl;
    stream << "# nPp: Number of pressure grid points in lookup arrays for data interpolation, >=1" << std::endl;
    stream << "<nPp>  " << "1" << std::endl;
    stream << "# iGrd: Flag for allocation of array of diffusition coefficients in DATACH structure (DCH file)" << std::endl;
    stream << "<iGrd>  " << "0" << std::endl;
    stream << "# fAalp: Flag for keeping specific surface areas of phases in DATABR structure (1) or ignoring them (0)" << std::endl;
    stream << "<fAalp>  " << "1" << std::endl;
    stream << "# mLook: Lookup mode: 0 interpolation over nTp*nPp grid; 1 data for T,P pairs, no interpolation" << std::endl;
    stream << "<mLook>  " << "0" << std::endl;

    stream << std::endl;
    stream << "<END_DIM>" << std::endl;
    stream << std::endl;

    stream << "## (4) DBR node recipe connection index lists" << std::endl;
    stream << "# xIC: DATACH access index list for ICs kept in the DATABR structure and in DBR files [nICb]" << std::endl;
    stream << "<xic>" << std::endl;
    stream << data.getICindexes() << std::endl;
    stream << "# xDC: DATACH access index list of DCs kept in the DATABR  structure and in DBR files [nDCb]" << std::endl;
    stream << "<xdc>" << std::endl;
    stream << data.getDCindexes() << std::endl;
    stream << "# xPH: DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]" << std::endl;
    stream << "<xph>" << std::endl;
    stream << data.getPhaseindexes() << std::endl;
    stream << std::endl;

    stream << "## (5) Independent Components and their properties" << std::endl;
    stream << "# ICNL: List of Independent Component names (<=4 characters per name) [nIC]" << std::endl;
    stream << "<ICNL>" << std::endl;
    stream << data.getICNames() << std::endl;
    stream << "# ccIC: Class codes of ICs (Independent Components) [nIC]" << std::endl;
    stream << "<ccIC>" << std::endl;
    stream << data.getICCodes() << std::endl;
    stream << "# ICmm: Atomic (molar) masses of ICs,  kg/mol [nIC]" << std::endl;
    stream << "<ICmm>" << std::endl;
    stream << data.getICAtomicMasses() << std::endl;
    stream << std::endl;

    stream << "## (6) Dependent Components and their codes" << std::endl;
    stream << "# DCNL: Name list of Dependent Components (<=16 characters per name) [nDC]" << std::endl;
    stream << "<DCNL>" << std::endl;
    stream << data.getDCNames() << std::endl;
    //'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'T' 'W' 'G' 'G' 'G' 'G' 'O' 'O' 'O' 'O' 'O'
    stream << "# ccDC: Class codes of DCs (Dependent Components) [nDC]" << std::endl;
    stream << "<ccDC>" << std::endl;
    stream << data.getDCCodes() << std::endl;
    stream << "# DCmm: Molar masses of DCs, kg/mol [nDC]" << std::endl;
    stream << "<DCmm>" << std::endl;
    stream << data.getDCAtomicMasses() << std::endl;
    stream << std::endl;

    stream << "## (7) Phases and their codes" << std::endl;
    stream << "# PHNL: List of Phase names (<=16 characters per name) [nPH]" << std::endl;
    stream << "<PHNL>" << std::endl;
    stream << data.getPhasesNames() << std::endl;
    stream << "# ccPH: Codes of phase aggregate state [nPH]" << std::endl;
    stream << "<ccPH>" << std::endl;
    stream << data.getPhasesCodes() << std::endl;
    stream << "# nDCinPH: Number of DCs included in each phase [nPH]" << std::endl;
    stream << "<nDCinPH>" << std::endl;
    stream << data.getCDinPhases() << std::endl;
    stream << std::endl;

    stream << "# (8) Data for Dependent Components" << std::endl;
    stream << "# A: Stoichiometry matrix A (expanded formulae) for DCs [nDC*nIC]" << std::endl;
    stream << "<A>" << std::endl;
    stream << data.getStohiometryMatrix() << std::endl;


    stream << "## (9) Thermodynamic data for Dependent Components" << std::endl;
    stream << "# Ttol: Tolerance for the temperature interpolation, K" << std::endl;
    stream << "<Ttol>  " << "0.1" << std::endl;

    stream << "# TKval: Temperature values, K for lookup arrays of thermodynamic data [nTp]" << std::endl;
    stream << "<TKval>" << std::endl;
    stream << data.T << std::endl;

    stream << "# Psat: Pressure Pa at saturated H2O vapour at given temperature [nTp]" << std::endl;
    stream << "<Psat>" << std::endl;
    stream << "1e-05" << std::endl;
    stream << std::endl;

    stream << "# Ptol: Tolerance for the pressure interpolation, Pa" << std::endl;
    stream << "<Ptol>  " << "50000" << std::endl;

    stream << "# Pval: Pressure values, Pa for lookup arrays of thermodynamic data [nPp]" << std::endl;
    stream << "<Pval>" << std::endl;
    stream << data.P << std::endl;
    stream << std::endl;

    stream << "# denW: Look-up array for the density of water-solvent, kg/m3, and its derivatives [5*nPp*nTp]" << std::endl;
    stream << "<denW>" << std::endl;
    stream << "998.356005405459" << std::endl;
    stream << "-0.351245995240374" << std::endl;
    stream << "-0.00744141355821083" << std::endl;
    stream << "0.04320907388247" << std::endl;
    stream << "0" << std::endl;
    stream << std::endl;

    stream << "# denWg: Look-up array for the density of water vapour, kg/m3, and its derivatives [5*nPp*nTp]" << std::endl;
    stream << "<denWg>" << std::endl;
    stream << "0" << std::endl;
    stream << "0" << std::endl;
    stream << "0" << std::endl;
    stream << "0" << std::endl;
    stream << "0" << std::endl;
    stream << std::endl;

    stream << "# epsW: Look-up array for the dielectric constant of water-solvent and its derivatives [5*nPp*nTp]" << std::endl;
    stream << "<epsW>" << std::endl;
    stream << "74.0844379009381" << std::endl;
    stream << "-0.33831198731748 " << std::endl;
    stream << "0.00141618917416284 " << std::endl;
    stream << "0.00385738431965382 " << std::endl;
    stream << "0.0" << std::endl;
    stream << std::endl;

    stream << "# epsWg: Look-up array for the dielectric constant of water vapour and its derivatives [5*nPp*nTp]" << std::endl;
    stream << "<epsWg>" << std::endl;
    stream << "0" << std::endl;
    stream << "0 " << std::endl;
    stream << "0 " << std::endl;
    stream << "0 " << std::endl;
    stream << "0" << std::endl;
    stream << std::endl;

    stream << "# V0: Look-up array for DC (standard) molar volumes, J/Pa [nDC*nPp*nTp]" << std::endl;
    stream << "<V0>" << std::endl;
    stream << data.getMolarVolumes() << std::endl;

    stream << "# G0: Look-up array for DC molar Gibbs energy function g(T,P), J/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<G0>" << std::endl;
    stream << data.getMolarGibbsEnery() << std::endl;

    stream << "# H0: Look-up array for DC molar enthalpy h(T,P), J/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<H0>" << std::endl;
    stream << data.getMolarEnthalpy() << std::endl;

    stream << "# S0: Look-up array for DC absolute entropy S(T,P), J/K/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<S0>" << std::endl;
    stream << data.getMolarEntropy() << std::endl;

    stream << "# Cp0: Look-up array for DC heat capacity Cp(T,P), J/K/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<Cp0>" << std::endl;
    stream << data.getMolarHeatCapacity() << std::endl;

    stream << "# A0: reserved: Look-up array for DC Helmholtz energy function, J/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<A0>" << std::endl;
    stream << data.getMolarHelmholzEnergy() << std::endl;

    stream << "# U0: reserved: Look-up array for DC internal energy function, J/mol [nDC*nPp*nTp]" << std::endl;
    stream << "<U0>" << std::endl;
    stream << data.getMolarInternalEnergy() << std::endl;

    stream << "# End of file" << std::endl;
    file << stream.str();
    file.close();
}

void Reader::generateIPMfile(SYSDATA & data, std::string& fileOut) {
    std::fstream file;
    std::string text = "";
    std::stringstream stream;

    file.open(fileOut.c_str(), std::fstream::out);

    stream.str("");
    stream.clear();

    stream << "# IPM text input file for the internal GEM IPM3 kernel data" << std::endl;
    stream << "# (should be read after the DCH file and before DBR files)" << std::endl;
    stream << std::endl;

    stream << "# ID key of the initial chemical system definition" << std::endl;
    stream << "<ID_key> " << "\"" << data.project << "\"" << std::endl;
    stream << std::endl;

    stream << "# PE: Flag for using electroneutrality condition in GEM IPM calculations (1 or 0)" << std::endl;
    stream << "<pa_PE>  " <<  "1" << std::endl;
    stream << std::endl;

    stream << "# PV: Flag for the volume balance constraint (on Vol IC) for indifferent equilibria at P_Sat (0 or 1)" << std::endl;
    stream << "<PV>  " <<  "1" << std::endl;
    stream << std::endl;


    stream << "# PSOL: Total number of DCs in liquid hydrocarbon phases (0; reserved)" << std::endl;
    stream << "<PSOL>  " <<  "0" << std::endl;
    stream << std::endl;

    stream << "# PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases" << std::endl;
    stream << "<PAalp>  " <<  "'+'" << std::endl;
    stream << std::endl;

    stream << "# PSigm: Flag for using (+) or ignoring (-) specific surface free energies" << std::endl;
    stream << "<PSigm>  " <<  "'+'" << std::endl;
    stream << std::endl;


    stream << "## (2) Dimensionalities that affect memory allocation" << std::endl;
    stream << "# Lads: Total number of Dependent Components in sorption phases included into this system" << std::endl;
    stream << "<Lads> " <<  "0" << std::endl;

    stream << "# FIa: Number of sorption phases included in this system (0 if no sorption phases )" << std::endl;
    stream << "<FIa>  " <<  "0" << std::endl;

    stream << "# FIat: Maximum number of surface types per adsorption phase (if FIa > 0, set FIat = 6)" << std::endl;
    stream << "<FIat>  " <<  "0" << std::endl;

    stream << std::endl;
    stream << "<END_DIM>" << std::endl;
    stream << std::endl;

    stream << "## (3) Numerical controls and tolerances of GEM IPM-3 kernel" << std::endl;
    stream << "#      - Need to be changed only in special cases (see gems3k_ipm.html)" << std::endl;
    stream << "# DB: Minimum amount of IC in the bulk composition, moles (except charge Zz) { 1e-17 }" << std::endl;
    stream << "<pa_DB>  " <<  "1e-17" << std::endl;
    stream << std::endl;

    stream << "# DHB: Maximum allowed relative mass balance residual for ICs { 1e-13 }" << std::endl;
    stream << "<pa_DHB>  " <<  "1e-13" << std::endl;
    stream << std::endl;

    stream << "# EPS: Tolerance of the SolveSimplex() balance residual for ICs { 1e-10 }" << std::endl;
    stream << "<pa_EPS>  " <<  "1e-10" << std::endl;
    stream << std::endl;

    stream << "# DK: Tolerance for the Dikin's criterion of IPM convergence { 1e-6 }" << std::endl;
    stream << "<pa_DK>  " <<  "1e-06" << std::endl;
    stream << std::endl;

    stream << "# DS: Cutoff minimum amount of stable phase in GEM IPM primal solution, moles { 1e-20 }" << std::endl;
    stream << "<pa_DS>  " <<  "1e-20" << std::endl;
    stream << std::endl;

    stream << "# DF: Tolerance DF of the stability criterion for a lost phase to be inserted to mass balance { 0.01 }" << std::endl;
    stream << "<pa_DF>  " <<  "0.01" << std::endl;
    stream << std::endl;

    stream << "# DFM: Tolerance for stability criterion for a phase to be eliminated from mass balance { 0.01 }" << std::endl;
    stream << "<pa_DFM>  " <<  "0.001" << std::endl;
    stream << std::endl;

    stream << "# DP: Maximal number of iterations in MassBalanceRefinement MBR() procedure { 130 }" << std::endl;
    stream << "<pa_DP>  " <<  "130" << std::endl;
    stream << std::endl;

    stream << "# IIM: Maximum allowed number of iterations in one main GEM IPM descent run { 7000 }" << std::endl;
    stream << "<pa_IIM>  " <<  "7000" << std::endl;
    stream << std::endl;

    stream << "# PD: Mode of calculation of DC activity coefficients ( 1 -IPM, 2 +MBR, 3 IPM ) { 2 }" << std::endl;
    stream << "<pa_PD>  " <<  "2" << std::endl;
    stream << std::endl;

    stream << "# PRD: Disable (0) or activate (-4 or less- max.dec.exp.for DC amount correction) SpeciationCleanup() { -5 }" << std::endl;
    stream << "<pa_PRD>  " <<  "-5" << std::endl;
    stream << std::endl;

    stream << "# AG: Smoothing parameter 1 for non-ideal primal chemical potential increments (-1 to +1) { 1.0 }" << std::endl;
    stream << "<pa_AG>  " <<  "1" << std::endl;
    stream << std::endl;

    stream << "# DGC: Smoothing parameter 2- exponent in smoothing function (-1 to +1) { 1 or 0.001 for adsorption }" << std::endl;
    stream << "<pa_DGC>  " <<  "0" << std::endl;
    stream << std::endl;

    stream << "# PSM: Level of diagnostic messages { 0- disabled (no ipmlog file); 1- default; 2-including warnings }" << std::endl;
    stream << "<pa_PSM>  " <<  "1" << std::endl;
    stream << std::endl;

    stream << "# GAR: Activity coefficient for major (M) species in solution phases at Simplex LP AIA { 1 }" << std::endl;
    stream << "<pa_GAR>  " <<  "1" << std::endl;
    stream << std::endl;

    stream << "# GAH: Activity coefficient for minor (J) species in solution phases at Simplex LP AIA { 1000 }" << std::endl;
    stream << "<pa_GAH>  " <<  "1000" << std::endl;
    stream << std::endl;

    stream << "# _Min: Cutoff amounts for elimination of: Xw - water-solvent { 1e-11 }; Sc - solid sorbent {1e-11}" << std::endl;
    stream << "#       Dc - solution- or surface species { 1e-30 }; Ph - non-electrolyte solution phase with all its components { 1e-20 }" << std::endl;
    stream << "# XwMin: Cutoff mole amount of water-solvent for aqueous phase elimination { 1e-13 }" << std::endl;
    stream << "<pa_XwMin>  " <<  "1e-13" << std::endl;
    stream << std::endl;

    stream << "# ScMin: Cutoff mole amount of solid sorbent for sorption phase elimination { 1e-13 }" << std::endl;
    stream << "<pa_ScMin>  " <<  "1e-13" << std::endl;
    stream << std::endl;

    stream << "# DcMin: Cutoff mole amount for elimination of DC (species) in multi-component phase { 1e-33 }" << std::endl;
    stream << "<pa_DcMin>  " <<  "1e-33" << std::endl;
    stream << std::endl;

    stream << "# PhMin: Cutoff mole amount for elimination of solution phases other than aqueous { 1e-20 }" << std::endl;
    stream << "<pa_PhMin>  " <<  "1e-20" << std::endl;
    stream << std::endl;

    stream << "# ICmin: Cutoff effective molal ionic strength for calculation of aqueous activity coefficients { 1e-5 }" << std::endl;
    stream << "<pa_ICmin>  " <<  "1e-05" << std::endl;
    stream << std::endl;

    stream << "# PC: Mode of Phase Selection: 1 old (Select-2), 2 new (PSSC), default { 2 }" << std::endl;
    stream << "<pa_PC>  " <<  "2" << std::endl;
    stream << "# DFY: Insertion mole amounts used after the LPP AIA and in PhaseSelection() algorithm" << std::endl;
    stream << "# DFYw: Insertion mole amount for water-solvent at Simplex()->MBR() bridge { 1e-5 }" << std::endl;
    stream << "<pa_DFYw>  " <<  "1e-05" << std::endl;
    stream << "# DFYaq: Insertion mole amount for aqueous species at Simplex()->MBR() bridge { 1e-5 }" << std::endl;
    stream << "<pa_DFYaq>  " <<  "1e-05" << std::endl;
    stream << std::endl;

    stream << "# DFYid: Insertion mole amount for DCs of ideal solution phases at Simplex()->MBR() bridge { 1e-5 }" << std::endl;
    stream << "<pa_DFYid>  " <<  "1e-05" << std::endl;
    stream << "# DFYr: Insertion mole amount for major DCs in solution phases at Simplex()->MBR()bridge { 1e-5 }" << std::endl;
    stream << "<pa_DFYr>  " <<  "1e-05" << std::endl;
    stream << "# DFYh: Insertion mole amount for junior DCs in solution phases Simplex()->MBR() bridge{ 1e-5 }" << std::endl;
    stream << "<pa_DFYh>  " <<  "1e-05" << std::endl;
    stream << "# DFYc: Insertion mole amount for single-component phase at Simplex()->MBR() bridge { 1e-5 }" << std::endl;
    stream << "<pa_DFYc>  " <<  "1e-05" << std::endl;
    stream << "# DFYs: Insertion mole amount for single-component phase in PSSC() algorithm { 1e-6 }" << std::endl;
    stream << "<pa_DFYs>  " <<  "1e-06" << std::endl;
    stream << std::endl;

    stream << "# Tolerances and controls of high-precision IPM-3 algorithm" << std::endl;
    stream << "# DW: Activate (1) or disable (0) error condition on maximum number of MBR() iterations DP { 1 }" << std::endl;
    stream << "<pa_DW>  " <<  "1" << std::endl;
    stream << "# DT: use DHB as relative maximum mass balance cutoff for all ICs (0), default, or for major ICs:" << std::endl;
    stream << "# decimal exponent (<-6) applied to DHB cutoff; (1) use DHB also as an absolute cutoff { 1 }" << std::endl;
    stream << "<pa_DT>  " <<  "0" << std::endl;
    stream << std::endl;

    stream << "# GAS: Threshold for primal-dual chemical potential difference used in SpeciationCleanup() { 0.0001 }" << std::endl;
    stream << "<pa_GAS>  " <<  "0.001" << std::endl;
    stream << "# Total number of moles used in internal re-scaling of the system (disabled if < 1e-4) { 1e3 }" << std::endl;
    stream << "<pa_DG>  " <<  "1000" << std::endl;
    stream << "# DNS: Standard surface number density, nm-2 for calculating activity of surface species { 12.05 }" << std::endl;
    stream << "<pa_DNS>  " <<  "12.05" << std::endl;
    stream << "# IEPS: Tolerance for calculation of surface activity coefficient terms for surface species { 1e-3 }" << std::endl;
    stream << "<pa_IEPS>  " <<  "0.001" << std::endl;
    stream << std::endl;

    stream << "# pKin: Flag for using metastability constraints on DC amounts in primal GEM solution { 1 }" << std::endl;
    stream << "<pKin>  " <<  "1" << std::endl;
    stream << "# DKIN: Tolerance for non-trivial metastability constraints on DC amounts, moles { 1e-10 }" << std::endl;
    stream << "<pa_DKIN>  " <<  "1e-10" << std::endl;
    stream << "# pa_PLLG: Tolerance for checking divergence in IPM dual solution, 1 to 30000 { 3000 }, 0 disables" << std::endl;
    stream << "<pa_PLLG>  " <<  "10000" << std::endl;
    stream << "# tMin: Type of thermodynamic potential to minimize (reserved)" << std::endl;
    stream << "<tMin>  " <<  "0" << std::endl;
    stream << std::endl;

    stream << "## (4) Initial data for multicomponent phases (see DCH file for dimension nPHs)" << std::endl;
    stream << std::endl;

    stream << "# sMod: Codes for TSolMod built-in  models of mixing in multicomponent phases [nPS*8]" << std::endl;
    stream << "<sMod>" << std::endl;
    stream << data.getMixModelCodes() << std::endl;
    stream << std::endl;

    stream << "# LsMod: Dimensions of TSolMod <IPxPH> and <PMc> data arrays [nPS*3]. In each row (for phase)" << std::endl;
    stream << "# [0] number of interaction parameters (rows in <IPx>); [1] max. parameter order (columns in <IPx>)" << std::endl;
    stream << "# [2] number of coefficients per interaction parameter in <PMc> array" << std::endl;
    stream << "<LsMod>" << std::endl;
    stream << data.getSolModDimentions() << std::endl;
    stream << std::endl;

    /*
    stream << "# PMc: Tables (in TSolMod convention) of interaction parameter coefficients  for non-ideal solutions" << std::endl;
    stream << "<PMc>" << std::endl;
    stream << data.getInteractionParameters() << std::endl;
    stream << std::endl;
    */

    stream << "# LsMdc: Dimensions of TSolMod <DMc> and <MoiSN> arrays [nPS*3]: In each row (for phase):" << std::endl;
    stream << "# [0] number of parameters per component; [1] 0; [2] 0. For multi-site (sublattice) models:" << std::endl;
    stream << "#   [1] number of sublattices nS; [2] total number of moieties nM acting in sublattice sites" << std::endl;
    stream << "<LsMdc>" << std::endl;
    stream << data.getDMcMoiSNDimentions() << std::endl;
    stream << std::endl;

    stream << "## (5) Data arrays which are provided neither in DCH nor in DBR files" << std::endl;
    stream << std::endl;

    stream << "# B: Full total bulk composition (vector b), moles [nIC] (will be partially re-written from DBR files)" << std::endl;
    stream << "<B>" << std::endl;
    stream << data.getICComposition() << std::endl;
    stream << std::endl;

    stream << "# Initial data for DCs - see DATACH file for dimensions nDC, nDCs" << std::endl;
    stream << "# Pparc: Partial pressures or fugacities of pure Dependent Components [nDC] (reserved)" << std::endl;
    stream << "<Pparc>" << std::endl;
    stream << data.getDCPartialPressures() << std::endl;
    stream << std::endl;

    stream << "# fDQF: DQF parameters of end members or pure gas fugacities, (J/mol/(RT) [nDC]" << std::endl;
    stream << "<fDQF>" << std::endl;
    stream << data.getDCFugacities() << std::endl;
    stream << std::endl;

    stream << "# lnGmf: Natural logarithms of DC activity coefficients used at Simplex LP approximation only [nDC]" << std::endl;
    stream << "<lnGmf>" << std::endl;
    stream << data.getDClnActivities() << std::endl;
    stream << std::endl;

    stream << "# (6) Metastability constraints on DC amounts from above (DUL) and below (DLL)" << std::endl;
    stream << std::endl;

    stream << "# RLC: Code of metastability constraints for DCs {L U B (default)} [nDC]" << std::endl;
    stream << "<RLC>" << std::endl;
    stream << data.getDCMetastabilityConstraints() << std::endl;
    stream << std::endl;

    stream << "# RSC: Units of metastability/kinetic constraints for DCs {M} moles [nDC]" << std::endl;
    stream << "<RSC>" << std::endl;
    stream << data.getDCMetastabilityUnits() << std::endl;
    stream << std::endl;

    stream << "# DLL: Lower metastability constraints on DC amounts <xDC>, moles [nDC] (default: 0)" << std::endl;
    stream << "<DLL>" << std::endl;
    stream << data.getDCLowerMetastConstr() << std::endl;
    stream << std::endl;

    stream << "# DUL: Upper metastability constraints on DC amounts <xDC>, moles [nDC] (default: 1e6)" << std::endl;
    stream << "<DUL>" << std::endl;
    stream << data.getDCUpperMetastConstr() << std::endl;
    stream << std::endl;


    stream << "# (7) Initial data for Phases" << std::endl;
    stream << std::endl;

    stream << "# Aalp: Specific surface areas of phases, m2/g [nPH]" << std::endl;
    stream << "<Aalp>" << std::endl;
    stream << data.getPHSpecifSurfaces() << std::endl;
    stream << std::endl;

    stream << "# Sigw: Specific surface free energy for phase-water interface, J/m2 [nPH] (reserved)" << std::endl;
    stream << "<Sigw>" << std::endl;
    stream << data.getPHSurfaceEnergyOfPHWater() << std::endl;
    stream << std::endl;

    stream << "# Sigg: Specific surface free energy for phase-gas interface, J/m2 (not yet used) [nPH]" << std::endl;
    stream << "<Sigg>" << std::endl;
    stream << data.getPHSurfaceEnergyOfPHGas() << std::endl;
    stream << std::endl;

    stream << "# YOF: Surface free energy parameter for phases in J/g (to accomodate for variable phase composition)  [nPH]" << std::endl;
    stream << "<YOF>" << std::endl;
    stream << data.getPHSurfaceEnergyParameter() << std::endl;
    stream << std::endl;

    stream << "# dcMod: Codes for PT corrections of DC thermodynamic data [nDC] (reserved)" << std::endl;
    stream << "<dcMod>" << std::endl;
    stream << data.getDCPTCorrectionCodes() << std::endl;
    stream << std::endl;

    stream << "# mui: IC indices in parent RMULTS IC list (not used in standalone GEMS3K)" << std::endl;
    stream << "<mui>" << std::endl;
    stream << data.getICindexes() << std::endl;
    stream << std::endl;

    stream << "# muk: Phase indices in parent RMULTS Phase list (not used in standalone GEMS3K)" << std::endl;
    stream << "<muk>" << std::endl;
    stream << data.getPhaseindexes() << std::endl;
    stream << std::endl;

    stream << "# muj: DC indices in parent RMULTS DC list (not used in standalone GEMS3K)" << std::endl;
    stream << "<muj>" << std::endl;
    stream << data.getDCindexes() << std::endl;
    stream << std::endl;

    stream << "# End of file" << std::endl;

    file << stream.str();
    file.close();

}

void Reader::generateDBRfile(SYSDATA & data, std::string& fileOut) {
    std::fstream file;
    std::string text = "";
    std::stringstream stream;

    file.open(fileOut.c_str(), std::fstream::out);

    stream.str("");
    stream.clear();

    stream << "#  GEMS3K v.3.3 r.1036 (rc) " << std::endl;
    stream << "# File: ../data/dbr.dat" << std::endl;
    stream << "# Comments can be marked with # $ ; as the first character in the line" << std::endl;
    stream << "# DBR text input file for node system recipe and speciation data" << std::endl;
    stream << "# (should be read only after the DCH and the IPM files)" << std::endl;
    stream << std::endl;

    stream << "# (1): Flags controlling GEM IPM-3 operation and data exchange" << std::endl;
    stream << "# NodeHandle: Node identification handle" << std::endl;
    stream << "<NodeHandle> " <<  "0" << std::endl;

    stream << "# NodeTypeHY:  Node type code (hydraulic), not used on TNode level; see typedef NODETYPE" << std::endl;
    stream << "<NodeTypeHY> " <<  "0" << std::endl;

    stream << "# NodeTypeMT:  Node type (mass transport), not used on TNode level; see typedef NODETYPE" << std::endl;
    stream << "<NodeTypeMT> " <<  "0" << std::endl;

    stream << "# NodeStatusFMT:  Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT" << std::endl;
    stream << "<NodeStatusFMT> " <<  "-1" << std::endl;

    stream << "# NodeStatusCH: Node status code and control in GEM input and output; see typedef NODECODECH" << std::endl;
    stream << "<NodeStatusCH> " <<  "1" << std::endl;
    stream << std::endl;

    stream << "# IterDone:  Number of iterations performed by GEM IPM in the last run (GEM output)" << std::endl;
    stream << "<IterDone>  1" << std::endl;

    stream << "## (2) Chemical scalar properies of the node system" << std::endl;
    stream << "# TK: Node temperature T, Kelvin. This value must always be provided (GEM input)" << std::endl;
    stream << "<TK> " <<  data.T << std::endl;

    stream << "# P:  Node Pressure P, Pa. This value must always be provided (GEM input)" << std::endl;
    stream << "<P> " <<  data.P << std::endl;
    stream << std::endl;


    stream << "### Arrays: for dimensions and index lists, see Section (2) of DCH file" << std::endl;
    stream << std::endl;

    stream << "## (4) Data for Independent Components" << std::endl;
    stream << "#  " << data.getICNames() << std::endl;
    stream << "# bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]" << std::endl;
    stream << "<bIC>" << std::endl;
    stream << data.getICComposition() << std::endl;
    stream << std::endl;

    stream << "# End of file" << std::endl;
    file << stream.str();
    file.close();
}

void Reader::setupData(SYSDATA & data) {
    generateDHC(data);
    generateIPM(data);
    //generateDBR(data);
}

void Reader::clean() {

}

void Reader::generateDHC(SYSDATA & data) {

    //(1) Dimensions for memory allocation" << std::endl;
    // nIC: Number of Independent Components (usually chemical elements and charge)
    data.inputDHC.nIC   = data.totElements.size();
    // nDC: Number of Dependent Components (chemical species made of Independent Components)
    data.inputDHC.nDC   = data.Compounds.size();
    // nPH: Number of phases (into which Dependent Components are grouped)
    data.inputDHC.nPH   = data.Phases.size();
    // nPS: Number of phases-solutions (multicomponent phases) <= nPH
    data.inputDHC.nPS   = data.getSolutionsCount();
    // nDCs: Number of Dependent Components in phases-solutions <= nDC
    data.inputDHC.nDCs  = data.getCoumpoundsInSolutions();

    // (2) Dimensions for DBR node recipe (memory allocation)
    // nICb: Number of ICs kept in the DBR file and DATABR memory structure (<= nIC)
    data.inputDHC.nICb  = data.totElements.size();
    // nDCb: Number of DCs kept in the DBR file and DATABR memory structure (<=nDC)
    data.inputDHC.nDCb  = data.Compounds.size();
    // nPHb: Number of phases kept in the DBR file and DATABR structure (<=nPH)
    data.inputDHC.nPHb  = data.Phases.size();
    // nPSb: Number of phases-solutions kept in the DBR file and DATABR structure (<=nPS)
    data.inputDHC.nPSb  = data.getSolutionsCount();

    // (3) Dimensions for thermodynamic data arrays
    // nTp: Number of temperature grid points in lookup arrays for data interpolation, >=1
    data.inputDHC.nTp   = 1;
    // nPp: Number of pressure grid points in lookup arrays for data interpolation, >=1
    data.inputDHC.nPp   = 1;
    // iGrd: Flag for allocation of array of diffusition coefficients in DATACH structure (DCH file)
    data.inputDHC.iGrd  = 0;
    // mLook: Lookup mode: 0 interpolation over nTp*nPp grid; 1 data for T,P pairs, no interpolation
    data.inputDHC.mLook = 0;

    allocateDHC(data);
    allocateDBR(data);

    // (4) DBR node recipe connection index lists
    // xIC: DATACH access index list for ICs kept in the DATABR structure and in DBR files [nICb]
    data.inputDHC.xic   = data.getICindexesArray();
    // xDC: DATACH access index list of DCs kept in the DATABR  structure and in DBR files [nDCb]
    data.inputDHC.xdc   = data.getDCindexesArray();
    // xPH: DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]
    data.inputDHC.xph   = data.getPhaseindexesArray();

    // (5) Independent Components and their properties
    // ICNL: List of Independent Component names (<=4 characters per name) [nIC]
    data.getICNamesArray(data.inputDHC.ICNL);
    // ccIC: Class codes of ICs (Independent Components) [nIC]
    data.inputDHC.ccIC = data.getICCodesArray();
    // ICmm: Atomic (molar) masses of ICs,  kg/mol [nIC]
    data.inputDHC.ICmm = data.getICAtomicMassesArray();


    // (6) Dependent Components and their codes
    // DCNL: Name list of Dependent Components (<=16 characters per name) [nDC]
    data.getDCNamesArray(data.inputDHC.DCNL);
    // ccDC: Class codes of DCs (Dependent Components) [nDC]
    data.inputDHC.ccDC = data.getDCCodesArray();

    // DCmm: Molar masses of DCs, kg/mol [nDC]
    data.inputDHC.DCmm = data.getDCAtomicMassesArray();

    // (7) Phases and their codes
    // PHNL: List of Phase names (<=16 characters per name) [nPH]
    data.getPhasesNamesArray(data.inputDHC.PHNL);
    // ccPH: Codes of phase aggregate state [nPH]
    data.inputDHC.ccPH = data.getPhasesCodesArray();
    // nDCinPH: Number of DCs included in each phase [nPH]
    data.inputDHC.nDCinPH = data.getCDinPhasesArray();

    // (8) Data for Dependent Components
    // A: Stoichiometry matrix A (expanded formulae) for DCs [nDC*nIC]
    data.inputDHC.A = data.getStohiometryMatrixArray();

    // (9) Thermodynamic data for Dependent Components
    // Ttol: Tolerance for the temperature interpolation, K
    data.inputDHC.Ttol = 0.1;
    // TKval: Temperature values, K for lookup arrays of thermodynamic data [nTp]
    data.inputDHC.TKval = &(data.T);
    // Psat: Pressure Pa at saturated H2O vapour at given temperature [nTp]
    data.inputDHC.Psat = new double{1.0e-5};
    // Ptol: Tolerance for the pressure interpolation, Pa
    data.inputDHC.Ptol = 50000;
    // Pval: Pressure values, Pa for lookup arrays of thermodynamic data [nPp]
    data.inputDHC.Pval = &data.P;

    // denW: Look-up array for the density of water-solvent, kg/m3, and its derivatives [5*nPp*nTp]
    data.inputDHC.denW = new double[5]{998.356005405459, -0.351245995240374, -0.00744141355821083, 0.04320907388247, 0.0};
    // denWg: Look-up array for the density of water vapour, kg/m3, and its derivatives [5*nPp*nTp]
    data.inputDHC.denWg = new double[5]{0.0, 0.0, 0.0, 0.0, 0.0};
    // epsW: Look-up array for the dielectric constant of water-solvent and its derivatives [5*nPp*nTp]
    data.inputDHC.epsW = new double[5]{74.0844379009381, -0.33831198731748, 0.00141618917416284, 0.00385738431965382, 0.0};
    // epsWg: Look-up array for the dielectric constant of water vapour and its derivatives [5*nPp*nTp]
    data.inputDHC.epsWg = new double[5]{0.0, 0.0, 0.0, 0.0, 0.0};
    // V0: Look-up array for DC (standard) molar volumes, J/Pa [nDC*nPp*nTp]
    data.inputDHC.V0 = data.getMolarVolumesArray();
    // G0: Look-up array for DC molar Gibbs energy function g(T,P), J/mol [nDC*nPp*nTp]
    data.inputDHC.G0 = data.getMolarGibbsEneryArray();
    // H0: Look-up array for DC molar enthalpy h(T,P), J/mol [nDC*nPp*nTp]
    data.inputDHC.H0 = data.getMolarEnthalpyArray();
    // S0: Look-up array for DC absolute entropy S(T,P), J/K/mol [nDC*nPp*nTp]
    data.inputDHC.S0 = data.getMolarEntropyArray();
    // Cp0: Look-up array for DC heat capacity Cp(T,P), J/K/mol [nDC*nPp*nTp]
    data.inputDHC.Cp0 = data.getMolarHeatCapacityArray();
    // A0: reserved: Look-up array for DC Helmholtz energy function, J/mol [nDC*nPp*nTp]
    data.inputDHC.A0 = data.getMolarHelmholzEnergyArray();
    // U0: reserved: Look-up array for DC internal energy function, J/mol [nDC*nPp*nTp]
    data.inputDHC.U0 = data.getMolarInternalEnergyArray();
}

void Reader::generateIPM(SYSDATA & data) {

    // IPM text input file for the internal GEM IPM3 kernel data
    // (should be read after the DCH file and before DBR files)

    // ID key of the initial chemical system definition
    strcpy( data.inputPM.stkey, data.project.c_str() );
    // PE: Flag for using electroneutrality condition in GEM IPM calculations (1 or 0)
    data.inputPM.E = 1.0;
    data.inputPA.p.PE = 1.0;
    // PV: Flag for the volume balance constraint (on Vol IC) for indifferent equilibria at P_Sat (0 or 1)
    data.inputPM.PV = 1.0;
    // PSOL: Total number of DCs in liquid hydrocarbon phases (0; reserved)
    data.inputPM.PSOL = 0.0;

    // PAalp: Flag for using (+) or ignoring (-) specific surface areas of phases
    //data.inputPM.PA = 0.0;
    //stream << "<PAalp>  " <<  "'-'

    // PSigm: Flag for using (+) or ignoring (-) specific surface free energies
    //stream << "<PSigm>  " <<  "'-'

    // (2) Dimensionalities that affect memory allocation
    // Lads: Total number of Dependent Components in sorption phases included into this system
    data.inputPM.Lads = 0.0;
    // FIa: Number of sorption phases included in this system (0 if no sorption phases )
    data.inputPM.FIa = 0.0;
    // FIat: Maximum number of surface types per adsorption phase (if FIa > 0, set FIat = 6)
    data.inputPM.FIat = 0.0;

    allocateIPM(data);

    // (3) Numerical controls and tolerances of GEM IPM-3 kernel
    //      - Need to be changed only in special cases (see gems3k_ipm.html)
    // DB: Minimum amount of IC in the bulk composition, moles (except charge Zz) { 1e-17 }
    data.inputPA.p.DB = 1.0e-17;
    // DHB: Maximum allowed relative mass balance residual for ICs { 1e-13 }
    data.inputPA.p.DHB = 1.0e-13;
    // EPS: Tolerance of the SolveSimplex() balance residual for ICs { 1e-10 }
    data.inputPA.p.EPS = 1.0e-10;
    // DK: Tolerance for the Dikin's criterion of IPM convergence { 1e-6 }
    data.inputPA.p.DK = 1.0e-6;
    // DS: Cutoff minimum amount of stable phase in GEM IPM primal solution, moles { 1e-20 }
    data.inputPA.p.DS = 1.0e-20;
    // DF: Tolerance DF of the stability criterion for a lost phase to be inserted to mass balance { 0.01 }
    data.inputPA.p.DF = 0.01;
    // DFM: Tolerance for stability criterion for a phase to be eliminated from mass balance { 0.01 }
    data.inputPA.p.DFM = 0.001;
    // DP: Maximal number of iterations in MassBalanceRefinement MBR() procedure { 130 }
    data.inputPA.p.DP = 130;
    // IIM: Maximum allowed number of iterations in one main GEM IPM descent run { 7000 }
    data.inputPA.p.IIM = 7000;
    // PD: Mode of calculation of DC activity coefficients ( 1 -IPM, 2 +MBR, 3 IPM ) { 2 }
    data.inputPA.p.PD = 2;
    // PRD: Disable (0) or activate (-4 or less- max.dec.exp.for DC amount correction) SpeciationCleanup() { -5 }
    data.inputPA.p.PRD = -5;
    // AG: Smoothing parameter 1 for non-ideal primal chemical potential increments (-1 to +1) { 1.0 }
    data.inputPA.p.AG = 1;
    // DGC: Smoothing parameter 2- exponent in smoothing function (-1 to +1) { 1 or 0.001 for adsorption }
    data.inputPA.p.DGC = 0;
    // PSM: Level of diagnostic messages { 0- disabled (no ipmlog file); 1- default; 2-including warnings }
    data.inputPA.p.PSM = 1;
    // GAR: Activity coefficient for major (M) species in solution phases at Simplex LP AIA { 1 }
    data.inputPA.p.GAR = 1;
    // GAH: Activity coefficient for minor (J) species in solution phases at Simplex LP AIA { 1000 }
    data.inputPA.p.GAH = 1000;
    // _Min: Cutoff amounts for elimination of: Xw - water-solvent { 1e-11 }; Sc - solid sorbent {1e-11}
    //       Dc - solution- or surface species { 1e-30 }; Ph - non-electrolyte solution phase with all its components { 1e-20 }
    // XwMin: Cutoff mole amount of water-solvent for aqueous phase elimination { 1e-13 }
    data.inputPA.p.XwMin = 1.0e-13;
    // ScMin: Cutoff mole amount of solid sorbent for sorption phase elimination { 1e-13 }
    data.inputPA.p.ScMin = 1.0e-13;
    // DcMin: Cutoff mole amount for elimination of DC (species) in multi-component phase { 1e-33 }
    data.inputPA.p.DcMin = 1.0e-33;
    // PhMin: Cutoff mole amount for elimination of solution phases other than aqueous { 1e-20 }
    data.inputPA.p.PhMin = 1.0e-20;
    // ICmin: Cutoff effective molal ionic strength for calculation of aqueous activity coefficients { 1e-5 }
    data.inputPA.p.ICmin = 1.0e-05;
    // PC: Mode of Phase Selection: 1 old (Select-2), 2 new (PSSC), default { 2 }
    data.inputPA.p.PC = 2;
    // DFY: Insertion mole amounts used after the LPP AIA and in PhaseSelection() algorithm
    // DFYw: Insertion mole amount for water-solvent at Simplex()->MBR() bridge { 1e-5 }
    data.inputPA.p.DFYw = 1.0e-05;
    // DFYaq: Insertion mole amount for aqueous species at Simplex()->MBR() bridge { 1e-5 }
    data.inputPA.p.DFYaq = 1.0e-05;
    // DFYid: Insertion mole amount for DCs of ideal solution phases at Simplex()->MBR() bridge { 1e-5 }
    data.inputPA.p.DFYid = 1.0e-05;
    // DFYr: Insertion mole amount for major DCs in solution phases at Simplex()->MBR()bridge { 1e-5 }
    data.inputPA.p.DFYr = 1.0e-05;
    // DFYh: Insertion mole amount for junior DCs in solution phases Simplex()->MBR() bridge{ 1e-5 }
    data.inputPA.p.DFYh = 1.0e-05;
    // DFYc: Insertion mole amount for single-component phase at Simplex()->MBR() bridge { 1e-5 }
    data.inputPA.p.DFYc = 1.0e-05;
    // DFYs: Insertion mole amount for single-component phase in PSSC() algorithm { 1e-6 }
    data.inputPA.p.DFYs = 1.0e-06;
    // Tolerances and controls of high-precision IPM-3 algorithm
    // DW: Activate (1) or disable (0) error condition on maximum number of MBR() iterations DP { 1 }
    data.inputPA.p.DW = 1;
    // DT: use DHB as relative maximum mass balance cutoff for all ICs (0), default, or for major ICs:
    // decimal exponent (<-6) applied to DHB cutoff; (1) use DHB also as an absolute cutoff { 1 }
    data.inputPA.p.DT = 0;
    // GAS: Threshold for primal-dual chemical potential difference used in SpeciationCleanup() { 0.0001 }
    data.inputPA.p.GAS = 0.001;
    // Total number of moles used in internal re-scaling of the system (disabled if < 1e-4) { 1e3 }
    data.inputPA.p.DG = 1000;
    // DNS: Standard surface number density, nm-2 for calculating activity of surface species { 12.05 }
    data.inputPA.p.DNS = 12.05;
    // IEPS: Tolerance for calculation of surface activity coefficient terms for surface species { 1e-3 }
    data.inputPA.p.IEPS = 0.001;
    // pKin: Flag for using metastability constraints on DC amounts in primal GEM solution { 1 }
    data.inputPM.PLIM = 1;
    // DKIN: Tolerance for non-trivial metastability constraints on DC amounts, moles { 1e-10 }
    data.inputPA.p.DKIN = 1.0e-10;
    // pa_PLLG: Tolerance for checking divergence in IPM dual solution, 1 to 30000 { 3000 }, 0 disables
    data.inputPA.p.PLLG = 10000;
    // tMin: Type of thermodynamic potential to minimize (reserved)
    data.inputPM.tMin = 0;

    // (4) Initial data for multicomponent phases (see DCH file for dimension nPHs)
    // sMod: Codes for TSolMod built-in  models of mixing in multicomponent phases [nPS*6]
    data.getMixModelCodesArray(data.inputPM.sMod);
    // LsMod: Dimensions of TSolMod <IPxPH> and <PMc> data arrays [nPS*3]. In each row (for phase)
    // [0] number of interaction parameters (rows in <IPx>); [1] max. parameter order (columns in <IPx>)
    // [2] number of coefficients per interaction parameter in <PMc> array
    data.inputPM.LsMod = data.getSolModDimentionsArray();

    // PMc: Tables (in TSolMod convention) of interaction parameter coefficients  for non-ideal solutions
    //stream << "<PMc>
    //stream << data.getInteractionParameters() << std::endl;
    //stream << std::endl;

    // LsMdc: Dimensions of TSolMod <DMc> and <MoiSN> arrays [nPS*3]: In each row (for phase):
    // [0] number of parameters per component; [1] 0; [2] 0. For multi-site (sublattice) models:
    //   [1] number of sublattices nS; [2] total number of moieties nM acting in sublattice sites
    //stream << "<LsMdc>
    //stream << data.getDMcMoiSNDimentions() << std::endl;
    //stream << std::endl;


    // (5) Data arrays which are provided neither in DCH nor in DBR files
    // B: Full total bulk composition (vector b), moles [nIC] (will be partially re-written from DBR files)
    data.inputPM.B = data.getICCompositionArray();

    // Initial data for DCs - see DATACH file for dimensions nDC, nDCs
    // Pparc: Partial pressures or fugacities of pure Dependent Components [nDC] (reserved)
    data.inputPM.Pparc = data.getDCPartialPressuresArray();
    // fDQF: DQF parameters of end members or pure gas fugacities, (J/mol/(RT) [nDC]
    data.inputPM.fDQF = data.getDCFugacitiesArray();
    // lnGmf: Natural logarithms of DC activity coefficients used at Simplex LP approximation only [nDC]
    data.inputPM.lnGmf = data.getDClnActivitiesArray();

    // (6) Metastability constraints on DC amounts from above (DUL) and below (DLL)
    // RLC: Code of metastability constraints for DCs {L U B (default)} [nDC]
    data.inputPM.RLC = data.getDCMetastabilityConstraintsArray();
    // RSC: Units of metastability/kinetic constraints for DCs {M} moles [nDC]
    data.inputPM.RSC = data.getDCMetastabilityUnitsArray();
    // DLL: Lower metastability constraints on DC amounts <xDC>, moles [nDC] (default: 0)
    data.inputPM.DLL = data.getDCLowerMetastConstrArray();
    // DUL: Upper metastability constraints on DC amounts <xDC>, moles [nDC] (default: 1e6)
    data.inputPM.DUL = data.getDCUpperMetastConstrArray();

    // (7) Initial data for Phases
    // Aalp: Specific surface areas of phases, m2/g [nPH]
    data.inputPM.Aalp = data.getPHSpecifSurfacesArray();
    // Sigw: Specific surface free energy for phase-water interface, J/m2 [nPH] (reserved)
    data.inputPM.Sigw = data.getPHSurfaceEnergyOfPHWaterArray();
    // Sigg: Specific surface free energy for phase-gas interface, J/m2 (not yet used) [nPH]
    data.inputPM.Sigg = data.getPHSurfaceEnergyOfPHGasArray();
    // YOF: Surface free energy parameter for phases in J/g (to accomodate for variable phase composition)  [nPH]
    data.inputPM.YOF = data.getPHSurfaceEnergyParameterArray();
    // dcMod: Codes for PT corrections of DC thermodynamic data [nDC] (reserved)
    data.getDCPTCorrectionCodesArray(data.inputPM.dcMod);
    // mui: IC indices in parent RMULTS IC list (not used in standalone GEMS3K)
    data.inputPM.mui = data.getICindexesArray();
    // muk: Phase indices in parent RMULTS Phase list (not used in standalone GEMS3K)
    data.inputPM.muk = data.getPhaseindexesArray();
    // muj: DC indices in parent RMULTS DC list (not used in standalone GEMS3K)
    data.inputPM.muk = data.getDCindexesArray();

    // End of file
}

void Reader::generateDBR(SYSDATA & data) {
    //  GEMS3K v.3.3 r.1036 (rc)
    // File: ../data/dbr.dat
    // Comments can be marked with # $ ; as the first character in the line
    // DBR text input file for node system recipe and speciation data
    // (should be read only after the DCH and the IPM files)

    // (1): Flags controlling GEM IPM-3 operation and data exchange
    // NodeHandle: Node identification handle
    data.inputBR.NodeHandle = 0;
    // NodeTypeHY:  Node type code (hydraulic), not used on TNode level; see typedef NODETYPE
    data.inputBR.NodeTypeHY = 0;
    // NodeTypeMT:  Node type (mass transport), not used on TNode level; see typedef NODETYPE
    data.inputBR.NodeTypeMT = 0;
    // NodeStatusFMT:  Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT
    data.inputBR.NodeStatusFMT = -1;
    // NodeStatusCH: Node status code and control in GEM input and output; see typedef NODECODECH
    data.inputBR.NodeStatusCH = 1;

    // IterDone:  Number of iterations performed by GEM IPM in the last run (GEM output)
    data.inputBR.IterDone = 1;
    // (2) Chemical scalar properies of the node system
    // TK: Node temperature T, Kelvin. This value must always be provided (GEM input)
    data.inputBR.TK = data.T;
    // P:  Node Pressure P, Pa. This value must always be provided (GEM input)
    data.inputBR.P = data.P;

    // Arrays: for dimensions and index lists, see Section (2) of DCH file

    // (4) Data for Independent Components
    // bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]
    data.inputBR.bIC = data.getICCompositionArray();
    // End of file
}

int Reader::gridTP(SYSDATA& data) {
    if( data.inputDHC.mLook == 1L )
      return data.inputDHC.nTp;
    else
      return (data.inputDHC.nPp * data.inputDHC.nTp);
}

void Reader::allocateDHC(SYSDATA& data) {
    if( data.inputDHC.mLook == 1 &&  (data.inputDHC.nPp != data.inputDHC.nTp) )
        Error( "No-interpolation mode",
           "Different number of points for temperature and pressure ");

     data.inputDHC.nDCinPH = new long int[data.inputDHC.nPH];

     if( data.inputDHC.nICb >0 )
       data.inputDHC.xic = new long int[data.inputDHC.nICb];
     else  data.inputDHC.xic = 0;
     if( data.inputDHC.nDCb >0 )
       data.inputDHC.xdc = new long int[data.inputDHC.nDCb];
     else  data.inputDHC.xdc = 0;
     if( data.inputDHC.nPHb >0 )
       data.inputDHC.xph = new long int[data.inputDHC.nPHb];
     else  data.inputDHC.xph = 0;

      data.inputDHC.A = new double[data.inputDHC.nIC*data.inputDHC.nDC];
      data.inputDHC.ICmm = new double[data.inputDHC.nIC];
      data.inputDHC.DCmm = new double[data.inputDHC.nDC];
    data.inputDHC.DCmm[0] = 0.0;   // Added by DK on 03.03.2007

      data.inputDHC.TKval = new double[data.inputDHC.nTp];
      data.inputDHC.Psat = new double[data.inputDHC.nTp];
      data.inputDHC.Pval = new double[data.inputDHC.nPp];

      data.inputDHC.denW = new double[ 5*gridTP(data)];
      data.inputDHC.denWg = new double[ 5*gridTP(data)];
      data.inputDHC.epsW = new double[ 5*gridTP(data)];
      data.inputDHC.epsWg = new double[ 5*gridTP(data)];

      data.inputDHC.G0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.V0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.H0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.S0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.Cp0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.A0 = new double[data.inputDHC.nDC*gridTP(data)];
      data.inputDHC.U0 = new double[data.inputDHC.nDC*gridTP(data)];

      if(  data.inputDHC.iGrd  )
           data.inputDHC.DD = new double[data.inputDHC.nDCs*gridTP(data)];
      else
           data.inputDHC.DD = 0;
      data.inputDHC.ICNL = new char[data.inputDHC.nIC][MaxICN];
      data.inputDHC.DCNL = new char[data.inputDHC.nDC][MaxDCN];
      data.inputDHC.PHNL = new char[data.inputDHC.nPH][MaxPHN];

      data.inputDHC.ccIC = new char[data.inputDHC.nIC];
      data.inputDHC.ccDC = new char[data.inputDHC.nDC];
      data.inputDHC.ccPH = new char[data.inputDHC.nPH];
}

void Reader::allocateDBR(SYSDATA& data) {
    long int j,k;
    data.inputBR.bIC = new double[data.inputDHC.nICb];
    data.inputBR.rMB = new double[data.inputDHC.nICb];
    data.inputBR.uIC = new double[data.inputDHC.nICb];
    data.inputBR.bSP = new double[data.inputDHC.nICb];

    for(  j=0; j<data.inputDHC.nICb; j++ ) {
        data.inputBR.rMB[j] = 0.;
        data.inputBR.uIC[j] = 0.;
        data.inputBR.bSP[j] = 0.;
     }

    data.inputBR.xDC = new double[data.inputDHC.nDCb];
    data.inputBR.gam = new double[data.inputDHC.nDCb];

    for(  j=0; j<data.inputDHC.nDCb; j++ ) {
      data.inputBR.xDC[j] = 0.;
      data.inputBR.gam[j] = 1.;
    }

    //  default assignment
    data.inputBR.dul = new double[data.inputDHC.nDCb];
    for(  j=0; j<data.inputDHC.nDCb; j++ )
        data.inputBR.dul[j] = 1.0e6;            // default assignment

    data.inputBR.dll = new double[data.inputDHC.nDCb];
    for(  j=0; j<data.inputDHC.nDCb; j++ )
        data.inputBR.dll[j] = 0.0;              // default assignment

    if( data.inputDHC.nAalp >0 )
    {
        data.inputBR.aPH = new double[data.inputDHC.nPHb];
        for(  k=0; k<data.inputDHC.nPHb; k++ )
            data.inputBR.aPH[k] = 0.0;       // default assignment
    }

    else
        data.inputBR.aPH = 0;

   data.inputBR.xPH = new double[data.inputDHC.nPHb];
   data.inputBR.omPH = new double[data.inputDHC.nPHb];

   for(  k=0; k<data.inputDHC.nPHb; k++ ) {
       data.inputBR.xPH[k] = 0.0;       // default assignment
       data.inputBR.omPH[k] = 0.0;
   }

   data.inputBR.vPS = new double[data.inputDHC.nPSb];
   data.inputBR.mPS = new double[data.inputDHC.nPSb];
   data.inputBR.bPS = new double[data.inputDHC.nPSb*data.inputDHC.nICb];
   data.inputBR.xPA = new double[data.inputDHC.nPSb];
   data.inputBR.amru = new double[data.inputDHC.nPSb];
   data.inputBR.amrl = new double[data.inputDHC.nPSb];

   for(  k=0; k<data.inputDHC.nPSb; k++ ) {
       data.inputBR.vPS[k] = 0.0;
       data.inputBR.mPS[k] = 0.0;
       data.inputBR.xPA[k] = 0.0;
       for(  j=0; j<data.inputDHC.nICb; j++ )
          data.inputBR.bPS[k*data.inputDHC.nICb+j] = 0.0;
       data.inputBR.amru[k] = 1.0e6;
       data.inputBR.amrl[k] = 0.0;
   }
}

void Reader::allocateIPM(SYSDATA & data) {
    long int ii, jj ;

    if( data.inputPM.N < 2 || data.inputPM.L < 2 || data.inputPM.FI < 1 )
        std::cout << "data.inputPM.N < 2 || data.inputPM.L < 2 || data.inputPM.FI < 1" << std::endl;

     // need  always to alloc vectors
     data.inputPM.L1 = new long int[data.inputPM.FI];
     data.inputPM.muk = new long int[data.inputPM.FI];
     for( ii=0; ii<data.inputPM.FI; ii++)
     {   data.inputPM.L1[ii] = 0;
         data.inputPM.muk[ii] = ii;
     }
     data.inputPM.mui = new long int[data.inputPM.N];
     for( ii=0; ii<data.inputPM.N; ii++)
       data.inputPM.mui[ii] = ii;
     data.inputPM.muj = new long int[data.inputPM.L];
     for( ii=0; ii<data.inputPM.L; ii++)
       data.inputPM.muj[ii] = ii;

     data.inputPM.DUL = new double[data.inputPM.L];
     data.inputPM.DLL = new double[data.inputPM.L];
     data.inputPM.Vol = new double[data.inputPM.L];
     data.inputPM.Pparc = new double[data.inputPM.L];
     data.inputPM.MM = new double[data.inputPM.L];
     data.inputPM.G = new double[data.inputPM.L];
     data.inputPM.G0 = new double[data.inputPM.L];
     data.inputPM.lnGam = new double[data.inputPM.L];
     data.inputPM.lnGmo = new double[data.inputPM.L];
     data.inputPM.X = new double[data.inputPM.L];
     data.inputPM.Y = new double[data.inputPM.L];
     data.inputPM.XY = new double[data.inputPM.L];
     data.inputPM.MU = new double[data.inputPM.L];
     data.inputPM.EMU = new double[data.inputPM.L];
     data.inputPM.NMU = new double[data.inputPM.L];
     data.inputPM.W = new double[data.inputPM.L];
     data.inputPM.F = new double[data.inputPM.L];
     data.inputPM.F0 = new double[data.inputPM.L];
     data.inputPM.RLC = new char[data.inputPM.L];
     data.inputPM.RSC = new char[data.inputPM.L];
     data.inputPM.DCC = new char[data.inputPM.L];
     data.inputPM.DCCW = new char[data.inputPM.L];
     data.inputPM.lnGmM = new double[data.inputPM.L];
     data.inputPM.fDQF = new double[data.inputPM.L]; //24
     for( ii=0; ii<data.inputPM.L; ii++ )
     {
         data.inputPM.DUL[ii] = 1e6;
         data.inputPM.DLL[ii] = 0.0;
         data.inputPM.Vol[ii] = 0.0;
         data.inputPM.Pparc[ii] = 1.;
         data.inputPM.MM[ii] = 0.0;
         data.inputPM.G[ii] = 0.0;
         data.inputPM.G0[ii] = 0.0;
         data.inputPM.lnGam[ii] = 0.0;
         data.inputPM.lnGmo[ii] = 0.0;
         data.inputPM.X[ii] = 0.0;
         data.inputPM.Y[ii] = 0.0;
         data.inputPM.XY[ii] = 0.0;
         data.inputPM.MU[ii] = 0.0;
         data.inputPM.EMU[ii] = 0.0;
         data.inputPM.NMU[ii] = 0.0;
         data.inputPM.W[ii] = 0.0;
         data.inputPM.F[ii] = 0.0;
         data.inputPM.F0[ii] = 0.0;
         data.inputPM.RLC[ii] = 'B';
         data.inputPM.RSC[ii] = 'M';
         data.inputPM.DCC[ii] = 0;
         data.inputPM.DCCW[ii] = 0;
         data.inputPM.lnGmM[ii] = 0.0;
         data.inputPM.fDQF[ii] = 0.0;
     }

     data.inputPM.A = new double[data.inputPM.N*data.inputPM.L];
     for( ii=0; ii<data.inputPM.N*data.inputPM.L; ii++ )
         data.inputPM.A[ii] = 0.0;

     data.inputPM.Awt = new double[data.inputPM.N];
     data.inputPM.B = new double[data.inputPM.N];
     data.inputPM.U = new double[data.inputPM.N];
     data.inputPM.U_r = new double[data.inputPM.N];
     data.inputPM.C = new double[data.inputPM.N];
     data.inputPM.ICC = new char[data.inputPM.N];  //6
     for( ii=0; ii<data.inputPM.N; ii++ )
     {
         data.inputPM.Awt[ii] = 0.0;
         data.inputPM.B[ii] = 0.0;
         data.inputPM.U[ii] = 0.0;
         data.inputPM.U_r[ii] = 0.0;
         data.inputPM.C[ii] = 0.0;
         data.inputPM.ICC[ii] = 0;
     }

     data.inputPM.XFs = new double[data.inputPM.FI];
     data.inputPM.Falps = new double[data.inputPM.FI];
     data.inputPM.XF = new double[data.inputPM.FI];
     data.inputPM.YF = new double[data.inputPM.FI];
     data.inputPM.Falp = new double[data.inputPM.FI];
     data.inputPM.YOF = new double[data.inputPM.FI];
     data.inputPM.PHC = new char[data.inputPM.FI];
     data.inputPM.FVOL = new double[data.inputPM.FI];
     data.inputPM.FWGT = new double[data.inputPM.FI]; //9
     for( ii=0; ii<data.inputPM.FI; ii++ )
     {
         data.inputPM.XFs[ii] = 0.0;
         data.inputPM.Falps[ii] = 0.0;
         data.inputPM.XF[ii] = 0.0;
         data.inputPM.YF[ii] = 0.0;
         data.inputPM.Falp[ii] = 0.0;
         data.inputPM.YOF[ii] = 0.0;
         data.inputPM.PHC[ii] = 0;
         data.inputPM.FVOL[ii] = 0.0;
         data.inputPM.FWGT[ii] = 0.0;
     }

      data.inputPM.SB = new char[data.inputPM.N][MAXICNAME+MAXSYMB];
      data.inputPM.SB1 = new char[data.inputPM.N][MAXICNAME];
      for( ii=0; ii<data.inputPM.N; ii++)
      {
          fillValue( data.inputPM.SB[ii], '\0', MAXICNAME+MAXSYMB);
          fillValue( data.inputPM.SB1[ii], '\0', MAXICNAME);
      }
      data.inputPM.SF = new char[data.inputPM.FI][MAXPHNAME+MAXSYMB];
      data.inputPM.SFs = new char[data.inputPM.FI][MAXPHNAME+MAXSYMB];
      for( ii=0; ii<data.inputPM.FI; ii++)
      {
          fillValue( data.inputPM.SF[ii], '\0', MAXPHNAME+MAXSYMB);
          fillValue( data.inputPM.SFs[ii], '\0',MAXPHNAME+MAXSYMB);
       }
      data.inputPM.SM = new char[data.inputPM.L][MAXDCNAME];
      for( ii=0; ii<data.inputPM.L; ii++)
          fillValue( data.inputPM.SM[ii], '\0', MAXDCNAME);
      data.inputPM.SM2 = new char[data.inputPM.Ls][MAXDCNAME];
      for( ii=0; ii<data.inputPM.Ls; ii++)
          fillValue( data.inputPM.SM2[ii], '\0', MAXDCNAME);
      data.inputPM.SF2 = new char[data.inputPM.FIs][MAXPHNAME+MAXSYMB];
      for( ii=0; ii<data.inputPM.FIs; ii++)
          fillValue( data.inputPM.SF2[ii], '\0', MAXPHNAME+MAXSYMB);
      data.inputPM.dcMod = new char[data.inputPM.L][6];

     if( data.inputPM.L > 0 )
     {
       data.inputPM.Y_la = new double[data.inputPM.L];
       data.inputPM.Y_w = new double[data.inputPM.L];
       data.inputPM.Fx = new double[data.inputPM.L];
       data.inputPM.Wx = new double[data.inputPM.L];
       data.inputPM.VL = new double[data.inputPM.L];
       data.inputPM.Gamma = new double[data.inputPM.L];
       data.inputPM.lnGmf = new double[data.inputPM.L]; //7
    data.inputPM.GamFs = new double[data.inputPM.L];
       for( ii=0; ii<data.inputPM.L; ii++ )
       {
           data.inputPM.Y_la[ii] = 0.0;
           data.inputPM.Y_w[ii] = 0.0;
           data.inputPM.Fx[ii] = 0.0;
           data.inputPM.Wx[ii] = 0.0;
           data.inputPM.VL[ii] = 0.0;
           data.inputPM.Gamma[ii] = 0.0;
           data.inputPM.lnGmf[ii] = 0.0;
           data.inputPM.GamFs[ii] = 0.0;
       }
       //   data.inputPM.D = new double[data.inputPM.L];
     }
     else
     {
       data.inputPM.Y_la = 0;
       data.inputPM.Y_w = 0;
       data.inputPM.Fx = 0;
       data.inputPM.Wx = 0;
       data.inputPM.VL = 0;
       data.inputPM.Gamma = 0;
       data.inputPM.lnGmf = 0;
    data.inputPM.GamFs = 0;
    //   data.inputPM.D = 0;
     }

       // Part 2  not always required arrays

     if( data.inputPM.FIs > 0 && data.inputPM.Ls > 0 )
     {
       data.inputPM.BF = new double[data.inputPM.FIs*data.inputPM.N];
       for( ii=0; ii<data.inputPM.FIs*data.inputPM.N; ii++ )
           data.inputPM.BF[ii] = 0.0;
       data.inputPM.BFC = new double[data.inputPM.N];
       for( ii=0; ii<data.inputPM.N; ii++ )
           data.inputPM.BFC[ii] = 0.0;

       data.inputPM.XFA = new double[data.inputPM.FIs];
       data.inputPM.YFA = new double[data.inputPM.FIs];
       data.inputPM.PUL = new double[data.inputPM.FIs];
       data.inputPM.PLL = new double[data.inputPM.FIs]; //5
       for( ii=0; ii<data.inputPM.FIs; ii++ )
       {
           data.inputPM.XFA[ii] = 0.0;
           data.inputPM.YFA[ii] = 0.0;
               data.inputPM.PUL[ii] = 1e6;
           data.inputPM.PLL[ii] = 0.0;
       }
       data.inputPM.RFLC = new char[data.inputPM.FIs];
       data.inputPM.RFSC = new char[data.inputPM.FIs];
       for( ii=0; ii<data.inputPM.FIs; ii++)
       {
          data.inputPM.RFLC[ii] = 0;
          data.inputPM.RFSC[ii] = 0;
       }

     }
     else
     {
       data.inputPM.BF = 0;
       data.inputPM.BFC = 0;
       data.inputPM.XFA = 0;
       data.inputPM.YFA = 0;
       data.inputPM.PUL = 0;
       data.inputPM.PLL = 0;
       data.inputPM.RFLC = 0;
       data.inputPM.RFSC = 0;
     }

     if( data.inputPM.LO > 1 )
     {
       data.inputPM.Y_m = new double[data.inputPM.L];
       for( ii=0; ii<data.inputPM.L; ii++ )
           data.inputPM.Y_m[ii] = 0.0;
       data.inputPM.IC_m = new double[data.inputPM.N];
       data.inputPM.IC_lm = new double[data.inputPM.N];
       data.inputPM.IC_wm = new double[data.inputPM.N];
       for( ii=0; ii<data.inputPM.N; ii++ )
       {
           data.inputPM.IC_m[ii] = 0.0;
           data.inputPM.IC_lm[ii] = 0.0;
           data.inputPM.IC_wm[ii] = 0.0;
       }
     }
     else
     {
       data.inputPM.Y_m = 0;
       data.inputPM.IC_m = 0;
       data.inputPM.IC_lm = 0;
       data.inputPM.IC_wm = 0;
     }

     // dispersion and sorption phases
     /*
     if( PAalp != S_OFF )
     {
       data.inputPM.Aalp = new double[data.inputPM.FI];
       for( ii=0; ii<data.inputPM.FI; ii++ )
           data.inputPM.Aalp[ii] = 0.0;
       data.inputPM.Xr0h0 = new double[data.inputPM.FI][2];
       for( ii=0; ii<data.inputPM.FI; ii++ )
          data.inputPM.Xr0h0[ii][0] =  data.inputPM.Xr0h0[ii][1] = 0.0;
     }
     else
     {*/
    data.inputPM.Aalp = 0;
    data.inputPM.Xr0h0 = 0;
     //}


     /*
     if( PSigm != S_OFF )
     {   data.inputPM.Sigw = new double[data.inputPM.FI];
         data.inputPM.Sigg = new double[data.inputPM.FI];
         for( ii=0; ii<data.inputPM.FI; ii++ )
         {
             data.inputPM.Sigw[ii] = 0.0;
             data.inputPM.Sigg[ii] = 0.0;
         }
     }
     else
     { */
    data.inputPM.Sigw = 0;
    data.inputPM.Sigg = 0;
    //}

     if( data.inputPM.E )
     {
        data.inputPM.EZ = new double[data.inputPM.L];
        for( ii=0; ii<data.inputPM.L; ii++ )
            data.inputPM.EZ[ii] = 0.0;
        data.inputPM.Xcond = new double[data.inputPM.FI];
        data.inputPM.Xeps = new double[data.inputPM.FI];
        for( ii=0; ii<data.inputPM.FI; ii++ )
        {
            data.inputPM.Xcond[ii] = 0.0;
            data.inputPM.Xeps[ii] = 0.0;
        }
     }
     else
     {
        data.inputPM.EZ = 0;
        data.inputPM.Xcond = 0;
        data.inputPM.Xeps = 0;
     }

     if( data.inputPM.FIat > 0 /*&& data.inputPM.Lads > 0*/ && data.inputPM.FIs > 0 )
     { // ADSORBTION AND ION IXCHANDG
       data.inputPM.SATX = new long int[data.inputPM.Lads][4];
       data.inputPM.MASDJ = new double[data.inputPM.Lads][DFCN];
       data.inputPM.lnSAC = new double[data.inputPM.Lads][4];
       for( ii=0; ii<data.inputPM.Lads; ii++ )
       {
           data.inputPM.SATX[ii][0] = data.inputPM.SATX[ii][1] = data.inputPM.SATX[ii][2] = data.inputPM.SATX[ii][3] = 0;
           data.inputPM.lnSAC[ii][0] = data.inputPM.lnSAC[ii][1] = data.inputPM.lnSAC[ii][2] = data.inputPM.lnSAC[ii][3] = 0.0;
          for( jj=0; jj<MST; jj++ )
              data.inputPM.MASDJ[ii][jj] = 0.0;
       }

       data.inputPM.SCM  = new char[data.inputPM.FIs][MST];
       data.inputPM.Nfsp = new double[data.inputPM.FIs][MST];
       data.inputPM.MASDT = new double[data.inputPM.FIs][MST];
       data.inputPM.XcapA = new double[data.inputPM.FIs][MST];
       data.inputPM.XcapB = new double[data.inputPM.FIs][MST];
       data.inputPM.XcapD = new double[data.inputPM.FIs][MST];
       data.inputPM.XcapF = new double[data.inputPM.FIs][MST];
       data.inputPM.XdlA = new double[data.inputPM.FIs][MST];
       data.inputPM.XdlB = new double[data.inputPM.FIs][MST];
       data.inputPM.XdlD = new double[data.inputPM.FIs][MST];
       data.inputPM.XpsiA = new double[data.inputPM.FIs][MST];
       data.inputPM.XpsiB = new double[data.inputPM.FIs][MST];
       data.inputPM.XpsiD = new double[data.inputPM.FIs][MST];
       data.inputPM.XlamA = new double[data.inputPM.FIs][MST];
       data.inputPM.Xetaf = new double[data.inputPM.FIs][MST];
       data.inputPM.XetaA = new double[data.inputPM.FIs][MST];
       data.inputPM.XetaB = new double[data.inputPM.FIs][MST];
       data.inputPM.XetaD = new double[data.inputPM.FIs][MST];
       data.inputPM.XFTS = new double[data.inputPM.FIs][MST];  //19
       for( ii=0; ii<data.inputPM.FIs; ii++ )
          for( jj=0; jj<MST; jj++ )
          {
              data.inputPM.SCM[ii][jj]  = 0;
              data.inputPM.Nfsp[ii][jj] = 0.0;
              data.inputPM.MASDT[ii][jj] = 0.0;
              data.inputPM.XcapA[ii][jj] = 0.0;
              data.inputPM.XcapB[ii][jj] = 0.0;
              data.inputPM.XcapD[ii][jj] = 0.0;
              data.inputPM.XcapF[ii][jj] = 0.0;
              data.inputPM.XdlA[ii][jj] = 0.0;
              data.inputPM.XdlB[ii][jj] = 0.0;
              data.inputPM.XdlD[ii][jj] = 0.0;
              data.inputPM.XpsiA[ii][jj] = 0.0;
              data.inputPM.XpsiB[ii][jj] = 0.0;
              data.inputPM.XpsiD[ii][jj] = 0.0;
              data.inputPM.XlamA[ii][jj] = 0.0;
              data.inputPM.Xetaf[ii][jj] = 0.0;
              data.inputPM.XetaA[ii][jj] = 0.0;
              data.inputPM.XetaB[ii][jj] = 0.0;
              data.inputPM.XetaD[ii][jj] = 0.0;
              data.inputPM.XFTS[ii][jj] = 0.0;
          }

      data.inputPM.SATT = new char[data.inputPM.Lads];
      data.inputPM.SM3 = new char[data.inputPM.Lads][MAXDCNAME];
      data.inputPM.DCC3 = new char[data.inputPM.Lads];
      for( ii=0; ii<data.inputPM.Lads; ii++)
      {
          fillValue( data.inputPM.SM3[ii], '\0', MAXDCNAME);
          data.inputPM.SATT[ii] = 0;
         data.inputPM.DCC3[ii] = 0;
      }

      data.inputPM.D = new double[MST][MST];
      for( ii=0; ii<MST; ii++ )
          for( jj=0; jj<MST; jj++ )
              data.inputPM.D[ii][jj] = 0.0;

     }
    else
     { // ADSORPTION AND ION EXCHANGE
       data.inputPM.SCM  = 0;
        data.inputPM.Nfsp = 0;
        data.inputPM.MASDT = 0;
        data.inputPM.XcapA = 0;
        data.inputPM.XcapB = 0;
        data.inputPM.XcapD = 0;
        data.inputPM.XcapF = 0;
        data.inputPM.XdlA = 0;
        data.inputPM.XdlB = 0;
        data.inputPM.XdlD = 0;
        data.inputPM.XpsiA = 0;
        data.inputPM.XpsiB = 0;
        data.inputPM.XpsiD = 0;
        data.inputPM.XlamA = 0;
        data.inputPM.Xetaf = 0;
        data.inputPM.XetaA = 0;
        data.inputPM.XetaB = 0;
        data.inputPM.XetaD = 0;
        data.inputPM.MASDJ = 0;
        data.inputPM.XFTS = 0;
        data.inputPM.lnSAC = 0;
        data.inputPM.SATT = 0;
        data.inputPM.SM3 = 0;
        data.inputPM.DCC3 = 0;
        data.inputPM.D = 0;
     }

     if( data.inputPM.PG > 0 )
     {
      data.inputPM.Fug = new double[data.inputPM.PG];
      data.inputPM.Fug_l = new double[data.inputPM.PG];
      data.inputPM.Ppg_l = new double[data.inputPM.PG];
      for( ii=0; ii<data.inputPM.PG; ii++ )
      {
          data.inputPM.Fug[ii] = 0.;
          data.inputPM.Fug_l[ii] = 0.;
          data.inputPM.Ppg_l[ii] = 0.;
      }
     }
    else
     {
      data.inputPM.Fug = 0;
      data.inputPM.Fug_l = 0;
      data.inputPM.Ppg_l = 0;
     }

       // Part 3
     if( data.inputPM.Ls > 1 && data.inputPM.FIs > 0 )
     {
        data.inputPM.Wb = new double[data.inputPM.Ls];
        data.inputPM.Wabs = new double[data.inputPM.Ls];
        data.inputPM.Rion = new double[data.inputPM.Ls];
        for( ii=0; ii<data.inputPM.Ls; ii++ )
        {
            data.inputPM.Wb[ii] = 0.;
            data.inputPM.Wabs[ii] = 0.;
            data.inputPM.Rion[ii] = 0.;
        }
        data.inputPM.Qp = new double[data.inputPM.FIs*QPSIZE];
        data.inputPM.Qd = new double[data.inputPM.FIs*QDSIZE];
        for( ii=0; ii<data.inputPM.FIs*QPSIZE; ii++ )
                data.inputPM.Qp[ii] = 0.;
        for( ii=0; ii<data.inputPM.FIs*QDSIZE; ii++ )
                data.inputPM.Qd[ii] = 0.;
     }
     else
     {
        data.inputPM.Wb = 0;
        data.inputPM.Wabs = 0;
        data.inputPM.Rion = 0;
        data.inputPM.Qp = 0;
        data.inputPM.Qd = 0;

     }

     // added SD 03/02/2009
     data.inputPM.XU = new double[data.inputPM.L];
     for( ii=0; ii<data.inputPM.L; ii++ )
          data.inputPM.XU[ii] = 0.;
     data.inputPM.Uc = new double[data.inputPM.N][2];
     data.inputPM.Uefd = new double[data.inputPM.N];
      for( ii=0; ii<data.inputPM.N; ii++ )
      {
          data.inputPM.Uc[ii][0] = 0.;
          data.inputPM.Uc[ii][1] = 0.;
          data.inputPM.Uefd[ii] = 0.;
      }

      data.inputPM.Cp0   = new double[data.inputPM.L];
      data.inputPM.H0    = new double[data.inputPM.L];
      data.inputPM.U0    = new double[data.inputPM.L];
      data.inputPM.S0    = new double[data.inputPM.L];
      data.inputPM.A0    = new double[data.inputPM.L];
      for( ii=0; ii<data.inputPM.L; ii++ )
      {
          data.inputPM.Cp0[ii]   = 0.;
          data.inputPM.H0[ii]    = 0.;
          data.inputPM.U0[ii]    = 0.;
          data.inputPM.S0[ii]    = 0.;
          data.inputPM.A0[ii]    = 0.;

      }
      data.inputPM.VPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.GPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.HPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.SPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.CPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.APh   = new double[data.inputPM.FIs][MIXPHPROPS];
      data.inputPM.UPh   = new double[data.inputPM.FIs][MIXPHPROPS];
      for( ii=0; ii<data.inputPM.FIs; ii++ )
        for( jj=0; jj<MIXPHPROPS; jj++ )
      {
          data.inputPM.VPh[ii][jj]  = 0.;
          data.inputPM.GPh[ii][jj]  = 0.;
          data.inputPM.HPh[ii][jj]  = 0.;
          data.inputPM.SPh[ii][jj]  = 0.;
          data.inputPM.CPh[ii][jj]  = 0.;
          data.inputPM.APh[ii][jj]  = 0.;
          data.inputPM.UPh[ii][jj]  = 0.;
      }

      // NEW phase definition

      if( data.inputPM.FIs > 0 && data.inputPM.Ls > 0 )
      {
        data.inputPM.IPx = 0;
        data.inputPM.PMc = 0;
        data.inputPM.DMc = 0;
        data.inputPM.MoiSN = 0;
        data.inputPM.SitFr = 0;
        data.inputPM.sMod = new char[data.inputPM.FIs][8];
        for( ii=0; ii<data.inputPM.FIs; ii++)
        {
           fillValue( data.inputPM.sMod[ii], '\0', 8);
        }
        data.inputPM.LsMod = new long int[data.inputPM.FIs*3];
        data.inputPM.LsMdc = new long int[data.inputPM.FIs*3];
        data.inputPM.LsMdc2 = new long int[data.inputPM.FIs*3];
        for( ii=0; ii<data.inputPM.FIs*3; ii++ )
        {
            data.inputPM.LsMod[ii] =0;
            data.inputPM.LsMdc[ii] = 0;
            data.inputPM.LsMdc2[ii] = 0;
        }
        data.inputPM.PhLin = 0;
        data.inputPM.lPhc  = 0;
        data.inputPM.LsPhl = new long int[data.inputPM.FI*2];
        for( ii=0; ii<data.inputPM.FI*2; ii++ )
            data.inputPM.LsPhl[ii] =0;
     // TSolMod stuff
        data.inputPM.lPhc   = 0;
        data.inputPM.DQFc   = 0;
    //    data.inputPM.rcpc   = 0;
        data.inputPM.lnDQFt   = new double[data.inputPM.Ls];
        data.inputPM.lnRcpt   = new double[data.inputPM.Ls];
        data.inputPM.lnExet   = new double[data.inputPM.Ls];
        data.inputPM.lnCnft   = new double[data.inputPM.Ls];
        for( ii=0; ii<data.inputPM.Ls; ii++ )
        {    data.inputPM.lnDQFt[ii] =0.;
            data.inputPM.lnRcpt[ii] =0.;
            data.inputPM.lnExet[ii] =0.;
            data.inputPM.lnCnft[ii] =0.;
         }
    //TSorpMod & TKinMet stuff
        data.inputPM.SorMc   = new double[data.inputPM.FIs*16];
        for( ii=0; ii<data.inputPM.FIs*16; ii++ )
            data.inputPM.SorMc[ii] =0.;
    // TSorpMod stuff
        data.inputPM.LsESmo   = new long int[data.inputPM.FIs*4];
        data.inputPM.LsISmo   = new long int[data.inputPM.FIs*4];
        for( ii=0; ii<data.inputPM.FIs*4; ii++ )
        {    data.inputPM.LsESmo[ii] =0;
            data.inputPM.LsISmo[ii] =0;
        }
        data.inputPM.xSMd   = 0;
        data.inputPM.EImc   = 0;
        data.inputPM.mCDc   = 0;
        data.inputPM.IsoPc   = 0;
        data.inputPM.IsoSc   = 0;
        data.inputPM.lnScalT   = new double[data.inputPM.Ls];
        data.inputPM.lnSACT   = new double[data.inputPM.Ls];
        data.inputPM.lnGammF   = new double[data.inputPM.Ls];
        data.inputPM.CTerms   = new double[data.inputPM.Ls];
        data.inputPM.IsoCt   = 0;
        for( ii=0; ii<data.inputPM.Ls; ii++ )
        {    data.inputPM.lnScalT[ii] =0.;
            data.inputPM.lnSACT[ii] =0.;
            data.inputPM.lnGammF[ii] =0.;
            data.inputPM.CTerms[ii] =0.;
         }
    // TKinMet stuff
        data.inputPM.LsKin   = new long int[data.inputPM.FI*6];
        for( ii=0; ii<data.inputPM.FI*6; ii++ )
            data.inputPM.LsKin[ii] =0;
        data.inputPM.LsUpt   = new long int[data.inputPM.FIs*2];
        for( ii=0; ii<data.inputPM.FIs*2; ii++ )
            data.inputPM.LsUpt[ii] =0;
        data.inputPM.xSKrC   = 0;
        data.inputPM.ocPRkC   = 0;
        data.inputPM.feSArC   = 0;
        data.inputPM.rpConC   = 0;
        data.inputPM.apConC   = 0;
        data.inputPM.AscpC   = 0;
        data.inputPM.UMpcC   = 0;
        data.inputPM.kMod   = new char[data.inputPM.FI][6];
        data.inputPM.PfFact  = new double[data.inputPM.FI];
        data.inputPM.PrT   = new double[data.inputPM.FI];
        data.inputPM.PkT   = new double[data.inputPM.FI];
        data.inputPM.PvT   = new double[data.inputPM.FI];
        for( ii=0; ii<data.inputPM.FI; ii++)
        {
           fillValue( data.inputPM.kMod[ii], 'N', 6);
           data.inputPM.PfFact[ii] =0.;
           data.inputPM.PrT[ii] =0.;
           data.inputPM.PkT[ii] =0.;
           data.inputPM.PvT[ii] =0.;
        }
        data.inputPM.emRd   = new double[data.inputPM.Ls];
        data.inputPM.emDf   = new double[data.inputPM.Ls];
        for( ii=0; ii<data.inputPM.Ls; ii++)
        {
           data.inputPM.emRd[ii] =0.;
           data.inputPM.emDf[ii] =0.;
        }
        data.inputPM.xICuC = 0;

      }
      else
      {
        data.inputPM.LsMod = 0;
        data.inputPM.LsMdc = 0;
        data.inputPM.PMc = 0;
        data.inputPM.DMc = 0;
        data.inputPM.MoiSN = 0;
        data.inputPM.SitFr = 0;
        data.inputPM.sMod = 0;

        data.inputPM.LsMdc2  = 0;
        data.inputPM.LsPhl   = 0;
        data.inputPM.PhLin   = 0;
    // TSolMod stuff
        data.inputPM.lPhc   = 0;
        data.inputPM.DQFc   = 0;
    //    data.inputPM.rcpc   = 0;
        data.inputPM.lnDQFt   = 0;
        data.inputPM.lnRcpt   = 0;
        data.inputPM.lnExet   = 0;
        data.inputPM.lnCnft   = 0;
    //TSorpMod & TKinMet stuff
        data.inputPM.SorMc   = 0;
    // TSorpMod stuff
        data.inputPM.LsESmo   = 0;
        data.inputPM.LsISmo   = 0;
        data.inputPM.xSMd   = 0;
        data.inputPM.EImc   = 0;
        data.inputPM.mCDc   = 0;
        data.inputPM.IsoPc   = 0;
        data.inputPM.IsoSc   = 0;
        data.inputPM.lnScalT   = 0;
        data.inputPM.lnSACT   = 0;
        data.inputPM.lnGammF   = 0;
        data.inputPM.CTerms   = 0;
        data.inputPM.IsoCt   = 0;
    // TKinMet stuff
        data.inputPM.LsKin   = 0;
        data.inputPM.LsUpt   = 0;
        data.inputPM.xSKrC   = 0;
        data.inputPM.ocPRkC   = 0;
        data.inputPM.feSArC   = 0;
        data.inputPM.rpConC   = 0;
        data.inputPM.apConC   = 0;
        data.inputPM.AscpC   = 0;
        data.inputPM.UMpcC   = 0;
        data.inputPM.kMod   = 0;
        // new
        data.inputPM.PfFact  = 0;
        data.inputPM.PrT   = 0;
        data.inputPM.PkT   = 0;
        data.inputPM.PvT   = 0;
        data.inputPM.emRd   = 0;
        data.inputPM.emDf   = 0;
        data.inputPM.xICuC = 0;
      }


     //Alloc_TSolMod( data.inputPM.FIs );
     //Alloc_TSorpMod( data.inputPM.FIs );
     //Alloc_TKinMet( data.inputPM.FI );
    }

    // free dynamic memory
    void Reader::freeDHC(SYSDATA& data)
    {
    if( data.inputDHC.nDCinPH )
    { delete[] data.inputDHC.nDCinPH;
    data.inputDHC.nDCinPH = 0;
    }
    if( data.inputDHC.xic )
    { delete[] data.inputDHC.xic;
    data.inputDHC.xic = 0;
    }
    if( data.inputDHC.xdc )
    { delete[] data.inputDHC.xdc;
    data.inputDHC.xdc = 0;
    }
    if( data.inputDHC.xph )
    { delete[] data.inputDHC.xph;
    data.inputDHC.xph = 0;
    }
    if( data.inputDHC.A )
    { delete[] data.inputDHC.A;
    data.inputDHC.A = 0;
    }
    if( data.inputDHC.ICmm )
    { delete[] data.inputDHC.ICmm;
    data.inputDHC.ICmm = 0;
    }
    if( data.inputDHC.DCmm )
    { delete[] data.inputDHC.DCmm;
    data.inputDHC.DCmm = 0;
    }

    if( data.inputDHC.TKval )
    { delete[] data.inputDHC.TKval;
    data.inputDHC.TKval = 0;
    }
    if( data.inputDHC.Psat )
    { delete[] data.inputDHC.Psat;
    data.inputDHC.Psat = 0;
    }
    if( data.inputDHC.Pval )
    { delete[] data.inputDHC.Pval;
    data.inputDHC.Pval = 0;
    }

    if( data.inputDHC.denW )
    { delete[] data.inputDHC.denW;
    data.inputDHC.denW = 0;
    }
    if( data.inputDHC.denWg )
    { delete[] data.inputDHC.denWg;
    data.inputDHC.denWg = 0;
    }
    if( data.inputDHC.epsW )
    { delete[] data.inputDHC.epsW;
    data.inputDHC.epsW = 0;
    }
    if( data.inputDHC.epsWg )
    { delete[] data.inputDHC.epsWg;
    data.inputDHC.epsWg = 0;
    }
    if( data.inputDHC.G0 )
    { delete[] data.inputDHC.G0;
    data.inputDHC.G0 = 0;
    }
    if( data.inputDHC.V0 )
    { delete[] data.inputDHC.V0;
    data.inputDHC.V0 = 0;
    }
    if( data.inputDHC.H0 )
    { delete[] data.inputDHC.H0;
    data.inputDHC.H0 = 0;
    }
    if( data.inputDHC.Cp0 )
    { delete[] data.inputDHC.Cp0;
    data.inputDHC.Cp0 = 0;
    }
    if( data.inputDHC.S0 )
    { delete[] data.inputDHC.S0;
     data.inputDHC.S0 = 0;
    }
    if( data.inputDHC.A0 )
    { delete[] data.inputDHC.A0;
     data.inputDHC.A0 = 0;
    }
    if( data.inputDHC.U0 )
    { delete[] data.inputDHC.U0;
     data.inputDHC.U0 = 0;
    }
    if( data.inputDHC.DD )
    { delete[] data.inputDHC.DD;
     data.inputDHC.DD = 0;
    }

    if( data.inputDHC.ICNL )
    { delete[] data.inputDHC.ICNL;
    data.inputDHC.ICNL = 0;
    }
    if( data.inputDHC.DCNL )
    { delete[] data.inputDHC.DCNL;
    data.inputDHC.DCNL = 0;
    }
    if( data.inputDHC.PHNL )
    { delete[] data.inputDHC.PHNL;
    data.inputDHC.PHNL = 0;
    }

    if( data.inputDHC.ccIC )
    { delete[] data.inputDHC.ccIC;
    data.inputDHC.ccIC = 0;
    }
    if( data.inputDHC.ccDC )
    { delete[] data.inputDHC.ccDC;
    data.inputDHC.ccDC = 0;
    }
    if( data.inputDHC.ccPH )
    { delete[] data.inputDHC.ccPH;
    data.inputDHC.ccPH = 0;
    }
    // delete[] CSD;
}

// free dynamic memory
void Reader::freeDBR(SYSDATA& data) {
 if( data.inputBR.bIC )
 { delete[] data.inputBR.bIC;
   data.inputBR.bIC = 0;
 }
 if( data.inputBR.rMB )
 { delete[] data.inputBR.rMB;
   data.inputBR.rMB = 0;
 }
 if( data.inputBR.uIC )
 { delete[] data.inputBR.uIC;
   data.inputBR.uIC = 0;
 }

 if( data.inputBR.xDC )
  { delete[] data.inputBR.xDC;
    data.inputBR.xDC = 0;
  }
 if( data.inputBR.gam )
  { delete[] data.inputBR.gam;
    data.inputBR.gam = 0;
  }
 if( data.inputBR.dul )
   { delete[] data.inputBR.dul;
     data.inputBR.dul = 0;
   }
 if( data.inputBR.dll )
   { delete[] data.inputBR.dll;
     data.inputBR.dll = 0;
   }

 if( data.inputBR.aPH )
 { delete[] data.inputBR.aPH;
   data.inputBR.aPH = 0;
 }
 if( data.inputBR.xPH )
  { delete[] data.inputBR.xPH;
    data.inputBR.xPH = 0;
  }
 if( data.inputBR.omPH )
  { delete[] data.inputBR.omPH;
    data.inputBR.omPH = 0;
  }
 if( data.inputBR.vPS )
  { delete[] data.inputBR.vPS;
    data.inputBR.vPS = 0;
  }
 if( data.inputBR.mPS )
  { delete[] data.inputBR.mPS;
    data.inputBR.mPS = 0;
  }
 if( data.inputBR.bPS )
  { delete[] data.inputBR.bPS;
    data.inputBR.bPS = 0;
  }
 if( data.inputBR.xPA )
  { delete[] data.inputBR.xPA;
    data.inputBR.xPA = 0;
  }

 if( data.inputBR.bSP )
 { delete[] data.inputBR.bSP;
   data.inputBR.bSP = 0;
 }
 if( data.inputBR.amru )
  { delete[] data.inputBR.amru;
    data.inputBR.amru = 0;
  }
 if( data.inputBR.amrl )
  { delete[] data.inputBR.amrl;
    data.inputBR.amrl = 0;
  }
}

void Reader::generateFiles(SYSDATA & data) {
    generateDHCfile(data, outFileDhc);
    generateIPMfile(data, outFileIpm);
    generateDBRfile(data, outFileDbr);
}

ValueUnit Reader::getValue(std::string& str) {

    ValueUnit res;

    std::string::iterator c;
    std::string val = "";
    std::string unit = "";


    for (c = str.begin(); c < str.end(); ++c) {

        if (isdigit(*c) || (*c) == '.' || (*c) == ',') {
            val += *c;
        } else if (isalpha(*c)) {
            unit += (*c);
        } else {
        }
    }


    res.scale = 1.0;
    res.shift = 0.0;

    if (unit.size() > 0) {
        res.unut  = unit;
        res.value = ToDouble(val);

        if (res.unut == "fs")
            res.scale = 1.0;
        if (res.unut == "ps")
            res.scale = 1.0e3;
        if (res.unut == "ns")
            res.scale = 1.0e6;

        if (res.unut == "A")
            res.scale = 1.0;
        if (res.unut == "nm")
            res.scale = 1.0e1;
        if (res.unut == "pm")
            res.scale = 1.0e-2;

        if (res.unut == "K")
            res.shift = 0.0;
        if (res.unut == "C")
            res.shift = +273.15;

        if (res.unut == "MPa")
            res.scale = 1.0;
        if (res.unut == "kPa")
            res.scale = 1.0e-3;
        if (res.unut == "GPa")
            res.scale = 1.0e3;
        if (res.unut == "Bar" || res.unut == "bar")
            res.scale = 1.0e-1;
    } else {
        res.scale = 1.0;
        res.shift = 0.0;
        res.value = ToDouble(val);
    }

    res.value = (res.value + res.shift) * res.scale;

    return res;

}

double Reader::ToDouble(const std::string& s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw std::runtime_error("convertToDouble(\"" + s + "\")");
    return x;
}
