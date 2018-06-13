#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include "data.h"
#include "node.h"

struct ValueUnit {
    std::string unut;
    double value, scale, shift;
};

struct groupStr
{
    std::string group;
    int sc;
    groupStr() {
        sc = 1;
        group = "";
    }

public:
    void Clear(){
        group = "";
        sc = 1;
    }
};

struct CompoundStr {
public:
    std::string name;
    std::string formula;
    std::string code;
    std::map<std::string, int> atoms;
    double initX;
    double V0;
    double G0;
    double S0;
    double H0;
    double Cp0;
    double A0;
    double U0;

    CompoundStr() {
        name = "";
        code = "";
        V0 = 0.0;
        G0 = 0.0;
        S0 = 0.0;
        H0 = 0.0;
        Cp0 = 0.0;
        A0 = 0.0;
        U0 = 0.0;
        initX = 0.0;
        atoms.clear();
    }
};

struct PhaseStr {
public:
    std::string name;
    std::string code;
    std::vector<std::string> compounds;
    PhaseStr(){
        name = "";
        code = "";
        compounds.clear();
    }
    void clear() {
        name = "";
        code = "";
        compounds.clear();
    }
};

struct SYSDATA  {
    MULTI inputPM;
    DATACH inputDHC;
    DATABR inputBR;
    SPP_SETTING inputPA;

    std::string project;
    double P, T;
    DataInfo dataInfo;
    std::vector<CompoundStr> Compounds;
    std::vector<double> Concentration, Ammount;
    std::map<std::string, int> totElements;
    std::vector<PhaseStr> Phases;

    SYSDATA() {
        project = "";
        P = 0.0;
        T = 298.15;
    }

    void Clear() {
        Ammount.clear();
        Concentration.clear();
        Compounds.clear();
        Phases.clear();
        project = "";
        P = 0.0;
        T = 298.15;
    }

    int getSolutionsCount() {
        int count = 0;
        for (auto& el : Phases)
            count += el.compounds.size() > 1 ? 1 : 0;
        return count;
    }

    int getCoumpoundsInSolutions() {
        int count = 0;
        for (auto& el : Phases)
            if (el.compounds.size() > 1)
                count += el.compounds.size();
        return count;
    }

    std::string getICComposition() {
        std::stringstream res;
        double conc;

        for (auto& el : totElements) {
            conc = 0.0;
            for (auto& c : Compounds) {
                conc += c.initX * c.atoms[el.first];
            }
            res << conc << "  " ;
        }
        return res.str();
    }

    double* getICCompositionArray() {
        std::stringstream res;
        double* buff = new double[200];
        double conc;
        int i = 0;
        for (auto& el : totElements) {
            conc = 0.0;
            for (auto& c : Compounds) {
                conc += c.initX * c.atoms[el.first];
            }
            buff[i] = conc;
        }
        return buff;
    }
    std::string getDCComposition() {
        std::stringstream res;

        for (auto& c : Compounds) {
            res << c.initX << " " ;
        }
        return res.str();
    }

    std::string getICindexes() {
        std::stringstream res;
        for (unsigned int i = 0; i < totElements.size(); i++)
            res << i << " ";
        return res.str();
    }

    long int* getICindexesArray() {
        long* buff = new long[200];
        for (unsigned int i = 0; i < totElements.size(); i++)
            buff[i] = i;
        return buff;
    }


    std::string getDCindexes() {
        std::stringstream res;
        for (unsigned int i = 0; i < Compounds.size(); i++)
            res << i << " ";
        return res.str();
    }

    long int* getDCindexesArray() {
        long* buff = new long[200];
        for (unsigned int i = 0; i < Compounds.size(); i++)
            buff[i] = i;
        return buff;
    }

    std::string getPhaseindexes() {
        std::stringstream res;
        for (unsigned int i = 0; i < Phases.size(); i++)
            res << i << " ";
        return res.str();
    }

    long int* getPhaseindexesArray() {
        long* buff = new long[200];
        for (unsigned int i = 0; i < Phases.size(); i++)
            buff[i] = i;
        return buff;
    }

    std::string getICNames() {
        std::stringstream res;
        for (auto& el : totElements)
            res << "'" << el.first << "' ";
        return res.str();
    }

    void getICNamesArray(char (*ICNL)[6]) {
        char buff[200][6];
        int i=0;
        for (auto& el : totElements) {
            for (unsigned int j = 0; j < el.first.size(); j++) {
                buff[i][j] = el.first.c_str()[j];
            }
            i++;
        }
        ICNL = buff;
    }


    std::string getICCodes() {
        std::stringstream res;
        for (auto& el : totElements)
            res << "'" << dataInfo.ClassCodes[el.first] << "' ";
        return res.str();
    }

    char* getICCodesArray() {
        char* buff = new char[200];
        int i=0;
        for (auto& el : totElements){
            buff[i] = *(dataInfo.ClassCodes[el.first].c_str());
            i++;
        }
        return buff;
    }

    std::string getICAtomicMasses() {
        std::stringstream res;
        for (auto& el : totElements)
            res << dataInfo.AtomicMasses[el.first] << " ";
        return res.str();
    }

    double* getICAtomicMassesArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& el : totElements) {
            buff[i] = dataInfo.AtomicMasses[el.first];
            i++;
        }
        return buff;
    }

    std::string getDCNames() {
        std::stringstream res;
        for (auto& p : Phases) {
            for (auto& c : p.compounds)
                res << "'" << c << "' ";
        }
        return res.str();
    }

    void getDCNamesArray(char (*ICNL)[16]) {
        char buff[200][16];
        int i=0;
        for (auto& p : Phases) {
            for (auto& c : p.compounds){
                for (unsigned int j = 0; j < c.size(); j++) {
                    buff[i][j] = c.c_str()[j];
                }
                i++;
            }
        }
        ICNL = buff;
    }

    std::string getDCCodes() {
        std::stringstream res;

        for (auto& p : Phases) {
            for (auto& c : p.compounds)
                res << "'" << p.code << "' ";
        }
        return res.str();
    }

    char* getDCCodesArray() {
        char* buff = new char[200];
        int i = 0;
        for (unsigned int p = 0; p < Phases.size(); p++) {
            for (auto& c : Phases[p].compounds) {
                buff[i] =  *(Phases[p].code.c_str());
                i++;
            }
        }
        return buff;
    }


    std::string getDCAtomicMasses() {
        std::stringstream res;
        double val;
        for (std::vector<CompoundStr>::const_iterator m = Compounds.begin();
             m != Compounds.end(); m++) {
            val = 0;
            for (auto& el : (*m).atoms) {
                val += el.second * dataInfo.AtomicMasses[el.first];
            }
            res << val << " ";
        }
        return res.str();
    }

    double* getDCAtomicMassesArray() {
        double* buff = new double[200];
        double val;
        int i = 0;
        for (std::vector<CompoundStr>::const_iterator m = Compounds.begin();
            m != Compounds.end(); m++) {
            val = 0;
            for (auto& el : (*m).atoms) {
                val += el.second * dataInfo.AtomicMasses[el.first];
            }
            buff[i] = val;
            i++;
        }
        return buff;
    }


    std::string getPhasesNames() {
        std::stringstream res;
        for (auto& el : Phases) {
            res << "'" << el.name << "' ";
        }
        return res.str();
    }

    void getPhasesNamesArray(char (*PHNL)[16]) {
        char buff[200][16];
        int i=0;
        for (auto& el : Phases) {
            for (unsigned int j = 0; j < el.name.size(); j++) {
                buff[i][j] = el.name.c_str()[j];
            }
            i++;
        }
        PHNL = buff;
    }

    std::string getPhasesCodes() {
        std::stringstream res;
        for (auto& el : Phases)
            res << "'" << el.code << "' ";
        return res.str();
    }

    char* getPhasesCodesArray() {
        char* buff = new char[200];
        int i=0;
        for (auto& el : Phases) {
            buff[i] = *(el.code.c_str());
            i++;
        }
        return buff;
    }

    std::string getCDinPhases() {
        std::stringstream res;
        for (auto& el : Phases)
            res << el.compounds.size() << " ";
        return res.str();
    }

    long* getCDinPhasesArray() {
        long* buff = new long[200];
        for (unsigned int i = 0; i < Phases.size(); i++)
            buff[i] = Phases[i].compounds.size();
        return buff;
    }

    std::string getStohiometryMatrix() {
        std::stringstream res;
        for (auto& c : Compounds) {
            for (auto& e : totElements) {
                res << c.atoms[e.first] << "  " ;
            }
            res << std::endl;
        }
        return res.str();
    }

    double* getStohiometryMatrixArray() {
        double* buff = new double[2000];
        int i = 0;
        for (auto& c : Compounds) {
            for (auto& e : totElements) {
                buff[i] = c.atoms[e.first];
                i++;
            }
        }
        return buff;
    }

    std::string getMolarVolumes() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.V0 << std::endl;
        return res.str();
    }

    double* getMolarVolumesArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.V0;
            i++;
        }
        return buff;
    }

    std::string getMolarGibbsEnery() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.G0 << std::endl;
        return res.str();
    }

    double* getMolarGibbsEneryArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.G0;
            i++;
        }
        return buff;
    }

    std::string getMolarEnthalpy() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.H0 << std::endl;
        return res.str();
    }

    double* getMolarEnthalpyArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.H0;
            i++;
        }
        return buff;
    }

    std::string getMolarEntropy() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.S0 << std::endl;
        return res.str();
    }

    double* getMolarEntropyArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.S0;
            i++;
        }
        return buff;
    }

    std::string getMolarHeatCapacity() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.Cp0 << std::endl;
        return res.str();
    }

    double* getMolarHeatCapacityArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.Cp0;
            i++;
        }
        return buff;
    }

    std::string getMolarHelmholzEnergy() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.A0 << std::endl;
        return res.str();
    }

    double* getMolarHelmholzEnergyArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.A0;
            i++;
        }
        return buff;
    }

    std::string getMolarInternalEnergy() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << c.U0 << std::endl;
        return res.str();
    }

    double* getMolarInternalEnergyArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = c.U0;
            i++;
        }
        return buff;
    }

    std::string getMixModelCodes() {
        std::stringstream res;
        for (auto& el : Phases)
            if (el.compounds.size() > 1)
                res << "'INNINNNN'" << " ";
        return res.str();
    }

    void getMixModelCodesArray(char (*Codes)[8]) {
        char buff[200][8];
        std::string str = "INNINNNN";
        int i=0;
        for (auto& el : Phases) {
            if (el.compounds.size() > 1) {
                for (unsigned int j = 0; j < str.size(); j++){
                    buff[i][j] = str.c_str()[j];
                }
                i++;
            }
        }
        Codes = buff;
    }

    std::string getSolModDimentions() {
        std::stringstream res;
        for (auto& el : Phases)
            if (el.compounds.size() > 1)
                res << "0 0 0" << std::endl;
        return res.str();
    }

    long* getSolModDimentionsArray() {
        long* buff = new long[200];
        int i = 0;
        for (auto& el : Phases)
            if (el.compounds.size() > 1)
                buff[i+0] = 0;
                buff[i+1] = 0;
                buff[i+2] = 0;
                i+=3;
        return buff;
    }

    std::string getInteractionParameters() {
        std::stringstream res;
        //for (auto& el : Phases)
            res << "0.0640000030398369 3.72000002861023 1 0 1 0 0 0";
            //if (el.compounds.size() > 1)
            //    res << "0" << " ";
        return res.str();
    }

    std::string getDMcMoiSNDimentions() {
        std::stringstream res;
        for (auto& el : Phases) {
            //if (el.code != "G") {
                if (el.compounds.size() > 1)
                    res << "0 0 0" << std::endl;
            //}
        }
        return res.str();
    }

    std::string getDCPartialPressures() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "1 ";
        return res.str();
    }

    double* getDCPartialPressuresArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 1.0;
            i++;
        }
        return buff;
    }

    std::string getDCFugacities() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "0 ";
        return res.str();
    }

    double* getDCFugacitiesArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getDClnActivities() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "0 ";
        return res.str();
    }

    double* getDClnActivitiesArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getDCMetastabilityConstraints() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "'B' ";
        return res.str();
    }

    char* getDCMetastabilityConstraintsArray() {
        char *buff = new char[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 'B';
            i++;
        }
        return buff;
    }

    std::string getDCMetastabilityUnits() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "'M' ";
        return res.str();
    }

    char* getDCMetastabilityUnitsArray() {
        char *buff = new char[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 'M';
            i++;
        }
        return buff;
    }

    std::string getDCLowerMetastConstr() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << 0 << " ";
        return res.str();
    }

    double* getDCLowerMetastConstrArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getDCUpperMetastConstr() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << 1.0e6 << " ";
        return res.str();
    }

    double* getDCUpperMetastConstrArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Compounds) {
            buff[i] = 1.0e6;
            i++;
        }
        return buff;
    }

    std::string getPHSpecifSurfaces() {
        std::stringstream res;
        for (auto& c : Phases)
            res << 0 << " ";
        return res.str();
    }

    double* getPHSpecifSurfacesArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Phases) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }


    std::string getPHSurfaceEnergyOfPHWater() {
        std::stringstream res;
        for (auto& c : Phases)
            res << 0 << " ";
        return res.str();
    }

    double* getPHSurfaceEnergyOfPHWaterArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Phases) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getPHSurfaceEnergyOfPHGas() {
        std::stringstream res;
        for (auto& c : Phases)
            res << 0 << " ";
        return res.str();
    }

    double* getPHSurfaceEnergyOfPHGasArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Phases) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getPHSurfaceEnergyParameter() {
        std::stringstream res;
        for (auto& c : Phases)
            res << 0 << " ";
        return res.str();
    }

    double* getPHSurfaceEnergyParameterArray() {
        double* buff = new double[200];
        int i = 0;
        for (auto& c : Phases) {
            buff[i] = 0.0;
            i++;
        }
        return buff;
    }

    std::string getDCPTCorrectionCodes() {
        std::stringstream res;
        for (auto& c : Compounds)
            res << "'HKF' ";
        return res.str();
    }

    void getDCPTCorrectionCodesArray(char (*ICNL)[6]) {
        char buff[200][6];
        int i=0;
        std::string str = "HKF";
        for (auto& c : Compounds){
            for (int j = 0; j < 3; j++) {
                buff[i][j] = str.c_str()[j];
            }
            i++;
        }
        ICNL = buff;
    }

};

#endif // DATATYPES_H
