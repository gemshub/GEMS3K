#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>
#include <vector>

struct ValueUnit {
    std::string unut;
    double value, scale, shift;
};

struct SYSDATA  {
    double P, T;
    std::vector<std::string> Compounds;
    std::vector<double> Concentration, Ammount;

    SYSDATA() {
        P = 0.0;
        T = 298.15;
    }

    void Clear() {
        Ammount.clear();
        Concentration.clear();
        Compounds.clear();
        P = 0.0;
        T = 298.15;
    }
};

#endif // DATATYPES_H
