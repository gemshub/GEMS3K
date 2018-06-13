#ifndef FORMULAPARSER_H
#define FORMULAPARSER_H

#include <string>
#include <vector>
#include <map>
#include "datatypes.h"

class FormulaParser {
public:
    std::vector<std::string> atoms;
    std::vector<int> num;

    void Clear();

    std::map<std::string, int> readFormula(std::string&);
    std::map<std::string, int> readGroup(std::string& formula, int group);
    std::vector<groupStr> getGroups(std::string& formula);
    int ToInt(const std::string&);
    int ToInt(const char& s);
};

#endif // FORMULAPARSER_H
