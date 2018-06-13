#ifndef FORMULAPARSER_H
#define FORMULAPARSER_H

#include <string>
#include <vector>

class FormulaParser {
public:
    std::vector<std::string> atoms;
    std::vector<int> num;

    void Clear();

    void readFormula(std::string&);
    int ToInt(const std::string&);
};

#endif // FORMULAPARSER_H
