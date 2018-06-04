#include "formulaparser.h"

void FormulaParser::Clear() {
    atoms.clear();
    num.clear();
}

void FormulaParser::readFormula(std::string & formula) {

    std::string::iterator c;
    std::string val;
    std::string at;

    at = "";
    val = "";

    for (c = formula.begin(); c < formula.end(); ++c) {

        if (isalpha(*c)) {
            if (isupper(*c) && at == "") {
                at += *c;
                val = "1";
            }

            if (isupper(*c) && at != "") {
                atoms.push_back(at);
                num.push_back(ToInt(val));
            }

            if (islower(*c)) {
                at += *c;
                atoms.push_back(at);
                num.push_back(ToInt(val));
                at = "";
                val = "";
            }
        }

        if (isdigit(*c) || (*c) == '.' || (*c) == ',') {
            val = *c;
            atoms.push_back(at);
            num.push_back(ToInt(val));
            at = "";
            val = "";
        }
    }

}


int FormulaParser::ToInt(const std::string& s){
    int r = 0;
    std::string::const_iterator p = s.begin();
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (p != s.end() && *p >= '0' && *p <= '9') {
        r = (r*10) + (*p - '0');
        ++p;
    }
    if (neg) {
        r = -r;
    }
    return r;
}
