#include "formulaparser.h"
#include <iostream>

void FormulaParser::Clear() {
    atoms.clear();
    num.clear();
}

std::map<std::string, int> FormulaParser::readGroup(std::string& formula, int group) {
    std::string::iterator c;
    std::string val;
    std::string at;
    std::map<std::string, int> el;

    at = "";
    val = "";

    for (c = formula.begin(); c < formula.end(); ++c) {

        if (isdigit(*c)) {
            val = *c;
        }

        if (islower(*c)) {
            at += *c;
        }

        if (isalpha(*c)) {
            if (isupper(*c) && at == "") {
                at += *c;
                continue;
            }

            if (isupper(*c) && at != "") {
                if (val == "") val = "1";
                el[at] += group*ToInt(val);
                at = *c;
                val = "";
                continue;
            }
        }

        if (c == (formula.end()-1) && at != "") {
            if (val == "") val = "1";
            el[at] += group*ToInt(val);
        }
    }

    return el;
}

std::vector<groupStr> FormulaParser::getGroups(std::string& formula) {
    std::vector<groupStr> groups;
    groups.clear();

    std::string::iterator c;
    std::string val;
    std::string at;
    std::string gro;
    groupStr group;
    int groupSc = 1;

    group.group = "";
    group.Clear();
    for (c = formula.begin(); c < formula.end(); ++c) {

        if (isalnum(*c)) {
            group.group += *c;
        }

        if ((*c) == '(' || (*c) == '[' || (*c) == '{') {
            if (group.group != "") {
                groups.push_back(group);
            }

            group.Clear();
            c++;
            while (*c != ')' && *c != ']' && *c != '}' && c != formula.end()) {
                group.group += *c;
                c++;
            }
            c++;
            group.sc = ToInt(*c);
            groups.push_back(group);
            group.Clear();
        }
    }
    if (group.group != "")
        groups.push_back(group);
    return groups;
}


std::map<std::string, int> FormulaParser::readFormula(std::string& formula) {
    std::map<std::string, int> el;
    std::map<std::string, int> res;

    res.clear();
    std::vector<groupStr> groups = getGroups(formula);

    for (auto& r : groups) {
        el = readGroup(r.group, r.sc);
        for(auto& it : el)
        {
            res[it.first] += it.second;
        }
    }
    return res;
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

int FormulaParser::ToInt(const char& s){
    int r = 0;
    r = (r*10) + (s - '0');
    return r;
}
