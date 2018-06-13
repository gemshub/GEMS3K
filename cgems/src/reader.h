#ifndef READER_H
#define READER_H

#include "datatypes.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include "formulaparser.h"
/*
#include <ctime>
#include <numeric>
#include <stdexcept>
*/

class Reader {
public:
    std::string outFileDbr = "data/dbr.dat";
    std::string outFileDhc = "data/dhc.dat";
    std::string outFileIpm = "data/ipm.dat";
    std::string inFiles    = "data/infile.dat";

    FormulaParser parser;

    void openInputFile(std::string&, SYSDATA&);
    void readPairs(std::vector<std::string>::iterator begin,
                   std::vector<std::string>::iterator end,
                   SYSDATA & dat);
    void generateFiles(SYSDATA&);
    void setupData(SYSDATA&);
    void generateDHCfile(SYSDATA&, std::string&);
    void generateIPMfile(SYSDATA&, std::string&);
    void generateDBRfile(SYSDATA&, std::string&);
    void generateDHC(SYSDATA&);
    void generateIPM(SYSDATA&);
    void generateDBR(SYSDATA&);
    void clearFiles();
    void clean();
    void allocateDHC(SYSDATA&);
    void allocateDBR(SYSDATA&);
    void allocateIPM(SYSDATA&);
    void freeDHC(SYSDATA&);
    void freeDBR(SYSDATA&);
    void freeIPM(SYSDATA&);
    void readFormula(std::string&, SYSDATA&);
    int gridTP(SYSDATA&);
    bool fileExists(const std::string& filename);
    std::vector<std::string> split(const std::string& str);
    ValueUnit getValue(std::string& str);
    double ToDouble(const std::string&);

    Reader();
    //~Reader();
};

#endif // READER_H
