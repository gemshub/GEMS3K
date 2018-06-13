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

    void openInputFile(std::string&, SYSDATA&);
    void readPairs(std::vector<std::string>::iterator begin,
                   std::vector<std::string>::iterator end,
                   SYSDATA & dat);
    void generateFiles(SYSDATA&);
    void clearFiles();

    bool fileExists(const std::string& filename);
    std::vector<std::string> split(const std::string& str);
    ValueUnit getValue(std::string& str);
    double ToDouble(const std::string&);

    Reader();
};

#endif // READER_H
