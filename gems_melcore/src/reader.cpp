#import "reader.h"
#import <vector>

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
    std::cout << "Done reading the input file" << std::endl;
}

void Reader::readPairs(std::vector<std::string>::iterator begin,
                       std::vector<std::string>::iterator end,
                       SYSDATA & data) {
    std::string unit;
    std::string value;
    ValueUnit ValUn;

    while (begin != end) {
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

        if (*begin == "Compound" || *begin == "Comp" || *begin=="comp" || *begin == "compound") {
            data.Compounds.clear();
            ++begin;
            while (begin < end){
                data.Compounds.push_back(*begin);
                ++begin;
            }
            continue;
        }

        if (*begin == "Concentration" || *begin == "Conc" || *begin=="conc" || *begin == "concentration") {
            data.Concentration.clear();
            ++begin;
            while (begin < end){
                data.Concentration.push_back(ToDouble(*begin));
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

void Reader::generateFiles(SYSDATA & data) {
    std::string text = "";

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
