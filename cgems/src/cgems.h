#ifndef CGEMS_H
#define CGEMS_H

#include "reader.h"
#include "node.h"

class cGems {
public:
    Reader reader;
    SYSDATA data;
    std::string path;
    TNode* node;
    DATABR* dBR;

    cGems() {}
    cGems(std::string file);
    void showSystem();
    void run();
};

#endif // CGEMS_H
