#ifndef GEMS3K_IMPEX_H
#define GEMS3K_IMPEX_H

#include "sstream"
#include "v_detail.h"


/// Descripton of data to generate MULTI, DATACH and DATABR files structure prepared from GEMS.
class GEMS3KImpexData
{

public:

    /// These are used io formats
    enum FileIOModes {
        f_key_value,
        f_binary,
        f_json
    };

    /// IPM work structure file path&name
    std::string ipmfiles_lst_name;

    /// Number of allocated nodes
    long int nIV = 1;

    /// Write IPM, DCH and DBR files in binary, txt or json mode)
    FileIOModes io_mode = f_json;

    /// Do not write data items that contain only default values
    bool brief_mode = false;

    /// Write files with comments for all data entries ( in text mode )
    bool with_comments = false;

    /// Print internal indices in RMULTS to IPM file for reading into Gems back
    bool add_mui = false;

    /// Prints files for separate coupled FMT-GEM programs that use GEMS3K module
    /// or if putNodT1 == true  as a break point for the running FMT calculation
    bool putNodT1 = false;

};



#endif // GEMS3K_IMPEX_H
