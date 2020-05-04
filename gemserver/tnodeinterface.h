#ifndef TNODEINTERFACE_H
#define TNODEINTERFACE_H

#include <string>
#include <vector>

extern const char* one_system_task;
extern const char* only_dbr_task;
extern const char* nodearray_task;


class NodeInterface
{

public:

    /// Destructor
    virtual ~NodeInterface()  {}

    /// Initialization of internal data from json/key-value strings
    /// Parameters:
    ///  @param dch_json -  DATACH - the Data for CHemistry data structure as a json string
    ///  @param ipm_json -  Multi structure as a json string
    ///  @param dbr_json -  DATABR - the data bridge structure as a json string
    ///  @return array of error type, error title and error message in case of error or empty array if all OK
    virtual std::vector<std::string> initData( const std::string& dch_json, const std::string& ipm_json, const std::string& dbr_json ) = 0;

    // Initialization of GEM IPM3 data structures in coupled programs
    // that use GEMS3K module. Also reads in the IPM, DCH and one or many DBR text input files.
    //virtual std::vector<std::string> initData( const char *ipmfiles_lst_name ) = 0;

    /// Run process of calculate equilibria into the GEMS3K/Reaktoro side
    /// Parameters:
    ///  @param dbr_json -  DATABR - the data bridge structure as a json string
    ///  @return array with strings contains:
    ///    - NodeStatus codes with respect to GEMIPM calculations
    ///    - the result DATABR - the data bridge structure as a json string
    ///    - calculation time in seconds elapsed during the last call of GEM_run (obsolete)
    ///    - number of IPM loops performed ( >1 up to 6 because of PSSC() ) (obsolete)
    ///    - Number of completed IA EFD iterations (obsolete)
    ///    - Number of completed GEM IPM iterations (obsolete)
    ///    or if NodeStatus is error
    ///    - error title
    ///    - error message
    virtual std::vector<std::string> calculateEquilibrium( const std::string& new_dbr ) = 0;

};

#endif // TNODEINTERFACE_H
