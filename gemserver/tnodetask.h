#ifndef TNODETASK_H
#define TNODETASK_H

#include <string>
#include <vector>
#include <memory>

class TNode;

class TNodeTask
{

public:

    /// Constructor
    TNodeTask();
    /// Copy constructor
    TNodeTask(const TNodeTask &other );
    /// Move constructor
    TNodeTask( TNodeTask &&other ) = default;
    /// Destructor
    virtual ~TNodeTask()  {}

    /// Copy assignment
    TNodeTask &operator =( const TNodeTask &other);
    /// Move assignment
    TNodeTask &operator =( TNodeTask &&other) = default;

    /// Initialization of GEMS3K internal data from json/key-value strings
    std::vector<std::string> initGEM( const std::string& dch_json, const std::string& ipm_json, const std::string& dbr_json );

    /// Initialization of GEM IPM3 data structures in coupled programs
    /// that use GEMS3K module. Also reads in the IPM, DCH and one or many DBR text input files.
    std::vector<std::string> initGEM( const char *ipmfiles_lst_name );

    /// Run process of calculate equilibria into the GEMS3K side
    std::vector<std::string> calculateEquilibrium( const std::string& new_dbr );


protected:

    /// Internal TNode structure instance accessible through the "node" pointer
    std::shared_ptr<TNode> current_task;

    /// Last loaded task name ( could be used to test )
    std::string task_name;

    /// copy data
    void copy(const TNodeTask &other);
};

#endif // TNODETASK_H
