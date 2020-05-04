#ifndef TNODETASK_H
#define TNODETASK_H

#include "tnodeinterface.h"
#include <memory>

class TNode;


class NodeGEMSTask: public NodeInterface
{

public:

    /// Constructor
    NodeGEMSTask();
    /// Copy constructor
    NodeGEMSTask(const NodeGEMSTask &other );
    /// Move constructor
    NodeGEMSTask( NodeGEMSTask &&other ) = default;
    /// Destructor
    ~NodeGEMSTask()  {}

    /// Copy assignment
    NodeGEMSTask &operator =( const NodeGEMSTask &other);
    /// Move assignment
    NodeGEMSTask &operator =( NodeGEMSTask &&other) = default;

    /// Initialization of GEMS3K internal data from json/key-value strings
    std::vector<std::string> initData( const std::string& dch_json, const std::string& ipm_json, const std::string& dbr_json ) override;

    /// Initialization of GEM IPM3 data structures in coupled programs
    /// that use GEMS3K module. Also reads in the IPM, DCH and one or many DBR text input files.
    virtual std::vector<std::string> initData( const char *ipmfiles_lst_name );

    /// Run process of calculate equilibria into the GEMS3K side
    std::vector<std::string> calculateEquilibrium( const std::string& new_dbr ) override;


protected:

    /// Internal TNode structure instance accessible through the "node" pointer
    std::shared_ptr<TNode> current_task;

    /// Last loaded task name ( could be used to test )
    std::string task_name;

    /// copy data
    void copy(const NodeGEMSTask &other);
};

#endif // TNODETASK_H
