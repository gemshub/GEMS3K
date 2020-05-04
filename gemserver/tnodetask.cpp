#include "tnodetask.h"
#include "node.h"

NodeGEMSTask::NodeGEMSTask():
    current_task(nullptr), task_name()
{ }

NodeGEMSTask::NodeGEMSTask(const NodeGEMSTask &other)
{
    copy(other);
}

NodeGEMSTask &NodeGEMSTask::operator =(const NodeGEMSTask &other)
{
    if ( this != &other)
    {
        copy(other);
    }
    return *this;
}

std::vector<std::string> NodeGEMSTask::initData(const std::string &dch_json, const std::string &ipm_json, const std::string &dbr_json)
{
    std::vector<std::string> ret_mess;

    current_task.reset(new TNode());
    if( current_task->GEM_init( dch_json, ipm_json, dbr_json ) )
    {
        cout << "error occured during deserialize the data" << endl;
        ret_mess.push_back(std::to_string(T_ERROR_GEM));
        ret_mess.push_back("Error occured during deserialize the data");
    }
    return ret_mess;
}

std::vector<std::string> NodeGEMSTask::initData(const char *ipmfiles_lst_name)
{
    std::vector<std::string> ret_mess;

    current_task.reset(new TNode());
    if( current_task->GEM_init( ipmfiles_lst_name ) )
    {
        cout << "error occured during deserialize the data" << endl;
        ret_mess.push_back(std::to_string(T_ERROR_GEM));
        ret_mess.push_back("Error occured during deserialize the data");
    }
    return ret_mess;
}

std::vector<std::string> NodeGEMSTask::calculateEquilibrium( const std::string& new_dbr )
{
    long int NodeStatusCH = T_ERROR_GEM;
    std::vector<std::string> ret_mess;
    try
    {
        // (1) Initialization of GEMS3K internal data from DataBR
        if( !new_dbr.empty() )
        {
            current_task->databr_from_string(new_dbr);
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        NodeStatusCH = current_task->GEM_run( true );
        long int pmm_K2,  pmm_ITF, pmm_ITG;
        auto time  = current_task->GEM_CalcTime(pmm_K2, pmm_ITF, pmm_ITG);

        // (3) Writing results in default DBR file
        ret_mess.push_back(std::to_string(NodeStatusCH));
        ret_mess.push_back(current_task->databr_to_string( false, false ));
        ret_mess.push_back(std::to_string(time));

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
        {
            // (4) Addition information
            ret_mess.push_back(std::to_string(pmm_K2));  // iter
            ret_mess.push_back(std::to_string(pmm_ITF));  // iter
            ret_mess.push_back(std::to_string(pmm_ITG));  // iter

            current_task->GEM_print_ipm( "GEMipmOK.txt" );   // possible debugging printout
        }
        else
        {
            // (5) possible return status analysis, error message
            ret_mess.push_back(current_task->code_error_IPM());
            ret_mess.push_back(current_task->description_error_IPM());

            current_task->GEM_print_ipm( "GEMipmError.txt" );   // possible debugging printout
        }
    }
    catch(TError& err)
    {
        ret_mess.push_back(std::to_string(NodeStatusCH));
        ret_mess.push_back(err.title);
        ret_mess.push_back(err.mess);
    }
    catch(...)
    {
        ret_mess.push_back(std::to_string(NodeStatusCH));
        ret_mess.push_back("Undefined error");
    }

    return ret_mess;
}

void NodeGEMSTask::copy(const NodeGEMSTask &other)
{
    task_name = other.task_name;
    TNode* new_task = new TNode(*other.current_task);
    current_task.reset( new_task );
}
