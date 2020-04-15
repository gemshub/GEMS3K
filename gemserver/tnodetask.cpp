#include "tnodetask.h"
#include "node.h"

const char* ipm_ok_head = "ipmOK";
const char* ipm_error_head = "ipmError";
const char* error_head = "error";

TNodeTask::TNodeTask():
    current_task(new TNode()), task_name()
{ }

TNodeTask::TNodeTask(const TNodeTask &other)
{
    copy(other);
}

TNodeTask &TNodeTask::operator =(const TNodeTask &other)
{
    if ( this != &other)
    {
        copy(other);
    }
    return *this;
}

std::vector<std::string> TNodeTask::initGEM(const std::string &dch_json, const std::string &ipm_json, const std::string &dbr_json)
{
    std::vector<std::string> ret_mess;

    current_task.reset(new TNode());
    if( current_task->GEM_init( dch_json, ipm_json, dbr_json ) )
    {
        cout << "error occured during deserialize the data" << endl;
        ret_mess.push_back(error_head);
        ret_mess.push_back("Error occured during deserialize the data");
    }
    return ret_mess;
}

std::vector<std::string> TNodeTask::initGEM(const char *ipmfiles_lst_name)
{
    std::vector<std::string> ret_mess;

    current_task.reset(new TNode());
    if( current_task->GEM_init( ipmfiles_lst_name ) )
    {
        cout << "error occured during deserialize the data" << endl;
        ret_mess.push_back(error_head);
        ret_mess.push_back("Error occured during deserialize the data");
    }
    return ret_mess;
}

std::vector<std::string> TNodeTask::calculateEquilibrium( const std::string& new_dbr )
{
    std::vector<std::string> ret_mess{ error_head };
    try
    {
        if( !new_dbr.empty() )
        {
            // update dbr
            current_task->databr_from_string(new_dbr);
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = current_task->GEM_run( true );
        long int pmm_K2,  pmm_ITF, pmm_ITG;
        auto time  = current_task->GEM_CalcTime(pmm_K2, pmm_ITF, pmm_ITG);

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
        {
            // (3) Writing results in default DBR file
            ret_mess[0] = ipm_ok_head;
            ret_mess.push_back(current_task->databr_to_string( false, false ));
            ret_mess.push_back(std::to_string(pmm_K2));  // iter
            ret_mess.push_back(std::to_string(pmm_ITF));  // iter
            ret_mess.push_back(std::to_string(pmm_ITG));  // iter
            ret_mess.push_back(std::to_string(time));

            current_task->GEM_print_ipm( "GEMipmOK.txt" );   // possible debugging printout
        }
        else
        {
            // (4) possible return status analysis, error message
            ret_mess[0] = ipm_error_head;
            ret_mess.push_back(current_task->databr_to_string( false, false ));
            ret_mess.push_back(current_task->code_error_IPM());
            ret_mess.push_back(current_task->description_error_IPM());
            ret_mess.push_back(std::to_string(time));

            current_task->GEM_print_ipm( "GEMipmError.txt" );   // possible debugging printout
        }
    }
    catch(TError& err)
    {
        ret_mess.push_back(err.title);
        ret_mess.push_back(err.mess);
    }
    catch(...)
    {
        ret_mess.push_back("Undefined error");
    }

    return ret_mess;
}

void TNodeTask::copy(const TNodeTask &other)
{
    task_name = other.task_name;
    TNode* new_task = new TNode(*other.current_task);
    current_task.reset( new_task );
}
