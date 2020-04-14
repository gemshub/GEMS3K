
#include <vector>
#include <memory>
#include <iostream>
#include "node.h"

// The method replace() returns a copy of the string
// in which the first occurrence of old have been replaced with new
std::string replace( const std::string& str, const char* old_part, const char* new_part)
{
    size_t pos = str.find( old_part );
    if( pos == gstring::npos )
      return str;

    std::string res = str.substr(0, pos);
            res += new_part;
            res += str.substr( pos+strlen(old_part));
    return res;
}

// Run process of calculate equilibria into the GEMS3K side (read files)
double  CalculateEquilibriumServer( const std::string& lst_f_name )
{
    long int pmm_K2,  pmm_ITF, pmm_ITG;
    double ret=0.;
    try
    {
        auto dbr_lst_f_name  = replace( lst_f_name,"-dat.lst","-dbr-out.dat");

        // Creates TNode structure instance accessible through the "node" pointer
        std::shared_ptr<TNode> node(new TNode());

        // (1) Initialization of GEMS3K internal data by reading  files
        //     whose names are given in the lst_f_name
        if( node->GEM_init( lst_f_name.c_str() ) )
        {
            // error occured during reading the files
            cout << "error occured during reading the files" << endl;
            return ret;
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = node->GEM_run( true );
        ret  = node->GEM_CalcTime(pmm_K2,  pmm_ITF, pmm_ITG);

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  ){
            // (3) Writing results in default DBR file
            node->GEM_write_dbr( nullptr/*dbr_lst_f_name.c_str()*/, false, false, false );
            node->GEM_print_ipm( "GEMipmOK.txt" );   // possible debugging printout
        }
        else {
            // (4) possible return status analysis, error message
            node->GEM_print_ipm( "GEMipmError.txt" );   // possible debugging printout
            return ret; // GEM IPM did not converge properly - error message needed
        }

    }catch(TError& err)
    {
    }
    catch(...)
    {
    }
    return ret;
}

// Run process of calculate equilibria into the GEMS3K side (read from strings)
std::vector<std::string>  CalculateEquilibriumServer( const std::vector<std::string>& dch_ipm_dbr )
{
    std::vector<std::string> ret_mess{ "error" };
    try
    {
        // Creates TNode structure instance accessible through the "node" pointer
        std::shared_ptr<TNode> node(new TNode());

        // (1) Initialization of GEMS3K internal data from json/key-value strings
        if( node->GEM_init( dch_ipm_dbr[1], dch_ipm_dbr[2], dch_ipm_dbr[3] ) )
        {
            cout << "error occured during deserialize the data" << endl;
            ret_mess.push_back("Error occured during deserialize the data");
            return ret_mess;
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = node->GEM_run( true );
        long int pmm_K2,  pmm_ITF, pmm_ITG;
        auto time  = node->GEM_CalcTime(pmm_K2, pmm_ITF, pmm_ITG);

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  ){

            // (3) Writing results in default DBR file
            ret_mess[0] = "ipmOK";
            ret_mess.push_back(node->databr_to_string( false, false ));
            ret_mess.push_back(std::to_string(pmm_K2));  // iter
            ret_mess.push_back(std::to_string(pmm_ITF));  // iter
            ret_mess.push_back(std::to_string(pmm_ITG));  // iter
            ret_mess.push_back(std::to_string(time));

            node->GEM_print_ipm( "GEMipmOK.txt" );   // possible debugging printout
        }
        else {
            // (4) possible return status analysis, error message
            ret_mess[0] = "ipmError";
            ret_mess.push_back(node->databr_to_string( false, false ));
            ret_mess.push_back(node->code_error_IPM());
            ret_mess.push_back(node->description_error_IPM());
            ret_mess.push_back(std::to_string(time));

            node->GEM_print_ipm( "GEMipmError.txt" );   // possible debugging printout
            // GEM IPM did not converge properly - error message needed
        }

    }catch(TError& err)
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


