
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
        ret  = node->GEM_CalcTime();

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
double  CalculateEquilibriumServer( const std::vector<std::string>& dch_ipm_dbr, std::string& dbr_result )
{
    double ret=0.;
    dbr_result = "";

    try
    {
        // Creates TNode structure instance accessible through the "node" pointer
        std::shared_ptr<TNode> node(new TNode());

        // (1) Initialization of GEMS3K internal data from json/key-value strings
        if( node->GEM_init( dch_ipm_dbr[1], dch_ipm_dbr[2], dch_ipm_dbr[3] ) )
        {
            cout << "error occured during deserialize the data" << endl;
            return ret;
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = node->GEM_run( true );
        ret  = node->GEM_CalcTime();

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  ){
            // (3) Writing results in default DBR file
            dbr_result = node->databr_to_string( false, false );
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


