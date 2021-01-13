//--------------------------------------------------------------------
// $Id: main.cpp 686 2012-06-01 14:10:22Z kulik $
//
// Demo test of usage of the TNode class for implementing a simple
// batch-like calculation of equilibria using text file input and
// GEM-IPM-3 numerical kernel.

// TNode class implements a simple C/C++ interface of GEMS3K code.
// It works with DATACH and work DATABR structures and respective
// DCH (chemical system definition) and DBR (recipe or data bridge)
// data files. In addition, the program reads an IPM input file which
// can be used for tuning up numerical controls of GEM IPM-3 algorithm
// and for setting up the parameters of non-ideal mixing models.
//
// Copyright (C) 2007-2012 D.Kulik, S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <memory>
#include "nodearray.h"
#include "v_detail.h"
#include "args_tool.h"


void show_usage( const std::string &name );
int extract_args( int argc, char* argv[], std::string& input_lst_path, std::string& dbr_lst_path, GEMS3KImpexData& export_data );


// -j -c -i solvus-in/series1-dat.lst -e solvus-out/series1-dat.lst
// -t -i solvus-in/series1-dat.lst -e solvus-kv/series1-dat.lst

// -t -c -i Kaolinite-in/pHtitr-dat.lst -e Kaolinite-out/pHtitr-dat.lst
// -j  -c -i Kaolinite-in/pHtitr-dat.lst -e Kaolinite-json/pHtitr-dat.lst
// -t  -c -i Kaolinite-json/pHtitr-dat.lst -e Kaolinite-out/pHtitr-dat.lst

// -t -c -i process/Kaolinite-dat.lst -e process-kv/Kaolinite-dat.lst -l process/Kaolinite-dbr.lst

//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{
    try{
        std::string input_lst_path;
        std::string  dbr_lst_path;
        GEMS3KImpexData export_data;

        if( extract_args( argc, argv, input_lst_path, dbr_lst_path, export_data ))
            return 1;

        // Creates TNodeArray structure instance accessible through the "node_arr" pointer
        std::shared_ptr<TNodeArray> node_arr  = std::make_shared<TNodeArray>( export_data.nIV );

        // (1) Initialization of GEMS3K internal data by reading  files
        //     whose names are given in the input_system_file_list_name
        if( node_arr->GEM_init( input_lst_path.c_str(), nullptr, nullptr, false ) )
        {
            // error occured during reading the files
            std::cout << "error occured during reading the files" << std::endl;
            return 1;
        }

        TestModeGEMParam calc_param;  // use default data
        FILE* diffile = fopen( "tools-ICdif-log.dat", "w+" );
        if( !diffile )
            return 1;

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        if( !node_arr->CalcIPM_List( calc_param, 0, export_data.nIV-1, diffile ) )
        {
            std::cout << "error occured during inital calculation" << std::endl;
            return 1;
        }

        // (3) Writing results in defined output format
        ProcessProgressFunction messageF = [](const std::string& , long ){
            //std::cout << "TProcess GEM3k output" <<  message.c_str() << point << std::endl;
            return false;
        };

        auto dbr_list =  node_arr->genGEMS3KInputFiles(  export_data.ipmfiles_lst_name, messageF, export_data.nIV,
                                                         export_data.io_mode, export_data.brief_mode,
                                                         export_data.with_comments, export_data.putNodT1, export_data.add_mui );

        // (4) Here a possible loop on more dbr files begins
        if( !dbr_lst_path.empty() )
        {
            GEMS3KGenerator input_data( input_lst_path );
            auto number_files = input_data.load_dbr_lst_file( dbr_lst_path );

            std::string export_dir, export_name, export_ext, tmp;
            u_splitpath( export_data.ipmfiles_lst_name, export_dir, export_name, export_ext );

            std::string dbr_name;
            u_splitpath( dbr_lst_path, tmp, dbr_name, export_ext );
            auto export_dbr_lst_name = u_makepath( export_dir, dbr_name, export_ext);
            std::ofstream export_dbr_fout( export_dbr_lst_name );

            for(size_t dbr_index=0; dbr_index < number_files; dbr_index++ )
            {
                auto dbr_file = input_data.get_next_dbr_file( dbr_index );
                if( dbr_file.empty() )
                    continue;

                // (5) Reading the next DBR file with different input composition or temperature
                node_arr->GEMS3k_read_dbr( 0, dbr_file, input_data.files_mode() );

                if( !node_arr->CalcIPM_One( calc_param, 0, diffile ) )
                {
                    std::cout << "error occured during calculation: " << dbr_file <<  std::endl;
                    return 1;
                }

                u_splitpath( dbr_file, tmp, export_name, export_ext );
                export_ext = GEMS3KGenerator::ext( export_data.io_mode );
                auto export_dbr_file = u_makepath( export_dir, export_name, export_ext );

                if( dbr_index )
                    export_dbr_fout << ",";
                export_dbr_fout <<  " \"" << export_name << "." << export_ext << "\"";

                node_arr->GEMS3k_write_dbr( export_dbr_file.c_str(), export_data.io_mode,
                                            export_data.with_comments,export_data.brief_mode );
            }
        }

        return 0;
    }
    catch(TError& err)
    {
        std::cout  << err.title << err.mess << std::endl;
    }
    catch(std::exception& e)
    {
        std::cout  << "std::exception: " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout  << "unknown exception" << std::endl;
    }
    return -1;
}

void show_usage( const std::string &name )
{
    std::cout << "Usage: " << name << " [ option(s) ] -i|--import-from LST_IMPORT -e|--export-to LST_EXPORT [ -l|--dbr-list DBR_LIST ]"
              << "\nRecalculate task and export to other mode\n"
              << "Options:\n"
              << "\t-h,\t--help\t\tshow this help message\n\n"
                 // file type
              << "\t-j,\t--json        \twrite IPM, DCH and DBR files in json mode (default) \n"
              << "\t-n,\t--nlohmanjson \twrite IPM, DCH and DBR files in json mode  \n"
              << "\t-t,\t--key-value   \twrite IPM, DCH and DBR files in txt mode \n"
              << "\t-b,\t--binary      \twrite IPM, DCH and DBR files in binary mode \n\n"
                 // method
              << "\t-d,\t--brife       \tdo not write data items that contain only default values (default false) \n"
              << "\t-c,\t--comments    \twrite files with comments for all data entries ( in text mode ) (default false) \n\n"
              << std::endl;
}


int extract_args( int argc, char* argv[], std::string& input_lst_path,
                  std::string& dbr_lst_path, GEMS3KImpexData& export_data )
{
    int i=0;
    for( i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_usage( "kva2json" );;
            return 1;
        }
        else if ((arg == "-j") || (arg == "--json"))
        {
            export_data.io_mode = GEMS3KGenerator::f_json;
        }
        else if ((arg == "-n") || (arg == "--nlohmanjson"))
        {
            export_data.io_mode = GEMS3KGenerator::f_json;
        }
        else if ((arg == "-t") || (arg == "--key-value"))
        {
            export_data.io_mode = GEMS3KGenerator::f_key_value;
        }
        else if ((arg == "-b") || (arg == "--binary"))
        {
            export_data.io_mode = GEMS3KGenerator::f_binary;
        }
        else if ((arg == "-d") || (arg == "--brife"))
        {
            export_data.brief_mode = true;
        }
        else if ((arg == "-c") || (arg == "--comments"))
        {
            export_data.with_comments = true;
        }
        else if ((arg == "-i") || (arg == "--import-from"))
        {
            if (i + 1 < argc) {
                input_lst_path = argv[++i];
            } else {
                std::cerr << "--import-from option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-l") || (arg == "--dbr-list"))
        {
            if (i + 1 < argc) {
                dbr_lst_path = argv[++i];
            } else {
                std::cerr << "--dbr-list option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-e") || (arg == "-e|--export-to"))
        {
            if (i + 1 < argc) {
                export_data.ipmfiles_lst_name = argv[++i];
            } else {
                std::cerr << "-e|--export-to option requires one argument." << std::endl;
                return 1;
            }
        }
    }
    if( input_lst_path.empty() || export_data.ipmfiles_lst_name.empty() )
        return 1;
    return 0;
}
