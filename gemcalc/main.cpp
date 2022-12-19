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

#include <ctime>
#include <cmath>
#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include "node.h"
#include "v_detail.h"


//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{
    try{

        // Analyzing command line arguments
        // Default arguments
        char input_system_file_list_name[256] = "system-dat.lst";

        // list of DCH, IPM and DBR input files for initializing GEMS3K
        if (argc >= 2 )
            strncpy( input_system_file_list_name, argv[1], 256);

        // Creates TNode structure instance accessible through the "node" pointer
        std::shared_ptr<TNode> node( new TNode() );

        // (1) Initialization of GEMS3K internal data by reading  files
        //     whose names are given in the input_system_file_list_name
        if( node->GEM_init( input_system_file_list_name ) )
        {
            // error occured during reading the files
            std::cout << "error occured during reading the files" << std::endl;
            return 1;
        }

        // Getting direct access to work node DATABR structure which exchanges the
        // data with GEM IPM3 (already filled out by reading the DBR input file)
        DATABR* dBR = node->pCNode();

        // test internal functions
        long int xCa = node->IC_name_to_xCH("Ca");
        long int xCa_ion = node->DC_name_to_xCH("Ca+2");
        long int xCal = node->DC_name_to_xCH("Cal");
        long int xaq = node->Ph_name_to_xCH("aq_gen");
        long int xCalcite = node->Ph_name_to_xCH("Calcite");
        long int xbCa = node->IC_xCH_to_xDB(xCa);
        long int xbCa_ion = node->DC_xCH_to_xDB(xCa_ion);
        long int xbCal = node->DC_xCH_to_xDB(xCal);
        long int xbaq = node->Ph_xCH_to_xDB(xaq);
        long int xbCalcite = node->Ph_xCH_to_xDB(xCalcite);

        std::cout << "          CH  BR" << std::endl;
        std::cout << " Ca       " << xCa << "   " << xbCa << std::endl;
        std::cout << " Ca+2     " << xCa_ion << "   " << xbCa_ion << std::endl;
        std::cout << " Cal     " << xCal << "  " << xbCal << std::endl;
        std::cout << " aq_gen   " << xaq << "   " << xbaq << std::endl;
        std::cout << " Calcite  " << xCalcite << "   " << xbCalcite << std::endl;
        std::cout << std::setprecision(7) << std::setw(10) << std::endl;

        // Asking GEM to run with automatic initial approximation
        dBR->NodeStatusCH = NEED_GEM_AIA;

        //node->Ph_Enthalpy(0);
        //node->Ph_Enthalpy(1);
        //node->Ph_Enthalpy(2);
        //node->Ph_Enthalpy(3);

        node->GEM_print_ipm( "BeforeCalcPhase.txt" );   // possible debugging printout

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = node->GEM_run( false );

        GEMS3KGenerator input_data(input_system_file_list_name);
        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
        {    // (3) Writing results in default DBR file
            node->GEM_write_dbr( nullptr, input_data.files_mode(), true, false );
            node->GEM_print_ipm( "AfterCalcPhase.txt" );   // possible debugging printout
        }
        else {
            // (4) possible return status analysis, error message
            node->GEM_print_ipm( nullptr );   // possible debugging printout
            return 5; // GEM IPM did not converge properly - error message needed
        }
        std::cout << "SatIndx aq" << node->Ph_SatInd(0) << " gas " << node->Ph_SatInd(1) << " s1 " << node->Ph_SatInd(2)
             << " s2 " << node->Ph_SatInd(3) << " s3 " << node->Ph_SatInd(4) << " s4 " << node->Ph_SatInd(5) << std::endl;
   
        // test internal functions
        std::cout << "Ph_Volume   Aq: " << node->Ph_Volume(xbaq) <<  " Calcite: " << node->Ph_Volume(xbCalcite) << std::endl;
        std::cout << "Ph_Mass     Aq: " << node->Ph_Mass(xbaq) <<  " Calcite: " << node->Ph_Mass(xbCalcite) << std::endl;
        std::cout << "Ph_SatInd   Aq: " << node->Ph_SatInd(xbaq) <<  " Calcite: " << node->Ph_SatInd(xbCalcite) << std::endl;

        std::cout << std::endl;
        std::cout << "Ca+2    Get_nDC  " << node->Get_nDC(xbCa_ion) <<  " DC_n  " << node->DC_n(xCa_ion) << std::endl;
        std::cout << "Cal     Get_nDC  " << node->Get_nDC(xbCal) <<  " DC_n  " << node->DC_n(xCal) << std::endl;
        std::cout << "Ca+2    Get_muDC " << node->Get_muDC(xbCa_ion) <<  " DC_mu " << node->DC_mu(xCa_ion) << std::endl;
        std::cout << "Cal     Get_muDC " << node->Get_muDC(xbCal) <<  " DC_mu " << node->DC_mu(xCal) << std::endl;
        std::cout << "Ca+2    Get_aDC  " << node->Get_aDC(xbCa_ion) <<  " DC_a  " << node->DC_a(xCa_ion) << std::endl;
        std::cout << "Cal     Get_aDC  " << node->Get_aDC(xbCal) <<  " DC_a  " << node->DC_a(xCal) << std::endl;
        std::cout << "Ca+2    Get_cDC  " << node->Get_cDC(xbCa_ion) <<  " DC_c  " << node->DC_c(xCa_ion) << std::endl;
        std::cout << "Cal     Get_cDC  " << node->Get_cDC(xbCal) <<  " DC_c  " << node->DC_c(xCal) << std::endl;
        std::cout << "Ca+2    Get_gDC  " << node->Get_gDC(xbCa_ion) <<  " DC_g  " << node->DC_g(xCa_ion) << std::endl;
        std::cout << "Cal     Get_gDC  " << node->Get_gDC(xbCal) <<  " DC_g  " << node->DC_g(xCal) << std::endl;
        std::cout << std::endl;
        std::cout << "G0   Ca+2: " << node->DC_G0( xCa_ion, node->cP(), node->cTK(), false ) <<  " Cal: " << node->DC_G0( xCal, node->cP(), node->cTK(), false ) << std::endl;
        std::cout << "V0   Ca+2: " << node->DC_V0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_V0( xCal, node->cP(), node->cTK() ) << std::endl;
        std::cout << "H0   Ca+2: " << node->DC_H0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_H0( xCal, node->cP(), node->cTK() ) << std::endl;
        std::cout << "S0   Ca+2: " << node->DC_S0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_S0( xCal, node->cP(), node->cTK() ) << std::endl;
        std::cout << "Cp0  Ca+2: " << node->DC_Cp0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_Cp0( xCal, node->cP(), node->cTK() ) << std::endl;
        std::cout << "A0   Ca+2: " << node->DC_A0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_A0( xCal, node->cP(), node->cTK() ) << std::endl;
        std::cout << "U0   Ca+2: " << node->DC_U0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_U0( xCal, node->cP(), node->cTK() ) << std::endl;

        // Here a possible loop on more input recipes begins
        if( argc >= 3 )
        {
            std::string  input_recipes_file_list_name = argv[2];
            std::string  NextRecipeFileName, NextRecipeOutFileName;
            // list of additional recipes (dbr.dat files) e.g. for simulation
            // of a titration or another irreversible process
            // Reading list of recipes names from file
            auto nRecipes = input_data.load_dbr_lst_file( input_recipes_file_list_name );

            for(size_t cRecipe=0; cRecipe < nRecipes; cRecipe++ )
            {
                NextRecipeFileName = input_data.get_next_dbr_file( cRecipe );
                if( NextRecipeFileName.empty() )
                    continue;

                // (5) Reading the next DBR file with different input composition or temperature
                if( node->GEM_read_dbr( NextRecipeFileName.c_str(), input_data.files_mode() ) )
                    continue;

                // Asking GEM IPM2 to run (faster) with smart initial approximation
                dBR->NodeStatusCH = NEED_GEM_SIA;

                NodeStatusCH = node->GEM_run( false );

                if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
                {
                    NextRecipeOutFileName = NextRecipeFileName + ".nc.out";
                    node->GEM_write_dbr( NextRecipeOutFileName.c_str(), input_data.files_mode(), false, false );
                    NextRecipeOutFileName = NextRecipeFileName + ".nc.Dump.out";
                    node->GEM_print_ipm( NextRecipeOutFileName.c_str() );
                }
                else {
                    // error message, debugging printout
                    NextRecipeOutFileName = NextRecipeFileName + ".Dump.out";
                    node->GEM_print_ipm( NextRecipeOutFileName.c_str() );
                    //              return 5; // GEM IPM did not converge properly - error message needed
                }
            }
        }
        // end of possible loop on input recipes

        // End of example
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
    return 1;
}

//---------------------------------------------------------------------------
// end of main.cpp for TNode class usage - GEMS3K single calculation example
