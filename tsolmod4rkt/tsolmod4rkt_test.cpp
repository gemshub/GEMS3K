//--------------------------------------------------------------------
//
/// Demo test of usage of the standalone SolModEngine.
//
// Copyright (C) S.Dmytriyeva, D.Kulik
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

// Can this ifdef be moved somewhere so that it is not needed in examples?
#ifdef OVERFLOW_EXCEPT
#ifdef __linux__
#include <cfenv>
#elif _MSC_VER
#include <float.h>
#else
#include <cfenv>
#endif
#endif

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <memory>
#include "jsonconfig.h"
#include "v_service.h"
#include "solmodfactory.h"

// Thermo-time-in/series1-dat.lst
//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{

// Can this ifdef be moved somewhere so that it is not needed in examples?
// It can be removed, but now better left check for different model cases or we make addition
// application to test other models posible with  check threading and memory leaks
#if  defined(OVERFLOW_EXCEPT)
#ifdef __linux__
feenableexcept (FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
#elif _MSC_VER
    _clearfp();
    _controlfp(_controlfp(0, 0) & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW),
               _MCW_EM);
#else

#endif
#endif

    // Read gems3k logger config file                   (is this needed in the example?)
    // this is example to config loggers
    gemsSettings().gems3k_update_loggers(true, "", 1);

    try{
        // Analyzing command line arguments
        // Defaults
        std::string after_reading = "_AfterReading.txt";
        std::string input_system_file_list_name = "Thermo-time-in/series1-dat.lst";

        // Get the list file name (of DCH, IPM and DBR input files) for a GEMS3K fileset
        if (argc >= 2 )
            input_system_file_list_name = argv[1];

        // Initialize SolModFactory from the GEMS3K file set
        SolModFactory task(input_system_file_list_name);

        // Optional: Check the data read for SolModFactory initialization
        std::string task_data_file_name = input_system_file_list_name + after_reading;
        task.to_text_file( task_data_file_name );

        // Getting the number of solution phases
        auto PhSolNumber = task.Get_SolPhasesNumber();

        // Getting the list of names of solution phases
        auto PhSolNames = task.Get_SolPhasesNames();

        std::cout << "Task: " << task_data_file_name << " PhSolNumber: " << PhSolNumber
                  << " PhSolNames:" << std::endl;
        for_each(PhSolNames.begin(), PhSolNames.end(), [](const auto& element)
        { //printing using " " separator
           std::cout << element << " ";
        });
        std::cout << std::endl;

        // Example using model calculations for task (Thermo-time-in/series1-dat.lst)
        // Two feldspar phases co-existing on binodal solvus
        // In future could be test input from dbr file or direct input
        std::vector<double> x1v = { 0.674, 7e-08, 0.326 };
        std::vector<double> x2v = { 0.187, 3.5e-09, 0.813 };
        std::vector<double> lnGamma1v(3, 0.);
        std::vector<double> lnGamma2v(3, 0.);

        // Getting SolModEngine for a feldspar phase 1 by name
        auto phase1 = task.SolPhase("Alkali feldspar");
        phase1.Set_MoleFractions(x1v.data());
        // Calculating activity coefficients of end members
        phase1.SolModActivityCoeffs();
        phase1.to_text_file("solmod_act_coef.txt", true);
        phase1.Get_lnActivityCoeffs(lnGamma1v.data());
        std::cout  << lnGamma1v[0] << " "<< lnGamma1v[1] << " " << lnGamma1v[2] << " " << std::endl;

        // Getting SolModEngine for a felsdpar phase 2 by index
        size_t p2index = 2;
        auto phase2 = task.Sol_Phase( p2index );
        phase2.Set_MoleFractions(x2v.data());

        phase2.SolModActivityCoeffs();
        phase2.to_text_file("solmod_act_coef.txt", true);
        phase2.Get_lnActivityCoeffs(lnGamma2v.data());
        std::cout  << lnGamma2v[0] << " "<< lnGamma2v[1] << " " << lnGamma2v[2] << " " << std::endl;

        auto map_ideal = phase1.SolModIdealProps();
        phase1.to_text_file("solmod_act_coef.txt", true);
        for(const auto& item: map_ideal ) {
           std::cout  << item.first << " "<< item.second  << std::endl;
        }

        auto map_excess = phase1.SolModExcessProps();
        phase1.to_text_file("solmod_act_coef.txt", true);
        for(const auto& item: map_excess ) {
           std::cout  << item.first << " "<< item.second  << std::endl;
        }

        // working with maps (dicts)
        auto& phase2m = task.SolPhase("Plagioclase");
        std::map<std::string, double> x2m = {
            {"Albite", 0.187},
            {"Anorthite", 3.5e-09},
            {"Sanidine", 0.813}};

        phase2m.SetMoleFractions(x2m);
        phase2m.SolModActivityCoeffs();
        auto ln_gamma = phase2.GetlnActivityCoeffs();
        for(const auto& item: ln_gamma ) {
           std::cout  << item.first << " "<< item.second  << std::endl;
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

