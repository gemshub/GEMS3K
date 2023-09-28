//--------------------------------------------------------------------
// $Id: main.cpp 686 2012-06-01 14:10:22Z kulik $
//
// Demo test of usage of the standalone TSolMod.

// TNode class implements a simple C/C++ interface of GEMS3K code.
// It works with DATACH and work DATABR structures and respective
// DCH (chemical system definition) and DBR (recipe or data bridge)
// data files. In addition, the program reads an IPM input file which
// can be used for tuning up numerical controls of GEM IPM-3 algorithm
// and for setting up the parameters of non-ideal mixing models.
//
// Copyright (C) D.Kulik, S.Dmytriyeva
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
#include "tsolmod_multi.h"

// solvus-in/series1-dat.lst
// ThermoTest-in/pHtitr-json.lst
// ThermoTest-in/pHtitr-fun.lst
// Thermo-time-json/series1-dat.lst
// Neutral-fun/Neutral-dat.lst


// Thermo-time-in/series1-dat.lst
//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{

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

    // Read config file
    gemsSettings().gems3k_update_loggers(true, "", 1);

    try{
        // Analyzing command line arguments
        // Default arguments
        char input_system_file_list_name[256] = "Thermo-time-in/series1-dat.lst";

        // list of DCH, IPM and DBR input files for initializing GEMS3K
        if (argc >= 2 )
            strncpy( input_system_file_list_name, argv[1], 256);

        // Creates TSolModMulti structure instance accessible through the "multi" pointer
        std::shared_ptr<TSolModMulti> multi(new TSolModMulti());

        // (1) Initialization of GEMS3K internal data by reading  files
        //     whose names are given in the input_system_file_list_name
        if( multi->GEM_init(input_system_file_list_name) )
        {
            std::cout << "error occured during reading the files" << std::endl;
            return 1;
        }

        // Trace data after input
        multi->to_text_file( "AfterRead.txt" );

        // Example using model calculations for task (Thermo-time-in/series1-dat.lst)
        // In future could be test input from dbr file
        std::vector<double> x = { 0.673936901807716, 6.97280484339258e-08, 0.326063028464235 };
        std::vector<double> lngam(3, 0.);

        auto phase = multi->get_phase("Alkali feldspar");
        phase.Set_MoleFractionsWx(x.data());

        phase.SolModActCoeff();
        phase.to_text_file("solmod_act_coef.txt", true);
        phase.Get_lnGamma(lngam.data());
        std::cout  << lngam[0] << " "<< lngam[1] << " " << lngam[2] << " " << std::endl;

        auto map_ideal = phase.SolModIdealProp();
        phase.to_text_file("solmod_act_coef.txt", true);
        for(const auto& item: map_ideal ) {
           std::cout  << item.first << " "<< item.second  << std::endl;
        }

        auto map_excess = phase.SolModExcessProp();
        phase.to_text_file("solmod_act_coef.txt", true);
        for(const auto& item: map_excess ) {
           std::cout  << item.first << " "<< item.second  << std::endl;
        }

        phase.SolModActCoeff();
        phase.to_text_file("solmod_act_coef.txt", true);
        phase.Get_lnGamma(lngam.data());
        std::cout  << lngam[0] << " "<< lngam[1] << " " << lngam[2] << " " << std::endl;

        auto& phase2 = multi->get_phase("Plagioclase");
        std::map<std::string, double> wx = {
            {"Albite", 0.186993363098213},
            {"Anorthite", 3.45294711467247e-09},
            {"Sanidine", 0.81300663344884}};

        phase2.SetMoleFractionsWx(wx);
        phase2.SolModActCoeff();
        auto ln_gamma = phase2.GetlnGamma();
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

