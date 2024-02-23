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

#include <iostream>
#include <vector>
#include "jsonconfig.h"
#include "solmodfactory.h"
#include "v_service.h"

void print_vector( const std::string& header, const std::vector<std::string>& data);

template < class T >
void print_vector( const std::string& header, const std::vector<T>& data) {
    std::cout << header << "\n";
    for(const auto& item: data)   {
        std::cout << item << " ";
    };
    std::cout << std::endl;
}

void print_vector( const std::string& header, const std::vector<std::string>& data) {
    std::cout << "\n" << header << "\n";
    for(const auto& item: data)   {
        std::cout << "'" << item << "' ";
    };
    std::cout << std::endl;
}

// Thermo-time-all/series1-dat.lst
//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{
    // This is example to config loggers
    // param show/hide logging to stdout
    //       add logging to rotating file name (hide if empty)
    //       set login level for all loggers (trace = 0, debug = 1, info = 2, warn = 3, err = 4, critical = 5, off = 6)
    gemsSettings().gems3k_update_loggers(true, "", 3);

    try{
        // Analyzing command line arguments  (with defaults)
        std::string input_system_file_list_name = "Thermo-time-all/series1-dat.lst";
        std::string after_reading;
        if (argc >= 2 ) {
            input_system_file_list_name = argv[1];
        }
        after_reading = input_system_file_list_name;
        replace(after_reading, "-dat.lst", "-AfterReading.txt");

        // Initialize SolModFactory from the GEMS3K file set
        SolModFactory task(input_system_file_list_name);

        // Check the data read for SolModFactory initialization
        task.to_text_file( after_reading );

        std::cout << "Task: " << input_system_file_list_name <<
                     "\n T(K): " << task.Get_Temperature() << " P(bar): " << task.Get_Pressure();

        std::cout << "\nNumber of Phases: " << task.Get_AllPhasesNumber();
        std::cout << ";  of SolPhases: " << task.Get_SolPhasesNumber();
        print_vector( "AllPhasesNames: ", task.Get_AllPhasesNames());

        for(size_t k=0; k<task.Get_AllPhasesNumber(); ++k) {
          auto phase = task.Sol_Phase(k);

          std::cout << "\nPhase: '" << phase.Get_SolPhaseName() << "'; mixing/activity model type: '"
                    << phase.Get_MixModelType() << "'; model code: '" <<  phase.Get_MixModelCode()
                    << "'; N endmembers: " << phase.Get_SpeciesNumber();
          print_vector( "SpeciesNames: ", phase.Get_SpeciesNames());

          // Calculate activity coefficients of endmembers (components, species)
          phase.SolModActivityCoeffs();

          // Calculate (a dict) of ideal properties of mixing in the phase2
          auto map_ideal = phase.SolModIdealProps();
          std::cout << "Ideal properties of mixing in phase:  \n";
          for(const auto& item: map_ideal ) {
              std::cout << "   '" << item.first << "': " << item.second << "; ";
          }

          print_vector( "\nMoleFractions: ", phase.Get_MoleFractions());
          print_vector( "Molalities: ", phase.Get_Molalities());
          print_vector( "lnActivities: ", phase.Get_lnActivities());
          print_vector( "lnActivityCoeffs: ", phase.Get_lnActivityCoeffs());
          print_vector( "lnConfTerms: ", phase.Get_lnConfTerms());
          print_vector( "lnRecipTerms: ", phase.Get_lnRecipTerms());
          print_vector( "lnExcessTerms: ", phase.Get_lnExcessTerms());
          print_vector( "lnDQFTerms: ", phase.Get_lnDQFTerms());
          print_vector( "G0Increments: ", phase.Get_G0Increments());
          print_vector( "MolarVolumes: ", phase.Get_MolarVolumes());
          print_vector( "PartialPressures: ", phase.Get_PartialPressures());
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

