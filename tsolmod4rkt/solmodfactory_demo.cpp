//--------------------------------------------------------------------
//
/// Demo test of usage of the of subset of TMulti class, configuration,
/// and related functions
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
        std::string input_system_file_list_name = "test01/gems3k-files/series1-dat.lst";
        std::string after_reading;
        if (argc >= 2 ) {
            input_system_file_list_name = argv[1];
        }
        after_reading = input_system_file_list_name;
        replace(after_reading, "-dat.lst", "-AfterReading.txt");

        // Initialize SolModFactory from the GEMS3K file set
        SolModFactory task(input_system_file_list_name);

        // Check the data read for SolModFactory initialization
        task.to_text_file(after_reading);

        std::cout << "Task: " << input_system_file_list_name <<
                     "\n T(K): " << task.Get_Temperature() << " P(bar): " << task.Get_Pressure();

        std::cout << "\nNumber of Elements: " << task.Get_AllElementsNumber();
        print_vector( "AllElementNames: ", task.Get_AllElementNames());
        print_vector( "ElementClassCodes: ", task.Get_ElementClassCodes());
        print_vector( "ElementMolarMasses: ", task.Get_ElementMolarMasses());
        print_vector( "ElementMoleAmounts: ", task.Get_ElementMoleAmounts());

        std::cout << "\nNumber of Phases: " << task.Get_AllPhasesNumber();
        std::cout << "\nNumber of SolPhases: " << task.Get_SolPhasesNumber();
        print_vector( "AllPhasesNames: ", task.Get_AllPhasesNames());
        print_vector( "PhasesAggrStateCodes: ", task.Get_PhasesAggrStateCodes());
        print_vector( "SpeciesInPhasesNumbers: ", task.Get_SpeciesInPhasesNumbers());

        print_vector( "SurfaceFreeEnergyParameter: ", task.Get_SurfaceFreeEnergyParameter());
        print_vector( "PhaseMoleAmounts: ", task.Get_PhaseMoleAmounts());

        print_vector( "SolPhasesNames: ", task.Get_SolPhasesNames());
        print_vector( "PhaseCarrierMoleAmounts: ", task.Get_PhaseCarrierMoleAmounts());

        std::cout << "\nNumber of Species: " << task.Get_AllSpeciesNumber();
        std::cout << "\nNumber of SolSpecies: " << task.Get_SolSpeciesNumber();
        std::cout << "\nIndex of water-solvent: " << task.Get_IndexOfH2O_solvent();
        std::cout << "\nGases Number: " << task.GetGasesNumber();
        print_vector( "AllSpeciesNames: ", task.Get_AllSpeciesNames());
        print_vector( "SpeciesClassCodes: ", task.Get_SpeciesClassCodes());
        print_vector( "SpeciesGenericClassCodes: ", task.Get_SpeciesGenericClassCodes());

        print_vector( "SpeciesUpperBounds: ", task.Get_SpeciesUpperBounds());
        print_vector( "SpeciesLowerBounds: ", task.Get_SpeciesLowerBounds());
        print_vector( "SpeciesBoundCodes: ", task.Get_SpeciesBoundCodes());
        print_vector( "SpeciesBoundUnitCodes: ", task.Get_SpeciesBoundUnitCodes());

        print_vector( "SpeciesMolarMasses: ", task.Get_SpeciesMolarMasses());
        print_vector( "SpeciesMoleAmounts: ", task.Get_SpeciesMoleAmounts());
        print_vector( "SpeciesMoleFractions: ", task.Get_SpeciesMoleFractions());

        print_vector( "IncrementsMolarG0: ", task.Get_IncrementsMolarG0());
        print_vector( "SpeciesMolarVolumes: ", task.Get_SpeciesMolarVolumes());
        print_vector( "MolarVolumeV0: ", task.Get_MolarVolumeV0());
        print_vector( "GibbsEnergyG0: ", task.Get_GibbsEnergyG0());
        print_vector( "MolarEnthalpyH0: ", task.Get_MolarEnthalpyH0());
        print_vector( "MolarEnropyS0: ", task.Get_MolarEnropyS0());
        print_vector( "HeatCapacityCp0: ", task.Get_HeatCapacityCp0());
        print_vector( "HelmholtzEnergyA0: ", task.Get_HelmholtzEnergyA0());
        print_vector( "InternalEnergyU0: ", task.Get_InternalEnergyU0());

        std::cout << "\nStoichiometryMatrix\n";
        for(const auto& line: task.Get_StoichiometryMatrix())   {
            for(const auto& item: line)   {
                std::cout << item << " ";
            };
            std::cout << std::endl;
        };
        std::cout << std::endl;

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

