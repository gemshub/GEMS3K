#include <sstream>
// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "jsonconfig.h"
#include "solmodfactory.h"

void exportSolMod(py::module& m)
{
    m.def("update_loggers", [](bool use_stdout, const std::string &logfile_name, size_t log_level)
           { gemsSettings().gems3k_update_loggers(use_stdout, logfile_name, log_level); });

    py::class_<SolModEngine>(m, "SolModEngine")
            .def(py::init<long int, long int, const std::string& >() )

            .def("SolModPTParams", &SolModEngine::SolModPTParams )
            .def("SolModActivityCoeffs", &SolModEngine::SolModActivityCoeffs )
            .def("SolModExcessProps", &SolModEngine::SolModExcessProps )
            .def("SolModIdealProps", &SolModEngine::SolModIdealProps )
            .def("SolModDarkenProps", &SolModEngine::SolModDarkenProps )
            .def("SolModStandProps", &SolModEngine::SolModStandProps )
            // Get functions
            .def("Get_SpeciesNumber", &SolModEngine::Get_SpeciesNumber )
            .def("Get_SpeciesNames", &SolModEngine::Get_SpeciesNames )
            .def("Get_SolPhaseName", &SolModEngine::Get_SolPhaseName )
            .def("Get_MixModelCode", &SolModEngine::Get_MixModelCode )
            .def("Get_MixModelType", &SolModEngine::Get_MixModelType )

            .def("GetMoleFractions", &SolModEngine::GetMoleFractions )
            .def("Get_MoleFractions", py::overload_cast<>(&SolModEngine::Get_MoleFractions) )
            .def("GetMolalities", &SolModEngine::GetMolalities )
            .def("Get_Molalities", py::overload_cast<>(&SolModEngine::Get_Molalities) )
            .def("GetlnActivities", &SolModEngine::GetlnActivities )
            .def("Get_lnActivities", py::overload_cast<>(&SolModEngine::Get_lnActivities) )
            .def("GetlnActivityCoeffs", &SolModEngine::GetlnActivityCoeffs )
            .def("Get_lnActivityCoeffs", py::overload_cast<>(&SolModEngine::Get_lnActivityCoeffs) )
            .def("GetlnConfTerms", &SolModEngine::GetlnConfTerms )
            .def("Get_lnConfTerms", py::overload_cast<>(&SolModEngine::Get_lnConfTerms) )
            .def("GetlnRecipTerms", &SolModEngine::GetlnRecipTerms )
            .def("Get_lnRecipTerms", py::overload_cast<>(&SolModEngine::Get_lnRecipTerms) )
            .def("GetlnExcessTerms", &SolModEngine::GetlnExcessTerms )
            .def("Get_lnExcessTerms", py::overload_cast<>(&SolModEngine::Get_lnExcessTerms) )
            .def("GetlnDQFTerms", &SolModEngine::GetlnDQFTerms )
            .def("Get_lnDQFTerms", py::overload_cast<>(&SolModEngine::Get_lnDQFTerms) )
            .def("GetG0Increments", &SolModEngine::GetG0Increments )
            .def("Get_G0Increments", py::overload_cast<>(&SolModEngine::Get_G0Increments) )
            .def("GetMolarVolumes", &SolModEngine::GetMolarVolumes )
            .def("Get_MolarVolumes", py::overload_cast<>(&SolModEngine::Get_MolarVolumes) )
            .def("GetPhaseVolume", &SolModEngine::GetPhaseVolume )
            .def("GetPartialPressures", &SolModEngine::GetPartialPressures )
            .def("Get_PartialPressures", py::overload_cast<>(&SolModEngine::Get_PartialPressures) )
            // Set functions
            .def("SetMoleFractions", &SolModEngine::SetMoleFractions, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesMolality", &SolModEngine::SetSpeciesMolality, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesAmounts", &SolModEngine::SetSpeciesAmounts, py::arg("val_map"), py::arg("def_val")=0. )
            .def("Set_PhaseMass", &SolModEngine::Set_PhaseMass )
            // Print to file
            .def("to_json_stream_short", &SolModEngine::to_json_stream_short )
            .def("to_json_file", &SolModEngine::to_json_file )
            .def("to_text_file", &SolModEngine::to_text_file,  py::arg("path"), py::arg("append")=false )
            .def("__repr__", [](const SolModEngine& self) { std::stringstream ss;  ss<<self; return ss.str(); })
            ;

    py::class_<SolModFactory>(m, "SolModFactory")
            .def(py::init<const std::string&>() )
            .def(py::init<const std::string&, const std::string&, const std::string&, const std::string&>() )

            .def("Get_Temperature", &SolModFactory::Get_Temperature )
            .def("Get_Pressure", &SolModFactory::Get_Pressure )
            .def("UpdateThermoData", &SolModFactory::UpdateThermoData )
            .def("SolPhase", py::overload_cast<const std::string&>(&SolModFactory::SolPhase))
            .def("Sol_Phase", py::overload_cast<std::size_t>(&SolModFactory::Sol_Phase))

            .def("Get_SolPhasesNumber", &SolModFactory::Get_SolPhasesNumber )
            .def("Get_SolPhasesNames", &SolModFactory::Get_SolPhasesNames )
            .def("Get_AllElementsNumber", &SolModFactory::Get_AllElementsNumber )
            .def("Get_AllSpeciesNumber", &SolModFactory::Get_AllSpeciesNumber )
            .def("Get_SolSpeciesNumber", &SolModFactory::Get_SolSpeciesNumber )
            .def("Get_IndexOfH2O_solvent", &SolModFactory::Get_IndexOfH2O_solvent )
            .def("GetGasesNumber", &SolModFactory::GetGasesNumber )
            .def("Get_AllPhasesNumber", &SolModFactory::Get_AllPhasesNumber )
            .def("Get_SolPhasesNumber", &SolModFactory::Get_SolPhasesNumber )

            .def("Get_AllElementNames", &SolModFactory::Get_AllElementNames )
            .def("Get_AllSpeciesNames", &SolModFactory::Get_AllSpeciesNames )
            .def("Get_AllPhasesNames", &SolModFactory::Get_AllPhasesNames )
            //.def("Get_SolPhasesNames", &SolModFactory::Get_SolPhasesNames )
            .def("Get_SpeciesInPhasesNumbers", &SolModFactory::Get_SpeciesInPhasesNumbers )

            .def("Get_StoichiometryMatrix", &SolModFactory::Get_StoichiometryMatrix )
            .def("Get_ElementMolarMasses", &SolModFactory::Get_ElementMolarMasses )
            .def("Get_IncrementsMolarG0", &SolModFactory::Get_IncrementsMolarG0 )
            .def("Get_SurfaceFreeEnergyParameter", &SolModFactory::Get_SurfaceFreeEnergyParameter )
            .def("Get_SpeciesMolarVolumes", &SolModFactory::Get_SpeciesMolarVolumes )

            .def("Get_GibbsEnergyG0", &SolModFactory::Get_GibbsEnergyG0 )
            .def("Get_MolarVolumeV0", &SolModFactory::Get_MolarVolumeV0 )
            .def("Get_MolarEnthalpyH0", &SolModFactory::Get_MolarEnthalpyH0 )
            .def("Get_MolarEnropyS0", &SolModFactory::Get_MolarEnropyS0 )
            .def("Get_HeatCapacityCp0", &SolModFactory::Get_HeatCapacityCp0 )
            .def("Get_HelmholtzEnergyA0", &SolModFactory::Get_HelmholtzEnergyA0 )
            .def("Get_InternalEnergyU0", &SolModFactory::Get_InternalEnergyU0 )

            .def("Get_SpeciesMolarMasses", &SolModFactory::Get_SpeciesMolarMasses )
            .def("Get_ElementMoleAmounts", &SolModFactory::Get_ElementMoleAmounts )
            .def("Get_SpeciesMoleAmounts", &SolModFactory::Get_SpeciesMoleAmounts )
            .def("Get_SpeciesMoleFractions", &SolModFactory::Get_SpeciesMoleFractions )
            .def("Get_PhaseMoleAmounts", &SolModFactory::Get_PhaseMoleAmounts )
            .def("Get_PhaseCarrierMoleAmounts", &SolModFactory::Get_PhaseCarrierMoleAmounts )

            .def("Get_ElementClassCodes", &SolModFactory::Get_ElementClassCodes )
            .def("Get_SpeciesClassCodes", &SolModFactory::Get_SpeciesClassCodes )
            .def("Get_SpeciesGenericClassCodes", &SolModFactory::Get_SpeciesGenericClassCodes )
            .def("Get_PhasesAggrStateCodes", &SolModFactory::Get_PhasesAggrStateCodes )

            .def("Get_SpeciesUpperBounds", &SolModFactory::Get_SpeciesUpperBounds )
            .def("Get_SpeciesLowerBounds", &SolModFactory::Get_SpeciesLowerBounds )
            .def("Get_SpeciesBoundCodes", &SolModFactory::Get_SpeciesBoundCodes )
            .def("Get_SpeciesBoundUnitCodes", &SolModFactory::Get_SpeciesBoundUnitCodes )
            .def("__repr__", [](const SolModFactory& self) { std::stringstream ss; ss << self; return ss.str(); })
            ;
}

