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
            .def("Get_MAmodel_code", &SolModEngine::Get_MAmodel_code )
            .def("Get_MAmodel_type", &SolModEngine::Get_MAmodel_type )

            .def("GetlnActivityCoeffs", &SolModEngine::GetlnActivityCoeffs )
            .def("GetlnConfTerms", &SolModEngine::GetlnConfTerms )
            .def("GetlnRecipTerms", &SolModEngine::GetlnRecipTerms )
            .def("GetlnExcessTerms", &SolModEngine::GetlnExcessTerms )
            .def("GetlnDQFTerms", &SolModEngine::GetlnDQFTerms )
            .def("GetG0Increments", &SolModEngine::GetG0Increments )
            .def("GetMolarVolumes", &SolModEngine::GetMolarVolumes )
            .def("GetPhaseVolume", &SolModEngine::GetPhaseVolume )
            .def("GetPartialPressures", &SolModEngine::GetPartialPressures )
            // Set functions
            .def("SetMoleFractions", &SolModEngine::SetMoleFractions, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesMolality", &SolModEngine::SetSpeciesMolality, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesAmounts", &SolModEngine::SetSpeciesAmounts, py::arg("val_map"), py::arg("def_val")=0. )
            .def("Set_PhaseMass", &SolModEngine::Set_PhaseMass )
            // Print to file
            .def("to_json_file", &SolModEngine::to_json_file )
            .def("to_text_file", &SolModEngine::to_text_file,  py::arg("path"), py::arg("append")=false )
            .def("__repr__", [](const SolModEngine& self) { std::stringstream ss; self.to_json(ss); return ss.str(); })
            ;

    py::class_<SolModFactory>(m, "SolModFactory")
            .def(py::init<const std::string&>() )
            .def(py::init<const std::string&, const std::string&, const std::string&, const std::string&>() )
            .def("Get_Temperature", &SolModFactory::Get_Temperature )
            .def("Get_Pressure", &SolModFactory::Get_Pressure )
            .def("UpdateThermoData", &SolModFactory::UpdateThermoData )
            .def("Get_SolPhasesNumber", &SolModFactory::Get_SolPhasesNumber )
            .def("Get_SolPhasesNames", &SolModFactory::Get_SolPhasesNames )
            .def("SolPhase", py::overload_cast<const std::string&>(&SolModFactory::SolPhase))
            .def("Sol_Phase", py::overload_cast<std::size_t>(&SolModFactory::Sol_Phase))
            //.def("__repr__", [](const TSolModMulti& self) { std::stringstream ss; /*ss << self.to_text_file();*/ return ss.str(); })
            ;
}
