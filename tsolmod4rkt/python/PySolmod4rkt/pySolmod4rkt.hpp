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
            .def("model_code", &SolModEngine::model_code )

            .def("SolModParPT", &SolModEngine::SolModParPT )
            .def("SolModActCoeff", &SolModEngine::SolModActCoeff )
            .def("SolModExcessProp", &SolModEngine::SolModExcessProp )
            .def("SolModIdealProp", &SolModEngine::SolModIdealProp )
            .def("SolModDarkenProp", &SolModEngine::SolModDarkenProp )
            .def("SolModStandProp", &SolModEngine::SolModStandProp )
            // Get functions
            .def("GetlnGamma", &SolModEngine::GetlnGamma )
            .def("GetlnGamConf", &SolModEngine::GetlnGamConf )
            .def("GetlnGamRecip", &SolModEngine::GetlnGamRecip )
            .def("GetlnGamEx", &SolModEngine::GetlnGamEx )
            .def("GetlnGamDQF", &SolModEngine::GetlnGamDQF )
            .def("GetIncrementstoG0", &SolModEngine::GetIncrementstoG0 )
            .def("GetMolarVolumes", &SolModEngine::GetMolarVolumes )
            .def("GetPhaseVolume", &SolModEngine::GetPhaseVolume )
            .def("GetPartialPressures", &SolModEngine::GetPartialPressures )
            // Set functions
            .def("SetMoleFractionsWx", &SolModEngine::SetMoleFractionsWx, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesMolality", &SolModEngine::SetSpeciesMolality, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetDCquantities", &SolModEngine::SetDCquantities, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetPhaseMasses", &SolModEngine::SetPhaseMasses )
            // Print to file
            .def("to_json_file", &SolModEngine::to_json_file )
            .def("to_text_file", &SolModEngine::to_text_file,  py::arg("path"), py::arg("append")=false )
            .def("__repr__", [](const SolModEngine& self) { std::stringstream ss; self.to_json(ss); return ss.str(); })
            ;

    py::class_<SolModFactory>(m, "SolModFactory")
            .def(py::init<const std::string&>() )
            .def(py::init<const std::string&, const std::string&, const std::string&, const std::string&>() )
            .def("UpdateThermodynamic", &SolModFactory::UpdateThermodynamic )
            .def("solution_phase", py::overload_cast<const std::string&>(&SolModFactory::solution_phase))
            .def("solution_phase", py::overload_cast<std::size_t>(&SolModFactory::solution_phase))
            //.def("__repr__", [](const TSolModMulti& self) { std::stringstream ss; /*ss << self.to_text_file();*/ return ss.str(); })
            ;
}
