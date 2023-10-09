#include <sstream>
// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "jsonconfig.h"
#include "tsolmod_multi.h"

void exportSolMod(py::module& m)
{
    m.def("update_loggers", [](bool use_stdout, const std::string &logfile_name, size_t log_level)
           { gemsSettings().gems3k_update_loggers(use_stdout, logfile_name, log_level); });

    py::class_<SolModCalc>(m, "SolModCalc")
            .def(py::init<long int, long int, const std::string& >() )
            .def("modCode", &SolModCalc::modCode )

            .def("SolModParPT", &SolModCalc::SolModParPT )
            .def("SolModActCoeff", &SolModCalc::SolModActCoeff )
            .def("SolModExcessProp", &SolModCalc::SolModExcessProp )
            .def("SolModIdealProp", &SolModCalc::SolModIdealProp )
            .def("SolModDarkenProp", &SolModCalc::SolModDarkenProp )
            .def("SolModStandProp", &SolModCalc::SolModStandProp )
            // Get functions
            .def("GetlnGamma", &SolModCalc::GetlnGamma )
            .def("GetlnGamConf", &SolModCalc::GetlnGamConf )
            .def("GetlnGamRecip", &SolModCalc::GetlnGamRecip )
            .def("GetlnGamEx", &SolModCalc::GetlnGamEx )
            .def("GetlnGamDQF", &SolModCalc::GetlnGamDQF )
            .def("GetIncrementstoG0", &SolModCalc::GetIncrementstoG0 )
            .def("GetMolarVolumes", &SolModCalc::GetMolarVolumes )
            .def("GetPhaseVolume", &SolModCalc::GetPhaseVolume )
            .def("GetPartialPressures", &SolModCalc::GetPartialPressures )
            // Set functions
            .def("SetMoleFractionsWx", &SolModCalc::SetMoleFractionsWx, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetSpeciesMolality", &SolModCalc::SetSpeciesMolality, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetDCquantities", &SolModCalc::SetDCquantities, py::arg("val_map"), py::arg("def_val")=0. )
            .def("SetPhaseMasses", &SolModCalc::SetPhaseMasses )
            // Print to file
            .def("to_json_file", &SolModCalc::to_json_file )
            .def("to_text_file", &SolModCalc::to_text_file,  py::arg("path"), py::arg("append")=false )
            .def("__repr__", [](const SolModCalc& self) { std::stringstream ss; self.to_json(ss); return ss.str(); })
            ;

    py::class_<TSolModMulti>(m, "TSolModMulti")
            .def(py::init<>() )
            .def("GEM_init", py::overload_cast<const std::string& >(&TSolModMulti::GEM_init))
            .def("GEM_init", py::overload_cast<const std::string&, const std::string&, const std::string&, const std::string&>(&TSolModMulti::GEM_init))
            .def("UpdateThermodynamic", &TSolModMulti::UpdateThermodynamic )
            .def("get_phase", py::overload_cast<const std::string&>(&TSolModMulti::get_phase))
            .def("get_phase", py::overload_cast<std::size_t>(&TSolModMulti::get_phase))
            //.def("__repr__", [](const TSolModMulti& self) { std::stringstream ss; /*ss << self.to_text_file();*/ return ss.str(); })
            ;
}
