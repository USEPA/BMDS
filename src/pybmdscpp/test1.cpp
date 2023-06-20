#include "../code_base/bmds_helper.h"
#include <pybind11/pybind11.h>


namespace py = pybind11;

void init_test1(py::module &m) {
    m.def("add2", &add2, "A function which adds two numbers", py::arg("i"), py::arg("j"));
    m.def("version", &version, "A function that returns version number");
//    m.def("testFun", &testFun, "A function to test python objs");
    m.def("testFun", &testFun, "A function to test python objs", py::arg("test_struct"));
    m.def("pythonBMDSDicho", &pythonBMDSDicho, "Entry point to run BMDS dichotomous models", py::arg("python_dichotomous_analysis"), py::arg("python_dichotomous_model_result"));
    m.def("pythonBMDSDichoMA", &pythonBMDSDichoMA, "Entry point to run BMDS dichotomous MA", py::arg("python_dichotomousMA_analysis"), py::arg("python_dichotomousMA_result"));
    m.def("pythonBMDSCont", &pythonBMDSCont, "Entry point to run BMDS continuous models", py::arg("python_continuous_analysis"), py::arg("python_continuous_model_result"));
    m.def("pythonBMDSMultitumor", &pythonBMDSMultitumor, "Entry point to run Multitumor analysis", py::arg("python_multitumor_analysis"), py::arg("python_multitumor_result"));
}
