#include "../code_base/bmds_helper.h"
#include <pybind11/pybind11.h>


namespace py = pybind11;

void init_test1(py::module &m) {
    m.def("add2", &add2, "A function which adds two numbers", py::arg("i"), py::arg("j"));
    m.def("version", &version, "A function that returns version number");
//    m.def("testFun", &testFun, "A function to test python objs");
    m.def("testFun", &testFun, "A function to test python objs", py::arg("test_struct"));
    m.def("pythonBMDSDicho", &pythonBMDSDicho, "Entry point to run BMDS dichotomous models", py::arg("python_dichotomous_analysis"), py::arg("python_dichotomous_model_result"), py::arg("dichotomous_GOF"), py::arg("BMDS_results"), py::arg("dicho_AOD"));
}
