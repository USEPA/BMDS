#include "../code_base/bmds_helper.h"
#include <pybind11/pybind11.h>


namespace py = pybind11;

void init_test1(py::module &m) {
    m.def("add2", &add2, "A function which adds two numbers", py::arg("i"), py::arg("j"));
    m.def("version", &version, "A function that returns version number");
}
