
#include <pybind11/pybind11.h>
#include "bmds_helper.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_test1(py::module &);


PYBIND11_MODULE(bmdscore, m) {
    m.doc() = "bmdscore c++ interface";
    py::class_<test_struct>(m, "test_struct")
       .def(py::init<>())
       .def_readwrite("BMD", &test_struct::BMD)
       .def_readwrite("n", &test_struct::n)
       .def_readwrite("validResult", &test_struct::validResult);
    // functions
    init_test1(m);

}
