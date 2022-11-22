#include <pybind11/pybind11.h>

namespace py = pybind11;

int sub(int i, int j) {
    return i - j;
}

void init_ex2(py::module_ &m) {
    m.def("sub", &sub, "A function which subtracts two numbers", py::arg("i"), py::arg("j"));
}
