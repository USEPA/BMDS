#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

void init_ex1(py::module_ &m) {
    m.def("add", &add, "A function which adds two numbers", py::arg("i"), py::arg("j"));
}
