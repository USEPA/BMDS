#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_ex1(py::module_ &);
void init_ex2(py::module_ &);


// test_binary_operators
enum Flags {
    Read = 4,
    Write = 2,
    Execute = 1
};


PYBIND11_MODULE(bleep, m) {
    m.doc() = "Simple c++ interface examples";

    // functions
    init_ex1(m);
    init_ex2(m);

    // enums
    py::enum_<Flags>(m, "Flags")
        .value("Read", Flags::Read)
        .value("Write", Flags::Write)
        .value("Execute", Flags::Execute)
        .export_values();

    // pass data from python at build-time
    #ifdef VERSION_INFO
        m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
    #else
        m.attr("__version__") = "dev";
    #endif
}
