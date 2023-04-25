
#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_test1(py::module &);


PYBIND11_MODULE(bmdscore, m) {
    m.doc() = "bmdscore c++ interface";

    // functions
    init_test1(m);

}
