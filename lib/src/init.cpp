#include "pybind11/pybind11.h"

namespace py = pybind11;

void init_stringsetpool(pybind11::module &);

PYBIND11_MODULE(lib_iranges, m) {
    m.doc() = "cpp implementations";

    init_stringsetpool(m);
}
