#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <stdexcept>

namespace py = pybind11;

// Equivalent to R's new_XStringSet_from_CHARACTER C-function
// Returns: (pool (bytes), starts (numpy array), widths (numpy array))
py::tuple create_dnastringset_pool(py::list py_seqs) {
    size_t n = py_seqs.size();

    py::array_t<int32_t> np_starts(n);
    py::array_t<int32_t> np_widths(n);

    int32_t* starts_ptr = np_starts.mutable_data();
    int32_t* widths_ptr = np_widths.mutable_data();

    std::stringstream pool_stream;
    int32_t current_start = 0;
    const std::string valid_chars = "ACGTRYSWKMBDHVN-";

    for (size_t i = 0; i < n; ++i) {
        std::string s = py_seqs[i].cast<std::string>();
        int32_t current_width = static_cast<int32_t>(s.length());
        starts_ptr[i] = current_start;
        widths_ptr[i] = current_width;

        for (char &c : s) {
            c = std::toupper(c);
            if (valid_chars.find(c) == std::string::npos) {
                throw std::invalid_argument(
                    "Sequence " + std::to_string(i) + " contains invalid DNA character: " + c
                );
            }
        }

        pool_stream.write(s.c_str(), current_width);
        current_start += current_width;
    }

    py::bytes pool = py::bytes(pool_stream.str());
    return py::make_tuple(pool, np_starts, np_widths);
}

void init_stringsetpool(pybind11::module &m) {
    m.doc() = "C++ extensions for biostrings";
    m.def(
        "create_dnastringset_pool",
        &create_dnastringset_pool,
        "Efficiently create the pool and ranges for a DnaStringset from a list of strings."
    );
}
