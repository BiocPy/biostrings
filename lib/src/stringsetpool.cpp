#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <sstream>

namespace py = pybind11;

py::object import_iranges() {
    py::module_ iranges_module = py::module_::import("iranges");
    return iranges_module.attr("IRanges");
}

// Equivalent to R's new_XStringSet_from_CHARACTER C-function
py::tuple create_dnastringset_pool(py::list py_seqs, py::list py_names) {
    std::stringstream pool_stream;
    size_t n = py_seqs.size();
    std::vector<int> starts(n);
    std::vector<int> widths(n);
    int current_start = 0;

    // A simple validation set
    const std::string valid_chars = "ACGTRYSWKMBDHVN-";

    for (size_t i = 0; i < n; ++i) {
        std::string s = py_seqs[i].cast<std::string>();
        int current_width = s.length();
        
        starts[i] = current_start;
        widths[i] = current_width;

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
    
    // 1. Create the py::bytes object for the pool
    py::bytes pool = py::bytes(pool_stream.str());

    // 2. Create the biocpy.IRanges object
    py::object IRanges = import_iranges();
    
    // Note: This depends on the IRanges constructor signature.
    // We pass numpy arrays or lists to it.
    py::object ranges = IRanges(
        py::array(starts.size(), starts.data()),
        py::array(widths.size(), widths.data()),
        py_names
    );

    return py::make_tuple(pool, ranges);
}

void init_coverage(pybind11::module &m) {
    m.doc() = "C++ extensions for biocpy.strings";
    m.def(
        "create_dnastringset_pool",
        &create_dnastringset_pool,
        "Efficiently create the pool and ranges for a DnaStringset from a list of strings."
    );
}