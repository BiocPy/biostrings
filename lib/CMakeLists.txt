cmake_minimum_required(VERSION 3.24)

project(biostring
    VERSION 1.0.0
    DESCRIPTION "Building the biostrings shared library"
    LANGUAGES CXX)

find_package(pybind11 CONFIG)

# pybind11 method:
pybind11_add_module(biostring
    src/stringsetpool.cpp
    src/init.cpp
)

set_property(TARGET biostring PROPERTY CXX_STANDARD 17)

target_link_libraries(biostring PRIVATE pybind11::pybind11)

set_target_properties(biostring PROPERTIES
    OUTPUT_NAME lib_biostrings
    PREFIX ""
)
