cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(PythonInterp 3 REQUIRED)

set(PYBIND11_CPP_STANDARD -std=c++11)

find_package(pybind11 REQUIRED)

message(STATUS "Found pybind11: ${pybind11_VERSION}")
include_directories(${PYBIND11_INCLUDE_DIR})

link_directories(../build)
pybind11_add_module(pyrandspg MODULE pyrandspg.cpp)
target_link_libraries(pyrandspg PRIVATE RandSpgLib)
