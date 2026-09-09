#pragma once

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Forward declare module builders in various source files
void init_planar_locator(py::module& m);
void init_mesh(py::module& m);
void init_geodesic(py::module& m);
void init_parallel_taskflow(py::module& m);
void init_meshio(py::module& m);