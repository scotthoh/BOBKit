// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer-ncsbuild.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;

void init_ca_ncsbuild(py::module &m)
{
  py::class_<Ca_ncsbuild>(m, "Ca_ncsbuild")
      .def(py::init<double, double, int>(), py::arg("reliability") = 0.5, py::arg("rmsd") = 1.0, py::arg("nmin") = 12)
      .def("__call__", &Ca_ncsbuild::operator());
}