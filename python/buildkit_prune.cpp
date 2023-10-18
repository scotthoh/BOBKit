// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer/buccaneer-prune.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;

void init_ca_prune(py::module &m)
{
  py::class_<Ca_prune>(m, "Ca_prune")
      .def(py::init<double>(), py::arg("rad") = 3.0)
      .def_static("prune", &Ca_prune::prune)
      .def("__call__", &Ca_prune::operator());
}