// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer-build.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;

void init_ca_build(py::module &m)
{
  py::class_<Ca_build>(m, "Ca_build")
      .def(py::init<clipper::String, bool>(), py::arg("newrestype") = "ALA", py::arg("flexible") = false)
      .def_static("build", &Ca_build::build)
      .def("__call__", &Ca_build::operator());
}