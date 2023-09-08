// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer-correct.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;

void init_ca_correct(py::module &m)
{
  py::class_<Ca_correct>(m, "Ca_correct")
      .def(py::init<int>(), py::arg("torsion_sampling") = 23)
      .def("__call__", &Ca_correct::operator())
      .def("num_corrected", &Ca_correct::num_corrected);
}