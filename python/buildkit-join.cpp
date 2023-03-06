// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer-join.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <future>
#include <thread>
#include <exception>

namespace py = pybind11;

void init_ca_join(py::module &m)
{
  // "Forward declaration" of python classes to avoid C++ signatures in docstrings
  py::class_<Ca_join> pyCaJoin(m, "Ca_join");

  pyCaJoin
      .def(py::init<double &, double &>(), py::arg("rad_merge") = 2.0, py::arg("rad_join") = 2.0)
      .def(
          "__call__", [](Ca_join &self, clipper::MiniMol &mol)
          { return self(mol); }, // problem segmentation fault
          py::arg("mol"),
          py::return_value_policy::reference_internal);
}
