// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-merge.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;
using namespace clipper;

void init_ca_merge(py::module &m)
{
  // "Forward declaration" of python class to avoid C++ signatures in docstrings
  py::class_<Ca_merge> camerge(m, "Ca_merge");
  camerge
      .def(py::init<double>(), py::arg("reliability") = 0.5)
      .def("__call__", &Ca_merge::operator())
      .def("merge_mr", &Ca_merge::merge_mr);
}