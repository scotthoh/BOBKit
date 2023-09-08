// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer-filter.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>

#include <future>
#include <thread>
#include <exception>

namespace py = pybind11;
using namespace clipper;

void init_ca_filter(py::module &m)
{
  // "Forward declaration" of python class to avoid C++ signatures in docstrings
  py::class_<Ca_filter> cafilter(m, "Ca_filter");
  cafilter
      .def(py::init<double>(), py::arg("sig_cut") = 3.0)
      .def("__call__", &Ca_filter::operator())
      .def_static("filter", static_cast<bool (*)(MiniMol &, const Xmap<float> &, double, bool)>(&Ca_filter::filter))
      .def_static("filter", static_cast<bool (*)(MiniMol &, double)>(&Ca_filter::filter));
}