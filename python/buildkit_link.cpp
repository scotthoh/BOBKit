// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer-link.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

void init_ca_link(py::module &m)
{
  py::class_<Ca_link>(m, "Ca_link")
      .def(py::init<clipper::ftype, int>(), py::arg("rad_link") = 5.0, py::arg("torsion_sampling") = 24)
      .def("__call__", &Ca_link::operator())
      .def("num_linked", &Ca_link::num_linked);
}