// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include <buccaneer/buccaneer-prune.h>

#include "type_conversions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_ca_prune(py::module &m) {
  py::class_<Ca_prune>(m, "Ca_prune")
      .def(py::init<double>(), py::arg("rad") = 3.0)
      .def_static("prune", &Ca_prune::prune, py::arg("mol"), py::arg("xmap"),
                  py::arg("rad") = 3.0)
      .def("__call__", &Ca_prune::operator(), py::arg("mol"), py::arg("xmap"))
      .def("__repr__",
           [](const Ca_prune &self) { return "<buccaneer.Ca_prune class.>"; });
}