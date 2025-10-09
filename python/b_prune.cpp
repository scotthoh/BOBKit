// Nanobind bindings for buccaneer-prune
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prune.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>


void add_ca_prune(nb::module_ &m) {
  nb::class_<Ca_prune>(m, "Ca_prune",
                       "Class for pruning clashing Ca chains using density.")
      .def(nb::init<double>(), nb::arg("rad") = 3.0, "Constructor with radius.")
      .def_static("prune", &Ca_prune::prune, nb::arg("mol"), nb::arg("xmap"),
                  nb::arg("rad") = 3.0,
                  "Static function to prune clashing Ca chains using density "
                  "and specified radius.")
      .def("__call__", &Ca_prune::operator(), nb::arg("mol"), nb::arg("xmap"),
           "Prune clashing Ca chains using density.")
      .def("__repr__",
           [](const Ca_prune &self) { return "<buccaneer.Ca_prune class.>"; });
}