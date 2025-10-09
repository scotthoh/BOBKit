// Nanobind bindings for buccaneer-ncsbuild
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-ncsbuild.h"
#include "commons.h"
#include <nanobind/operators.h>

void add_ca_ncsbuild(nb::module_ &m) {
  nb::class_<Ca_ncsbuild>(m, "Ca_ncsbuild")
      .def(nb::init<double, double, int>(), nb::arg("reliability") = 0.5,
           nb::arg("rmsd") = 1.0, nb::arg("nmin") = 12,
           "Constructor for Ca_ncsbuild class.")
      .def("__call__", &Ca_ncsbuild::operator(), nb::arg("mol"),
           nb::arg("xmap"), nb::arg("llktarget"), nb::arg("seq"),
           "Build NCS related chains.")
      .def("__repr__",
           [](const Ca_ncsbuild &self) {
             return "<buccaneer.Ca_ncsbuild class.>";
           })
      .doc() = "Class for build NCS related chains using density.";
}