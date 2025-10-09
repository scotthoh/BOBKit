// Nanobind bindings for buccaneer-build
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-build.h"
#include "commons.h"
#include <nanobind/operators.h>

void add_ca_build(nb::module_ &m) {
  nb::class_<Ca_build>(m, "Ca_build")
      .def(nb::init<clipper::String, bool>(), nb::arg("newrestype") = "ALA",
           nb::arg("flexible") = false, "Constructor for Ca_build.")
      .def_static("build", &Ca_build::build, nb::arg("mol"), nb::arg("xmap"),
                  nb::arg("newrestype") = "ALA", nb::arg("flexible") = false,
                  "Build Ca chains using density, static function.")
      .def("__call__", &Ca_build::operator(), nb::arg("mol"), nb::arg("xmap"),
           "Build Ca chains using density.")
      .def("__repr__",
           [](const Ca_build &self) {
             return "<buccaneer.Ca_build class>";
           })
      .doc() = "Class for building Ca chains using density.";
}