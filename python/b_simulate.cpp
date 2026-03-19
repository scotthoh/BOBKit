// Nanobind bindings for simulate-lib
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/simulate-lib.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>

void add_map_simulate(nb::module_ &m) {
  nb::class_<MapSimulate>(m, "MapSimulate")
      .def(nb::init<int, int>(), nb::arg("nresbins") = 100,
           nb::arg("binmin") = 20,
           "Constructor from number of resolution bins and minimum number of "
           "bins.")
      .def("__call__", &MapSimulate::operator(), nb::arg("sim_f"),
           nb::arg("sim_hl"), nb::arg("ref_f"), nb::arg("ref_hl"),
           nb::arg("wrk_f"), nb::arg("wrk_hl"), "Simulate map.")
      .def("__repr__",
           [](const MapSimulate &self) {
             return "<buccaneer.MapSimulate class.>";
           })
      .doc() =
      "Map simulation class.\n"
      "This class simulates a set of HL coeffs for a reference structure "
      "(A & B only) matching the properties of the coefficients of a work "
      "structure.";
}