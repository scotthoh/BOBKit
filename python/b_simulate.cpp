// Wrapper for simulate-lib
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/simulate-lib.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;

void init_map_simulate(py::module &m) {
  py::class_<MapSimulate>(m, "MapSimulate")
      .def(py::init<int, int>(), py::arg("nresbins") = 100,
           py::arg("binmin") = 20,
           "Constructor from number of resolution bins and minimum number of "
           "bins.")
      .def("__call__", &MapSimulate::operator(), py::arg("sim_f"),
           py::arg("sim_hl"), py::arg("ref_f"), py::arg("ref_hl"),
           py::arg("wrk_f"), py::arg("wrk_hl"), "Simulate map.")
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