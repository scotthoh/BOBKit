// Wrapper for buccaneer-ncsbuild
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-ncsbuild.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"

namespace py = pybind11;

void init_ca_ncsbuild(py::module &m) {
  py::class_<Ca_ncsbuild>(m, "Ca_ncsbuild")
      .def(py::init<double, double, int>(), py::arg("reliability") = 0.5,
           py::arg("rmsd") = 1.0, py::arg("nmin") = 12,
           "Constructor for Ca_ncsbuild class.")
      .def("__call__", &Ca_ncsbuild::operator(), py::arg("mol"),
           py::arg("xmap"), py::arg("llktarget"), py::arg("seq"),
           "Build NCS related chains.")
      .def("__repr__",
           [](const Ca_ncsbuild &self) {
             return "<buccaneer.Ca_ncsbuild class.>";
           })
      .doc() = "Class for build NCS related chains using density.";
}