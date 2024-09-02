// Wrapper for buccaneer-correct
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-correct.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"

namespace py = pybind11;

void init_ca_correct(py::module &m) {
  py::class_<Ca_correct>(m, "Ca_correct")
      .def(py::init<int>(), py::arg("torsion_sampling") = 12,
           "Constructor for Ca correct class.")
      .def("__call__", &Ca_correct::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktargets"), py::arg("seq"),
           "Rebuild chain to fix insertion/deletions.")
      .def_property_readonly("num_corrected", &Ca_correct::num_corrected,
                             "Get number of C-alphas corrected.")
      .def("__repr__",
           [](const Ca_correct &self) {
             std::stringstream stream;
             stream << "<buccaneer.Ca_correct with ";
             stream << self.num_corrected() << " C-alphas corrected.>";
             return stream.str();
           })
      .doc() = "Class for correcting Ca chains using density.";
}