// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer/buccaneer-correct.h"

#include "type_conversions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_ca_correct(py::module &m) {
  py::class_<Ca_correct>(m, "Ca_correct")
      .def(py::init<int>(), py::arg("torsion_sampling") = 12)
      .def("__call__", &Ca_correct::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktarget"), py::arg("seq"))
      .def_property_readonly("num_corrected", &Ca_correct::num_corrected)
      .def("__repr__", [](const Ca_correct &self) {
        std::stringstream stream;
        stream << "<buccaneer.Ca_correct with ";
        stream << self.num_corrected() << " C-alphas corrected.>";
        return stream.str();
      });
}