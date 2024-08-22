// Wrapper for buccaneer-link
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-link.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void init_ca_link(py::module &m) {
  py::class_<Ca_link>(m, "Ca_link")
      .def(py::init<ftype, int>(), py::arg("rad_link") = 5.0,
           py::arg("torsion_sampling") = 24, "Constructor for Ca_link class.")
      .def("__call__", &Ca_link::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktarget"), "Merge overlapped Ca chains.")
      .def_property_readonly("num_linked", &Ca_link::num_linked,
                             "Get number of Ca alphas linked.")
      .def("__repr__",
           [](const Ca_link &self) {
             return "<buccaneer.Ca_link with " + String(self.num_linked()) +
                    " C-alphas linked.>";
           })
      .doc() =
      "Class for merging overlapped Ca chains and grouping by symmetry.";
}