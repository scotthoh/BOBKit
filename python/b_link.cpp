// Nanobind bindings for buccaneer-lib
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-link.h"
#include "commons.h"
#include <nanobind/operators.h>

using namespace clipper;

void add_ca_link(nb::module_ &m) {
  nb::class_<Ca_link>(m, "Ca_link")
      .def(nb::init<ftype, int>(), nb::arg("rad_link") = 5.0,
           nb::arg("torsion_sampling") = 24, "Constructor for Ca_link class.")
      .def("__call__", &Ca_link::operator(), nb::arg("mol"), nb::arg("xmap"),
           nb::arg("llktarget"), "Merge overlapped Ca chains.")
      .def_prop_ro("num_linked", &Ca_link::num_linked,
                             "Get number of Ca alphas linked.")
      .def("__repr__",
           [](const Ca_link &self) {
             return "<buccaneer.Ca_link with " + String(self.num_linked()) +
                    " C-alphas linked.>";
           })
      .doc() =
      "Class for merging overlapped Ca chains and grouping by symmetry.";
}