// Nanobind bindings for buccaneer-correct
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-correct.h"
#include "commons.h"
#include <nanobind/operators.h>

using namespace clipper;

void add_ca_correct(nb::module_ &m)
{
  nb::class_<Ca_correct>(m, "Ca_correct")
      .def(nb::init<int>(), nb::arg("torsion_sampling") = 12,
           "Constructor for Ca correct class.")
      .def("__call__", &Ca_correct::operator(), nb::arg("mol"), nb::arg("xmap"),
           nb::arg("llktargets"), nb::arg("seq"),
           "Rebuild chain to fix insertion/deletions.")
      .def_prop_ro("num_corrected", &Ca_correct::num_corrected,
                   "Get number of C-alphas corrected.")
      .def("__repr__",
           [](const Ca_correct &self)
           {
             return "<buccaneer.Ca_correct with " + String(self.num_corrected()) + " C-alphas corrected.>";
           })
      .doc() = "Class for correcting Ca chains using density.";
}