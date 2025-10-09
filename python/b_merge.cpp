// Nanobind bindings for buccaneer-merge
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-merge.h"
#include "commons.h"
#include <nanobind/operators.h>

using namespace clipper;

void add_ca_merge(nb::module_ &m) {

  nb::class_<Ca_merge> camerge(m, "Ca_merge");
  camerge.def(nb::init<double>(), nb::arg("reliability") = 0.5)
      .def("__call__", &Ca_merge::operator(), nb::arg("mol"), nb::arg("xmap"),
           nb::arg("llktarget"), nb::arg("seq"), "Merge model.")
      .def_static("merge_mr", &Ca_merge::merge_mr, nb::arg("mol"),
                  nb::arg("mol_mr"), nb::arg("sigcut"), nb::arg("nseed"),
                  nb::arg("mr_filter"), nb::arg("mr_seed"),
                  "Merge with molecular replacement model.")
      .def("__repr__",
           [](const Ca_merge &self) { return "<buccaneer.Ca_merge class.>"; })
      .doc() = "Class for augmenting model with MERGE model.";
}