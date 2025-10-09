// Nanobind bindings for buccaneer-filter
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-filter.h"
#include "commons.h"
#include <nanobind/operators.h>


using namespace clipper;

void add_ca_filter(nb::module_ &m) {
  // "Forward declaration" of python class to avoid C++ signatures in docstrings
  nb::class_<Ca_filter> cafilter(m, "Ca_filter");
  cafilter
      .def(nb::init<double>(), nb::arg("sig_cut") = 3.0,
           "Constructor for Ca filter class.")
      .def("__call__", &Ca_filter::operator(), nb::arg("mol"), nb::arg("xmap"),
           "Merge overlapped Ca chains and grouping by symmetry.")
      .def_static(
          "filter",
          nb::overload_cast<MiniMol &, const Xmap<float> &, double, bool>(&Ca_filter::filter),
          nb::arg("mol"), nb::arg("xmap"), nb::arg("sigcut"),
          nb::arg("keep") = true)
      .def_static("filter",
                  nb::overload_cast<MiniMol &, double>(&Ca_filter::filter),
                  nb::arg("mol"), nb::arg("sigcut"))
      .def("__repr__",
           [](const Ca_filter &self) {
             return "<buccaneer.Ca_filter class.>";
           })
      .doc() =
      "Class for merging overlapped Ca chains and grouping by symmetry.";
}