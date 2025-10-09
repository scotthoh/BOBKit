// Nanobind bindings for buccaneer-grow
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-grow.h"
#include "commons.h"
#include <nanobind/stl/vector.h>

using namespace clipper;

void declare_ca_grow(nb::module_ &m) {
  nb::class_<Ca_grow>(m, "Ca_grow")
      .def(nb::init<int>(), nb::arg("n_grow") = 25,
           "Constructor for Ca_grow class.")
      .def("__call__", &Ca_grow::operator(), nb::arg("mol"), nb::arg("xmap"),
           nb::arg("llktarget"), "Grow the chain given a map and target.")
      .def_static("grow", &Ca_grow::grow, nb::arg("chain"), nb::arg("xmap"),
                  nb::arg("llktarget"), nb::arg("rama1"), nb::arg("rama2"),
                  nb::arg("cutoff"), nb::arg("ngrow"),
                  "Grow a chain at both ends.")
      .def_static("next_ca_group", &Ca_grow::next_ca_group, nb::arg("chain"),
                  nb::arg("xmap"), nb::arg("llktarget"), nb::arg("rama1"),
                  nb::arg("rama2"),
                  "Add a new Ca-group to the C-terminus of a chain by best fit "
                  "to density.")
      .def_static("prev_ca_group", &Ca_grow::prev_ca_group, nb::arg("chain"),
                  nb::arg("xmap"), nb::arg("llktarget"), nb::arg("rama1"),
                  nb::arg("rama2"),
                  "Add a new Ca-group to the N-terminus of a chain by best fit "
                  "to density.")
      .def_static("set_cpus", &Ca_grow::set_cpus, nb::arg("ncpus"),
                  "Set number of cpu threads to use.")
      .def("__repr__",
           [](const Ca_grow &self) {
             return "<buccaneer.Ca_grow class>";
           })
      .doc() = "Class for growing Ca chains using density.";
}

void declare_grow_threaded(nb::module_ &m) {
  nb::class_<Grow_threaded>(m, "Grow_threaded")
      .def(nb::init<>())
      .def(nb::init<const std::vector<Ca_chain> &, const Xmap<float> &,
                    const LLK_map_target &, const double &, const int &>(),
           nb::arg("chains"), nb::arg("xmap"), nb::arg("llktarget"),
           nb::arg("cutoff"), nb::arg("n_grow"))
      .def("grow", &Grow_threaded::grow, nb::arg("chn"), "Grow chain.")
      .def_prop_ro("result", &Grow_threaded::result,
                             "Accessor for list of chains.")
      .def("__call__", &Grow_threaded::operator(), nb::arg("nthread") = 0,
           "Run single or multithreaded.")
      .def("merge", &Grow_threaded::merge, nb::arg("other"),
           "Merge results from multiple threads.")
      // inherited function/property
      .def_prop_ro("id", &Grow_threaded::id, "Accessor for id.")
      .def("__repr__",
           [](const Grow_threaded &self) {
             return "<buccaneer.Grow_threaded class, id=" + String(self.id()) + ">";
           })
      .doc() = "Class for growing Ca groups.";
}

// Target_fn_refine_n_terminal_build defined in b_simplex.cpp
// to be within same scope as Target_fn_zero_order trampoline definition

void add_ca_grow(nb::module_ &m) {
  declare_ca_grow(m);
  declare_grow_threaded(m);
}