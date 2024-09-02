// Wrapper for buccaneer-grow
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-grow.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_ca_grow(py::module &m) {
  py::class_<Ca_grow>(m, "Ca_grow")
      .def(py::init<int>(), py::arg("n_grow") = 25,
           "Constructor for Ca_grow class.")
      .def("__call__", &Ca_grow::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktarget"), "Grow the chain given a map and target.")
      .def_static("grow", &Ca_grow::grow, py::arg("chain"), py::arg("xmap"),
                  py::arg("llktarget"), py::arg("rama1"), py::arg("rama2"),
                  py::arg("cutoff"), py::arg("ngrow"),
                  "Grow a chain at both ends.")
      .def_static("next_ca_group", &Ca_grow::next_ca_group, py::arg("chain"),
                  py::arg("xmap"), py::arg("llktarget"), py::arg("rama1"),
                  py::arg("rama2"),
                  "Add a new Ca-group to the C-terminus of a chain by best fit "
                  "to density.")
      .def_static("prev_ca_group", &Ca_grow::prev_ca_group, py::arg("chain"),
                  py::arg("xmap"), py::arg("llktarget"), py::arg("rama1"),
                  py::arg("rama2"),
                  "Add a new Ca-group to the N-terminus of a chain by best fit "
                  "to density.")
      .def_static("set_cpus", &Ca_grow::set_cpus, py::arg("ncpus"),
                  "Set number of cpu threads to use.")
      .def("__repr__",
           [](const Ca_grow &self) {
             std::stringstream stream;
             stream << "<buccaneer.Ca_grow class>";
             return stream.str();
           })
      .doc() = "Class for growing Ca chains using density.";
}

void declare_grow_threaded(py::module &m) {
  py::class_<Grow_threaded>(m, "Grow_threaded")
      .def(py::init<>())
      .def(py::init<const std::vector<Ca_chain> &, const Xmap<float> &,
                    const LLK_map_target &, const double &, const int &>(),
           py::arg("chains"), py::arg("xmap"), py::arg("llktarget"),
           py::arg("cutoff"), py::arg("n_grow"))
      .def("grow", &Grow_threaded::grow, py::arg("chn"), "Grow chain.")
      .def_property_readonly("result", &Grow_threaded::result,
                             "Accessor for list of chains.")
      .def("__call__", &Grow_threaded::operator(), py::arg("nthread") = 0,
           "Run single or multithreaded.")
      .def("merge", &Grow_threaded::merge, py::arg("other"),
           "Merge results from multiple threads.")
      // inherited function/property
      .def_property_readonly("id", &Grow_threaded::id, "Accessor for id.")
      .def("__repr__",
           [](const Grow_threaded &self) {
             std::stringstream stream;
             stream << "<buccaneer.Grow_threaded class, id=";
             stream << self.id() << ">";
             return stream.str();
           })
      .doc() = "Class for growing Ca groups.";
}

// Target_fn_refine_n_terminal_build defined in b_simplex.cpp
// to be within same scope as Target_fn_zero_order trampoline definition

void init_ca_grow(py::module &m) {
  declare_ca_grow(m);
  declare_grow_threaded(m);
}