// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prep.h"

#include "helper_functions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void init_ca_prep(py::module &m) {
  py::class_<Ca_prep> ca_prep(m, "Ca_prep");

  using Class = Ca_prep::Rama_flt;
  py::class_<Class>(ca_prep, "Rama_flt")
      .def(py::init<>())
      .def(py::init<double, double, double>(), py::arg("phi"), py::arg("psi"),
           py::arg("rad"))
      .def_readwrite("phi", &Class::phi)
      .def_readwrite("psi", &Class::psi)
      .def_readwrite("rad", &Class::rad)
      .def("__repr__", [](const Class &self) {
        return "<buccaneer.Ca_prep.Rama_flt phi = " +
               clipper::String(self.phi, 6, 6) +
               ", psi = " + clipper::String(self.psi, 6, 6) +
               ", rad = " + clipper::String(self.rad, 6, 6) + ">";
      });

  ca_prep
      .def(py::init<double, double, Ca_prep::Rama_flt, bool, bool, bool>(),
           py::arg("main_tgt_rad"), py::arg("side_tgt_rad"),
           py::arg("rama_flt"), py::arg("correl"), py::arg("seqnc"),
           py::arg("debug") = false)
      .def("__call__", &Ca_prep::operator())
      .def_property_readonly_static(
          "rama_flt_all", [](py::object) { return Ca_prep::rama_flt_all; })
      .def_property_readonly_static(
          "rama_flt_helix", [](py::object) { return Ca_prep::rama_flt_helix; })
      .def_property_readonly_static(
          "rama_flt_strand",
          [](py::object) { return Ca_prep::rama_flt_strand; })
      .def_property_readonly_static(
          "rama_flt_nonhelix",
          [](py::object) { return Ca_prep::rama_flt_nonhelix; })
      .def_static("set_cpus", &Ca_prep::set_cpus)
      .def("__repr__",
           [](const Ca_prep &self) { return "<buccaneer.Ca_prep class>"; });

  py::class_<Prep_threaded>(m, "Prep_threaded")
      .def(py::init<>())
      .def(py::init<std::vector<LLK_map_target> &, const clipper::Xmap<float> &,
                    const std::vector<std::vector<clipper::RTop_orth>> &>())
      .def("prep", &Prep_threaded::prep)
      .def("__call__", &Prep_threaded::operator())
      .def("merge", &Prep_threaded::merge)
      .def("__repr__", [](const Prep_threaded &self) {
        return "<buccaneer.Prep_threaded class>";
      });
}
