// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prep.h"

#include "helper_functions.h"
#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// ssclass LLK_map_targetList : public std::vector<LLK_map_target> {
// sspublic:
// ss   LLK_map_targetList() {}
// ss};
// ss
// ssvoid declare_llkmaptgtlist(py::module &m) {
// ss   py::class_<LLK_map_targetList>(m, "LLK_map_targetList")
// ss       .def(py::init<>())
// ss       .def("__len__",
// ss            [](const LLK_map_targetList &self) { return self.size(); })
// ss       .def("__getitem__",
// ss            [](const LLK_map_targetList &self, const int &i) ->
// LLK_map_target
//            {
// ss              return self.at(i);
// ss            })
// ss       .def("__setitem__",
// ss            [](LLK_map_targetList &self, const int &i,
// ss               const LLK_map_target &llktgt) { self.at(i) = llktgt; })
// ss       .def(
// ss           "__iter__",
// ss           [](LLK_map_targetList &self) {
// ss             return py::make_iterator(&self[0], &self[self.size()]);
// ss           },
// ss           py::keep_alive<0, 1>());
// ss}

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
      //.def("__call__", &Ca_prep::operator())
      .def("check_params", &Ca_prep::check_params)
      .def("__call__",
           [](Ca_prep &self, LLK_map_target &llktgt,
              std::vector<LLK_map_target> &llkcls, const clipper::MiniMol &mol,
              const clipper::Xmap<float> &xmap) {
             // std::vector<LLK_map_target> v;
             std::vector<LLK_map_target> llkclsa(20);
             // v = llkcls;
             //  if (llkcls.size() != 0) {
             //  for (int i = 0; i < llkcls.size(); i++)
             //  v->at(i) = llkcls[i];
             //  }

             // auto buf = llkcls.request();
             // auto ptr = (LLK_map_target *)buf.ptr;
             // for (int i = 0; i < n; ++i)
             //   v[i] = ptr[i];
             //  fill_array_1d<std::vector<LLK_map_target>, LLK_map_target>(v,
             //  0,
             //                                                             llkcls);
             self(llktgt, llkclsa, mol, xmap);
             return llkclsa;
           })
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
