// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-lib.h>
#include <buccaneer/buccaneer-prep.h>

#include "helper_functions.h"
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//! Makeup class to hold a vector of LLK_map_target used in cbuccaneer.
/*! Simple hack to pass vectors as reference to C++ side and gets updated.
 */
class LLK_TargetList {
public:
  //! Null constructor
  LLK_TargetList() {}
  //! Constructor: takes std::vector<LLK_map_target>
  LLK_TargetList(const std::vector<LLK_map_target> &l) {
    llkcls_ = l;
    max = llkcls_.size();
  };
  //! Constructor: takes int size of vector to initialize
  LLK_TargetList(const int &num) { init(num); };
  void init(const int &n) {
    max = n;
    llkcls_.clear();
    llkcls_.resize(max);
  }
  // access methods to LLK_map_target stored in vector
  //! get LLK_map_target by index
  const LLK_map_target &operator[](const int &i) const { return llkcls_[i]; }
  //! set LLK_map_target by index
  LLK_map_target &operator[](const int &i) { return llkcls_[i]; }
  //! get vector<LLK_map_target>
  std::vector<LLK_map_target> &get_vector() { return llkcls_; }
  //! get size
  int size() const { return llkcls_.size(); }

private:
  // members
  std::vector<LLK_map_target> llkcls_;
  int max;
};

void declare_llktargetlist(py::module &m) {

  py::class_<LLK_TargetList>(m, "LLK_TargetList")
      .def(py::init<>())
      .def(py::init<const std::vector<LLK_map_target> &>())
      .def(py::init<const int &>())
      .def("init", &LLK_TargetList::init)
      .def("__len__", &LLK_TargetList::size)
      .def("size", &LLK_TargetList::size)
      .def(
          "__getitem__",
          [](const LLK_TargetList &self, const int &i) { return self[i]; },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](LLK_TargetList &self, const int &i,
              const LLK_map_target &llktgt) { self[i] = llktgt; })
      .def("get_vector", &LLK_TargetList::get_vector)
      .def("__repr__", [](const LLK_TargetList &self) {
        return "< List containing " + clipper::String(self.size()) +
               " LLK_map_target(s).>";
      });
}

void declare_ca_prep(py::module &m) {
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
      //.def("check_params", &Ca_prep::check_params)
      .def("__call__",
           [](Ca_prep &self, LLK_map_target &llktgt, LLK_TargetList &llkcls,
              const clipper::MiniMol &mol, const clipper::Xmap<float> &xmap) {
             self(llktgt, llkcls.get_vector(), mol, xmap);
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

void init_ca_prep(py::module &m) {
  declare_llktargetlist(m);
  declare_ca_prep(m);
}