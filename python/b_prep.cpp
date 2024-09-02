// Wrapper for buccaneer-prep
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-lib.h"
#include "buccaneer/buccaneer-prep.h"
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"

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
  //! Constructor: takes int size of vector to initialise
  LLK_TargetList(const int &num) { init(num); };
  //! Initiliaser: takes int size of vector to initilise
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
      .def(py::init<const std::vector<LLK_map_target> &>(),
           "Constructor from a list of LLK_map_targets.")
      .def(py::init<const int &>(),
           "Constructor with maximum size for the list.")
      .def("init", &LLK_TargetList::init,
           "Initialiser with maximum size for the list.")
      .def("__len__", &LLK_TargetList::size)
      .def("size", &LLK_TargetList::size, "Get size of the list.")
      .def(
          "__getitem__",
          [](const LLK_TargetList &self, const int &i) { return self[i]; },
          py::return_value_policy::reference_internal,
          "Get LLK_map_target at given index.")
      .def(
          "__setitem__",
          [](LLK_TargetList &self, const int &i, const LLK_map_target &llktgt) {
            self[i] = llktgt;
          },
          "Set LLK_map_target at given index.")
      .def("get_vector", &LLK_TargetList::get_vector,
           "Return the list of LLK_map_targets.")
      .def("__repr__",
           [](const LLK_TargetList &self) {
             return "< List containing " + clipper::String(self.size()) +
                    " LLK_map_target(s).>";
           })
      .doc() = "Class containing a list of LLK_map_target classes.";
}

void declare_ca_prep(py::module &m) {
  py::class_<Ca_prep> ca_prep(m, "Ca_prep");

  using Class = Ca_prep::Rama_flt;
  py::class_<Class>(ca_prep, "Rama_flt")
      .def(py::init<>())
      .def(py::init<double, double, double>(), py::arg("phi"), py::arg("psi"),
           py::arg("rad"),
           "Construct Ramachandran filter with phi, psi, radius.")
      .def_readwrite("phi", &Class::phi, "Accessor for phi.")
      .def_readwrite("psi", &Class::psi, "Accessor for psi.")
      .def_readwrite("rad", &Class::rad, "Accessor for radius.")
      .def("__repr__",
           [](const Class &self) {
             return "<buccaneer.Ca_prep.Rama_flt phi = " +
                    clipper::String(self.phi, 6, 6) +
                    ", psi = " + clipper::String(self.psi, 6, 6) +
                    ", rad = " + clipper::String(self.rad, 6, 6) + ">";
           })
      .doc() = "Struct for Ramachandran filter data.";

  ca_prep
      .def(py::init<double, double, Ca_prep::Rama_flt, bool, bool, bool>(),
           py::arg("main_tgt_rad"), py::arg("side_tgt_rad"),
           py::arg("rama_flt"), py::arg("correl"), py::arg("seqnc"),
           py::arg("debug") = false, "Constructor for Ca_prep class.")
      .def(
          "__call__",
          [](Ca_prep &self, LLK_map_target &llktgt, LLK_TargetList &llkcls,
             const clipper::MiniMol &mol, const clipper::Xmap<float> &xmap) {
            self(llktgt, llkcls.get_vector(), mol, xmap);
          },
          "Prepare LLK targets.")
      .def_property_readonly_static(
          "rama_flt_all", [](py::object) { return Ca_prep::rama_flt_all; },
          "Ramachandran filter data.")
      .def_property_readonly_static(
          "rama_flt_helix", [](py::object) { return Ca_prep::rama_flt_helix; },
          "Ramachandran filter data.")
      .def_property_readonly_static(
          "rama_flt_strand",
          [](py::object) { return Ca_prep::rama_flt_strand; },
          "Ramachandran filter data.")
      .def_property_readonly_static(
          "rama_flt_nonhelix",
          [](py::object) { return Ca_prep::rama_flt_nonhelix; },
          "Ramachandran filter data.")
      .def_static("set_cpus", &Ca_prep::set_cpus,
                  "Set number of threads to use.")
      .def("__repr__",
           [](const Ca_prep &self) { return "<buccaneer.Ca_prep class>"; })
      .doc() = "Class to prepare log likelihood targets.";

  py::class_<Prep_threaded>(m, "Prep_threaded")
      .def(py::init<>())
      .def(py::init<std::vector<LLK_map_target> &, const clipper::Xmap<float> &,
                    const std::vector<std::vector<clipper::RTop_orth>> &>(),
           "Constructor for Prep_threaded class.")
      .def("prep", &Prep_threaded::prep, "Prepare log likelihood targets.")
      .def("__call__", &Prep_threaded::operator(),
           "Run single or multi-threaded.")
      .def("merge", &Prep_threaded::merge,
           "Merge results from mutiple threads.")
      .def("__repr__",
           [](const Prep_threaded &self) {
             return "<buccaneer.Prep_threaded class>";
           })
      .doc() = "Class to prepare log likelihood targets.";
}

void init_ca_prep(py::module &m) {
  declare_llktargetlist(m);
  declare_ca_prep(m);
}