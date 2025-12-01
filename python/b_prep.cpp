// Nanobind bindings for buccaneer-prep
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-lib.h"
#include "buccaneer/buccaneer-prep.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

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
    max_ = llkcls_.size();
  };
  //! Constructor: takes int size of vector to initialise
  LLK_TargetList(const int &num) { init(num); };
  //! Initialiser: takes int size of vector to initialise
  void init(const int &n) {
    max_ = n;
    llkcls_.clear();
    llkcls_.resize(max_);
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
  int max_;
};

void declare_llktargetlist(nb::module_ &m) {

  nb::class_<LLK_TargetList>(m, "LLK_TargetList")
      .def(nb::init<>())
      .def(nb::init<const std::vector<LLK_map_target> &>(),
           "Constructor from a list of LLK_map_targets.")
      .def(nb::init<const int &>(),
           "Constructor with maximum size for the list.")
      .def("init", &LLK_TargetList::init,
           "Initialiser with maximum size for the list.")
      .def("__len__", &LLK_TargetList::size)
      .def("size", &LLK_TargetList::size, "Get size of the list.")
      .def(
          "__getitem__",
          [](const LLK_TargetList &self, const int &i) { return self[i]; },
          nb::rv_policy::reference_internal,
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

void declare_ca_prep(nb::module_ &m) {
  nb::class_<Ca_prep> ca_prep(m, "Ca_prep");

  using Class = Ca_prep::Rama_flt;
  nb::class_<Class>(ca_prep, "Rama_flt")
      .def(nb::init<>())
      .def(nb::init<double, double, double>(), nb::arg("phi"), nb::arg("psi"),
           nb::arg("rad"),
           "Construct Ramachandran filter with phi, psi, radius.")
      .def_rw("phi", &Class::phi, "Accessor for phi.")
      .def_rw("psi", &Class::psi, "Accessor for psi.")
      .def_rw("rad", &Class::rad, "Accessor for radius.")
      .def("__repr__",
           [](const Class &self) {
             return "<buccaneer.Ca_prep.Rama_flt phi = " +
                    clipper::String(self.phi, 6, 6) +
                    ", psi = " + clipper::String(self.psi, 6, 6) +
                    ", rad = " + clipper::String(self.rad, 6, 6) + ">";
           })
      .doc() = "Struct for Ramachandran filter data.";

  ca_prep
      .def(nb::init<double, double, Ca_prep::Rama_flt, bool, bool, bool>(),
           nb::arg("main_tgt_rad"), nb::arg("side_tgt_rad"),
           nb::arg("rama_flt"), nb::arg("correl"), nb::arg("seqnc"),
           nb::arg("debug") = false, "Constructor for Ca_prep class.")
      .def(
          "__call__",
          [](Ca_prep &self, LLK_map_target &llktgt, LLK_TargetList &llkcls,
             const clipper::MiniMol &mol, const clipper::Xmap<float> &xmap) {
            self(llktgt, llkcls.get_vector(), mol, xmap);
          },
          "Prepare LLK targets.")
      .def_prop_ro_static(
          "rama_flt_all", [](nb::object) { return Ca_prep::rama_flt_all; },
          "Ramachandran filter data.")
      .def_prop_ro_static(
          "rama_flt_helix", [](nb::object) { return Ca_prep::rama_flt_helix; },
          "Ramachandran filter data.")
      .def_prop_ro_static(
          "rama_flt_strand",
          [](nb::object) { return Ca_prep::rama_flt_strand; },
          "Ramachandran filter data.")
      .def_prop_ro_static(
          "rama_flt_nonhelix",
          [](nb::object) { return Ca_prep::rama_flt_nonhelix; },
          "Ramachandran filter data.")
      .def_static("set_cpus", &Ca_prep::set_cpus,
                  "Set number of threads to use.")
      .def("__repr__",
           [](const Ca_prep &self) { return "<buccaneer.Ca_prep class>"; })
      .doc() = "Class to prepare log likelihood targets.";

  nb::class_<Prep_threaded>(m, "Prep_threaded")
      .def(nb::init<>())
      .def(nb::init<std::vector<LLK_map_target> &, const clipper::Xmap<float> &,
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

void add_ca_prep(nb::module_ &m) {
  declare_llktargetlist(m);
  declare_ca_prep(m);
}