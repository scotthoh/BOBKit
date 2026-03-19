// Nanobind bindings for buccaneer-pro
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prot.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/make_iterator.h>


using namespace clipper;

void declare_ca_group(nb::module_ &m) {
  nb::class_<Ca_group>(m, "Ca_group")
      .def(nb::init<>())
      .def(nb::init<const Coord_orth &, const Coord_orth &,
                    const Coord_orth &>(),
           nb::arg("n"), nb::arg("ca"), nb::arg("c"),
           "Constructor from atom coordinates.")
      // constructor from list/array of coordinates N, CA, C
      .def("__init__", [](Ca_group *Cag, const std::array<ftype, 3> &n,
                       const std::array<ftype, 3> &ca,
                       const std::array<ftype, 3> &c) {
             Coord_orth N(n[0], n[1], n[2]);
             Coord_orth CA(ca[0], ca[1], ca[2]);
             Coord_orth C(c[0], c[1], c[2]);
             new ( Cag ) Ca_group(N, CA, C);
           },
           nb::arg("n"), nb::arg("ca"), nb::arg("c"),
           "Constructor from atom coordinates.")
      .def(nb::init<const MResidue &>(), "Initialise from monomer.")
      .def("is_null", &Ca_group::is_null, "Return true if group is null.")
      .def_prop_ro("coord_n", &Ca_group::coord_n,
                             "Get N atom coordinate.")
      .def_prop_ro("coord_ca", &Ca_group::coord_ca,
                             "Get C-alpha atom coordinate.")
      .def_prop_ro("coord_c", &Ca_group::coord_c,
                             "Get C atom coordinate.")
      .def_prop_ro("coord_cb", &Ca_group::coord_cb,
                             "Get C-beta atom coordinate.")
      .def("rtop_from_std_ori", &Ca_group::rtop_from_std_ori,
           "Get operator generating this group from standard orientation.")
      .def("rtop_beta_carbon", &Ca_group::rtop_beta_carbon,
           "Get operator centering C-beta standard orientation.")
      .def("next_ca_group", &Ca_group::next_ca_group, nb::arg("psi"),
           nb::arg("phi"), "Generate next Ca_group using Ramachandran angles.")
      .def("prev_ca_group", &Ca_group::prev_ca_group, nb::arg("phi"),
           nb::arg("psi"),
           "Generate previous Ca_group using Ramachandran angles")
      .def_static("std_coord_ca", &Ca_group::std_coord_ca,
                  "Return standard C-alpha atom coordinate")
      .def_static("std_coord_c", &Ca_group::std_coord_c,
                  "Return standard C atom coordinate")
      .def_static("std_coord_n", &Ca_group::std_coord_n,
                  "Return standard N atom coordinate")
      .def_static("std_coord_cb", &Ca_group::std_coord_cb,
                  "Return standard C-beta atom coordinate")
      .def_static("null", &Ca_group::null, "Return null Ca_group.")
      .def("__repr__",
           [](const Ca_group &self) {
             return "<buccaneer.Ca_group containing N,C-alpha,C atom "
                    "coordinates>";
           })
      .doc() =
      "C-alpha group.\n"
      "The Ca-group class represents a residue by the alpha Carbon and "
      "its neighbouring N and C main-chain atoms. It has methods to return "
      "the operator to generate the group from a standard orientation, and "
      "to generate the next or previous residue given two Ramachandran "
      "angles (one of this residue and one of the new one). ";
  // nb::bind_vector<std::deque<Ca_group>>(m, "DequeCagroup");
}

void declare_ca_chain(nb::module_ &m) {
  // Have to manually bind some of the deque member functions.
  nb::class_<Ca_chain>(m, "Ca_chain")
      .def(nb::init<>())
      // from buccaneer
      .def("ramachandran_phi", &Ca_chain::ramachandran_phi,
           "Return Ramachandran phi for any residue except first in a chain.")
      .def("ramachandran_psi", &Ca_chain::ramachandran_psi,
           "Return Ramachandran psi for any residue except last in a chain.")
      // defining functions from std::deque
      .def("__len__", [](const Ca_chain &self) { return self.size(); })
      .def(
          "size", [](const Ca_chain &self) { return self.size(); },
          "Get size of list.")
      .def(
          "empty", [](const Ca_chain &self) { return self.empty(); },
          "Return true if class is empty.")
      // element access
      .def(
          "__getitem__",
          [](Ca_chain &self, const int index) { return self.at(index); },
          "Get Ca_group at index.")
      .def(
          "__setitem__",
          [](Ca_chain &self, const int index, Ca_group &c) {
            self.at(index) = c;
          },
          "Set Ca_group at index.")
      // need to check if these are reference or copy
      .def(
          "front", [](Ca_chain &self) { return self.front(); },
          "Return front of list.")
      .def(
          "back", [](Ca_chain &self) { return self.back(); },
          "Return back of list.")
      // iterators
      .def(
          "__iter__",
          [](Ca_chain &self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<Ca_chain>(), "iterator", self.begin(), self.end());
          },
          nb::keep_alive<0, 1>())
      .def(
          "__reversed__",
          [](Ca_chain &self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<Ca_chain>(), "iterator", self.rbegin(), self.rend());
          },
          nb::keep_alive<0, 1>())
      // modifiers
      .def(
          "append", [](Ca_chain &self, Ca_group &c) { self.push_back(c); },
          "Append Ca_group to the end.")
      .def(
          "appendleft", [](Ca_chain &self, Ca_group &c) { self.push_front(c); },
          "Append Ca_group to the beginning.")
      .def(
          "remove",
          [](Ca_chain &self, const int &pos) { delitem_at_index(self, pos); },
          "Delete Ca_group at the given index.")
      .def(
          "remove",
          [](Ca_chain &self, const int &start, const int &end) {
            delitem_range(self, start, end);
          },
          "Delete Ca_group at range of indices given [start, end).")
      .def(
          "insert",
          [](Ca_chain &self, const int &pos, const Ca_group &c) {
            add_item(self, c, pos);
          },
          "Insert Ca_group at the given index.")
      .def(
          "pop",
          [](Ca_chain &self) {
            auto ret = self.back();
            self.pop_back();
            // delitem_at_index(self, pos);
            return ret;
          },
          "Remove the last Ca_group in the container.")
      .def(
          "popleft",
          [](Ca_chain &self) {
            auto ret = self.front();
            self.pop_front();
            return ret;
          },
          "Remove the first Ca_group in the container.")
      .def(
          "swap", [](Ca_chain &self, Ca_chain &other) { self.swap(other); },
          "Exchanges content of the container by the content of the given "
          "Ca_chain.")
      .def(
          "copy", [](const Ca_chain &self) { return self; },
          "Return a copy of the object.")
      .def("__repr__",
           [](const Ca_chain &self) {
             return "<buccaneer.Ca_chain with " + String(int(self.size())) +
                    " Ca_group(s).>";
           })
      .doc() =
      "Chain of Ca-groups.\n"
      "A Ca-chain is a std::deque (double ended queue) of Ca_group. In "
      "addition it has methods to return the Ramachandran angles of any "
      "residue.";
}

void declare_pr_group(nb::module_ &m) {
  nb::class_<Pr_group> pr_group(m, "Pr_group");

  nb::enum_<Pr_group::TYPE>(pr_group, "TYPE", "Atom types used in constructor.")
      .value("CaCN", Pr_group::TYPE::CaCN)
      .value("CaCO", Pr_group::TYPE::CaCO)
      .export_values();

  pr_group.def(nb::init<>())
      .def(nb::init<const Coord_orth &, const Coord_orth &, const Coord_orth &,
                    const Pr_group::TYPE &>(),
           nb::arg("ca"), nb::arg("c"), nb::arg("other"), nb::arg("type"),
           "Constructor from atom coordinates (Ca, C, N[+1] or Ca, C, O).")
      // constructor from list/array of coordinates
      .def("__init__", [](Pr_group *Prg, const std::array<ftype, 3> &ca,
                       const std::array<ftype, 3> &c,
                       const std::array<ftype, 3> &thirdatom,
                       const Pr_group::TYPE &type) {
             Coord_orth CA(ca[0], ca[1], ca[2]);
             Coord_orth C(c[0], c[1], c[2]);
             Coord_orth other(thirdatom[0], thirdatom[1], thirdatom[2]);
             new (Prg) Pr_group(CA, C, other, type);
           },
           nb::arg("ca"), nb::arg("c"), nb::arg("other"), nb::arg("type"),
           "Constructor from atom coordinates (Ca, C, N[+1] or Ca, C, O).")
      .def_prop_ro("coord_ca", &Pr_group::coord_ca,
                             "Get C-alpha atom coordinate.")
      .def_prop_ro("coord_c", &Pr_group::coord_c,
                             "Get C atom coordinate.")
      .def_prop_ro("coord_n_next", &Pr_group::coord_n_next,
                             "Get next N atom coordinate.")
      .def("coord_o", &Pr_group::coord_o, "Generate O atom coordinate.")
      .def("coord_ca_next", &Pr_group::coord_ca_next,
           "Get next C-alpha atom coordinate.")
      .def("rtop_from_std_ori", &Pr_group::rtop_from_std_ori,
           "Get operator generating this group from standard orientation.")
      .def("next_pr_group", &Pr_group::next_pr_group, nb::arg("phi"),
           nb::arg("psi"), "Generate next Pr_group using Ramachandran angles.")
      .def("prev_pr_group", &Pr_group::prev_pr_group, nb::arg("psi"),
           nb::arg("phi"),
           "Generate previous Pr_group using Ramachandran angles.")
      .def("__repr__",
           [](const Pr_group &self) {
             return "<buccaneer.Pr_group, Planar-residue-group class.>";
           })
      .doc() =
      "Planar residue group.\n"
      "The Planar-residue-group class represents the planar atoms "
      "surrounding a peptide bond. It is represented by the alpha Carbon, C "
      "atom, and N atom of the next residue. It may also be described within a "
      "single residue by referring to the O atom instead of the next N.  It "
      "has methods to return the operator to generate the group from a "
      "standard orientation, and to generate the next or previous residue "
      "given two Ramachandran angles.";
}

void add_buccaneer_prot(nb::module_ &m) {
  declare_ca_group(m);
  declare_ca_chain(m);
  declare_pr_group(m);
}