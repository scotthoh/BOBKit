// Wrapper for clipper coord's atom, atomlist
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <gemmi/math.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

void init_atomlist(py::module &m) {
  py::class_<Atom> atom(m, "Atom");
  atom.def(py::init<>())
      .def(py::init<const Atom &>(), "Constructor from atom-like object.")
      .def_property("element", &Atom::element, &Atom::set_element,
                    "Get/set atom element name: e.g. \'C\', \'N\', \'Zn2+\'.")
      .def_property("pos", &Atom::coord_orth, &Atom::set_coord_orth,
                    "Get/set atom orthogonal (Angstrom) coordinate.")
      .def_property("occupancy", &Atom::occupancy, &Atom::set_occupancy,
                    "Get/set atom occupancy.")
      .def_property("u_iso", &Atom::u_iso, &Atom::set_u_iso,
                    "Get/set atom orthogonal isotropic U value.")
      .def_property(
          "u_aniso_orth", &Atom::u_aniso_orth,
          [](Atom &self, py::array_t<ftype> u) {
            std::vector<ftype> uval(6);
            fill_array_1d<std::vector<ftype>, ftype>(uval, 6, u);
            auto uaniso = U_aniso_orth(uval[0], uval[1], uval[2], uval[3],
                                       uval[4], uval[5]);
            self.set_u_aniso_orth(uaniso);
          },
          "Get set atom orthogonal anisotropic U value.")
      .def_property(
          "b_iso",
          [](const Atom &self) { return clipper::Util::u2b(self.u_iso()); },
          [](Atom &self, const ftype &value) {
            self.set_u_iso(clipper::Util::b2u(value));
          },
          "Get/set atom isotropic B value.")
      .def(
          "x", [](const Atom &self) { return self.coord_orth().x(); },
          "Get atom x coordinate.")
      .def(
          "y", [](const Atom &self) { return self.coord_orth().y(); },
          "Get atom y coordinate.")
      .def(
          "z", [](const Atom &self) { return self.coord_orth().z(); },
          "Get atom z coordinate.")
      .def("transform", &Atom::transform, py::arg("rtop"),
           "Apply a rotation-translation operator (RTop) to the atom.")
      .def("is_null", &Atom::is_null,
           "Test for null atom: atom is null if coord is null.")
      .def_static("null", &Atom::null, "Return null atom.")
      .def(
          "copy", [](const Atom &self) { return self; },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .def("__repr__",
           [](const Atom &self) {
             return "<clipper.Atom " + self.element().trim() + " " +
                    self.coord_orth().format() + ">";
           })
      .doc() = "Atom class.\nThis class defines a minimal atom object "
               "providing only those properties required for an electron "
               "density calculation. A template constructor allows it to be "
               "constructed from any other object with appropriate properties.";

  py::class_<Atom_list> atomlist(m, "Atom_list");
  atomlist.def(py::init<>())
      .def(py::init<const std::vector<Atom> &>(), "Constructor from list[Atom]")
      //.def(py::init<const T &>())
      .def(
          "clear", [](Atom_list &self) { self.clear(); }, "Clear Atom_list.")
      .def(
          "pop_back", [](Atom_list &self) { self.pop_back(); },
          "Remove last atom in list. No data is returned. If data is needed, "
          "it should be retrieved before calling pop_back.")
      .def(
          "push_back",
          [](Atom_list &self, Atom &atom) { self.push_back(atom); },
          "Add atom to the end of the list.")
      .def("__repr__",
           [](const Atom_list &self) {
             std::stringstream stream;
             stream << "<clipper.Atom_list containing " << self.size()
                    << " atom(s).>";
             return stream.str();
           })
      .def("__len__", [](const Atom_list &self) { return self.size(); })
      .def(
          "size", [](const Atom_list &self) { return self.size(); },
          "Return size of the list.")
      .def(
          "__getitem__",
          [](const Atom_list &self, const int &i) -> Atom {
            // return self.at(i);
            return self.at(normalise_index(i, self.size()));
          },
          "Returns a copy of Atom at given index.")
      .def(
          "__setitem__",
          [](Atom_list &self, const int &i, const Atom &atom) {
            self.at(normalise_index(i, self.size())) = atom;
          },
          "Replaces Atom at given index.")
      .def(
          "delete_atom",
          [](Atom_list &self, const int &pos) { delitem_at_index(self, pos); },
          py::arg("index"), "Delete atom at given index.")
      .def(
          "delete_atoms",
          [](Atom_list &self, py::slice slice) { delitem_slice(self, slice); },
          py::arg("slice"), "Delete atoms in range from slice object.")
      .def(
          "insert_atom",
          [](Atom_list &self, const Atom &atom, const int &pos) {
            add_item(self, atom, pos);
          },
          py::arg("index"), py::arg("Atom"),
          "Insert atom to the position of given index.")
      .def(
          "__iter__",
          [](Atom_list &self) {
            return py::make_iterator(&self[0], &self[self.size()]);
          },
          py::keep_alive<0, 1>())
      .def(
          "insert_list",
          [](Atom_list &self, const std::vector<Atom> &a, const int &pos) {
            if (pos == -1)
              self.insert(self.end(), a.begin(), a.end());
            else
              self.insert(self.begin() + pos, a.begin(), a.end());
          },
          py::arg("list"), py::arg("pos") = -1,
          "Concatenate another Atom_list at the given position of this "
          "Atom_list. Default is at the end.")
      .def(
          "copy", [](const Atom_list &self) { return self; },
          "Return a copy of the object.")
      .doc() = "Atom list class.\nThis class defines a minimal atom list "
               "object providing only those properties required for an "
               "electron density calculation. It is a trivial derivation "
               "from std::vector<Atom>. In addition a template constructor "
               "allows it to be constructed from any other object with "
               "appropriate properties.";
}
