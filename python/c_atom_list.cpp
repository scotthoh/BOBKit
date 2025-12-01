// Nanobind bindings for clipper atom and atomlist
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "c_atom_list.h"
#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <clipper/clipper-gemmi.h>
#include <gemmi/model.hpp>
using namespace clipper;

void add_atomlist( nb::module_ &m ) {
  nb::class_<Atom> atom( m, "Atom" );
  atom.def( nb::init<>() )
      .def( nb::init<const Atom &>(), "Constructor from atom-like object." )
      .def_prop_rw( "element", &Atom::element, &Atom::set_element,
                    nb::for_getter( "Get atom element name: e.g. \'C\', \'N\', \'Zn2+\'." ),
                    nb::for_setter( "Set atom element name: e.g. \'C\', \'N\', \'Zn2+\'." ) )
      .def_prop_rw( "pos", &Atom::coord_orth, &Atom::set_coord_orth, "Get/set atom orthogonal (Angstrom) coordinate." )
      .def_prop_rw( "occupancy", &Atom::occupancy, &Atom::set_occupancy, "Get/set atom occupancy." )
      .def_prop_rw( "u_iso", &Atom::u_iso, &Atom::set_u_iso, "Get/set atom orthogonal isotropic U value." )
      .def_prop_rw( "u_aniso_orth", &Atom::u_aniso_orth, &Atom::set_u_aniso_orth, "Get/set atom anisotropic U values." )
      //.def_prop_rw( "u_aniso_orth", &Atom::u_aniso_orth, [](Atom &self, nb::array_t<ftype> u) {
      //              std::vector<ftype> uval(6);
      //              fill_array_1d<std::vector<ftype>, ftype>(uval, 6, u);
      //              auto uaniso = U_aniso_orth(uval[0], uval[1], uval[2], uval[3], uval[4], uval[5]);
      //              self.set_u_aniso_orth(uaniso); },
      //    "Get set atom orthogonal anisotropic U value.")
      .def_prop_rw(
          "b_iso", []( const Atom &self ) { return clipper::Util::u2b( self.u_iso() ); },
          []( Atom &self, const ftype &value ) { self.set_u_iso( clipper::Util::b2u( value ) ); },
          nb::for_getter( nb::sig( "def b_iso(self, /) -> float" ) ),
          nb::for_setter( nb::sig( "def b_iso(self, value: float, /) -> None" ) ),
          nb::for_getter( "Get atom isotropic B value." ), nb::for_setter( "Set atom isotropic B value." ) )
      .def_prop_ro(
          "x", []( const Atom &self ) { return self.coord_orth().x(); }, "Get atom x coordinate." )
      .def_prop_ro(
          "y", []( const Atom &self ) { return self.coord_orth().y(); }, "Get atom y coordinate." )
      .def_prop_ro(
          "z", []( const Atom &self ) { return self.coord_orth().z(); }, "Get atom z coordinate." )
      .def( "transform", &Atom::transform, nb::arg( "rtop" ),
            "Apply a rotation-translation operator (RTop) to the atom." )
      .def( "is_null", &Atom::is_null, "Test for null atom: atom is null if coord is null." )
      .def_static( "null", &Atom::null, "Return null atom." )
      .def(
          "clone", []( const Atom &self ) { return new Atom( self ); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .def( "__repr__",
            []( const Atom &self ) {
              return "<clipper.Atom " + self.element().trim() + " " + self.coord_orth().format() + ">";
            } )
      .def( "__getstate__", [](const Atom &a) {
        return nb::make_tuple(a.element(), a.coord_orth(), a.occupancy(), a.u_iso(), a.u_aniso_orth());
      })
      .def( "__setstate__", [](Atom &a, const nb::tuple &t) {
        if ( t.size() != 5 ) throw std::runtime_error("Invalide size!");
        new (&a) Atom();
        a.set_element(nb::cast<String>(t[0]));
        a.set_coord_orth(nb::cast<Coord_orth>(t[1]));
        a.set_occupancy(nb::cast<ftype>(t[2]));
        a.set_u_iso(nb::cast<ftype>(t[3]));
        a.set_u_aniso_orth(nb::cast<U_aniso_orth>(t[4]));
        //new (&a) Atom(atm);
      })
      .doc() = "Atom class.\nThis class defines a minimal atom object "
               "providing only those properties required for an electron "
               "density calculation. A template constructor allows it to be "
               "constructed from any other object with appropriate properties.";

  //nb::bind_vector<Atom_list, rv_ri>( m, "Atom_list",
  //                                   "Vector of atoms, do not use this class, use Atom_list instead." )
  //    // doesn't work
  //    .def( "__repr__", [](const Atom_list &self) {
  //          return "<clipper.Atom_list containing " + std::to_string(self.size()) +
  //          " atom(s).>";
  //          }, nb::sig("def __repr__(self, /) -> str"))
  //    //q.def( "__getitem__", [](const Atom_list &self, size_t i) {
  //    //q  if ( i >= self.size() ) throw std::out_of_range("Index out of range");
  //    //q  return self[i];
  //    //q}, nb::rv_policy::reference_internal)
  //    //q.def( "__setitem__", [](Atom_list &self, size_t i, const Atom &atm) {
  //    //q  if ( i >= self.size() ) throw std::out_of_range("Index out of range");
  //    //q  self[i] = atm;
  //    //q})
  //    .def( "__len__", [](const Atom_list &self) {
  //      return self.size();
  //    })
  //    .def( "__getstate__", [](const Atom_list &a) {
  //      return nb::make_tuple(static_cast<const std::vector<Atom>&>(a));
  //    })
  //    .def( "__setstate__", [](Atom_list &a, nb::tuple state) {
  //      const auto &vec = nb::cast<std::vector<Atom>>(state[0]);
  //      a.clear();
  //      new (&a) Atom_list(vec);
  //      //a.insert(a.end(), vec.begin(), vec.end());
  //    } )
  //    .doc() = "Atom list class.\nThis class defines a minimal atom list "
  //             "object providing only those properties required for an "
  //             "electron density calculation. It is a trivial derivation "
  //             "from std::vector<Atom>. In addition a template constructor "
  //             "allows it to be constructed from any other object with "
  //             "appropriate properties.";
      
  nb::bind_vector<Atom_list, rv_ri>( m, "Atom_list")
  /*nb::class_<Atom_list, std::vector<Atom>> atomlist( m, "Atom_list" );
  atomlist.def( nb::init<>() )
      .def( nb::init<const std::vector<Atom> &>(), "Constructor from list[Atom]" )
      .def( "__repr__",
            []( const Atom_list &self ) {
              return "<clipper.Atom_list containing " + std::to_string( self.size() ) + " atom(s).>";
            } )
      //.def( "__getitem__", [](const Atom_list &self, int i) {
      //  //if ( i >= self.size() ) throw std::out_of_range("Index out of range");
      //  return self[normalise_index(i, self.size())];
      //  //return self[i];
      //}, nb::rv_policy::reference_internal)
      //.def( "__setitem__", [](Atom_list &self, int i, const Atom &atm) {
      //  //if ( i >= self.size() ) throw std::out_of_range("Index out of range");
      //  self[normalise_index(i, self.size())] = atm;
      //  //self[i] = atm;
      //})
      .def( "__len__", [](const Atom_list &self) {
        return self.size();
      })*/
      //.def( "__init__", [](Atom_list *alist, const std::vector<Atom>&v) {
      //  new (alist) Atom_list(v);
      //})
      .def_static( "as_atomlist", []( gemmi::Model &mdl ) {
        GemmiAtom_list ga(mdl.all());
        return new Atom_list(ga);
        //std::cout << mdl.chains.size() << std::endl;
      })
      .def( "__getstate__", [](const Atom_list &a) {
        return nb::make_tuple(static_cast<const std::vector<Atom>&>(a));
      })
      .def( "__setstate__", [](Atom_list &a, nb::tuple state) {
        // this works fine
        new (&a) Atom_list(nb::cast<Atom_list>(state[0]));
      } )
      .doc() = "Atom list class.\nThis class defines a minimal atom list "
               "object providing only those properties required for an "
               "electron density calculation. It is a trivial derivation "
               "from std::vector<Atom>. In addition a template constructor "
               "allows it to be constructed from any other object with "
               "appropriate properties."; //*/
  //.def( nb::init<>() )
  //    .def( nb::init<const std::vector<Atom> &>(), "Constructor from list[Atom]" )
  //.def(nb::init<const T &>())
  //.def(
  //    "clear", [](Atom_list &self) { self.clear(); }, "Clear Atom_list.")
  //.def(
  //    "pop_back", [](Atom_list &self) { self.pop_back(); },
  //    "Remove last atom in list. No data is returned. If data is needed, "
  //    "it should be retrieved before calling pop_back.")
  //.def(
  //    "push_back",
  //    [](Atom_list &self, Atom &atom) { self.push_back(atom); },
  //    "Add atom to the end of the list.")
  //.def( "__len__", [](const Atom_list &self) { return self.size(); } )
  //.def( "size", [](const Atom_list &self) { return self.size(); },
  //      "Return size of the list." )
  //.def(
  //    "__getitem__",
  //    [](const Atom_list &self, const int &i) -> Atom {
  //      // return self.at(i);
  //      return self.at(normalise_index(i, self.size()));
  //    },
  //    "Returns a copy of Atom at given index.")
  //.def(
  //    "__setitem__",
  //    [](Atom_list &self, const int &i, const Atom &atom) {
  //      self.at(normalise_index(i, self.size())) = atom;
  //    },
  //    "Replaces Atom at given index.")
  //.def(
  //    "delete_atom",
  //    [](Atom_list &self, const int &pos) { delitem_at_index(self, pos); },
  //    nb::arg("index"), "Delete atom at given index.")
  //.def(
  //    "delete_atoms",
  //    [](Atom_list &self, nb::slice slice) { delitem_slice(self, slice); },
  //    nb::arg("slice"), "Delete atoms in range from slice object.")
  //.def(
  //    "insert_atom",
  //    [](Atom_list &self, const Atom &atom, const int &pos) {
  //      add_item(self, atom, pos);
  //    },
  //    nb::arg("index"), nb::arg("Atom"),
  //    "Insert atom to the position of given index.")
  //.def( "__iter__", [](Atom_list &self) {
  //      return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<Atom_list>(), "iterator", &self[0],
  //      &self[self.size()] ); }, nb::keep_alive<0, 1>() )
  //.def(
  //    "insert_list",
  //    [](Atom_list &self, const std::vector<Atom> &a, const int &pos) {
  //      if (pos == -1)
  //        self.insert(self.end(), a.begin(), a.end());
  //      else
  //        self.insert(self.begin() + pos, a.begin(), a.end());
  //    },
  //    nb::arg("list"), nb::arg("pos") = -1,
  //    "Concatenate another Atom_list at the given position of this "
  //    "Atom_list. Default is at the end.")
  //.def(
  //    "copy", [](const Atom_list &self) { return self; },
  //    "Return a copy of the object.")
  //
}
