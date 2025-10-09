// Nanobind bindings for clipper minimol
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-minimol-gemmi.h>
#include <clipper/minimol/minimol_utils.h>
#include <nanobind/operators.h>

using namespace clipper;

void add_matomindex( nb::module_ &m ) {
  nb::class_<MAtomIndex>( m, "MAtomIndex" )
      .def( nb::init<>() )
      .def( nb::init<const int &, const int &, const int &>(), nb::arg( "polymer" ), nb::arg( "monomer" ),
            nb::arg( "atom" ), "Constructor from polymer, monomer and atom numbers." )
      .def( "is_null", &MAtomIndex::is_null, "Test if object has been initialised." )
      .def_prop_ro( "polymer", &MAtomIndex::polymer, "Return polymer index." )
      .def_prop_ro( "monomer", &MAtomIndex::monomer, "Return monomer index." )
      .def_prop_ro( "atom", &MAtomIndex::atom, "Return atom index." )
      .def( nb::self < nb::self )
      //.def(
      //    "__lt__", []( MAtomIndex self, MAtomIndex other ) { return self < other; }, "Less than operator." )
      .doc() = "Class for holding the indices of an atom within a MiniMol molecule "
               "class.\nThe indices remain valid only while no changes are made to "
               "the MiniMol object.";

  nb::class_<MAtomIndexSymmetry, MAtomIndex>( m, "MAtomIndexSymmetry" )
      .def( nb::init<>() )
      .def( nb::init<const int &, const int &, const int &, const int &>(), nb::arg( "polymer" ), nb::arg( "monomer" ),
            nb::arg( "atom" ), nb::arg( "sym" ),
            "Constructor from polymer, monomer, atom numbers, and symmetry "
            "index." )
      .def( "symmetry", &MAtomIndexSymmetry::symmetry, "Return symmetry index." )
      .doc() = "Class for holding the indices of an atom within a MiniMol molecule "
               "class.\nThe indices remain valid only while no changes are made to the "
               "MiniMol object. This class can also hold the number of a symmetry "
               "operator.";

  nb::class_<MAtomNonBond>( m, "MAtomNonBond" )
      .def( nb::init<>(), "Null constructor." )
      .def( nb::init<const MiniMol &, double>(), "Constructor from MiniMol and grid radius." )
      .def( "atoms_near", &MAtomNonBond::atoms_near, nb::arg( "coord" ), nb::arg( "rad" ),
            "Get a list of atoms in the rough vicinity of a coordinate." )
      .def( "__call__", &MAtomNonBond::operator(), nb::arg( "coord" ), nb::arg( "rad" ),
            "Get a list of atoms in the rough vicinity of a coordinate." )
      .def( "debug", &MAtomNonBond::debug, "Output debugging details." )
      .doc() = "Find atoms in the vicinity of some coordinate in real space.\n"
               "Uses a fast non-bonded atom search.";
}