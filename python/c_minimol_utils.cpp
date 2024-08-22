// Wrapper for clipper minimol utils
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper-minimol.h>
#include <clipper/minimol/minimol_utils.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void init_matomindex(py::module &m) {
  py::class_<MAtomIndex>(m, "MAtomIndex")
      .def(py::init<>())
      .def(py::init<const int &, const int &, const int &>(),
           py::arg("polymer"), py::arg("monomer"), py::arg("atom"),
           "Constructor from polymer, monomer and atom numbers.")
      .def("is_null", &MAtomIndex::is_null,
           "Test if object has been initialised.")
      .def_property_readonly("polymer", &MAtomIndex::polymer,
                             "Return polymer index.")
      .def_property_readonly("monomer", &MAtomIndex::monomer,
                             "Return monomer index.")
      .def_property_readonly("atom", &MAtomIndex::atom, "Return atom index.")
      .def(
          "__lt__",
          [](MAtomIndex self, MAtomIndex other) { return self < other; },
          "Less than operator.")
      .doc() =
      "Class for holding the indices of an atom within a MiniMol molecule "
      "class.\nThe indices remain valid only while no changes are made to "
      "the MiniMol object.";

  py::class_<MAtomIndexSymmetry, MAtomIndex>(m, "MAtomIndexSymmetry")
      .def(py::init<>())
      .def(py::init<const int &, const int &, const int &, const int &>(),
           py::arg("polymer"), py::arg("monomer"), py::arg("atom"),
           py::arg("sym"),
           "Constructor from polymer, monomer, atom numbers, and symmetry "
           "index.")
      .def("symmetry", &MAtomIndexSymmetry::symmetry, "Return symmetry index.")
      .doc() =
      "Class for holding the indices of an atom within a MiniMol molecule "
      "class.\nThe indices remain valid only while no changes are made to the "
      "MiniMol object. This class can also hold the number of a symmetry "
      "operator.";

  py::class_<MAtomNonBond>(m, "MAtomNonBond")
      .def(py::init<>())
      .def(py::init<const MiniMol &, double>(),
           "Constructor from MiniMol and grid radius.")
      .def("atoms_near", &MAtomNonBond::atoms_near, py::arg("coord"),
           py::arg("rad"),
           "Get a list of atoms in the rough vicinity of a coordinate.")
      .def("__call__", &MAtomNonBond::operator(), py::arg("coord"),
           py::arg("rad"),
           "Get a list of atoms in the rough vicinity of a coordinate.")
      .def("debug", &MAtomNonBond::debug, "Output debugging details.")
      .doc() = "Find atoms in the vicinity of some coordinate in real space "
               "Uses a fast non-bonded atom search.";
}