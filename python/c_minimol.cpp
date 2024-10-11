// Wrapper for clipper minimol
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// ##include "buccaneer/buccaneer-util.h"
#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper-minimol.h>
#include <clipper/clipper.h>
#include <gemmi/model.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace clipper;

// refer to gemmi mol.cpp to wrap the add_child/set_child etc

// all these need to replace .begin() with something else like [0], and
// children ... there is no access to children (private)
// template <typename Item>
// void delitem_at_index(std::vector<Item> &items, pybind11::ssize_t index)
//{
//  items.erase(items.begin() + index);
//}

void init_minimol(py::module &m) {
  py::class_<MAtom, Atom, PropertyManager> pyAtom(m, "MAtom");
  py::class_<MResidue, PropertyManager> pyResidue(m, "MResidue");
  py::class_<MChain, PropertyManager> pyChain(m, "MChain");
  py::class_<MModel, PropertyManager> pyModel(m, "MModel");

  py::enum_<MM::MODE>(m, "MODE", "Search modes.")
      .value("UNIQUE", MM::MODE::UNIQUE)
      .value("ANY", MM::MODE::ANY);

  py::enum_<MM::COPY>(m, "COPY", "Copy modes.")
      .value("COPY_NONE", MM::COPY::COPY_NONE)
      .value("COPY_M", MM::COPY::COPY_M)
      .value("COPY_P", MM::COPY::COPY_P)
      .value("COPY_MP", MM::COPY::COPY_MP)
      .value("COPY_C", MM::COPY::COPY_C)
      .value("COPY_MC", MM::COPY::COPY_MC)
      .value("COPY_PC", MM::COPY::COPY_PC)
      .value("COPY_MPC", MM::COPY::COPY_MPC)
      .value("MEMBERS", MM::COPY::MEMBERS)
      .value("PROPERTIES", MM::COPY::PROPERTIES)
      .value("CHILDREN", MM::COPY::CHILDREN);

  py::enum_<MResidue::TYPE>(m, "TYPE", "Rotamer library type.")
      .value("Default", MResidue::TYPE::Default)
      .value("Dunbrack", MResidue::TYPE::Dunbrack)
      .value("Richardson", MResidue::TYPE::Richardson)
      .export_values();

  pyAtom.def(py::init<>())
      .def(py::init<const Atom &>(), py::arg("atom"),
           "Constructor from clipper::Atom.")
      // also inherit methods from Atom class e.g. element, coord_orth...
      .def_property("id", &MAtom::id, &MAtom::set_id,
                    "Get/set atom ID, e.g. \" N  \"")
      .def_property("name", &MAtom::name, &MAtom::set_name,
                    "Get/set atom name, i.e. the ID omitting any alternate "
                    "conformation code.")
      .def("__repr__",
           [](const MAtom &self) {
             std::stringstream stream;
             auto coord = self.coord_orth();
             stream << "<clipper.MAtom " << self.name().trim();
             stream << " at (" << coord.x() << ", " << coord.y() << ", "
                    << coord.z();
             stream << ")>";
             return stream.str();
           })
      .def(
          "atom", [](const MAtom &self) -> const Atom & { return self.atom(); },
          py::return_value_policy::reference_internal, "Get atom.")
      .def(
          "atom", [](MAtom &self, Atom atm) { self.atom() = atm; },
          py::return_value_policy::reference_internal, "Set atom.")
      .def("copy_from", &MAtom::copy, py::arg("other"),
           py::arg("mode") = MM::COPY::COPY_C, "Configurable copy function.")
      .def(
          "copy", [](const MAtom &self) { return self; },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .def_static("id_tidy", &MAtom::id_tidy, py::arg("id"),
                  py::arg("is_gemmi") = false, "Convert id to standard format.")
      .def_static("id_match", &MAtom::id_match, py::arg("id1"), py::arg("id2"),
                  py::arg("mode"), "Compare IDs.")
      .doc() =
      "MiniMol atom object.\nThe MiniMol atom is derived "
      "from the basic clipper::Atom, with the addition of an 'id', "
      "which is a unique identifier within a monomer in accordance "
      "with the mmCIF definition.\nIn addition, it is a "
      "clipper::PropertyManager, which means you can add labelled "
      "properties of any type to the object. These may be simple "
      "strings, or complex objects such as maps, function objects, "
      "or whatever.\nThe most commonly used properties are: "
      "`-` \"CID\" The original CID of this atom in an MMDB heirarchy. "
      "The id() is the unique key which identifies an atom.";

  pyResidue.def(py::init<>())
      .def_property("id", &MResidue::id, &MResidue::set_id,
                    "Get/set monomer ID.")
      .def_property("type", &MResidue::type, &MResidue::set_type,
                    "Get/set monomer type, e.g. LYS, VAL, G.")
      //.def("set_type", &MResidue::set_type)
      .def("seqnum", &MResidue::seqnum, "Get monomer sequence number.")
      .def("set_seqnum", &MResidue::set_seqnum, py::arg("s"),
           py::arg("inscode") = "", "Set full sequence id.")
      .def("atom_list", &MResidue::atom_list,
           "Return a list of contained atoms.")
      .def("transform", &MResidue::transform, py::arg("rtop"),
           "Apply transformation to object.")
      .def("size", &MResidue::size, "Return number of atoms in monomer.")
      .def("__len__", &MResidue::size)
      .def("__repr__",
           [](const MResidue &self) {
             std::stringstream stream;
             stream << "<clipper.MResidue ";
             stream << self.id().trim() << "(" << self.type()
                    << ") containing ";
             stream << self.size() << " atom(s)>";
             return stream.str();
           })
      .def(
          "__getitem__",
          [](MResidue &self, const int i) -> const MAtom & {
            return self[normalise_index(i, self.size())];
          },
          py::arg("i"), py::return_value_policy::reference_internal,
          "Get atom.")
      .def(
          "__getitem__",
          [](MResidue &self, const std::string &n) -> const MAtom & {
            return self.find(n);
          },
          py::arg("id"), py::return_value_policy::reference_internal,
          "Set atom.")
      .def(
          "find",
          [](const MResidue &self, const std::string &n, const MM::MODE mode)
              -> const MAtom & { return self.find(n, mode); },
          py::arg("id"), py::arg("mode") = MM::MODE::UNIQUE,
          py::return_value_policy::reference_internal,
          "Lookup by id. If mode=UNIQUE, the alternate conformation code must "
          "match, otherwise the first atom with the same name is returned.")
      .def(
          "find",
          [](MResidue &self, const std::string &n, const MAtom &atm,
             const MM::MODE mode) -> MAtom & {
            return self.find(n, mode) = atm;
          },
          "Set atom by looking up id. If mode=UNIQUE, the alternate "
          "conformation code must match, otherwise the first atom with "
          "the same name is returned.")
      .def(
          "__setitem__",
          [](MResidue &self, const int i, const MAtom atm) {
            self[normalise_index(i, self.size())] = atm;
          },
          py::arg("i"), py::arg("atom"),
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](MResidue &self, const std::string &n, const MAtom atm) {
            self.find(n) = atm;
          },
          py::arg("id"), py::arg("atom"),
          py::return_value_policy::reference_internal)
      .def("select", &MResidue::select, py::arg("sel"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Creates a copy of this monomet containing only atoms described by "
           "selection string. The atom selection must contain an atom ID or "
           "a comma separated list of atom IDs, or \"*\" to select all atoms.")
      .def("select_index", &MResidue::select_index, py::arg("sel"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Creates a list of inidices of children matching the given "
           "selection string.")
      .def("insert", &MResidue::insert, py::arg("add"), py::arg("pos"),
           "Add atom to given position.")
      .def("lookup", &MResidue::lookup, py::arg("id"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Lookup atom by ID and return the index position.")
      .def(py::self & py::self)
      .def(py::self | py::self)
      .def("copy_from", &MResidue::copy, py::arg("other"),
           py::arg("mode") = MM::COPY::COPY_C, "Configurable copy function.")
      .def(
          "copy", [](const MResidue &self) { return self; },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .def_static("id_match", &MResidue::id_match, py::arg("id1"),
                  py::arg("id2"), py::arg("mode"), "Compare two IDs.")
      .def_static("id_tidy", &MResidue::id_tidy, py::arg("id"),
                  "Convert ID to standard format.")
      // UTILITY
      .def("build_carbonyl_oxygen",
           (void(MResidue::*)(const MResidue &)) &
               MResidue::protein_mainchain_build_carbonyl_oxygen,
           py::arg("next"),
           "Build carbonyl oxygen, given next residue in chain.")
      .def("build_carbonyl_oxygen",
           (void(MResidue::*)()) &
               MResidue::protein_mainchain_build_carbonyl_oxygen,
           "Build carbonyl oxygen, without next residue in chain.")
      .def("number_of_rotamers",
           (int(MResidue::*)(MResidue::TYPE) const) &
               MResidue::protein_sidechain_number_of_rotamers,
           py::arg("t"),
           "Get number of rotamers for protein sidechain, given a rotamer "
           "library type.")
      .def("number_of_rotamers",
           (int(MResidue::*)() const) &
               MResidue::protein_sidechain_number_of_rotomers,
           "Get number of rotamers for protein sidechain from Richardson "
           "rotamer library.")
      .def("build_sidechain_numbered_rotamer",
           (ftype(MResidue::*)(const int &, MResidue::TYPE)) &
               MResidue::protein_sidechain_build_rotamer,
           py::arg("n"), py::arg("t"),
           "Build numbered rotamer for protein sidechain.")
      .def("build_sidechain_numbered_rotamer",
           (ftype(MResidue::*)(const int &)) &
               MResidue::protein_sidechain_build_rotomer,
           py::arg("n"), "Build numbered rotamer for protein sidechain.")
      .def_static("protein_peptide_bond", &MResidue::protein_peptide_bond,
                  py::arg("m1"), py::arg("m2"), py::arg("r") = 1.5,
                  "Test if two peptide are adjacent.")
      .def_static("protein_ramachandran_phi",
                  &MResidue::protein_ramachandran_phi, py::arg("m1"),
                  py::arg("m2"),
                  "Return Ramachandran phi, or NaN if atoms missing.")
      .def_static("protein_ramachandran_psi",
                  &MResidue::protein_ramachandran_psi, py::arg("m1"),
                  py::arg("m2"),
                  "Return Ramachandran psi, or NaN if atoms missing.")
      .def_static("default_type", &MResidue::default_type,
                  "Return default rotamer library type.")
      .doc() = "MiniMol monomer (e.g. residue) object.\nThe MiniMol "
               "monomer object contains a list of clipper::MAtom. "
               "It has two properties: a sequence number and a type. "
               "The sequence number need not reflect the order in which "
               "the monomers are stored in a polymer. MResidue is an alias "
               "for MMonomer. In addition, it is a clipper::PropertyManager, "
               "refer documented details in MAtom class.";

  pyChain.def(py::init<>())
      .def_property("id", &MChain::id, &MChain::set_id, "Get/set id.")
      .def("atom_list", &MChain::atom_list, "Return list of contained atoms.")
      .def("transform", &MChain::transform, py::arg("rtop"),
           "Apply transformation to object.")
      .def("size", &MChain::size, "Return number of monomers in polymer.")
      .def("__len__", &MChain::size)
      .def("__repr__",
           [](const MChain &self) {
             std::stringstream stream;
             stream << "<clipper.MChain ";
             stream << self.id() << " containing ";
             stream << self.size() << " residue(s)>";
             return stream.str();
           })
      .def(
          "__getitem__",
          [](MChain &self, const int i) -> const MResidue & {
            return self[normalise_index(i, self.size())];
          },
          py::arg("i"), py::return_value_policy::reference_internal,
          "Get residue.")
      .def(
          "__getitem__",
          [](MChain &self, const std::string &n) -> const MResidue & {
            return self.find(n);
          },
          py::arg("id"), py::return_value_policy::reference_internal,
          "Get residue.")
      .def(
          "find",
          [](const MChain &self, const std::string &n, const MM::MODE mode)
              -> const MResidue & { return self.find(n, mode); },
          py::arg("id"), py::arg("mode") = MM::MODE::UNIQUE,
          py::return_value_policy::reference_internal, "Find residue by id.")
      .def(
          "find",
          [](MChain &self, const std::string &n, const MResidue &res,
             const MM::MODE mode) { self.find(n, mode) = res; },
          py::arg("id"), py::arg("res"), py::arg("mode"),
          "Find and set residue by id.")
      .def(
          "__setitem__",
          [](MChain &self, const int i, const MResidue res) {
            self[normalise_index(i, self.size())] = res;
          },
          py::arg("i"), py::arg("res"),
          py::return_value_policy::reference_internal, "Set residue.")
      .def(
          "__setitem__",
          [](MChain &self, const std::string &n, const MResidue res) {
            self.find(n) = res;
          },
          py::arg("id"), py::arg("res"),
          py::return_value_policy::reference_internal, "Set residue.")
      .def(
          "__iter__",
          [](MChain &self) {
            return py::make_iterator(&self[0], &self[self.size()]);
          },
          py::keep_alive<0, 1>())
      .def("insert", &MChain::insert, py::arg("add"), py::arg("pos"),
           "Add residue to given position.")
      .def("select", &MChain::select, py::arg("sel"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Create a copy of this polymer containing only the monomers and "
           "atoms described by the selection string.")
      .def("select_index", &MChain::select_index, py::arg("sel"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Get child indices mathcing the selection criteria.")
      .def(py::self & py::self, "and operator")
      .def(py::self | py::self, "or operator")
      .def("copy_from", &MChain::copy, py::arg("other"),
           py::arg("mode") = MM::COPY::COPY_C, "Configurable copy function.")
      .def(
          "copy", [](const MChain &self) { return self; },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .def_static("id_tidy", &MChain::id_tidy, py::arg("id"),
                  "Convert ID to standard format.")
      .def_static("id_match", &MChain::id_match, py::arg("id1"), py::arg("id2"),
                  py::arg("mode") = MM::MODE::UNIQUE, "Compare two ids.")
      .doc() = "MiniMol polymer (e.g. chain) object.\nThe MiniMol "
               "polymer object has one property: an identifying name. "
               "It contains a list of clipper::MMonomer. In addition, "
               "it is a clipper::PropertyManager, refer documented "
               "details in MAtom class.";

  pyModel.def(py::init<>())
      .def("atom_list", &MModel::atom_list, "Return list of contained atoms")
      .def("transform", &MModel::transform,
           "Apply transformation to object.") // maybe can use array/matrix?
      .def("size", &MModel::size, "Return number of polymers in model.")
      .def("__len__", &MModel::size)
      .def("__repr__",
           [](const MModel &self) {
             return "<clipper.MModel containing " +
                    std::to_string(self.size()) + " chain(s)>";
           })
      .def(
          "__getitem__",
          [](MModel &self, const int i) -> const MChain & {
            return self[normalise_index(i, self.size())];
          },
          py::arg("i")) // py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](MModel &self, const std::string &n) -> const MChain & {
            return self.find(n);
          },
          py::arg("n")) // py::return_value_policy::reference_internal)
      .def(
          "find",
          [](const MModel &self, const std::string &n, const MM::MODE mode)
              -> const MPolymer & { return self.find(n, mode); },
          py::arg("n"), py::arg("mode") = MM::MODE::UNIQUE,
          py::return_value_policy::reference_internal, "Find polymer by id.")
      .def(
          "find",
          [](MModel &self, const std::string &n, const MChain &chn,
             const MM::MODE mode) { self.find(n, mode) = chn; },
          py::arg("n"), py::arg("chain"), py::arg("mode") = MM::MODE::UNIQUE,
          "Find and set polymer by id.")
      .def(
          "__setitem__",
          [](MModel &self, const int i, const MChain chn) {
            self[normalise_index(i, self.size())] = chn;
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](MModel &self, const std::string &n, const MChain chn) {
            self.find(n) = chn;
          },
          py::return_value_policy::reference_internal)
      .def(
          "__iter__",
          [](MModel &self) {
            return py::make_iterator(&self[0], &self[self.size()]);
          },
          py::keep_alive<0, 1>())
      .def("select", &MModel::select, py::arg("selection"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Creates copy of this model containing only the polymers, "
           "monomers and atoms described by the selection string. The "
           "selection string must be of the form \'X/Y/Z\' where X is a "
           "polymer selection, Y is a monomer selection described "
           "under MMonomer::select(), and Z is an atom selection "
           "described under MAtom::select(). The polymer selection must "
           "contain a polymer ID or a comma separated list of "
           "polymer IDs, or \"*\" to select all polymers.")
      .def("select_index", &MModel::select_index, py::arg("selection"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Creates a list of indices of children matching the given "
           "selection string.")
      .def("lookup", &MModel::lookup, py::arg("id"), py::arg("mode"),
           "Lookup polymer by id.")
      .def("insert", &MModel::insert, py::arg("add"), py::arg("pos") = -1,
           "Add polymer to given position.")
      .def(py::self & py::self, "and operator.")
      .def(py::self | py::self, "or operator.")
      .def("select_atom_index", &MModel::select_atom_index, py::arg("sel"),
           py::arg("mode") = MM::MODE::UNIQUE,
           "Select and return a list of MAtomIndex.")
      .def("copy_from", &MModel::copy, py::arg("other"),
           py::arg("mode") = MM::COPY::COPY_C, "Configurable copy function.")
      .def(
          "copy", [](const MModel &self) { return self; },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .doc() = "MiniMol model object.\nThe MiniMol model object contains "
               "a list of clipper::MPolymer. It is a clipper::PropertyManager, "
               "refer documented details in MAtom class.";

  py::class_<MiniMol, MModel> minimol(m, "MiniMol");
  minimol.def(py::init<>())
      .def(py::init<const Spacegroup &, const Cell &>(), py::arg("spacegroup"),
           py::arg("cell"), "Constructor from spacegroup and cell.")
      .def("init", &MiniMol::init, py::arg("spacegroup"), py::arg("cell"),
           "Initialiser from spacegroup and cell.")
      //.def("__len__", [](const MiniMol &self) { return
      // self.model().size(); })
      .def("__repr__",
           [](const MiniMol &self) {
             std::stringstream stream;
             stream << "<clipper.MiniMol containing model with ";
             stream << self.model().size() << " chain(s)>";
             return stream.str();
           })
      .def_property_readonly("cell", &MiniMol::cell, "Get cell.")
      .def_property_readonly("spacegroup", &MiniMol::spacegroup,
                             "Get spacegroup.")
      .def(
          "model",
          [](const MiniMol &self) -> const MModel & { return self.model(); },
          py::return_value_policy::reference_internal, "Get model.")
      .def(
          "model", [](MiniMol &self, MModel mol) { self.model() = mol; },
          py::return_value_policy::reference_internal, "Set model.")
      .def(
          "copy",
          [](const MiniMol &self) { return self; }, // new MiniMol(self); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy.")
      .def("is_null", &MiniMol::is_null,
           "Test for null model (Uninitialised). ")
      .def(
          "is_empty", [](const MiniMol &self) { return (self.size() == 0); },
          "Test if MiniMol object is empty.")
      .def("symmetry_atom", &MiniMol::symmetry_atom, py::arg("index"),
           "Return symmetry atom by MAtomIndexSymmetry.");
  //.def("__copy__", [](const MiniMol &self) { return MiniMol(self); })
  //.def(
  //    "__deepcopy__",
  //    [](const MiniMol &self, py::dict memo) { return MiniMol(self); },
  //    py::arg("memo"));
  // need to bind MAtomIndexSymmetry from minimol_util.h

  //.def_property_readonly("model",
  // py::overload_cast<MModel>(&MiniMol::model,
  // py::const_)); .def("model", py::overload_cast<>(&MiniMol::model));
  //.def("model", (MModel)&MiniMol::model);
  // py::keep_alive<0, 1>());
  //.def("model", (std::vector<MChain>(MModel::*)()) & MiniMol::model,
  // py::keep_alive<0, 1>());

  //  py::class_<bk::PyCMiniMol>(m, "PyCMiniMol")
  //      .def(py::init<std::string &, bool>(),
  //      py::arg("filepath_to_structure") = "undefined",
  //      py::arg("enable_messages") = true) //; .def("get_mmol",
  //      &bk::PyCMiniMol::get_mmol, py::arg("filepath_to_structure") =
  //      "undefined");
}