// Nanobind bindings for clipper minimol
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-minimol-gemmi.h>
#include <clipper/clipper.h>
#include <nanobind/make_iterator.h>
#include <nanobind/operators.h>

// #include <gemmi/model.hpp>

using namespace clipper;

// refer to gemmi mol.cpp to wrap the add_child/set_child etc

// all these need to replace .begin() with something else like [0], and
// children ... there is no access to children (private)
// template <typename Item>
// void delitem_at_index(std::vector<Item> &items, pybind11::ssize_t index)
//{
//  items.erase(items.begin() + index);
//}

void init_minimol( nb::module_ &m, nb::module_ &mm ) {
  // since nanobind doesn't support multiple inheritance, we leave out PropertyManager but need to bind the methods
  // individually
  nb::class_<MAtom, Atom> pyAtom( m, "MAtom" );
  nb::class_<MResidue, PropertyManager> pyResidue( m, "MResidue" );
  nb::class_<MChain, PropertyManager> pyChain( m, "MChain" );
  nb::class_<MModel, PropertyManager> pyModel( m, "MModel" );

  nb::enum_<MM::MODE>( mm, "MODE", "Search modes." )
      .value( "UNIQUE", MM::MODE::UNIQUE )
      .value( "ANY", MM::MODE::ANY )
      .export_values();

  nb::enum_<MM::COPY>( mm, "COPY", "Copy modes." )
      .value( "COPY_NONE", MM::COPY::COPY_NONE )
      .value( "COPY_M", MM::COPY::COPY_M )
      .value( "COPY_P", MM::COPY::COPY_P )
      .value( "COPY_MP", MM::COPY::COPY_MP )
      .value( "COPY_C", MM::COPY::COPY_C )
      .value( "COPY_MC", MM::COPY::COPY_MC )
      .value( "COPY_PC", MM::COPY::COPY_PC )
      .value( "COPY_MPC", MM::COPY::COPY_MPC )
      .value( "MEMBERS", MM::COPY::MEMBERS )
      .value( "PROPERTIES", MM::COPY::PROPERTIES )
      .value( "CHILDREN", MM::COPY::CHILDREN )
      .export_values();

  nb::enum_<MResidue::TYPE>( pyResidue, "TYPE", "Rotamer library type." )
      .value( "Default", MResidue::TYPE::Default )
      .value( "Dunbrack", MResidue::TYPE::Dunbrack )
      .value( "Richardson", MResidue::TYPE::Richardson )
      .export_values();

  pyAtom.def( nb::init<>() )
      .def( nb::init<const Atom &>(), nb::arg( "atom" ), "Constructor from clipper::Atom." )
      // also inherit methods from Atom class e.g. element, coord_orth...
      .def_prop_ro(
          "unpadded_id", []( const MAtom &self ) { return self.id().trim(); },
          "Get unpadded atom id, e.g. \"N\", \"CA\", \"CG1\", \"CA:A\"." )
      .def_prop_rw(
          "id", &MAtom::id, []( MAtom &self, const String &id ) { self.set_id( id ); },
          nb::for_getter( "Get atom id, e.g. \" N  \", \" CA \", \" CG1\", \" CA :A\"." ),
          nb::for_setter( "Set atom id." ) )
      .def_prop_rw(
          "name", &MAtom::name, []( MAtom &self, const String &name ) { self.set_name( name ); },
          nb::for_getter( "Get atom name, i.e. the ID, omitting any alternate conformation code." ),
          nb::for_setter( "Set atom name, include \":\" for alternate conformation i.e. CA:A." ) )
      //.def("set_id", &MAtom::set_id, nb::arg("id"), nb::arg("is_gemmi") = false,
      //     "Set atom id.")
      .def( "set_name", &MAtom::set_name, nb::arg( "name" ), nb::arg( "altconf" ) = "", nb::arg( "is_gemmi" ) = false,
            "Set atom name." )
      .def( "__repr__",
            []( const MAtom &self ) {
              return ( "<clipper.MAtom " + self.name().trim() + " at (" + String( self.coord_orth().x() ) + ", " +
                       String( self.coord_orth().y() ) + ", " + String( self.coord_orth().z() ) + ")>" );
            } )
      .def_prop_rw(
          "atom", []( const MAtom &self ) -> const Atom & { return self.atom(); },
          []( MAtom &self, Atom atm ) { self.atom() = atm; }, nb::for_getter( "Get atom." ),
          nb::for_setter( "Set atom." ), nb::rv_policy::reference_internal )
      //.def("set_atom", [](MAtom& self, Atom atm) { self.atom() = atm; },
      //    nb::rv_policy::reference_internal, "Set atom.")
      .def( "copy_from", &MAtom::copy, nb::arg( "other" ), nb::arg( "mode" ) = MM::COPY::COPY_C,
            "Configurable copy function." )
      .def(
          "clone", []( const MAtom &self ) { return new MAtom( self ); },
          "Return a clone of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .def_static( "id_tidy", &MAtom::id_tidy, nb::arg( "id" ), nb::arg( "is_gemmi" ) = false,
                   "Convert id to standard format." )
      .def_static( "id_match", &MAtom::id_match, nb::arg( "id1" ), nb::arg( "id2" ), nb::arg( "mode" ), "Compare IDs." )
      // property manager bindings
      //.def("copy_property", []( MAtom& self, const PropertyManager& mgr ) -> PropertyManager& { return
      // self.PropertyManager::copy(mgr); }, nb::arg("property_manager"),
      //     "Copy manager.")
      .def( "get_property", &MAtom::get_property, nb::arg( "label" ), nb::rv_policy::reference_internal,
            "Get a labelled property from the list." )
      .def( "set_property", &MAtom::set_property, nb::arg( "label" ), nb::arg( "property" ),
            "Add a labelled property to the list." )
      .def( "exists_property", &MAtom::exists_property, nb::arg( "label" ), "Test if property exists." )
      .def( "delete_property", &MAtom::delete_property, nb::arg( "label" ), "Delete property." )
      .doc() = "MiniMol atom object.\nThe MiniMol atom is derived "
               "from the basic clipper::Atom, with the addition of an 'id', "
               "which is a unique identifier within a monomer in accordance "
               "with the mmCIF definition.\nIn addition, it is a "
               "clipper::PropertyManager, which means you can add labelled "
               "properties of any type to the object. These may be simple "
               "strings, or complex objects such as maps, function objects, "
               "or whatever.\nThe most commonly used properties are: "
               "`-` \"CID\" The original CID of this atom in an MMDB heirarchy. "
               "The id() is the unique key which identifies an atom.";

  pyResidue.def( nb::init<>() )
      .def_prop_rw( "id", &MResidue::id, &MResidue::set_id, "Get/set monomer ID." )
      .def_prop_rw( "type", &MResidue::type, &MResidue::set_type, "Get/set monomer type, e.g. LYS, VAL, G." )
      //.def("set_type", &MResidue::set_type)
      .def_prop_ro( "seqnum", &MResidue::seqnum, "Get monomer sequence number." )
      .def( "set_seqnum", &MResidue::set_seqnum, nb::arg( "num" ), nb::arg( "inscode" ) = "", "Set full sequence id." )
      .def( "atom_list", &MResidue::atom_list, "Return a list of contained atoms." )
      .def( "transform", &MResidue::transform, nb::arg( "rtop" ), "Apply transformation to object." )
      .def( "size", &MResidue::size, "Return number of atoms in monomer." )
      .def( "__len__", &MResidue::size )
      .def( "__repr__",
            []( const MResidue &self ) {
              return "<clipper.MResidue " + self.id().trim() + "(" + self.type() + ") containing " +
                     String( self.size() ) + " atom(s)>";
            } )
      .def(
          "__getitem__",
          []( MResidue &self, const int i ) -> const MAtom & { return self[normalise_index( i, self.size() )]; },
          nb::arg( "i" ), nb::rv_policy::reference_internal, "Get atom." )
      .def(
          "__getitem__", []( MResidue &self, const std::string &n ) -> const MAtom & { return self.find( n ); },
          nb::arg( "id" ), nb::rv_policy::reference_internal, "Get atom." )
      .def(
          "find",
          []( const MResidue &self, const std::string &n, const MM::MODE mode ) -> const MAtom & {
            return self.find( n, mode );
          },
          nb::arg( "id" ), nb::arg( "mode" ) = MM::MODE::UNIQUE, nb::rv_policy::reference_internal,
          "Lookup by id. If mode=UNIQUE, the alternate conformation "
          "code must match, otherwise the first atom with the same name is returned." )
      .def(
          "find",
          []( MResidue &self, const std::string &n, const MAtom &atm, const MM::MODE mode ) -> MAtom & {
            return self.find( n, mode ) = atm;
          },
          nb::arg( "id" ), nb::arg( "atom" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
          "Set atom by looking up id. If mode=UNIQUE, the alternate "
          "conformation code must match, otherwise the first atom with the same name is returned." )
      .def(
          "__setitem__",
          []( MResidue &self, const int i, const MAtom atm ) { self[normalise_index( i, self.size() )] = atm; },
          nb::arg( "i" ), nb::arg( "atom" ), nb::rv_policy::reference_internal )
      .def(
          "__setitem__", []( MResidue &self, const std::string &n, const MAtom atm ) { self.find( n ) = atm; },
          nb::arg( "id" ), nb::arg( "atom" ), nb::rv_policy::reference_internal )
      .def( "select", &MResidue::select, nb::arg( "sel" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Creates a copy of this monomet containing only atoms described by "
            "selection string. The atom selection must contain an atom ID or "
            "a comma separated list of atom IDs, or \"*\" to select all atoms." )
      .def( "select_index", &MResidue::select_index, nb::arg( "sel" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Creates a list of inidices of children matching the given "
            "selection string." )
      .def( "insert", &MResidue::insert, nb::arg( "add" ), nb::arg( "pos" ) = -1, "Add atom to given position." )
      .def( "lookup", &MResidue::lookup, nb::arg( "id" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Lookup atom by ID and return the index position." )
      .def( nb::self & nb::self )
      .def( nb::self | nb::self )
      .def( "copy_from", &MResidue::copy, nb::arg( "other" ), nb::arg( "mode" ) = MM::COPY::COPY_C,
            "Configurable copy function." )
      .def(
          "clone", []( const MResidue &self ) { return new MResidue( self ); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .def_static( "id_match", &MResidue::id_match, nb::arg( "id1" ), nb::arg( "id2" ), nb::arg( "mode" ),
                   "Compare two IDs." )
      .def_static( "id_tidy", &MResidue::id_tidy, nb::arg( "id" ), "Convert ID to standard format." )
      // UTILITY
      .def( "build_carbonyl_oxygen",
            ( void ( MResidue::* )( const MResidue & ) )&MResidue::protein_mainchain_build_carbonyl_oxygen,
            nb::arg( "next" ), "Build carbonyl oxygen, given next residue in chain." )
      .def( "build_carbonyl_oxygen", ( void ( MResidue::* )() )&MResidue::protein_mainchain_build_carbonyl_oxygen,
            "Build carbonyl oxygen, without next residue in chain." )
      .def( "number_of_rotamers",
            ( int ( MResidue::* )( MResidue::TYPE ) const ) & MResidue::protein_sidechain_number_of_rotamers,
            nb::arg( "t" ),
            "Get number of rotamers for protein sidechain, given a rotamer "
            "library type." )
      .def( "number_of_rotamers", ( int ( MResidue::* )() const ) & MResidue::protein_sidechain_number_of_rotomers,
            "Get number of rotamers for protein sidechain from Richardson "
            "rotamer library." )
      .def( "build_sidechain_numbered_rotamer",
            ( ftype ( MResidue::* )( const int &, MResidue::TYPE ) )&MResidue::protein_sidechain_build_rotamer,
            nb::arg( "n" ), nb::arg( "t" ), "Build numbered rotamer for protein sidechain." )
      .def( "build_sidechain_numbered_rotamer",
            ( ftype ( MResidue::* )( const int & ) )&MResidue::protein_sidechain_build_rotomer, nb::arg( "n" ),
            "Build numbered rotamer for protein sidechain." )
      .def_static( "protein_peptide_bond", &MResidue::protein_peptide_bond, nb::arg( "m1" ), nb::arg( "m2" ),
                   nb::arg( "r" ) = 1.5, "Test if two peptide are adjacent." )
      .def_static( "protein_ramachandran_phi", &MResidue::protein_ramachandran_phi, nb::arg( "m1" ), nb::arg( "m2" ),
                   "Return Ramachandran phi, or NaN if atoms missing." )
      .def_static( "protein_ramachandran_psi", &MResidue::protein_ramachandran_psi, nb::arg( "m1" ), nb::arg( "m2" ),
                   "Return Ramachandran psi, or NaN if atoms missing." )
      .def_static( "default_type", &MResidue::default_type, "Return default rotamer library type." )
      .doc() = "MiniMol monomer (e.g. residue) object.\nThe MiniMol "
               "monomer object contains a list of clipper::MAtom. "
               "It has two properties: a sequence number and a type. "
               "The sequence number need not reflect the order in which "
               "the monomers are stored in a polymer. MResidue is an alias "
               "for MMonomer. In addition, it is a clipper::PropertyManager, "
               "refer documented details in MAtom class.";

  pyChain.def( nb::init<>() )
      .def_prop_rw( "id", &MChain::id, &MChain::set_id, "Get/set id." )
      .def( "atom_list", &MChain::atom_list, "Return list of contained atoms." )
      .def( "transform", &MChain::transform, nb::arg( "rtop" ), "Apply transformation to object." )
      .def( "size", &MChain::size, "Return number of monomers in polymer." )
      .def( "__len__", &MChain::size )
      .def( "__repr__",
            []( const MChain &self ) {
              return "<clipper.MChain " + self.id().trim() + " containing " + String( self.size() ) + " residue(s)>";
            } )
      .def(
          "__getitem__",
          []( MChain &self, const int i ) -> const MResidue & { return self[normalise_index( i, self.size() )]; },
          nb::arg( "i" ), nb::rv_policy::reference_internal, "Get residue." )
      .def(
          "__getitem__", []( MChain &self, const std::string &n ) -> const MResidue & { return self.find( n ); },
          nb::arg( "id" ), nb::rv_policy::reference_internal, "Get residue." )
      .def(
          "find",
          []( const MChain &self, const std::string &n, const MM::MODE mode ) -> const MResidue & {
            return self.find( n, mode );
          },
          nb::arg( "id" ), nb::arg( "mode" ) = MM::MODE::UNIQUE, nb::rv_policy::reference_internal,
          "Find residue by id." )
      .def(
          "find",
          []( MChain &self, const std::string &n, const MResidue &res, const MM::MODE mode ) {
            self.find( n, mode ) = res;
          },
          nb::arg( "id" ), nb::arg( "res" ), nb::arg( "mode" ), "Find and set residue by id." )
      .def(
          "__setitem__",
          []( MChain &self, const int i, const MResidue res ) { self[normalise_index( i, self.size() )] = res; },
          nb::arg( "i" ), nb::arg( "res" ), nb::rv_policy::reference_internal, "Set residue." )
      .def(
          "__setitem__", []( MChain &self, const std::string &n, const MResidue res ) { self.find( n ) = res; },
          nb::arg( "id" ), nb::arg( "res" ), nb::rv_policy::reference_internal, "Set residue." )
      .def(
          "__iter__",
          []( MChain &self ) {
            return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<MChain>(), "iterator", &self[0],
                                                                         &self[self.size()] );
          },
          nb::keep_alive<0, 1>() )
      .def( "insert", &MChain::insert, nb::arg( "add" ), nb::arg( "pos" ), "Add residue to given position." )
      .def( "select", &MChain::select, nb::arg( "sel" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Create a copy of this polymer containing only the monomers and "
            "atoms described by the selection string." )
      .def( "select_index", &MChain::select_index, nb::arg( "sel" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Get child indices matching the selection criteria." )
      .def( nb::self & nb::self, "and operator" )
      .def( nb::self | nb::self, "or operator" )
      .def( "copy_from", &MChain::copy, nb::arg( "other" ), nb::arg( "mode" ) = MM::COPY::COPY_C,
            "Configurable copy function." )
      .def(
          "clone", []( const MChain &self ) { return new MChain( self ); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .def_static( "id_tidy", &MChain::id_tidy, nb::arg( "id" ), "Convert ID to standard format." )
      .def_static( "id_match", &MChain::id_match, nb::arg( "id1" ), nb::arg( "id2" ),
                   nb::arg( "mode" ) = MM::MODE::UNIQUE, "Compare two ids." )
      .doc() = "MiniMol polymer (e.g. chain) object.\nThe MiniMol "
               "polymer object has one property: an identifying name. "
               "It contains a list of clipper::MMonomer. In addition, "
               "it is a clipper::PropertyManager, refer documented "
               "details in MAtom class.";

  pyModel.def( nb::init<>() )
      .def( "atom_list", &MModel::atom_list, "Return list of contained atoms" )
      .def( "transform", &MModel::transform,
            "Apply transformation to object." ) // maybe can use array/matrix?
      .def( "size", &MModel::size, "Return number of polymers in model." )
      .def( "__len__", &MModel::size )
      .def( "__repr__",
            []( const MModel &self ) {
              return "<clipper.MModel containing " + std::to_string( self.size() ) + " chain(s)>";
            } )
      .def(
          "__getitem__",
          []( MModel &self, const int i ) -> const MChain & { return self[normalise_index( i, self.size() )]; },
          nb::arg( "i" ), nb::rv_policy::reference_internal )
      .def(
          "__getitem__", []( MModel &self, const std::string &n ) -> const MChain & { return self.find( n ); },
          nb::arg( "n" ), nb::rv_policy::reference_internal )
      .def(
          "find",
          []( const MModel &self, const std::string &n, const MM::MODE mode ) -> const MPolymer & {
            return self.find( n, mode );
          },
          nb::arg( "n" ), nb::arg( "mode" ) = MM::MODE::UNIQUE, nb::rv_policy::reference_internal,
          "Find polymer by id." )
      .def(
          "find",
          []( MModel &self, const std::string &n, const MChain &chn, const MM::MODE mode ) {
            self.find( n, mode ) = chn;
          },
          nb::arg( "n" ), nb::arg( "chain" ), nb::arg( "mode" ) = MM::MODE::UNIQUE, "Find and set polymer by id." )
      .def(
          "__setitem__",
          []( MModel &self, const int i, const MChain chn ) { self[normalise_index( i, self.size() )] = chn; },
          nb::rv_policy::reference_internal )
      .def(
          "__setitem__", []( MModel &self, const std::string &n, const MChain chn ) { self.find( n ) = chn; },
          nb::rv_policy::reference_internal )
      .def(
          "__iter__",
          []( MModel &self ) {
            return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<MModel>(), "iterator", &self[0],
                                                                         &self[self.size()] );
          },
          nb::keep_alive<0, 1>() )
      .def( "select", &MModel::select, nb::arg( "selection" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Creates copy of this model containing only the polymers, "
            "monomers and atoms described by the selection string. The "
            "selection string must be of the form \'X/Y/Z\' where X is a "
            "polymer selection, Y is a monomer selection described "
            "under MMonomer::select(), and Z is an atom selection "
            "described under MAtom::select(). The polymer selection must "
            "contain a polymer ID or a comma separated list of "
            "polymer IDs, or \"*\" to select all polymers." )
      .def( "select_index", &MModel::select_index, nb::arg( "selection" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Creates a list of indices of children matching the given "
            "selection string." )
      .def( "lookup", &MModel::lookup, nb::arg( "id" ), nb::arg( "mode" ), "Lookup polymer by id." )
      .def( "insert", &MModel::insert, nb::arg( "add" ), nb::arg( "pos" ) = -1, "Add polymer to given position." )
      .def( nb::self & nb::self, "and operator." )
      .def( nb::self | nb::self, "or operator." )
      .def( "select_atom_index", &MModel::select_atom_index, nb::arg( "sel" ), nb::arg( "mode" ) = MM::MODE::UNIQUE,
            "Select and return a list of MAtomIndex." )
      .def( "copy_from", &MModel::copy, nb::arg( "other" ), nb::arg( "mode" ) = MM::COPY::COPY_C,
            "Configurable copy function." )
      .def(
          "clone", []( const MModel &self ) { return new MModel( self ); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .doc() = "MiniMol model object.\nThe MiniMol model object contains "
               "a list of clipper::MPolymer. It is a clipper::PropertyManager, "
               "refer documented details in MAtom class.";

  nb::class_<MiniMol, MModel> minimol( m, "MiniMol" );
  minimol.def( nb::init<>() )
      .def( nb::init<const Spacegroup &, const Cell &>(), nb::arg( "spacegroup" ), nb::arg( "cell" ),
            "Constructor from spacegroup and cell." )
      .def( "init", &MiniMol::init, nb::arg( "spacegroup" ), nb::arg( "cell" ),
            "Initialiser from spacegroup and cell." )
      //.def("__len__", [](const MiniMol &self) { return
      // self.model().size(); })
      .def( "__repr__",
            []( const MiniMol &self ) {
              return "<clipper.MiniMol containing model with " + std::to_string( self.model().size() ) + " chain(s)>";
            } )
      .def_prop_ro( "cell", &MiniMol::cell, "Get cell." )
      .def_prop_ro( "spacegroup", &MiniMol::spacegroup, "Get spacegroup." )
      .def_prop_rw(
          "model", []( const MiniMol &self ) -> const MModel & { return self.model(); },
          []( MiniMol &self, MModel mol ) { self.model() = mol; }, nb::for_getter( nb::sig("def get_model(self, /) -> MModel" )),
          nb::for_setter( nb::sig("def set_model(self, mol: MModel, /) -> None" )), nb::for_getter( "Get model." ),
          nb::for_setter( "Set model." ), nb::rv_policy::reference_internal, "Get model." )
      //.def( "set_model", [](MiniMol& self, MModel mol) { self.model() = mol; },
      //      nb::rv_policy::reference_internal, "Set model." )
      .def(
          "clone", []( const MiniMol &self ) { return new MiniMol( self ); },
          "Return a copy of object. Use this to make copy because "
          "assignment operator in Python only create bindings not copy." )
      .def( "is_null", &MiniMol::is_null, "Test for null model (Uninitialised). " )
      .def(
          "is_empty", []( const MiniMol &self ) { return ( self.size() == 0 ); }, "Test if MiniMol object is empty." )
      .def( "symmetry_atom", &MiniMol::symmetry_atom, nb::arg( "index" ),
            "Return symmetry atom by MAtomIndexSymmetry." );
      //add_convert_structure2mmol(m, minimol);
  //.def("__copy__", [](const MiniMol &self) { return MiniMol(self); })
  //.def(
  //    "__deepcopy__",
  //    [](const MiniMol &self, nb::dict memo) { return MiniMol(self); },
  //    nb::arg("memo"));
  // need to bind MAtomIndexSymmetry from minimol_util.h

  //.def_prop_ro("model",
  // nb::overload_cast<MModel>(&MiniMol::model,
  // nb::const_)); .def("model", nb::overload_cast<>(&MiniMol::model));
  //.def("model", (MModel)&MiniMol::model);
  // nb::keep_alive<0, 1>());
  //.def("model", (std::vector<MChain>(MModel::*)()) & MiniMol::model,
  // nb::keep_alive<0, 1>());

  //  nb::class_<bk::PyCMiniMol>(m, "PyCMiniMol")
  //      .def(nb::init<std::string &, bool>(),
  //      nb::arg("filepath_to_structure") = "undefined",
  //      nb::arg("enable_messages") = true) //; .def("get_mmol",
  //      &bk::PyCMiniMol::get_mmol, nb::arg("filepath_to_structure") =
  //      "undefined");
}