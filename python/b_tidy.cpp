// Nanobind bindings for buccaneer-tidy
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-tidy.h"
#include "commons.h"
#include "arrays.h"

using namespace clipper;

void add_model_tidy(nb::module_ &m) {
  nb::class_<ModelTidy> model_tidy(m, "ModelTidy", "Class for tidying model.");
  model_tidy
      .def( nb::init<double, int, String, bool>(), nb::arg( "rmsd" ) = 1.0, nb::arg( "nmin" ) = 12,
            nb::arg( "newrestype" ) = "ALA", nb::arg( "verbose" ) = false )
      .def( "tidy", &ModelTidy::tidy, nb::arg( "mol" ), nb::arg( "mol_mr" ), nb::arg( "seq" ), "Tidy model." )
      .def_static( "tidy_model",
            []( clipper::MiniMol &mol, const clipper::MiniMol &mol_mr, const clipper::MMoleculeSequence &seq,
                double &rmsd, int &nmin, String &newrestype, bool &verbose ) -> bool {
              ModelTidy mtidy( rmsd, nmin, newrestype, verbose );
              bool chk = mtidy.tidy( mol, mol_mr, seq );
              // if ( !chk ) std::cout << "ModelTidy error" << std::endl; // can't happen
              return chk;
            }, nb::arg( "mol" ), nb::arg( "mol_mr" ), nb::arg( "seq" ), nb::arg( "rmsd" ) = 1.0, nb::arg( "nmin" ) = 12,
            nb::arg( "newrestype" ) = "ALA", nb::arg( "verbose" ) = false, "Static function to tidy model" )
      .def_static( "chain_renumber",
                   nb::overload_cast<MiniMol &, const MMoleculeSequence &>( &ModelTidy::chain_renumber ),
                   nb::arg( "mol" ), nb::arg( "seq" ), "Renumber residues in chain." )
      //.def_static("chain_renumber",
      //            static_cast<std::vector<int> (*)(MiniMol &,
      //                                             const MMoleculeSequence &)>(
      //                &ModelTidy::chain_renumber),
      //            nb::arg("mol"), nb::arg("seq"), "Renumber residues in chain.")
      .def_static( "chain_assign", &ModelTidy::chain_assign, nb::arg( "mol" ), nb::arg( "mol_mr" ),
                   nb::arg( "seq_nums" ), nb::arg( "rmsd" ), nb::arg( "nmin" ), "Assign chain fragments to MR model." )
      .def_static( "chain_move", &ModelTidy::chain_move, nb::arg( "mol" ), nb::arg( "mol_mr" ), nb::arg( "chnnums" ) )
      .def_static( "sequence_correct", &ModelTidy::sequence_correct, nb::arg( "mol" ), nb::arg( "seq" ),
                   nb::arg( "seqnums" ), nb::arg( "newrestype" ), "Correct the sequence." )
      .def_static( "sequence_count", &ModelTidy::sequence_count, nb::arg( "mol" ), "Count sequence types." )
      // this returns clipper::array2d<int>
      .def_static( "sequence_flags", &ModelTidy::sequence_flags, nb::arg( "mol" ),
                   "Make flags of used residue numbers." )
      .def_static( "trim", &ModelTidy::trim, nb::arg( "mol" ), nb::arg( "seq" ),
                   "Trim ends of chains which go beyond the end of the sequnce." )
      .def( "__repr__", []( const ModelTidy &self ) { return "<buccaneer.ModelTidy class.>"; } );
}