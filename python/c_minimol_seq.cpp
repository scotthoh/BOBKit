// Nanobind bindings for clipper minimol sequence
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-minimol-gemmi.h>
#include <nanobind/operators.h>

using namespace clipper;

void declare_mpolymer_sequence( nb::module_ &m ) {
  nb::class_<MPolymerSequence>( m, "MPolymerSequence" )
      .def( nb::init<>() )
      .def_prop_rw( "id", &MPolymerSequence::id, &MPolymerSequence::set_id,
                    nb::for_getter( nb::sig( "def id(self, /) -> str" ) ),
                    nb::for_setter( nb::sig( "def id(self, arg: str, /) -> None" ) ),
                    nb::for_getter( "Get sequence id." ), nb::for_setter( "Set sequence id." ) )
      .def_prop_rw( "sequence", &MPolymerSequence::sequence, &MPolymerSequence::set_sequence, "Get/set sequence." )
      .def_static( "id_tidy", &MPolymerSequence::id_tidy, nb::arg( "id" ), "Convert id to standard format." )
      .def_static( "id_match", &MPolymerSequence::id_match, nb::arg( "id1" ), nb::arg( "id2" ), nb::arg( "mode" ),
                   "Compare two ids." )
      .def( "__repr__",
            []( const MPolymerSequence self ) {
              return "<clipper.MPolymerSequence " + self.id() + self.sequence() + ">";
            } )
      .def( "__str__", []( const MPolymerSequence self ) { return self.id() + self.sequence(); } )
      .doc() = "Polymer sequence object.\nThe polymer sequence object "
               "represents the named sequence of a single chain.";
}

void declare_mmolecule_sequence( nb::module_ &m ) {
  nb::class_<MMoleculeSequence>( m, "MMoleculeSequence" )
      .def( nb::init<>() )
      .def( "size", &MMoleculeSequence::size, "Return number of polymer sequences in model." )
      .def( "__len__", &MMoleculeSequence::size )
      .def( "__repr__",
            []( const MMoleculeSequence &self ) {
              return "<clipper.MMoleculeSequence of length " + String( self.size() ) + " >";
            } )
      .def( "__str__",
            []( const MMoleculeSequence &self ) {
              std::string s = "";
              for ( int i = 0; i < self.size(); i++ ) {
                s += self[i].sequence();
                s += "\n";
              }
              return s;
            } )
      .def(
          "__getitem__", []( const MMoleculeSequence &self, const int &i ) { return self[i]; },
          "Get polymer sequence." )
      .def(
          "__setitem__", []( MMoleculeSequence &self, const int &i, MPolymerSequence &value ) { self[i] = value; },
          "Set polymer sequence." )
      .def_prop_rw( "find",
                    ( const MPolymerSequence &( MMoleculeSequence::* )( const String &, const MM::MODE ) const ) &
                        MMoleculeSequence::find,
                    ( MPolymerSequence & ( MMoleculeSequence::* )( const String &, const MM::MODE ) ) &
                        MMoleculeSequence::find,
                    "Get/set polymer sequence by id." )
      .def( "lookup", &MMoleculeSequence::lookup, nb::arg( "id" ), nb::arg( "mode" ), "Lookup polymer sequence by id." )
      //.def("insert", &MMoleculeSequence::insert, nb::arg("seq"),
      //     nb::arg("pos") = -1)
      .def(
          "insert",
          []( MMoleculeSequence &self, const MPolymerSequence &id, const int &pos ) { self.insert( id, pos ); },
          nb::arg( "seq" ), nb::arg( "pos" ) = -1, "Add polymer sequence." )
      .def(
          "insert",
          []( MMoleculeSequence &self, const String &id, const String &seq, const int &pos ) {
            MPolymerSequence mpseq;
            mpseq.set_id( id );
            mpseq.set_sequence( seq );
            self.insert( mpseq, pos );
          },
          nb::arg( "id" ), nb::arg( "seq" ), nb::arg( "pos" ) = -1, "Add polymer sequence." )
      .def( "is_null", &MMoleculeSequence::is_null, "Test for null model." )
      .doc() = "Molecule sequence object.\nThe molecule sequence "
               "object is a list of polymer sequence objects representing "
               "the named sequences of all the chains in a molecule.";
}

void declare_msequence_align( nb::module_ &m ) {
  nb::class_<MSequenceAlign> mseqalign( m, "MSequenceAlign" );
  nb::enum_<MSequenceAlign::TYPE>( mseqalign, "TYPE", "Alignment method." )
      .value( "GLOBAL", MSequenceAlign::TYPE::GLOBAL )
      .value( "LOCAL", MSequenceAlign::TYPE::LOCAL )
      .export_values();

  mseqalign
      .def( nb::init<MSequenceAlign::TYPE, ftype, ftype, ftype>(), nb::arg( "type" ) = MSequenceAlign::TYPE::GLOBAL,
            nb::arg( "match_score" ) = 1.0, nb::arg( "miss_score" ) = -0.5, nb::arg( "gap_score" ) = -1.0,
            "Constructor from alignment method, match, miss and gap scores." )
      .def( "__call__", &MSequenceAlign::operator(), nb::arg( "seq1" ), nb::arg( "seq2" ), "Align sequences." )
      .def( "__repr__", []( const MSequenceAlign &self ) { return "<clipper.MSequenceAlign class.>"; } )
      .doc() = "Sequence alignment object.\nProvides methods to "
               "find an optimal alignment between two sequences.";
}

void declare_seqfile( nb::module_ &m ) {
  nb::class_<SEQfile>( m, "SEQfile", "SEQ file object for MiniMol sequence i/o." )
      .def( nb::init<>() )
      .def( "read_file", &SEQfile::read_file, nb::arg( "file" ), "Load SEQ data from file." )
      .def( "import_polymer_sequence", &SEQfile::import_polymer_sequence, nb::arg( "target" ),
            "Read a single sequence from SEQ file." )
      .def( "import_molecule_sequence", &SEQfile::import_molecule_sequence, nb::arg( "target" ),
            "Read a molecule from the SEQ file." );
}

void init_minimol_seq( nb::module_ &m ) {
  declare_mpolymer_sequence( m );
  declare_mmolecule_sequence( m );
  declare_msequence_align( m );
  declare_seqfile( m );
}