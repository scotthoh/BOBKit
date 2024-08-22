// Wrapper for clipper minimol sequence
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper-minimol.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_mpolymer_sequence(py::module &m) {
  py::class_<MPolymerSequence>(m, "MPolymerSequence")
      .def(py::init<>())
      .def_property("id", &MPolymerSequence::id, &MPolymerSequence::set_id,
                    "Get/set sequence id.")
      .def_property("sequence", &MPolymerSequence::sequence,
                    &MPolymerSequence::set_sequence, "Get/set sequence.")
      .def_static("id_tidy", &MPolymerSequence::id_tidy, py::arg("id"),
                  "Convert id to standard format.")
      .def_static("id_match", &MPolymerSequence::id_match, py::arg("id1"),
                  py::arg("id2"), py::arg("mode"), "Compare two ids.")
      .def("__repr__",
           [](const MPolymerSequence self) {
             return "<clipper.MPolymerSequence " + self.id() + self.sequence() +
                    ">";
           })
      .def("__str__",
           [](const MPolymerSequence self) {
             return self.id() + self.sequence();
           })
      .doc() = "Polymer sequence object.\nThe polymer sequence object "
               "represents the named sequence of a single chain.";
}

void declare_mmolecule_sequence(py::module &m) {
  py::class_<MMoleculeSequence>(m, "MMoleculeSequence")
      .def(py::init<>())
      .def("size", &MMoleculeSequence::size,
           "Return number of polymer sequences in model.")
      .def("__len__", &MMoleculeSequence::size)
      .def("__repr__",
           [](const MMoleculeSequence &self) {
             return "<clipper.MMoleculeSequence of length " +
                    String(self.size()) + " >";
           })
      .def("__str__",
           [](const MMoleculeSequence &self) {
             std::string s = "";
             for (int i = 0; i < self.size(); i++) {
               s += self[i].sequence();
               s += "\n";
             }
             return s;
           })
      .def(
          "__getitem__",
          [](const MMoleculeSequence &self, const int &i) { return self[i]; },
          "Get polymer sequence.")
      .def(
          "__setitem__",
          [](MMoleculeSequence &self, const int &i, MPolymerSequence &value) {
            self[i] = value;
          },
          "Set polymer sequence.")
      .def_property(
          "find",
          (const MPolymerSequence &(
              MMoleculeSequence::*)(const String &, const MM::MODE) const) &
              MMoleculeSequence::find,
          (MPolymerSequence &
           (MMoleculeSequence::*)(const String &, const MM::MODE)) &
              MMoleculeSequence::find,
          "Get/set polymer sequence by id.")
      .def("lookup", &MMoleculeSequence::lookup, py::arg("id"), py::arg("mode"),
           "Lookup polymer sequence by id.")
      //.def("insert", &MMoleculeSequence::insert, py::arg("seq"),
      //     py::arg("pos") = -1)
      .def(
          "insert",
          [](MMoleculeSequence &self, const MPolymerSequence &id,
             const int &pos) { self.insert(id, pos); },
          py::arg("seq"), py::arg("pos") = -1, "Add polymer sequence.")
      .def(
          "insert",
          [](MMoleculeSequence &self, const String &id, const String &seq,
             const int &pos) {
            MPolymerSequence mpseq;
            mpseq.set_id(id);
            mpseq.set_sequence(seq);
            self.insert(mpseq, pos);
          },
          py::arg("id"), py::arg("seq"), py::arg("pos") = -1,
          "Add polymer sequence.")
      .def("is_null", &MMoleculeSequence::is_null, "Test for null model.")
      .doc() = "Molecule sequence object.\nThe molecule sequence "
               "object is a list of polymer sequence objects representing "
               "the named sequences of all the chains in a molecule.";
}

void declare_msequence_align(py::module &m) {
  py::class_<MSequenceAlign> mseqalign(m, "MSequenceAlign");
  py::enum_<MSequenceAlign::TYPE>(mseqalign, "TYPE", "Alignment method.")
      .value("GLOBAL", MSequenceAlign::TYPE::GLOBAL)
      .value("LOCAL", MSequenceAlign::TYPE::LOCAL)
      .export_values();

  mseqalign
      .def(py::init<MSequenceAlign::TYPE, ftype, ftype, ftype>(),
           py::arg("type") = MSequenceAlign::TYPE::GLOBAL,
           py::arg("match_score") = 1.0, py::arg("miss_score") = -0.5,
           py::arg("gap_score") = -1.0,
           "Constructor from alignment method, match, miss and gap scores.")
      .def("__call__", &MSequenceAlign::operator(), py::arg("seq1"),
           py::arg("seq2"), "Align sequences.")
      .def("__repr__",
           [](const MSequenceAlign &self) {
             return "<clipper.MSequenceAlign class.>";
           })
      .doc() = "Sequence alignment object.\nProvides methods to "
               "find an optimal alignment between two sequences.";
}

void declare_seqfile(py::module &m) {
  py::class_<SEQfile>(m, "SEQfile", "SEQ file object for MiniMol sequence i/o.")
      .def(py::init<>())
      .def("read_file", &SEQfile::read_file, py::arg("file"),
           "Load SEQ daat from file.")
      .def("import_polymer_sequence", &SEQfile::import_polymer_sequence,
           py::arg("target"), "Read a single sequence from SEQ file.")
      .def("import_molecule_sequence", &SEQfile::import_molecule_sequence,
           py::arg("target"), "Read a molecule from the SEQ file.");
}

void init_minimol_seq(py::module &m) {
  declare_mpolymer_sequence(m);
  declare_mmolecule_sequence(m);
  declare_msequence_align(m);
  declare_seqfile(m);
}