#include "type_conversions.h"
#include <clipper/clipper-minimol.h>
#include <clipper/minimol/minimol_seq.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_mpolymer_sequence(py::module &m) {
  py::class_<MPolymerSequence>(m, "MPolymerSequence")
      .def(py::init<>())
      .def_property("id", &MPolymerSequence::id, &MPolymerSequence::set_id)
      .def_property("sequence", &MPolymerSequence::sequence,
                    &MPolymerSequence::set_sequence)
      .def_static("id_tidy", &MPolymerSequence::id_tidy, py::arg("id"))
      .def_static("id_match", &MPolymerSequence::id_match, py::arg("id1"),
                  py::arg("id2"), py::arg("mode"))
      .def("__repr__",
           [](const MPolymerSequence self) {
             return "<clipper.MPolymerSequence " + self.id() + self.sequence() +
                    ">";
           })
      .def("__str__", [](const MPolymerSequence self) {
        return self.id() + self.sequence();
      });
}

void declare_mmolecule_sequence(py::module &m) {
  py::class_<MMoleculeSequence>(m, "MMoleculeSequence")
      .def(py::init<>())
      .def("size", &MMoleculeSequence::size)
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
      .def("__getitem__",
           [](const MMoleculeSequence &self, const int &i) { return self[i]; })
      .def("__setitem__", [](MMoleculeSequence &self, const int &i,
                             MPolymerSequence &value) { self[i] = value; })
      //  .def("find_sequence",
      //      [](const MMoleculeSequence &self, const String &n, const MM::MODE
      //      mode){
      // return self.find(n, mode);
      //      })
      //  .def("set_sequence", [](MMolecule &self, const String &n , const
      //  MPolymerSequence &seq, const MM::MODE mode){
      // self.find(n, mode) = seq;
      //  })

      .def_property(
          "find",
          (const MPolymerSequence &(
              MMoleculeSequence::*)(const String &, const MM::MODE) const) &
              MMoleculeSequence::find,
          (MPolymerSequence &
           (MMoleculeSequence::*)(const String &, const MM::MODE)) &
              MMoleculeSequence::find)
      .def("lookup", &MMoleculeSequence::lookup, py::arg("id"), py::arg("mode"))
      //.def("insert", &MMoleculeSequence::insert, py::arg("seq"),
      //     py::arg("pos") = -1)
      .def(
          "insert",
          [](MMoleculeSequence &self, const MPolymerSequence &id,
             const int &pos) { self.insert(id, pos); },
          py::arg("seq"), py::arg("pos") = -1)
      .def(
          "insert",
          [](MMoleculeSequence &self, const String &id, const String &seq,
             const int &pos) {
            MPolymerSequence mpseq;
            mpseq.set_id(id);
            mpseq.set_sequence(seq);
            self.insert(mpseq, pos);
          },
          py::arg("id"), py::arg("seq"), py::arg("pos") = -1)
      .def("is_null", &MMoleculeSequence::is_null);
}

void declare_msequence_align(py::module &m) {
  py::class_<MSequenceAlign> mseqalign(m, "MSequenceAlign");
  py::enum_<MSequenceAlign::TYPE>(mseqalign, "TYPE")
      .value("GLOBAL", MSequenceAlign::TYPE::GLOBAL)
      .value("LOCAL", MSequenceAlign::TYPE::LOCAL);

  mseqalign
      .def(py::init<MSequenceAlign::TYPE, ftype, ftype, ftype>(),
           py::arg("type") = MSequenceAlign::TYPE::GLOBAL,
           py::arg("match_score") = 1.0, py::arg("miss_score") = -0.5,
           py::arg("gap_score") = -1.0)
      .def("__call__", &MSequenceAlign::operator(), py::arg("seq1"),
           py::arg("seq2"))
      .def("__repr__", [](const MSequenceAlign &self) {
        return "<clipper.MSequenceAlign class.>";
      });
}

void declare_seqfile(py::module &m) {
  py::class_<SEQfile>(m, "SEQfile")
      .def(py::init<>())
      .def("read_file", &SEQfile::read_file, py::arg("file"))
      .def("import_polymer_sequence", &SEQfile::import_polymer_sequence,
           py::arg("target"))
      .def("import_molecule_sequence", &SEQfile::import_molecule_sequence,
           py::arg("target"));
}

void init_minimol_seq(py::module &m) {
  declare_mpolymer_sequence(m);
  declare_mmolecule_sequence(m);
  declare_msequence_align(m);
  declare_seqfile(m);
}