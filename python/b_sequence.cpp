// Wrapper for buccaneer-sequence
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-sequence.h>
#include <clipper/clipper.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_ca_sequence(py::module &m) {
  py::class_<Ca_sequence> casequence(
      m, "Ca_sequence", "Class for sequence Ca chains using density");
  casequence
      .def(py::init<double>(), py::arg("reliability") = 0.5,
           "Constructor with sequence reliability.")
      .def("__call__", &Ca_sequence::operator(), py::arg("mol"),
           py::arg("xmap"), py::arg("llktargets"), py::arg("seq"),
           "Sequence Ca chains using density.")
      .def("num_sequenced", &Ca_sequence::num_sequenced,
           "Return number of C-alphas sequenced.")
      .def("format", &Ca_sequence::format,
           "Return sequencing results with up to top 5 scores for each chain "
           "as string.")
      // need to test these
      .def_static("phi_approx", &Ca_sequence::phi_approx, py::arg("z"),
                  "Approximate cumulative normal distribution function.")
      .def_static("prepare_score", &Ca_sequence::prepare_score, py::arg("mm"),
                  py::arg("xmap"), py::arg("llksample"),
                  "Cache scores in residue properties.")
      .def_static("prepare_scores", &Ca_sequence::prepare_scores, py::arg("mp"),
                  py::arg("xmap"), py::arg("llksample"),
                  "Cache scores in residue properties.")
      .def_static(
          "sequence_overlap", &Ca_sequence::sequence_overlap, py::arg("seq1"),
          py::arg("seq2"),
          "Return fraction of first sequence which overlaps the second.")
      .def_static("sequence_similarity", &Ca_sequence::sequence_similarity,
                  py::arg("seq1"), py::arg("seq2"),
                  "Return fraction of sequenced residues which match")
      .def_static("sequence_combine", &Ca_sequence::sequence_combine,
                  py::arg("seq"), py::arg("reliability"),
                  "Combine multiple non-conflicting sequence alignments.")
      .def_static("sequence_score", &Ca_sequence::sequence_score,
                  py::arg("scores"), py::arg("subseq"),
                  "Return highest scoring subsequence matching the given "
                  "sequence to the supplied LLK scores.")
      .def_static("sequence_align", &Ca_sequence::sequence_align,
                  py::arg("scores"), py::arg("seq"),
                  "Perform sequence alignment between the given chain scores "
                  "and the given sequence.")
      .def_static("sequence_match", &Ca_sequence::sequence_match,
                  py::arg("scores"), py::arg("seq"),
                  "Return a scored list of sequence matches between a set of "
                  "LLK scores and available sequence ranges.")
      .def_static("sequence_chain", &Ca_sequence::sequence_chain,
                  py::arg("chain"), py::arg("seq"),
                  "Sequence a chain based on the map LLK target, and available "
                  "sequence ranges.")
      .def_static("sequence_apply", &Ca_sequence::sequence_apply,
                  py::arg("chain"), py::arg("seq"), py::arg("flags"),
                  "Apply matching sequence to chain, taking account any "
                  "existing sequence.")
      .def_static("sequence", &Ca_sequence::sequence, py::arg("chain"),
                  py::arg("seq"), py::arg("reliability"),
                  "Run sequencing, combine sequence, and apply matching "
                  "sequence to chain.")

      .def_static("set_semet", &Ca_sequence::set_semet, py::arg("semet"),
                  "Set flag to translate MET to MSE.")
      .def_static("set_prior_model", &Ca_sequence::set_prior_model,
                  py::arg("mol"), "Set prior model.")
      .def_static("set_cpus", &Ca_sequence::set_cpus, py::arg("ncpu"),
                  "Set number of cpu threads to use.")
      .def("__repr__", [](const Ca_sequence &self) {
        return "<buccaneer.Ca_sequence class.>";
      });

  using Class = Ca_sequence::Sequence_data;
  py::class_<Class>(casequence, "Sequence_data",
                    "Class to hold sequence data with scores.")
      .def(py::init<>())
      .def(py::init<const Ca_group &, const std::vector<double> &>(),
           py::arg("Ca_group"), py::arg("data"),
           "Constructor with Ca_group and list of scores.")
      .def_readonly("ca", &Class::ca, "Get Ca_group.")
      .def_readonly("data", &Class::data, "Get list of scores.")
      .def("__repr__", [](const Class &self) {
        return "<buccaneer.Sequence_data class.>";
      });
}

void declare_sequence_score_threaded(py::module &m) {
  // inheritance public thread base
  py::class_<Sequence_score_threaded>(
      m, "Sequence_score_threaded",
      "Class with threads methods for sequencing for Ca groups.")
      .def(py::init<>())
      .def(py::init<MPolymer &, const Xmap<float> &,
                    const std::vector<LLK_map_target::Sampled> &>(),
           py::arg("mp"), py::arg("xmap"), py::arg("llksample"),
           "Constructor with chain, density map and list of LLK sampled.")
      .def("sequence_score", &Sequence_score_threaded::sequence_score,
           py::arg("chn"), "Runs prepare score.")
      .def("result", &Sequence_score_threaded::result,
           "Return sequenced chain.")
      .def("__call__", &Sequence_score_threaded::operator(),
           py::arg("nthread") = 0, "Run single or multi-threaded.")
      .def("merge", &Sequence_score_threaded::merge, py::arg("other"),
           "Merge results from multiple threads.")
      .def("__repr__", [](const Sequence_score_threaded &self) {
        return "<buccaneer.Sequence_score_threaded class.>";
      });
}
void declare_sequence_threaded(py::module &m) {

  py::class_<Sequence_threaded>(
      m, "Sequence_threaded",
      "Class with threads methods for sequencing for Ca groups.")
      .def(py::init<>())
      .def(py::init<const MiniMol &, const MMoleculeSequence &,
                    const double &>(),
           py::arg("mol"), py::arg("seq"), py::arg("reliability"),
           "Constructor with model, sequence and sequence reliability.")
      .def("sequence", &Sequence_threaded::sequence, py::arg("ichn"),
           "Run sequencing on chain.")
      .def("result", &Sequence_threaded::result,
           "Return result sequenced model.")
      .def("history", &Sequence_threaded::history,
           "Return a list of Score_list.")
      .def("__call__", &Sequence_threaded::operator(), py::arg("nthread") = 0,
           "Run single of multi-threaded.")
      .def("merge", &Sequence_threaded::merge, py::arg("other"),
           "Merge results from multiple threads.")
      .def("__repr__", [](const Sequence_threaded &self) {
        return "<buccaneer.Sequence_threaded class.>";
      });
}
void init_ca_sequence(py::module &m) {
  declare_ca_sequence(m);
  declare_sequence_score_threaded(m);
  declare_sequence_threaded(m);
}