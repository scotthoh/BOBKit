// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <buccaneer/buccaneer-sequence.h>
#include <clipper/clipper.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void declare_ca_sequence(py::module &m) {
  py::class_<Ca_sequence> casequence(m, "Ca_sequence");
  casequence.def(py::init<double>(), py::arg("reliability") = 0.5)
      .def("__call__", &Ca_sequence::operator(), py::arg("mol"),
           py::arg("xmap"), py::arg("llktargets"), py::arg("seq"))
      .def("num_sequenced", &Ca_sequence::num_sequenced)
      .def("format", &Ca_sequence::format)
      // need to test these
      .def_static("phi_approx", &Ca_sequence::phi_approx, py::arg("z"))
      .def_static("prepare_score", &Ca_sequence::prepare_score, py::arg("mm"),
                  py::arg("xmap"), py::arg("llksample"))
      .def_static("prepare_scores", &Ca_sequence::prepare_scores, py::arg("mp"),
                  py::arg("xmap"), py::arg("llksample"))
      .def_static("sequence_overlap", &Ca_sequence::sequence_overlap,
                  py::arg("seq1"), py::arg("seq2"))
      .def_static("sequence_similarity", &Ca_sequence::sequence_similarity,
                  py::arg("seq1"), py::arg("seq2"))
      .def_static("sequence_combine", &Ca_sequence::sequence_combine,
                  py::arg("seq"), py::arg("reliability"))
      .def_static("sequence_score", &Ca_sequence::sequence_score,
                  py::arg("scores"), py::arg("subseq"))
      .def_static("sequence_align", &Ca_sequence::sequence_align,
                  py::arg("scores"), py::arg("seq"))
      .def_static("sequence_match", &Ca_sequence::sequence_match,
                  py::arg("scores"), py::arg("seq"))
      .def_static("sequence_chain", &Ca_sequence::sequence_chain,
                  py::arg("chain"), py::arg("seq"))
      .def_static("sequence_apply", &Ca_sequence::sequence_apply,
                  py::arg("chain"), py::arg("seq"), py::arg("flags"))
      .def_static("sequence", &Ca_sequence::sequence, py::arg("chain"),
                  py::arg("seq"), py::arg("reliability"))

      .def_static("set_semet", &Ca_sequence::set_semet, py::arg("semet"))
      .def_static("set_prior_model", &Ca_sequence::set_prior_model,
                  py::arg("mol"))
      .def_static("set_cpus", &Ca_sequence::set_cpus, py::arg("ncpu"))
      .def("__repr__", [](const Ca_sequence &self) {
        return "<buccaneer.Ca_sequence class.>";
      });

  using Class = Ca_sequence::Sequence_data;
  py::class_<Class>(casequence, "Sequence_data")
      .def(py::init<>())
      .def(py::init<const Ca_group &, const std::vector<double> &>(),
           py::arg("Ca_group"), py::arg("data"))
      .def_readonly("ca", &Class::ca)
      .def_readonly("data", &Class::data)
      .def("__repr__", [](const Class &self) {
        return "<buccaneer.Sequence_data class.>";
      });
}

void declare_sequence_score_threaded(py::module &m) {
  // inheritance public thread base
  py::class_<Sequence_score_threaded>(m, "Sequence_score_threaded")
      .def(py::init<>())
      .def(py::init<clipper::MPolymer &, const clipper::Xmap<float> &,
                    const std::vector<LLK_map_target::Sampled> &>(),
           py::arg("mp"), py::arg("xmap"), py::arg("llksample"))
      .def("sequence_score", &Sequence_score_threaded::sequence_score,
           py::arg("chn"))
      .def("result", &Sequence_score_threaded::result)
      .def("__call__", &Sequence_score_threaded::operator(),
           py::arg("nthread") = 0)
      .def("merge", &Sequence_score_threaded::merge, py::arg("other"))
      .def("__repr__", [](const Sequence_score_threaded &self) {
        return "<buccaneer.Sequence_score_threaded class.>";
      });
}
void declare_sequence_threaded(py::module &m) {

  py::class_<Sequence_threaded>(m, "Sequence_threaded")
      .def(py::init<>())
      .def(py::init<const clipper::MiniMol &,
                    const clipper::MMoleculeSequence &, const double &>(),
           py::arg("mol"), py::arg("seq"), py::arg("reliability"))
      .def("sequence", &Sequence_threaded::sequence, py::arg("chn"))
      .def("result", &Sequence_threaded::result)
      .def("history", &Sequence_threaded::history)
      .def("__call__", &Sequence_threaded::operator(), py::arg("nthread") = 0)
      .def("merge", &Sequence_threaded::merge, py::arg("other"))
      .def("__repr__", [](const Sequence_threaded &self) {
        return "<buccaneer.Sequence_threaded class.>";
      });
}
void init_ca_sequence(py::module &m) {
  declare_ca_sequence(m);
  declare_sequence_score_threaded(m);
  declare_sequence_threaded(m);
}