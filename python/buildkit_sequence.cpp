// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-sequence.h"
#include <clipper/clipper.h>
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;

void declare_ca_sequence(py::module &m)
{
  py::class_<Ca_sequence> casequence(m, "Ca_sequence");
  casequence
      .def(py::init<double>(), py::arg("reliability") = 0.5)
      .def("__call__", &Ca_sequence::operator())
      .def("num_sequenced", &Ca_sequence::num_sequenced)
      .def("format", &Ca_sequence::format)
      // need to test these
      .def_static("phi_approx", &Ca_sequence::phi_approx)
      .def_static("prepare_score", &Ca_sequence::prepare_score)
      .def_static("prepare_scores", &Ca_sequence::prepare_scores)
      .def_static("sequence_overlap", &Ca_sequence::sequence_overlap)
      .def_static("sequence_similarity", &Ca_sequence::sequence_similarity)
      .def_static("sequence_combine", &Ca_sequence::sequence_combine)
      .def_static("sequence_score", &Ca_sequence::sequence_score)
      .def_static("sequence_align", &Ca_sequence::sequence_align)
      .def_static("sequence_match", &Ca_sequence::sequence_match)
      .def_static("sequence_chain", &Ca_sequence::sequence_chain)
      .def_static("sequence_apply", &Ca_sequence::sequence_apply)
      .def_static("sequence", &Ca_sequence::sequence)

      .def_static("set_semet", &Ca_sequence::set_semet)
      .def_static("set_prior_model", &Ca_sequence::set_prior_model)
      .def_static("set_cpus", &Ca_sequence::set_cpus);

  using Class = Ca_sequence::Sequence_data;
  py::class_<Class>(casequence, "Sequence_data")
      .def(py::init<>())
      .def(py::init<const Ca_group &, const std::vector<double> &>())
      .def_readonly("ca", &Class::ca)
      // changes done with .append() in python side does not reflect
      .def_readonly("data", &Class::data);
}

void declare_sequence_score_threaded(py::module &m)
{
  // inheritance public thread base
  py::class_<Sequence_score_threaded>(m, "Sequence_score_threaded")
      .def(py::init<>())
      .def(py::init<clipper::MPolymer &, const clipper::Xmap<float> &, const std::vector<LLK_map_target::Sampled> &>())
      .def("sequence_score", &Sequence_score_threaded::sequence_score)
      .def("result", &Sequence_score_threaded::result)
      .def("__call__", &Sequence_score_threaded::operator())
      .def("merge", &Sequence_score_threaded::merge);
}
void declare_sequence_threaded(py::module &m)
{

  py::class_<Sequence_threaded>(m, "Sequence_threaded")
      .def(py::init<>())
      .def(py::init<const clipper::MiniMol &, const clipper::MMoleculeSequence &, const double &>())
      .def("sequence", &Sequence_threaded::sequence)
      .def("result", &Sequence_threaded::result)
      .def("history", &Sequence_threaded::history)
      .def("__call__", &Sequence_threaded::operator(), py::arg("nthread") = 0)
      .def("merge", &Sequence_threaded::merge);
}
void init_ca_sequence(py::module &m)
{
  declare_ca_sequence(m);
  declare_sequence_score_threaded(m);
  declare_sequence_threaded(m);
}