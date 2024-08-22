// Wrapper for buccaneer-tidy
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-tidy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void init_model_tidy(py::module &m) {
  py::class_<ModelTidy> model_tidy(m, "ModelTidy", "Class for tidying model.");
  model_tidy
      .def(py::init<double, int, String, bool>(), py::arg("rmsd") = 1.0,
           py::arg("nmin") = 12, py::arg("newrestype") = "ALA",
           py::arg("verbose") = false)
      .def("tidy", &ModelTidy::tidy, py::arg("mol"), py::arg("mol_mr"),
           py::arg("seq"), "Tidy model.")
      .def_static("chain_renumber",
                  static_cast<std::vector<int> (*)(MiniMol &,
                                                   const MMoleculeSequence &)>(
                      &ModelTidy::chain_renumber),
                  py::arg("mol"), py::arg("seq"), "Renumber residues in chain.")
      .def_static("chain_assign", &ModelTidy::chain_assign, py::arg("mol"),
                  py::arg("mol_mr"), py::arg("seq_nums"), py::arg("rmsd"),
                  py::arg("nmin"), "Assign chain fragments to MR model.")
      .def_static("chain_move", &ModelTidy::chain_move, py::arg("mol"),
                  py::arg("mol_mr"), py::arg("chnnums"))
      .def_static("sequence_correct", &ModelTidy::sequence_correct,
                  py::arg("mol"), py::arg("seq"), py::arg("seqnums"),
                  py::arg("newrestype"), "Correct the sequence.")
      .def_static("sequence_count", &ModelTidy::sequence_count, py::arg("mol"),
                  "Count sequence types.")
      // this returns clipper::array2d<int>
      .def_static("sequence_flags", &ModelTidy::sequence_flags, py::arg("mol"),
                  "Make flags of used residue numbers.")
      .def_static("trim", &ModelTidy::trim, py::arg("mol"), py::arg("seq"),
                  "Trim ends of chains which go beyind the end of the sequnce.")
      .def("__repr__", [](const ModelTidy &self) {
        return "<buccaneer.ModelTidy class.>";
      });
}