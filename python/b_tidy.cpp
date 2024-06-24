// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-tidy.h"

#include "type_conversions.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_model_tidy(py::module &m) {
  py::class_<ModelTidy> model_tidy(m, "ModelTidy");
  model_tidy
      .def(py::init<double, int, clipper::String, bool>(),
           py::arg("rmsd") = 1.0, py::arg("nmin") = 12,
           py::arg("newrestype") = "ALA", py::arg("verbose") = false)
      .def("tidy", &ModelTidy::tidy, py::arg("mol"), py::arg("mol_mr"),
           py::arg("seq"))
      .def_static("chain_renumber",
                  static_cast<std::vector<int> (*)(
                      clipper::MiniMol &, const clipper::MMoleculeSequence &)>(
                      &ModelTidy::chain_renumber),
                  py::arg("mol"), py::arg("seq"))
      .def_static("chain_assign", &ModelTidy::chain_assign, py::arg("mol"),
                  py::arg("mol_mr"), py::arg("seq_nums"), py::arg("rmsd"),
                  py::arg("nmin"))
      .def_static("chain_move", &ModelTidy::chain_move, py::arg("mol"),
                  py::arg("mol_mr"), py::arg("chnnums"))
      .def_static("sequence_correct", &ModelTidy::sequence_correct,
                  py::arg("mol"), py::arg("seq"), py::arg("seqnums"),
                  py::arg("newrestype"))
      .def_static("sequence_count", &ModelTidy::sequence_count, py::arg("mol"))
      // this returns clipper::array2d<int>
      .def_static("sequence_flags", &ModelTidy::sequence_flags, py::arg("mol"))
      .def_static("trim", &ModelTidy::trim, py::arg("mol"), py::arg("seq"))
      .def("__repr__", [](const ModelTidy &self) {
        return "<buccaneer.ModelTidy class.>";
      });
}