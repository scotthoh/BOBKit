// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer-tidy.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "type_conversions.h"

namespace py = pybind11;

void init_model_tidy(py::module &m)
{
  py::class_<ModelTidy> model_tidy(m, "ModelTidy");
  model_tidy
      .def(py::init<double, int, clipper::String, bool>(),
           py::arg("rmsd") = 1.0, py::arg("nmin") = 12, py::arg("newrestype") = "ALA", py::arg("verbose") = false)
      .def("tidy", &ModelTidy::tidy)
      .def_static("chain_renumber", static_cast<std::vector<int> (*)(clipper::MiniMol &, const clipper::MMoleculeSequence &)>(&ModelTidy::chain_renumber))
      .def_static("chain_assign", &ModelTidy::chain_assign)
      .def_static("chain_move", &ModelTidy::chain_move)
      .def_static("sequence_correct", &ModelTidy::sequence_correct)
      .def_static("sequence_count", &ModelTidy::sequence_count)
      // this returns clipper::array2d<int>
      .def_static("sequence_flags", &ModelTidy::sequence_flags)
      .def_static("trim", &ModelTidy::trim);
}