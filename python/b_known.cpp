// Wrapper for buccaneer-known
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-known.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>

#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void init_knownstructure(py::module &m) {
  py::class_<KnownStructure>(m, "KnownStructure",
                             "Class for augmenting model with known model.")
      .def(py::init<const MiniMol &,
                    const std::vector<std::pair<String, double>> &, double>(),
           "Constructor from model and arguments.")
      .def("copy_to", &KnownStructure::copy_to, py::arg("mol"),
           py::arg("includeAll") = true,
           "Add known structure to existing structure.")
      .def("clash", &KnownStructure::clash, py::arg("coord"),
           "Check for clashes against known model")
      .def("prune", &KnownStructure::prune, py::arg("mol"),
           "Prune model where is clashes with known model.")
      .def_static("parse", &KnownStructure::parse, py::arg("arg"),
                  "Parse and store and input argument.")
      .def("debug", &KnownStructure::debug, "Print out selected atoms.");
}