// Nanobind bindings for buccaneer-known
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-known.h"
#include "commons.h"
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

using namespace clipper;

void add_knownstructure(nb::module_ &m) {
  nb::class_<KnownStructure>(m, "KnownStructure",
                             "Class for augmenting model with known model.")
      .def(nb::init<const MiniMol &,
                    const std::vector<std::pair<String, double>> &, double>(),
           "Constructor from model and arguments.")
      .def("copy_to", &KnownStructure::copy_to, nb::arg("mol"),
           nb::arg("includeAll") = true,
           "Add known structure to existing structure.")
      .def("clash", &KnownStructure::clash, nb::arg("coord"),
           "Check for clashes against known model")
      .def("prune", &KnownStructure::prune, nb::arg("mol"),
           "Prune model where is clashes with known model.")
      .def_static("parse", &KnownStructure::parse, nb::arg("arg"),
                  "Parse and store and input argument.")
      .def("debug", &KnownStructure::debug, "Print out selected atoms.");
}