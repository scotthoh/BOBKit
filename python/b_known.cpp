#include "buccaneer/buccaneer-known.h"
#include "type_conversions.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>

namespace py = pybind11;
using namespace clipper;

void init_knownstructure(py::module &m) {
  py::class_<KnownStructure>(m, "KnownStructure")
      //.def(py::init<const MiniMol &,
      //              const std::vector<std::pair<String, double>> &, double>())
      .def(py::init([](const MiniMol &mol, const py::list &l, double nprad) {
        std::vector<std::pair<String, double>> ids;
        for (auto it = l.begin(); it != l.end(); it += 2) {
          ids.emplace_back(std::pair<String, double>(it->cast<String>(),
                                                     (++it)->cast<double>()));
        }

        return std::unique_ptr<KnownStructure>(
            new KnownStructure(mol, ids, nprad));
      }))
      .def("copy_to", &KnownStructure::copy_to, py::arg("mol"),
           py::arg("includeAll") = true)
      .def("clash", &KnownStructure::clash, py::arg("coord"))
      .def("prune", &KnownStructure::prune, py::arg("mol"))
      .def_static("parse", &KnownStructure::parse, py::arg("arg"))
      .def("debug", &KnownStructure::debug);
}