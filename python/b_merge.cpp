// Wrapper for buccaneer-merge
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-merge.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void init_ca_merge(py::module &m) {

  py::class_<Ca_merge> camerge(m, "Ca_merge");
  camerge.def(py::init<double>(), py::arg("reliability") = 0.5)
      .def("__call__", &Ca_merge::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktarget"), py::arg("seq"), "Merge model.")
      .def_static("merge_mr", &Ca_merge::merge_mr, py::arg("mol"),
                  py::arg("mol_mr"), py::arg("sigcut"), py::arg("nseed"),
                  py::arg("mr_filter"), py::arg("mr_seed"),
                  "Merge with molecular replacement model.")
      .def("__repr__",
           [](const Ca_merge &self) { return "<buccaneer.Ca_merge class.>"; })
      .doc() = "Class for augmenting model with MERGE model.";
}