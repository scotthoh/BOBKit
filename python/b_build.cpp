// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

// #include "buildkit-join.h"
#include "buccaneer/buccaneer-build.h"

#include "type_conversions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_ca_build(py::module &m) {
  py::class_<Ca_build>(m, "Ca_build")
      .def(py::init<clipper::String, bool>(), py::arg("newrestype") = "ALA",
           py::arg("flexible") = false)
      .def_static("build", &Ca_build::build, py::arg("mol"), py::arg("xmap"),
                  py::arg("newrestype") = "ALA", py::arg("flexible") = false)
      .def("__call__", &Ca_build::operator(), py::arg("mol"), py::arg("xmap"))
      .def("__repr__", [](const Ca_build &self) {
        std::stringstream stream;
        stream << "<buccaneer.Ca_build class>";
        return stream.str();
      });
}