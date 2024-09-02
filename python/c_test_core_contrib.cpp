// Wrapper for clipper tests for core and contrib
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <clipper/contrib/test_contrib.h>
#include <clipper/core/test_core.h>
#include <clipper/minimol/test_minimol_gemmi.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

// Self test for clipper_core and clipper_contrib,
// useful to check if fftw3 compiled correctly
void init_clipper_tests(py::module &m) {
  py::class_<Test_core>(m, "Test_core", "Class test clipper core methods.")
      .def(py::init<>())
      .def("__call__", &Test_core::operator())
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_core::set_stream)
      .def("__repr__",
           [](const Test_core &self) { return "<clipper.Test_core class.>"; });

  py::class_<Test_contrib>(m, "Test_contrib",
                           "Class test clipper contrib methods.")
      .def(py::init<>())
      .def("__call__", &Test_contrib::operator())
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_contrib::set_stream)
      .def("__repr__", [](const Test_contrib &self) {
        return "<clipper.Test_contrib class.>";
      });

  py::class_<Test_minimol_gemmi>(m, "Test_minimol_gemmi",
                                 "Class test clipper minimol-gemmi methods.")
      .def(py::init<>())
      .def("run", &Test_minimol_gemmi::run, py::arg("input_file"))
      .def("__repr__", [](const Test_minimol_gemmi &self) {
        return "<clipper.Test_minimol_gemmi class.>";
      });
}