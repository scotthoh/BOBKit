#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <clipper/contrib/test_contrib.h>
#include <clipper/core/test_core.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

// Self test for clipper_core and clipper_contrib,
// useful to check if fftw3 compiled correctly
void init_clipper_tests(py::module &m) {
  py::class_<Test_core>(m, "Test_core")
      .def(py::init<>())
      .def("__call__", &Test_core::operator())
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_core::set_stream)
      .def("__repr__",
           [](const Test_core &self) { return "<clipper.Test_core class.>"; });

  py::class_<Test_contrib>(m, "Test_contrib")
      .def(py::init<>())
      .def("__call__", &Test_contrib::operator())
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_contrib::set_stream)
      .def("__repr__", [](const Test_contrib &self) {
        return "<clipper.Test_contrib class.>";
      });
  ;
}