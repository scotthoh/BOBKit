#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void init_clipper_util(py::module &m) {
  py::class_<Util>(m, "Util")
      .def(py::init<>())
      .def_static("u2b", Util::u2b, py::arg("x"))
      .def_static("b2u", Util::b2u, py::arg("x"))
      .def_static("is_nan", static_cast<bool (*)(const ftype32)>(&Util::is_nan))
      .def_static("is_nan",
                  static_cast<bool (*)(const ftype64)>(&Util::is_nan));
}