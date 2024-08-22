#include "type_conversions.h"
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void init_clipper_util(py::module &m) {
  py::class_<Util>(m, "Util", "Utility class.")
      .def(py::init<>())
      .def_static("u2b", &Util::u2b, py::arg("x"),
                  "Convert isotropic U-value to B-factor.")
      .def_static("b2u", &Util::b2u, py::arg("x"),
                  "Convert isotropic B-factor to U-value.")
      .def_static("is_nan", static_cast<bool (*)(const ftype32)>(&Util::is_nan),
                  py::arg("f"),
                  "Fast Util::nan() test. Used for missing entries: THIS DOES "
                  "NOT DISTINGUISH BETWEEN NAN & INF")
      .def_static("is_nan", static_cast<bool (*)(const ftype64)>(&Util::is_nan),
                  py::arg("f"),
                  "Fast Util::nan() test. Used for missing entries: THIS DOES "
                  "NOT DISTINGUISH BETWEEN NAN & INF")
      .def_static("d2rad", &Util::d2rad, py::arg("x"),
                  "degree-to-radian conversion.")
      .def_static("rad2d", &Util::rad2d, py::arg("x"),
                  "radian-to-degree conversion.")
      .def_static("pi", &Util::pi, "Return pi.")
      .def_static("twopi", &Util::twopi, "Return 2 pi.")
      .def_static("twopi2", &Util::twopi2, "Return two pi squared.")
      .def_static("eightpi2", &Util::eightpi2, "Return 8 pi squared.");
}