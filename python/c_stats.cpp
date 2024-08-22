// Wrapper for clipper stats
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

template <class T> void declare_range(py::module &m, const std::string &name) {
  using Class = Range<T>;
  std::string PyClass = std::string("Range_") + name;
  py::class_<Class>(m, PyClass.c_str(),
                    "Range upper and lower bounds of some type.")
      .def(py::init<>())
      .def(py::init<const T &, const T &>(), py::arg("min"), py::arg("max"),
           "Constructor from min and max.")
      .def_property_readonly("min", &Class::min, "Return minimum value.")
      .def_property_readonly("max", &Class::max, "Return maximum value.")
      .def_property_readonly("range", &Class::range, "Return range = max-min.")
      .def("include", &Class::include, py::arg("datum"),
           "Update limits to include new datum.")
      .def("contains", &Class::include, py::arg("datum"),
           "Test if data is within limits (min<=datum<=max).")
      .def("truncate", &Class::truncate, py::arg("datum"),
           "Truncate data to be within range.")
      .def("__repr__", [=](const Class &self) {
        return "<clipper.Range_" + name + String(ftype(self.min()), 10, 4) +
               " - " + String(ftype(self.max()), 10, 4) + ">";
      });
}

void declare_range_sampling(py::module &m) {
  py::class_<Range_sampling, Range<ftype>>(
      m, "Range_sampling", "Range sampling: discrete sampling of a real range.")
      .def(py::init<>())
      .def(py::init<const int &>(), "Constructor from a number of samplings.")
      .def(py::init<const Range<ftype> &, const int &>(),
           "Constructor from range and number of samplings.")
      .def("indexf", &Range_sampling::indexf, py::arg("pos"),
           "Return fractional posn in counting range from x-value (0..n).")
      .def("x",
           (ftype(Range_sampling::*)(const ftype &) const) & Range_sampling::x,
           py::arg("frac_pos"),
           "Return x-value (0..n) from fractional posn in counting range.")
      .def("index", &Range_sampling::index, py::arg("pos"),
           "Return nearest index to particular x-value.")
      .def("index_bounded", &Range_sampling::index_bounded, py::arg("pos"),
           "Return nearest index to particular x-value (bounded 0...n-1)")
      .def("x",
           (ftype(Range_sampling::*)(const int &) const) & Range_sampling::x,
           py::arg("pos"),
           "Return x-value corresponding to centre of i'th range.")
      .def("x_min", &Range_sampling::x_min,
           "Return x-value corresponding to bottom of i'th range.")
      .def("x_max", &Range_sampling::x_max,
           "Return x-value corresponding to top of i'th range.")
      .def("size", &Range_sampling::size,
           "Return number of samplings in range.")
      .def("__repr__", [](const Range_sampling &self) {
        return "<clipper.Range_sampling class.>";
      });
}

void declare_histogram(py::module m) {
  py::class_<Histogram, Range_sampling>(m, "Histogram")
      .def(py::init<>())
      .def(py::init<const Range<ftype> &, const int &>(), py::arg("range"),
           py::arg("n"))
      .def("accumulate",
           (void(Histogram::*)(const ftype &)) & Histogram::accumulate,
           py::arg("x"))
      .def("accumulate",
           (void(Histogram::*)(const ftype &, const ftype &)) &
               Histogram::accumulate,
           py::arg("x"), py::arg("value"))
      .def("sum", &Histogram::sum)
      .def("y", (const ftype &(Histogram::*)(const int &) const) & Histogram::y,
           py::arg("index"))
      .def("y", (ftype(Histogram::*)(const ftype &) const) & Histogram::y,
           py::arg("x"))
      .def(
          "__iadd__",
          [](Histogram &self, const Histogram &other) { return self += other; },
          py::is_operator())
      .def("__repr__",
           [](const Histogram &self) {
             return "<clipper.Histogram with sample size " +
                    String(self.size());
           })
      .doc() =
      "General histogram class.\nThis class is used to accumulate and "
      "access a histogram of values spread over a specified range. On "
      "storing data or retrieving by interpolation the range is checked.";
}

void init_clipper_stats(py::module &m) {
  declare_range<int>(m, "int");
  declare_range<float>(m, "float");
  declare_range<double>(m, "double");
  declare_range_sampling(m);
  declare_histogram(m);
}