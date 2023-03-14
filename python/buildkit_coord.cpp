#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <clipper/clipper.h>
#include <type_conversions.h>
#include <helper_functions.h>

namespace py = pybind11;
using namespace clipper;
void init_coord_orth(py::module &m)
{
  py::class_<Coord_orth>(m, "Coord_orth")
      .def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &>(), py::arg("x"), py::arg("y"), py::arg("z"))
      .def_property_readonly("x", &Coord_orth::x)
      .def_property_readonly("y", &Coord_orth::y)
      .def_property_readonly("z", &Coord_orth::z)
      .def("__setitem__", [](Coord_orth &self, const int i, ftype coord)
           { return self[i] = coord; })
      .def_property_readonly("lengthsq", &Coord_orth::lengthsq)
      .def("__iter__", [](Coord_orth &self)
           { return py::make_iterator(&self[0], &self[2]); })
      .def("__str__", &Coord_orth::format)
      .def(
          "to_numpy", [](const Coord_orth &self)
          { return to_numpy_1d<Coord_orth, ftype>(self, 3); },
          py::return_value_policy::reference_internal);
}
