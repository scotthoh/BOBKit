#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <clipper/core/cell.h>

namespace py = pybind11;
using namespace clipper;

void init_cell(py::module &m)
{
  py::class_<Cell_descr>(m, "Cell_descr")
      .def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &, const ftype &, const ftype &>(),
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("alpha"), py::arg("beta"), py::arg("gamma"))
      .def_property_readonly("a", &Cell_descr::a)
      .def_property_readonly("b", &Cell_descr::b)
      .def_property_readonly("c", &Cell_descr::c)
      .def_property_readonly("alpha", &Cell_descr::alpha)
      .def_property_readonly("beta", &Cell_descr::alpha)
      .def_property_readonly("gamma", &Cell_descr::alpha)
      .def_property_readonly("alpha_deg", &Cell_descr::alpha_deg)
      .def_property_readonly("beta_deg", &Cell_descr::beta_deg)
      .def_property_readonly("gamma_deg", &Cell_descr::gamma_deg);

  py::class_<Cell>(m, "Cell")
      .def(py::init<>())
      .def(py::init<const Cell_descr &>(), py::arg("Cell_description"))
      .def("init", &Cell::init, py::arg("Cell_description"))
      .def_property_readonly("a_star", &Cell::a_star)
      .def_property_readonly("b_star", &Cell::a_star)
      .def_property_readonly("c_star", &Cell::a_star)
      .def_property_readonly("alpha_star", &Cell::alpha_star)
      .def_property_readonly("beta_star", &Cell::beta_star)
      .def_property_readonly("gamma_star", &Cell::gamma_star)
      .def_property_readonly("description", &Cell::descr)
      .def_property_readonly("volume", &Cell::volume)
      .def_property_readonly("debug", &Cell::debug)
      .def_property_readonly("is_null", &MiniMol::is_null)
      .def_property_reaconly("descr", &Cell::descr)
      .def_property_readonly("volume", &Cell::&volume);
}