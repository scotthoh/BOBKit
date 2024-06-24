
// #include <Python.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>
#include <string>
// #include <clipper/core/cell.h>
#include "helper_functions.h"
#include "type_conversions.h"
// #include <clipper/clipper.h>
#include <clipper/core/cell.h>
namespace py = pybind11;
using namespace clipper;

void declare_metric_tensor(py::module &m) {
  py::class_<Metric_tensor> metric_tensor(m, "Metric_tensor");
  metric_tensor.def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("alpha"),
           py::arg("beta"), py::arg("gamma"))
      .def("lengthsq", (ftype(Metric_tensor::*)(const Vec3<> &) const) &
                           Metric_tensor::lengthsq)
      .def("lengthsq", (ftype(Metric_tensor::*)(const Vec3<int> &) const) &
                           Metric_tensor::lengthsq)
      .def("format", &Metric_tensor::format)
      .def("__str__", [](const Metric_tensor &self) { return self.format(); });
}
void declare_cell_descr(py::module &m) {
  py::class_<Cell_descr> cell_descr(m, "Cell_descr");
  cell_descr.def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("alpha"),
           py::arg("beta"), py::arg("gamma"))
      .def(py::init([](const py::array_t<ftype> &params) {
        return std::unique_ptr<Cell_descr>(
            new Cell_descr(params.at(0), params.at(1), params.at(2),
                           params.at(3), params.at(4), params.at(5)));
      })) // 14Jun2024
      .def_property_readonly("a", &Cell_descr::a)
      .def_property_readonly("b", &Cell_descr::b)
      .def_property_readonly("c", &Cell_descr::c)
      .def_property_readonly("alpha", &Cell_descr::alpha)
      .def_property_readonly("beta", &Cell_descr::beta)
      .def_property_readonly("gamma", &Cell_descr::gamma)
      .def_property_readonly("alpha_deg", &Cell_descr::alpha_deg)
      .def_property_readonly("beta_deg", &Cell_descr::beta_deg)
      .def_property_readonly("gamma_deg", &Cell_descr::gamma_deg)
      .def("format", &Cell_descr::format)
      .def("__str__", [](const Cell_descr &self) { return self.format(); })
      .def("__repr__", [](const Cell_descr &self) {
        std::stringstream stream;
        stream << "<Cell(" << self.a() << ", " << self.b() << ", ";
        stream << self.c() << ", " << self.alpha_deg() << ", ";
        stream << self.beta_deg() << ", " << self.gamma_deg() << ")>";
        return stream.str();
      });
}

void declare_cell(py::module &m) {
  py::class_<Cell, Cell_descr> cell(m, "Cell");
  cell.def(py::init<>())
      .def(py::init<const Cell_descr &>(), py::arg("Cell_description"))
      .def("init", &Cell::init, py::arg("Cell_description"))
      .def_property_readonly("a_star", &Cell::a_star)
      .def_property_readonly("b_star", &Cell::b_star)
      .def_property_readonly("c_star", &Cell::c_star)
      .def_property_readonly("alpha_star", &Cell::alpha_star)
      .def_property_readonly("beta_star", &Cell::beta_star)
      .def_property_readonly("gamma_star", &Cell::gamma_star)
      .def_property_readonly("description", &Cell::descr)
      .def_property_readonly("volume", &Cell::volume)
      .def_property_readonly("debug", &Cell::debug)
      .def("equals", &Cell::equals, py::arg("cell"), py::arg("tol") = 1.0)
      .def_property_readonly("matrix_orth", &Cell::matrix_orth)
      .def_property_readonly("matrix_frac", &Cell::matrix_frac)
      .def_property_readonly("metric_real", &Cell::metric_real)
      .def_property_readonly("metric_reci", &Cell::metric_reci);
}

void init_cell(py::module &m) {
  declare_metric_tensor(m);
  declare_cell_descr(m);
  declare_cell(m);
}