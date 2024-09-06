// Wrapper for clipper cell
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <sstream>
#include <string>

#include <gemmi/unitcell.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <clipper/clipper-gemmi.h>
#include <clipper/core/cell.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_metric_tensor(py::module &m) {
  py::class_<Metric_tensor> metric_tensor(m, "Metric_tensor");
  metric_tensor.def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("alpha"),
           py::arg("beta"), py::arg("gamma"),
           "Constructor: takes parameters of normal or inverse cell.")
      .def("lengthsq",
           (ftype(Metric_tensor::*)(const Vec3<> &) const) &
               Metric_tensor::lengthsq,
           "Apply metric to vector.")
      .def("lengthsq",
           (ftype(Metric_tensor::*)(const Vec3<int> &) const) &
               Metric_tensor::lengthsq,
           "Apply metric to int vector.")
      .def("format", &Metric_tensor::format,
           "Return formatted string representation.")
      .def("__str__", [](const Metric_tensor &self) { return self.format(); })
      .doc() =
      "The metric tensor is used to determine a distance in real or reciprocal"
      "space using fraction coordinates or Miller indices. It is symmetrical, "
      "so only the upper triangle is stored with the off-diagonal elements "
      "doubled.";
}
void declare_cell_descr(py::module &m) {
  py::class_<Cell_descr> cell_descr(m, "Cell_descr");
  cell_descr.def(py::init<>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("a"), py::arg("b"), py::arg("c"), py::arg("alpha"),
           py::arg("beta"), py::arg("gamma"),
           "Constructor from cell parameters.")
      .def(py::init([](const py::array_t<ftype> &params) {
             return std::unique_ptr<Cell_descr>(
                 new Cell_descr(params.at(0), params.at(1), params.at(2),
                                params.at(3), params.at(4), params.at(5)));
           }),
           "Constructor from cell parameters (list/numpy array).")
      .def_property_readonly("a", &Cell_descr::a, "Get a.")
      .def_property_readonly("b", &Cell_descr::b, "Get b.")
      .def_property_readonly("c", &Cell_descr::c, "Get c.")
      .def_property_readonly("alpha", &Cell_descr::alpha, "Get alpha.")
      .def_property_readonly("beta", &Cell_descr::beta, "Get beta.")
      .def_property_readonly("gamma", &Cell_descr::gamma, "Get gamma.")
      .def_property_readonly("alpha_deg", &Cell_descr::alpha_deg,
                             "Get alpha in degrees.")
      .def_property_readonly("beta_deg", &Cell_descr::beta_deg,
                             "Get beta in degrees.")
      .def_property_readonly("gamma_deg", &Cell_descr::gamma_deg,
                             "Get gamma in degrees.")
      .def("format", &Cell_descr::format,
           "Return formatted string representation.")
      .def("__str__", [](const Cell_descr &self) { return self.format(); })
      .def("__repr__",
           [](const Cell_descr &self) {
             std::stringstream stream;
             stream << "<Cell(" << self.a() << ", " << self.b() << ", ";
             stream << self.c() << ", " << self.alpha_deg() << ", ";
             stream << self.beta_deg() << ", " << self.gamma_deg() << ")>";
             return stream.str();
           })
      .doc() = "Cell description (automatically converts to radians)"
               "The cell description is a compact description of a cell, "
               "containing just the cell parameters. It is usually used to "
               "construct a full Cell object, which provides the expected "
               "functionality.";
}

void declare_cell(py::module &m) {
  py::class_<Cell, Cell_descr> cell(m, "Cell");
  cell.def(py::init<>(), "Null constructor, must initialise later.")
      .def(py::init<const Cell_descr &>(), py::arg("Cell_description"),
           "Constructor from Cell descriptor")
      //.def("init", 8&Cell::init, py::arg("Cell_description"))
      .def(
          "init", [](Cell &self, const Cell_descr &cd) { self.init(cd); },
          py::arg("Cell_description"), "Initialise with Cell descriptor.")
      .def(
          "init",
          [](Cell &self, const gemmi::UnitCell &c) {
            self.init(GEMMI::cell(c).descr());
          },
          py::arg("cell"), "Initialise with GEMMI UnitCell.")
      .def(
          "init", [](Cell &self, const Cell &c) { self.init(c.descr()); },
          py::arg("cell"), "Initialise with Cell descriptor.")

      .def_static(
          "to_gemmi_cell", [](const Cell &c) { return GEMMI::cell(c); },
          "Convert clipper Cell to GEMMI Unitcell.")
      .def_static(
          "from_gemmi_cell",
          [](const gemmi::UnitCell &c) { return GEMMI::cell(c); },
          "Convert GEMMI Unitcell from clipper Cell.")
      .def_property_readonly("a_star", &Cell::a_star, "Get a*")
      .def_property_readonly("b_star", &Cell::b_star, "Get b*")
      .def_property_readonly("c_star", &Cell::c_star, "Get c*")
      .def_property_readonly("alpha_star", &Cell::alpha_star, "Get alpha*")
      .def_property_readonly("beta_star", &Cell::beta_star, "Get beta*")
      .def_property_readonly("gamma_star", &Cell::gamma_star, "Get gamma*")
      .def_property_readonly("description", &Cell::descr,
                             "Return cell dimensions.")
      .def_property_readonly("volume", &Cell::volume, "Retuen cell volume.")
      .def_property_readonly("debug", &Cell::debug, "Output class details.")
      .def("equals", &Cell::equals, py::arg("cell"), py::arg("tol") = 1.0,
           "Test equality with another cell.")
      .def_property_readonly("matrix_orth", &Cell::matrix_orth,
                             "Return orthogonalisation matrix")
      .def_property_readonly("matrix_frac", &Cell::matrix_frac,
                             "Return fractionalisation matrix")
      .def_property_readonly("metric_real", &Cell::metric_real,
                             "Return real space metric tensor.")
      .def_property_readonly("metric_reci", &Cell::metric_reci,
                             "Return reciprocal space metric tensor.")
      .doc() = "Cell object.\n"
               "The Cell class is the fully functional description of the unit "
               "cell. In addition to the cell parameters, it stores derived "
               "information including the cell volume, orthogonalising and "
               "fractionalising matrices, and the metric tensors.";
}

void init_cell(py::module &m) {
  declare_metric_tensor(m);
  declare_cell_descr(m);
  declare_cell(m);
}