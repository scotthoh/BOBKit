#include "commons.h"
#include "arrays.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

template <class T> void declare_grad_classes(nb::module_ &m, const std::string &name ) {
  using GradOrth = Grad_orth<T>;
  using VecClass = Vec3<T>;
  std::string PyClass = std::string( "Grad_orth_" ) + name;
  nb::class_<GradOrth, VecClass>( m, PyClass.c_str() )
    .def( nb::init<>(), "Null constructor" )
    .def( nb::init<const VecClass& >(), nb::arg("vector"), "Constructor: copy/convert")
    .def( nb::init<const T&, const T&, const T&>(), nb::arg("dx"), nb::arg("dy"), nb::arg("dz"), "Constructor: from d/dx, d/dy, d/dz")
    .def( "dx", &GradOrth::dx, "Get d/dx" )
    .def( "dy", &GradOrth::dx, "Get d/dy" )
    .def( "dz", &GradOrth::dx, "Get d/dz" )
    .def( "grad_frac", [](const GradOrth& self, const Cell& cell) {
      return self.grad_frac(cell);
    }, nb::arg("cell"), "Orthogonal-fractional derivative conversion")
    .def( "grad_frac", [](const GradOrth& self, const float& a, const float& b, const float& c, const float& alpha, const float& beta, const float& gamma ) {
      Cell cell(Cell_descr(a, b, c, alpha, beta, gamma));
      return self.grad_frac(cell);
    }, nb::arg("a"),nb::arg("b"),nb::arg("c"),nb::arg("alpha"),nb::arg("beta"),nb::arg("gamma"), "Orthogonal-fractional derivative conversion")
    .def( "format", &GradOrth::format, "Return formatted string representation.")
    .doc() = "Orthogonal (Angstom) gradient, with respect to orthogonal x,y,z";

  using GradFrac = Grad_frac<T>;
  PyClass = std::string("Grad_frac_") + name;
  nb::class_<GradFrac, VecClass>( m, PyClass.c_str() )
    .def( nb::init<>(), "Null constructor" )
    .def( nb::init< const VecClass&>(), nb::arg("vector"), "Constructor: copy/convert")
    .def( nb::init< const T&, const T&, const T&>(), nb::arg("du"), nb::arg("dv"), nb::arg("dw"), "Constructor: from d/du, d/dv, d/dw")
    .def( "du", &GradFrac::du, "Get d/du" )
    .def( "dv", &GradFrac::dv, "Get d/dv" )
    .def( "dw", &GradFrac::dw, "Get d/dw" )
    .def( "grad_orth", [](const GradFrac& self, const Cell& cell) {
      return self.grad_orth(cell);
    }, nb::arg("cell"), "Fractional-orthogonal derivative conversion")
    .def( "grad_orth", [](const GradFrac& self, const float& a, const float& b, const float& c, const float& alpha, const float& beta, const float& gamma ) {
      Cell cell(Cell_descr(a, b, c, alpha, beta, gamma));
      return self.grad_orth(cell);
    }, nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("alpha"), nb::arg("beta"), nb::arg("gamma"), "Fractional-orthogonal derivative conversion")
    .def( "grad_map", [](const GradFrac& self, const Grid& grid) {
      return self.grad_map(grid);
    }, nb::arg("cell"), "Fractional-grid derivative conversion")
    .def( "grad_map", [](const GradFrac& self, const int& nu, const int& nv, const int& nw ) {
      return self.grad_map(Grid(nu, nv, nw));
    }, nb::arg("nu"), nb::arg("nv"), nb::arg("nw"), "Fractional-grid derivative conversion")
    .def( "format", &GradFrac::format, "Return formatted string representation.")
    .doc() = "Fractional (cell) gradient, with respect to fractional u,v,w";
  
    using GradMap = Grad_map<T>;
    PyClass = std::string("Grad_map_") + name;
    nb::class_<GradMap, VecClass>( m, PyClass.c_str() )
      .def( nb::init<>(), "Null constructor" )
      .def( nb::init< const VecClass& >(), nb::arg("vector"), "Constructor: copy/convert")
      .def( nb::init< const T&, const T&, const T&>(), nb::arg("du"), nb::arg("dv"), nb::arg("dw"), "Constructor: from d/du, d/dv, d/dw")
      .def( "du", &GradMap::du, "Get d/du" )
      .def( "dv", &GradMap::dv, "Get d/dv" )
      .def( "dw", &GradMap::dw, "Get d/dw" )
      .def( "grad_frac", [](const GradMap& self, const Grid& grid) {
        return self.grad_frac(grid);
      }, nb::arg("grid"), "Grid-fractional derivative conversion")
      .def( "grad_frac", [](const GradMap& self, const int& nu, const int& nv, const int& nw) {
        return self.grad_frac(Grid(nu, nv, nw));
      }, nb::arg("nu"), nb::arg("nv"), nb::arg("nw"), "Grid-fractional derivative conversion")
      .def( "format", &GradMap::format, "Return formatted string representation")
      .doc() = "Map coordinate gradient, with respect to grid u,v,w";
}

template <class T> void declare_curv_classes(nb::module_ &m, const std::string &name ) {
  using CurvOrth = Curv_orth<T>;
  using MatClass = Mat33<T>;
  std::string PyClass = std::string( "Curv_orth_" ) + name;
  nb::class_<CurvOrth, MatClass>( m, PyClass.c_str() )
    .def( nb::init<>(), "Null constructor")
    .def( nb::init< const MatClass&>(), nb::arg("matrix"), "Constructor: copy/convert" )
    .def( "curv_frac", [](const CurvOrth& self, const Cell& cell) {
      return self.curv_frac(cell);
    }, nb::arg("cell"), "Orthogonal-fractional derivative conversion")
    .def( "curv_frac", [](const CurvOrth& self, const float& a, const float& b, const float& c, const float& alpha, const float& beta, const float& gamma) {
      Cell cell(Cell_descr(a, b, c, alpha, beta, gamma));
      return self.curv_frac(cell);
    }, nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("alpha"), nb::arg("beta"), nb::arg("gamma"), "Orthogonal-fractional derivative conversion")
    .doc() = "Orthogonal (Angstrom) curvatures, with respect to orthogonal x,y,z";
  
  using CurvFrac = Curv_frac<T>;
  PyClass = std::string( "Curv_frac_" ) + name;
  nb::class_<CurvFrac, MatClass>( m, PyClass.c_str() )
    .def( nb::init<>(), "Null constructor" )
    .def( nb::init<const MatClass&>(), nb::arg("matrix"), "Constructor: copy/convert" )
    .def( "curv_orth", [](const CurvFrac& self, const Cell& cell) {
      return self.curv_orth(cell);
    }, nb::arg("cell"), "Fractional-orthogonal derivative conversion")
    .def( "curv_orth", [](const CurvFrac &self, const float& a, const float& b, const float& c, const float& alpha, const float& beta, const float& gamma) {
      Cell cell(Cell_descr(a, b, c, alpha, beta, gamma));
      return self.curv_orth(cell);
    }, nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("alpha"), nb::arg("beta"), nb::arg("gamma"), "Fractional-orthogonal derivative conversion")
    .def( "curv_map", [](const CurvFrac &self, const Grid& grid) {
      return self.curv_map(grid);
    }, nb::arg("grid"), "Fractional-grid derivative conversion")
    .def( "curv_map", [](const CurvFrac &self, const int& nu, const int& nv, const int& nw) {
      return self.curv_map(Grid(nu, nv, nw));
    }, nb::arg("nu"), nb::arg("nv"), nb::arg("nw"), "Fractional-grid derivative conversion")
    .doc() = "Fractional (cell) curvatures, with respect to fractional u,v,w";

  using CurvMap = Curv_map<T>;
  PyClass = std::string( "Curv_map_" ) + name;
  nb::class_<CurvMap, MatClass>( m, PyClass.c_str() )
    .def( nb::init<>(), "Null constructor" )
    .def( nb::init<const MatClass&>(), nb::arg("martrix"), "Constructor: copy/convert")
    .def( "curv_frac", [](const CurvMap& self, const Grid& grid) {
      return self.curv_frac(grid);
    }, nb::arg("grid"), "Grid-fractional derivative conversion")
    .def( "curv_frac", [](const CurvMap& self, const int& nu, const int& nv, const int& nw) {
      return self.curv_frac(Grid(nu, nv, nw));
    }, nb::arg("nu"), nb::arg("nv"), nb::arg("nw"), "Grid-fractional derivative conversion")
    .doc() = "Map coordinate curvatures, with respect to grid u,v,w";
    
}

void add_derivs(nb::module_ &m) {
  declare_grad_classes<ftype32>(m, "float");
  declare_grad_classes<ftype64>(m, "double");
  declare_curv_classes<ftype32>(m, "float");
  declare_curv_classes<ftype64>(m, "double");
}