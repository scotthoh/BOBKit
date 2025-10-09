// Nanobind bindings for clipper cell.h
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
//#include <clipper/clipper-gemmi.h>
#include <clipper/core/cell.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>

using namespace clipper;

// std::string six_params(double a, double b, double c, double alp, double bet, double gam) {
//     char buf[128];
//
// }

void declare_metric_tensor( nb::module_ &m ) {
  nb::class_<Metric_tensor> metric_tensor( m, "Metric_tensor" );
  metric_tensor.def( nb::init<>() )
      .def( nb::init<const double &, const double &, const double &, const double &, const double &, const double &>(),
            nb::arg( "a" ), nb::arg( "b" ), nb::arg( "c" ), nb::arg( "alpha" ), nb::arg( "beta" ), nb::arg( "gamma" ),
            "Constructor: takes parameters of normal or inverse cell." )
      .def( "lengthsq", ( double ( Metric_tensor::* )( const Vec3<> & ) const ) & Metric_tensor::lengthsq,
            "Apply metric to vector." )
      .def( "lengthsq", ( double ( Metric_tensor::* )( const Vec3<int> & ) const ) & Metric_tensor::lengthsq,
            "Apply metric to int vector." )
      .def( "format", &Metric_tensor::format, "Return formatted string representation." )
      .def( "__str__", []( const Metric_tensor &self ) { return self.format(); } )
      .def( "__repr__", []( const Metric_tensor &self ) { return "< clipper.Metric_tensor " + self.format() + " >"; } )
      .doc() = "The metric tensor is used to determine a distance in real or reciprocal"
               "space using fraction coordinates or Miller indices. It is symmetrical, "
               "so only the upper triangle is stored with the off-diagonal elements "
               "doubled.";
}

void declare_cell_descr( nb::module_ &m ) {
  nb::class_<Cell_descr> cell_descr( m, "Cell_descr" );
  cell_descr.def( nb::init<>() )
      .def( nb::init<const double &, const double &, const double &, const double &, const double &, const double &>(),
            nb::arg( "a" ), nb::arg( "b" ), nb::arg( "c" ), nb::arg( "alpha" ) = 90., nb::arg( "beta" ) = 90.,
            nb::arg( "gamma" ) = 90., "Constructor from cell parameters." )
      .def(
          "__init__",
          []( Cell_descr *cdes, std::array<double, 6> &arr ) {
            new ( cdes ) Cell_descr( arr[0], arr[1], arr[2], arr[3], arr[4], arr[5] );
          },
          "Constructor from cell parameters (list/numpy array)." )
      .def_prop_ro( "a", &Cell_descr::a, "Get a." )
      .def_prop_ro( "b", &Cell_descr::b, "Get b." )
      .def_prop_ro( "c", &Cell_descr::c, "Get c." )
      .def_prop_ro( "alpha", &Cell_descr::alpha, "Get alpha." )
      .def_prop_ro( "beta", &Cell_descr::beta, "Get beta." )
      .def_prop_ro( "gamma", &Cell_descr::gamma, "Get gamma." )
      .def_prop_ro( "alpha_deg", &Cell_descr::alpha_deg, "Get alpha in degrees." )
      .def_prop_ro( "beta_deg", &Cell_descr::beta_deg, "Get beta in degrees." )
      .def_prop_ro( "gamma_deg", &Cell_descr::gamma_deg, "Get gamma in degrees." )
      .def( "format", &Cell_descr::format, "Return formatted string representation." )
      .def( "__str__", []( const Cell_descr &self ) { return self.format(); } )
      .def( "__repr__", []( const Cell_descr &self ) { return "< clipper." + self.format().trim() + " >"; } )
      .doc() = "Cell description (automatically converts to radians)"
               "The cell description is a compact description of a cell, "
               "containing just the cell parameters. It is usually used to "
               "construct a full Cell object, which provides the expected "
               "functionality.";
}

void declare_cell( nb::module_ &m ) {
  nb::class_<Cell, Cell_descr> cell( m, "Cell" );
  cell.def( nb::init<>(), "Null constructor, must initialise later." )
      .def( nb::init<const Cell_descr &>(), nb::arg( "Cell_descr" ), "Constructor from Cell descriptor" )
      //.def("init", &Cell::init, nb::arg("Cell_description"))
      .def( "__init__",
          []( Cell *c, std::array<double, 6> &params ) {
              new ( c ) Cell( Cell_descr( params[0], params[1], params[2], params[3], params[4], params[5] ) );
            } )
      .def(
          "init", []( Cell &self, const Cell_descr &cdes ) { self.init( cdes ); }, nb::arg( "Cell_descr" ),
          "Initialise with Cell descriptor." )
      //.def(
      //    "init", []( Cell &self, const gemmi::UnitCell &c ) { self.init( GEMMI::cell( c ).descr() ); },
      //    nb::arg( "cell" ), "Initialise with GEMMI UnitCell." )
      .def(
          "init", []( Cell &self, const Cell &c ) { self.init( c.descr() ); }, nb::arg( "cell" ),
          "Initialise with Cell descriptor." )
      // from clipper-gemmi
      //.def_static(
      //    "to_gemmi_cell", []( const Cell &c ) { return GEMMI::cell( c ); }, "Convert CLIPPER cell to GEMMI Unitcell." )
      //.def_static(
      //    "from_gemmi_cell", []( const gemmi::UnitCell &c ) { return GEMMI::cell( c ); },
      //    "Convert GEMMI Unitcell to CLIPPER Cell." )
      .def_prop_ro( "a_star", &Cell::a_star, "Get a*" )
      .def_prop_ro( "b_star", &Cell::b_star, "Get b*" )
      .def_prop_ro( "c_star", &Cell::c_star, "Get c*" )
      .def_prop_ro( "alpha_star", &Cell::alpha_star, "Get alpha*" )
      .def_prop_ro( "beta_star", &Cell::beta_star, "Get beta*" )
      .def_prop_ro( "gamma_star", &Cell::gamma_star, "Get gamma*" )
      .def_prop_ro( "description", &Cell::descr, "Return cell dimensions." )
      .def_prop_ro( "volume", &Cell::volume, "Return cell volume." )
      .def( "debug", &Cell::debug, "Output class details." )
      .def( "equals", &Cell::equals, nb::arg( "cell" ), nb::arg( "tol" ) = 1.0, "Test equality with another cell." )
      .def_prop_ro( "matrix_orth", &Cell::matrix_orth, "Return orthogonalisation matrix" )
      .def_prop_ro( "matrix_frac", &Cell::matrix_frac, "Return fractionalisation matrix" )
      .def_prop_ro( "metric_real", &Cell::metric_real, "Return real space metric tensor." )
      .def_prop_ro( "metric_reci", &Cell::metric_reci, "Return reciprocal space metric tensor." )
      .def_prop_ro(
          "parameters",
          []( const Cell &c ) {
            return nb::make_tuple( c.a(), c.b(), c.c(), c.alpha_deg(), c.beta_deg(), c.gamma_deg() );
          },
          "Return a tuple of Cell parameters." )
      .doc() = "Cell object.\n"
               "The Cell class is the fully functional description of the unit "
               "cell. In addition to the cell parameters, it stores derived "
               "information including the cell volume, orthogonalising and "
               "fractionalising matrices, and the metric tensors.";
}

void add_cell( nb::module_ &m ) {
  declare_metric_tensor( m );
  declare_cell_descr( m );
  declare_cell( m );
}
