// Nanobind bindings for clipper containers
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>

using namespace clipper;

void add_clipper_util( nb::module_ &m ) {
  nb::class_<Util>( m, "Util", "Utility class." )
      .def( nb::init<>() )
      .def_static( "u2b", &Util::u2b, nb::arg( "x" ), "Convert isotropic U-value to B-factor." )
      .def_static( "b2u", &Util::b2u, nb::arg( "x" ), "Convert isotropic B-factor to U-value." )
      .def_static( "is_nan", static_cast<bool ( * )( const ftype32 )>( &Util::is_nan ), nb::arg( "f" ),
                   "Fast Util::nan() test. Used for missing entries: THIS DOES "
                   "NOT DISTINGUISH BETWEEN NAN & INF" )
      .def_static( "is_nan", static_cast<bool ( * )( const ftype64 )>( &Util::is_nan ), nb::arg( "f" ),
                   "Fast Util::nan() test. Used for missing entries: THIS DOES "
                   "NOT DISTINGUISH BETWEEN NAN & INF" )
      .def_static( "d2rad", &Util::d2rad, nb::arg( "x" ), "degree-to-radian conversion." )
      .def_static( "rad2d", &Util::rad2d, nb::arg( "x" ), "radian-to-degree conversion." )
      .def_static( "pi", &Util::pi, "Return pi." )
      .def_static( "twopi", &Util::twopi, "Return 2 pi." )
      .def_static( "twopi2", &Util::twopi2, "Return two pi squared." )
      .def_static( "eightpi2", &Util::eightpi2, "Return 8 pi squared." );
}