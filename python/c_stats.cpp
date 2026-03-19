// Nanobind bindings for clipper stats
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

template <class T> void declare_range( nb::module_ &m, const std::string &name ) {
  using Class = Range<T>;
  nb::class_<Class>( m, ( "Range_" + name ).c_str(), "Range upper and lower bounds of some type." )
      .def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "min" ), nb::arg( "max" ), "Constructor from min and max." )
      .def_prop_ro( "min", &Class::min, "Return minimum value." )
      .def_prop_ro( "max", &Class::max, "Return maximum value." )
      .def( "range", &Class::range, "Return range = max-min." )
      .def( "include", &Class::include, nb::arg( "datum" ), "Update limits to include new datum." )
      .def( "contains", &Class::include, nb::arg( "datum" ), "Test if data is within limits (min<=datum<=max)." )
      .def( "truncate", &Class::truncate, nb::arg( "datum" ), "Truncate data to be within range." )
      .def( "__repr__", [=]( const Class &self ) {
        return "<clipper.Range_" + name + String( ftype( self.min() ), 10, 4 ) + " - " +
               String( ftype( self.max() ), 10, 4 ) + ">";
      } );
}

void declare_range_sampling( nb::module_ &m ) {
  nb::class_<Range_sampling, Range<ftype>>( m, "Range_sampling", "Range sampling: discrete sampling of a real range." )
      .def( nb::init<>() )
      .def( nb::init<const int &>(), "Constructor from a number of samplings." )
      .def( nb::init<const Range<ftype> &, const int &>(), "Constructor from range and number of samplings." )
      .def( "indexf", &Range_sampling::indexf, nb::arg( "pos" ),
            "Return fractional posn in counting range from x-value (0..n)." )
      .def( "x", ( ftype ( Range_sampling::* )( const ftype & ) const ) & Range_sampling::x, nb::arg( "frac_pos" ),
            "Return x-value (0..n) from fractional posn in counting range." )
      .def( "index", &Range_sampling::index, nb::arg( "pos" ), "Return nearest index to particular x-value." )
      .def( "index_bounded", &Range_sampling::index_bounded, nb::arg( "pos" ),
            "Return nearest index to particular x-value (bounded 0...n-1)" )
      .def( "x", ( ftype ( Range_sampling::* )( const int & ) const ) & Range_sampling::x, nb::arg( "pos" ),
            "Return x-value corresponding to centre of i'th range." )
      .def( "x_min", &Range_sampling::x_min, "Return x-value corresponding to bottom of i'th range." )
      .def( "x_max", &Range_sampling::x_max, "Return x-value corresponding to top of i'th range." )
      .def_prop_ro( "size", &Range_sampling::size, "Return number of samplings in range." )
      .def( "__repr__", []( const Range_sampling &self ) { return "<clipper.Range_sampling class.>"; } );
}

void declare_histogram( nb::module_ &m ) {
  nb::class_<Histogram, Range_sampling>( m, "Histogram" )
      .def( nb::init<>() )
      .def( nb::init<const Range<ftype> &, const int &>(), nb::arg( "range" ), nb::arg( "n" ) )
      .def( "accumulate", ( void ( Histogram::* )( const ftype & ) )&Histogram::accumulate, nb::arg( "x" ) )
      .def( "accumulate", ( void ( Histogram::* )( const ftype &, const ftype & ) )&Histogram::accumulate,
            nb::arg( "x" ), nb::arg( "value" ) )
      .def( "sum", &Histogram::sum )
      .def( "y", ( const ftype &( Histogram::* )( const int & ) const ) & Histogram::y, nb::arg( "index" ) )
      .def( "y", ( ftype ( Histogram::* )( const ftype & ) const ) & Histogram::y, nb::arg( "x" ) )
	    .def( nb::self += nb::self, nb::rv_policy::none )
      .def( "__repr__",
            []( const Histogram &self ) { return "<clipper.Histogram with sample size " + String( self.size() ); } )
      .doc() = "General histogram class.\nThis class is used to accumulate and "
               "access a histogram of values spread over a specified range. On "
               "storing data or retrieving by interpolation the range is checked.";
}

void declare_generic_ordinal( nb::module_ &m ) {
  nb::class_<Generic_ordinal>( m, "Generic_ordinal" )
      .def( nb::init<>() )
      .def( nb::init<const Range<ftype> &, const int &>(), nb::arg( "range" ), nb::arg( "n" ),
            "Constructor: from range and sampling" )
      .def(
          "init",
          []( Generic_ordinal &self, Range<ftype> &range, const int num_ranges ) { self.init( range, num_ranges ); },
          nb::arg( "range" ), nb::arg( "num_ranges" ) = 1000 )
      .def(
          "init",
          []( Generic_ordinal &self, std::vector<ftype> &values, const int num_ranges ) {
            self.init( values, num_ranges );
          },
          nb::arg( "values" ), nb::arg( "num_ranges" ) = 1000 )
      .def( "ordinal", &Generic_ordinal::ordinal, nb::arg( "value" ) )
      .def(
          "accumulate", []( Generic_ordinal &self, const ftype &value ) { self.accumulate( value ); },
          nb::arg( "value" ) )
      .def(
          "accumulate",
          []( Generic_ordinal &self, const ftype &value, const ftype &weight ) { self.accumulate( value, weight ); },
          nb::arg( "value" ), nb::arg( "weight" ) )
      .def( "prep_ordinal", &Generic_ordinal::prep_ordinal )
      .def( "invert", &Generic_ordinal::invert );
}

void init_clipper_stats( nb::module_ &m ) {
  declare_range<int>( m, "int" );
  declare_range<float>( m, "float" );
  declare_range<double>( m, "double" );
  declare_range_sampling( m );
  declare_histogram( m );
  declare_generic_ordinal( m );
}