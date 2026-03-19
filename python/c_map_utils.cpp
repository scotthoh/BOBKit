// Nanobind bindings for clipper map_utils
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

void declare_map_stats( nb::module_ &m ) {
  nb::class_<Map_stats>( m, "Map_stats" )
      .def( nb::init<>() )
      .def( nb::init<const Xmap<ftype32> &>(), "Constructor from Xmap" )
      .def( nb::init<const Xmap<ftype64> &>(), "Constructor from Xmap" )
      .def( nb::init<const NXmap<ftype32> &>(), "Constructor from NXmap" )
      .def( nb::init<const NXmap<ftype64> &>(), "Constructor from NXmap" )
      .def_prop_ro( "mean", &Map_stats::mean, "Return mean of map." )
      .def_prop_ro( "std_dev", &Map_stats::std_dev, "Return standard deviation of map." )
      .def_prop_ro( "min", &Map_stats::min, "Return minimum of map." )
      .def_prop_ro( "max", &Map_stats::max, "Return maximum of map." )
      .def_prop_ro( "range", &Map_stats::range, "Return range." )
      .def( "__repr__", []( const Map_stats &self ) { return "<clipper.Map_stats class.>"; } )
      .def(
          "format",
          []( const Map_stats &self ) {
            return "Min: " + String( self.min(), 10, 4 ) + ", Max: " + String( self.max(), 10, 4 ) +
                   ", Mean: " + String( self.mean(), 10, 4 ) + ", SD: " + String( self.std_dev(), 10, 4 );
          },
          "Return formatted string representation." )
      .def( "__str__",
            []( const Map_stats &self ) {
              return String( self.min(), 10, 4 ) + ", " + String( self.max(), 10, 4 ) + ", " +
                     String( self.mean(), 10, 4 ) + ", " + String( self.std_dev(), 10, 4 );
            } )
      .doc() = "Generic map statistics class.\nThis class is used "
               "to calculate and store the mean and standard deviation "
               "of a generic map object of scalar types (e.g. Xmap, "
               "NXmap). If the map contains NaN values, those points "
               "are excluded for the calculation. In the case of an Xmap, "
               "the proper multiplicty corrections are applied to give "
               "statistics for a whole unit cell.\n\n";
}

template <class T> void apply_sort_methods( nb::class_<Map_index_sort> &pyclass ) {
  pyclass
      .def_static( "sort_increasing",
                   static_cast<void ( * )( const Xmap<T> &, std::vector<int> & )>( &Map_index_sort::sort_increasing ) )
      .def_static( "sort_decreasing",
                   static_cast<void ( * )( const Xmap<T> &, std::vector<int> & )>( &Map_index_sort::sort_decreasing ) );
}

void declare_map_index_sort( nb::module_ &m ) {
  nb::class_<Map_index_sort> misort( m, "Map_index_sort" );
  // these are templates declared in clipper/core/map_utils.cpp
  apply_sort_methods<ftype32>( misort );
  apply_sort_methods<ftype64>( misort );
  apply_sort_methods<int>( misort );
  apply_sort_methods<unsigned int>( misort );
}

void add_map_utils( nb::module_ &m ) {
  declare_map_stats( m );
  declare_map_index_sort( m );
}