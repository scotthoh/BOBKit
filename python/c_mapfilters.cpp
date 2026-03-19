// Nanobind bindings for clipper mapfilter
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <nanobind/operators.h>
#include <nanobind/trampoline.h>

using namespace clipper;

template <class T> class PyMapFilter_base : public MapFilter_base<T> {
public:
  NB_TRAMPOLINE( MapFilter_base<T>, 1 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( Xmap<T> &result, const Xmap<T> &xmap ) const override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), result, xmap );
  }
};

class PyMapFilterFn_base : public MapFilterFn_base {
public:
  NB_TRAMPOLINE( MapFilterFn_base, 1 );

  /* Trampoline (need one for each virtual function) */
  ftype operator()( const ftype &radius ) const override { NB_OVERRIDE_PURE_NAME( "__call__", operator(), radius ); }
};

template <class T> void declare_mapfilter_base( nb::module_ &m, const std::string &name ) {
  nb::class_<MapFilter_base<T>, PyMapFilter_base<T>>( m, ( "_MapFilter_base_" + name ).c_str() )
      .def( "__call__", [](const MapFilter_base<T> &self, Xmap<T> &result, const Xmap<T> &xmap) {
            return self(result, xmap);
      })
      //.def( "__call__", &MapFilter_base<T>::operator(), nb::arg( "result" ), nb::arg( "xmap_in" ) )
      .doc() = "Base class for map filters.";
}

void declare_mapfilterfn_base( nb::module_ &m ) {
  nb::class_<MapFilterFn_base, PyMapFilterFn_base>( m, "_MapFilterFn_base" )
      .def( "__call__", [](const MapFilterFn_base &self, const ftype &radius) {
            return self(radius);
      })
      //.def( "__call__", &MapFilterFn_base::operator(), nb::arg( "radius" ) )
      .doc() = "Base class for radial map filter.";
}

template <class T> void declare_mapfilters( nb::module_ &m, const std::string &name ) {
  using MFSlow = MapFilter_slow<T>;
  nb::class_<MFSlow, MapFilter_base<T>> mfltrslow( m, ( "MapFilter_slow_" + name ).c_str() );
  nb::enum_<typename MFSlow::TYPE>( mfltrslow, "TYPE", "MapFilter scaling modes" )
      .value( "NONE", MFSlow::TYPE::NONE )
      .value( "Absolute", MFSlow::TYPE::Absolute )
      .value( "Relative", MFSlow::TYPE::Relative )
      .export_values();

  mfltrslow
      .def( nb::init<const MapFilterFn_base &, const ftype, const typename MFSlow::TYPE>(), nb::arg( "fltr" ),
            nb::arg( "scale" ) = 1.0, nb::arg( "scale_mode" ) = MFSlow::TYPE::NONE,
            "Constructor: takes MapFilterFn, scale factor and scale mode." )
      .def( nb::init< Xmap<T> &, const Xmap<T> &, MapFilterFn_base &, const ftype, const typename MFSlow::TYPE>(),
            nb::arg( "result" ), nb::arg( "xmap_in" ), nb::arg( "fltr" ), nb::arg( "scale" ) = 1.0,
            nb::arg( "scale_mode" ) = MFSlow::NONE, "Constructor: Shorthand for constructor+operator." )
      .def( "__call__", [](const MFSlow &self, Xmap<T> &result,const Xmap<T> &xmap) {
            return self(result, xmap);
      })
      //.def( "__call__", &MFSlow::operator(), nb::arg( "result" ), nb::arg( "xmap_in" ) )
      .doc() = "Simple slow convolution-based radial mapfiltering implementation.\n"
               "This version is of course very slow, and is mainly provided so "
               "you can test the precision of the fft version.\n\n"
               ".. code-block:: python\n\n"
               " # make squared map\n"
               " xmap2 = Xmap_float( xmap )\n"
               " ix = xmap2.first()\n"
               " while not ix.last():\n"
               "   xmap2[ix] = pow( xmap2[ix], 2.0 )\n"
               "   ix.next()\n\n"
               " # now calculate local mean, local mean squared\n"
               " fn = MapFilterFn_step( filter_radius )\n"
               " fltr = MapFilter_slow_float( fn, 1.0, MapFilter_slow_float.Relative )\n"
               " lmean, lsigm = Xmap_float(), Xmap_float()\n"
               " fltr( lmean, xmap )\n"
               " fltr( lsigm, xmap2 )\n\n"
               " # calculate std deviation\n"
               " ix = lmean.first()\n"
               " while not ix.last():\n"
               "   lsigm[ix] = sqrt( lsigm[ix] - pow( lmean[ix], 2.0 ) )\n"
               "   ix.next()\n"
               "This would be a useful step in solvent mask determination, for example.";

  using MFFft = MapFilter_fft<T>;
  nb::class_<MFFft, MapFilter_base<T>> mfltrfft( m, ( "MapFilter_fft_" + name ).c_str() );
  nb::enum_<typename MFFft::TYPE>( mfltrfft, "TYPE", "MapFilter scaling modes" )
      .value( "NONE", MFFft::TYPE::NONE )
      .value( "Absolute", MFFft::TYPE::Absolute )
      .value( "Relative", MFFft::TYPE::Relative )
      .export_values();

  mfltrfft
      .def( nb::init<const MapFilterFn_base &, const ftype, const typename MFFft::TYPE>(), nb::arg( "fltr" ),
            nb::arg( "scale" ) = 1.0, nb::arg( "scale_mode" ) = MFFft::NONE,
            "Constructor: takes MapFilterFn, scale factor and scale mode." )
      .def( nb::init< Xmap<T> &, const Xmap<T> &, MapFilterFn_base &, const ftype, const typename MFFft::TYPE>(),
            nb::arg( "result" ), nb::arg( "xmap_in" ), nb::arg( "fltr" ), nb::arg( "scale" ) = 1.0,
            nb::arg( "scale_mode" ) = MFFft::NONE, "Constructor: Shorthand for constructor+operator." )
      .def( "__call__", ( bool ( MFFft::* )( Xmap<T> &, const Xmap<T> & ) const ) & MFFft::operator(),
            nb::arg( "result" ), nb::arg( "xmap_in" ) )
      .def( "__call__", ( bool ( MFFft::* )( NXmap<T> &, const NXmap<T> & ) const ) & MFFft::operator(),
            nb::arg( "result" ), nb::arg( "nxmap_in" ) )
      .doc() = "Simple fft-based radial mapfiltering implementation.\n"
               "The FFT method is fast, and also gives good precision. "
               "The following example demonstrates how to use the MapFilter to "
               "calculate the local mean and local deviation of an electron "
               "density map, in 'xmap':";
}

void declare_mapfilterfns( nb::module_ &m ) {
  using Step = MapFilterFn_step;
  nb::class_<Step, MapFilterFn_base>( m, "MapFilterFn_step" )
      .def( nb::init<const ftype &>(), nb::arg( "radius" ), "Constructor: takes radius for step function cutoff." )
      .def( "__call__", &Step::operator(), nb::arg( "radius" ),
            "Evaluate radial step function: 1.0 if inside or 0.0 outside." )
      .doc() = "Step-function radial map filter.";

  using Linear = MapFilterFn_linear;
  nb::class_<Linear, MapFilterFn_base>( m, "MapFilterFn_linear" )
      .def( nb::init<const ftype &>(), nb::arg( "radius" ), "Constructor: takes radius for linear function cutoff." )
      .def( "__call__", &Linear::operator(), nb::arg( "radius" ),
            "Evaluate radial linear function: 1-r/r0 if inside or 0.0 outside." )
      .doc() = "Linear-function radial map filter.";

  using Quad = MapFilterFn_quadratic;
  nb::class_<Quad, MapFilterFn_base>( m, "MapFilterFn_quadratic" )
      .def( nb::init<const ftype &>(), nb::arg( "radius" ), "Constructor: takes radius for quadratic function cutoff." )
      .def( "__call__", &Quad::operator(), nb::arg( "radius" ),
            "Evaluate radial quadratic function: 1-r^2/r0^2 if inside or 0.0 outside." )
      .doc() = "Quadratic-function radial map filter.";
}

void add_mapfilters( nb::module_ &m ) {
  declare_mapfilterfn_base( m );
  declare_mapfilterfns( m );
  declare_mapfilter_base<ftype32>( m, "float" );
  declare_mapfilter_base<ftype64>( m, "double" );
  declare_mapfilters<ftype32>( m, "float" );
  declare_mapfilters<ftype64>( m, "double" );
}
