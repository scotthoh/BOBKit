// Nanobind bindings for clipper fffear
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <nanobind/trampoline.h>

using namespace clipper;

// Trampoline Class, from function_object_bases.h
template <class T> class PyFFFear_base : public FFFear_base<T> {
public:
  NB_TRAMPOLINE( FFFear_base<T>, 1 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
                   const NX_operator &nxop ) const override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), result, srchval, srchwgt, nxop );
  }
};

template <class T> void declare_fffear_base( nb::module_ &m, const std::string &name ) {
  using Fbase = FFFear_base<T>;
  nb::class_<Fbase, PyFFFear_base<T> >( m, ( "_FFFear_base_" + name ).c_str() )
      .def(
          "__call__",
          []( Fbase &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
              const NX_operator &nxop ) { return self( result, srchval, srchwgt, nxop ); },
          nb::arg( "result" ), nb::arg( "srchval" ), nb::arg( "srchwgt" ), nb::arg( "nxop" ),
          "Perform actual fffear calculation. Function definition." )
      .doc() = "Base class for Fast Fourier Feature recognitio (FFFEAR) methods.";
}

template <class T> void declare_fffear_fft( nb::module_ &m, const std::string &name ) {
  using Fbase = FFFear_base<T>;
  using FFFear = FFFear_fft<T>;
  nb::class_<FFFear, Fbase> fffear( m, ( "FFFear_fft_" + name ).c_str() );
  nb::enum_<typename FFFear::FFTtype>( fffear, "FFTtype", "FFT backend selection" )
      .value( "Default", FFFear::FFTtype::Default )
      .value( "Normal", FFFear::FFTtype::Normal )
      .value( "Sparse", FFFear::FFTtype::Sparse )
      .export_values();

  fffear.def( nb::init<>(), "Null constructor" )
      .def( nb::init<const Xmap<T> &>(), nb::arg( "xmap" ), "Constructor with Xmap" )
      .def( nb::init< Xmap<T> &, const NXmap<T> &, const NXmap<T> &, const Xmap<T> &, const NX_operator &>(),
            nb::arg( "result" ), nb::arg( "srchval" ), nb::arg( "srchwgt" ), nb::arg( "xmap" ), nb::arg( "nxop" ),
            "Constructor: shorthand for constructor+operator." )
      .def(
          "init", []( FFFear &self, const Xmap<T> xmap ) { return self.init( xmap ); },
          "Initialiser: initiliase with the given target Xmap" )
      .def( "set_fft_type", &FFFear::set_fft_type, nb::arg( "type" ), "Set FFT backend" )
      .def(
          "set_resolution", []( FFFear &self, Resolution reso ) { self.set_resolution( reso ); },
          nb::arg( "resolution" ), "Set resolution cutoff" )
      .def(
          "set_resolution", []( FFFear &self, ftype reso ) { self.set_resolution( Resolution( reso ) ); },
          nb::arg( "resolution" ), "Set resolution cutoff" )
      .def(
          "__call__",
          []( FFFear &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
              const NX_operator &nxop ) { return self( result, srchval, srchwgt, nxop ); },
          nb::arg( "result" ), nb::arg( "srchval" ), nb::arg( "srchwgt" ), nb::arg( "nxop" ),
          "Search for given target." )
      .def(
          "__call__",
          []( FFFear &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt, const RTop_orth &rtop ) {
            return self( result, srchval, srchwgt, rtop );
          },
          nb::arg( "result" ), nb::arg( "srchval" ), nb::arg( "srchwgt" ), nb::arg( "rtop" ),
          "Search for given target." )
      .doc() = "FFT-based fffear implementation.\nThis implementation is currently unoptimised, but much faster than "
               "the simple implementation.";
}

void add_fffear( nb::module_ &m ) {
  declare_fffear_base<float>( m, "float" );
  declare_fffear_base<double>( m, "double" );
  declare_fffear_fft<float>( m, "float" );
  declare_fffear_fft<double>( m, "double" );
}