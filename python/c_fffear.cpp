// Wrapper for clipper fffear
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/clipper-contrib.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;


// trampoline class for property base
template <class T> class PyFFFear_base : public FFFear_base<T> {
public:
  /* Inherit the constructors. */
  using FFFear_base<T>::FFFear_base;
  /* Trampoline (need one for each virtual function) */
  bool operator()(Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
                   const NX_operator &nxop) const override {
      PYBIND11_OVERRIDE_PURE_NAME(bool, FFFear_base<T>,
                                  "__call__", operator(),
                                  (result, srchval, srchwgt, nxop) // arguments
      );
    }
};


template <class T> void declare_fffear_base( py::module &m, const std::string &name ) {
  using Fbase = FFFear_base<T>;
  py::class_<Fbase, PyFFFear_base<T> >( m, ( "_FFFear_base_" + name ).c_str() )
      .def(
          "__call__",
          []( Fbase &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
              const NX_operator &nxop ) { return self( result, srchval, srchwgt, nxop ); },
          py::arg( "result" ), py::arg( "srchval" ), py::arg( "srchwgt" ), py::arg( "nxop" ),
          "Perform actual fffear calculation. Function definition." )
      .doc() = "Base class for Fast Fourier Feature recognitio (FFFEAR) methods.";
}

template <class T> void declare_fffear_fft( py::module &m, const std::string &name ) {
  using Fbase = FFFear_base<T>;
  using FFFear = FFFear_fft<T>;
  py::class_<FFFear, Fbase> fffear( m, ( "FFFear_fft_" + name ).c_str() );
  py::enum_<typename FFFear::FFTtype>( fffear, "FFTtype", "FFT backend selection" )
      .value( "Default", FFFear::FFTtype::Default )
      .value( "Normal", FFFear::FFTtype::Normal )
      .value( "Sparse", FFFear::FFTtype::Sparse )
      .export_values();

  fffear.def( py::init<>(), "Null constructor" )
      .def( py::init<const Xmap<T> &>(), py::arg( "xmap" ), "Constructor with Xmap" )
      .def( py::init< Xmap<T> &, const NXmap<T> &, const NXmap<T> &, const Xmap<T> &, const NX_operator &>(),
            py::arg( "result" ), py::arg( "srchval" ), py::arg( "srchwgt" ), py::arg( "xmap" ), py::arg( "nxop" ),
            "Constructor: shorthand for constructor+operator." )
      .def(
          "init", []( FFFear &self, const Xmap<T> xmap ) { return self.init( xmap ); },
          "Initialiser: initiliase with the given target Xmap" )
      .def( "set_fft_type", &FFFear::set_fft_type, py::arg( "type" ), "Set FFT backend" )
      .def(
          "set_resolution", []( FFFear &self, Resolution reso ) { self.set_resolution( reso ); },
          py::arg( "resolution" ), "Set resolution cutoff" )
      .def(
          "set_resolution", []( FFFear &self, ftype reso ) { self.set_resolution( Resolution( reso ) ); },
          py::arg( "resolution" ), "Set resolution cutoff" )
      .def(
          "__call__",
          []( FFFear &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt,
              const NX_operator &nxop ) { return self( result, srchval, srchwgt, nxop ); },
          py::arg( "result" ), py::arg( "srchval" ), py::arg( "srchwgt" ), py::arg( "nxop" ),
          "Search for given target." )
      .def(
          "__call__",
          []( FFFear &self, Xmap<T> &result, const NXmap<T> &srchval, const NXmap<T> &srchwgt, const RTop_orth &rtop ) {
            return self( result, srchval, srchwgt, rtop );
          },
          py::arg( "result" ), py::arg( "srchval" ), py::arg( "srchwgt" ), py::arg( "rtop" ),
          "Search for given target." )
      .doc() = "FFT-based fffear implementation.\nThis implementation is currently unoptimised, but much faster than "
               "the simple implementation.";
}

void init_fffear( py::module &m ) {
  declare_fffear_base<ftype32>( m, "float" );
  declare_fffear_base<ftype64>( m, "double" );
  declare_fffear_fft<ftype32>( m, "float" );
  declare_fffear_fft<ftype64>( m, "double" );
}