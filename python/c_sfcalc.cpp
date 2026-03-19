// Nanobind bindings for clipper contrib sfcalc
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper sfcalc by Tristan Croll

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <nanobind/trampoline.h>

using namespace clipper;
using namespace clipper::datatypes;

// from function_object_bases.h
template <class T> class PySFcalc_base : public SFcalc_base<T> {
public:
  NB_TRAMPOLINE( SFcalc_base<T>, 1 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( HKL_data<F_phi<T> > &fphidata, const Atom_list &atoms ) const override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), fphidata, atoms );
  }
};

template <class T> void declare_sfcalc_base( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_base<T>;
  auto pyclass_name = std::string( "_SFcalc_base_" ) + dtype;
  nb::class_<Class, PySFcalc_base<T>>( m, pyclass_name.c_str(), "Base class for structure factor calculation methods" )
      .def( "__call__", []( const Class &self, HKL_data<F_phi<T>> &fphidata, const Atom_list &atoms ) {
        return self( fphidata, atoms );
      } );
} // declare_sfcalc_base

template <class T> void declare_sfcalc_iso_sum( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_iso_sum<T>;
  auto pyclass_name = std::string( "SFcalc_iso_sum_" ) + dtype;
  nb::class_<Class, SFcalc_base<T>>( m, pyclass_name.c_str() )
      .def( nb::init<>() )
      .def( nb::init<HKL_data<F_phi<T>> &, const Atom_list &>() );
} // declare_sfcalc_iso_sum

template <class T> void declare_sfcalc_aniso_sum( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_aniso_sum<T>;
  auto pyclass_name = std::string( "SFcalc_aniso_sum_" ) + dtype;
  nb::class_<Class, SFcalc_base<T>>( m, pyclass_name.c_str() )
      .def( nb::init<>() )
      .def( nb::init<HKL_data<F_phi<T>> &, const Atom_list &>() );
} // declare_sfcalc_aniso_sum

template <class T> void declare_sfcalc_iso_fft( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_iso_fft<T>;
  auto pyclass_name = std::string( "SFcalc_iso_fft_" ) + dtype;
  nb::class_<Class, SFcalc_base<T>>( m, pyclass_name.c_str() )
      .def( nb::init<const ftype, const ftype, const ftype>(), nb::arg( "radius" ) = 2.5, nb::arg( "rate" ) = 1.5,
            nb::arg( "uadd" ) = 0.0 )
      .def( nb::init<HKL_data<F_phi<T>> &, const Atom_list &, const ftype, const ftype, const ftype>(),
            nb::arg( "fphidata_out" ), nb::arg( "atoms" ), nb::arg( "radius" ) = 2.5, nb::arg( "rate" ) = 1.5,
            nb::arg( "uadd" ) = 0.0 );
} // declare_sfcalc_iso_fft

template <class T> void declare_sfcalc_aniso_fft( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_aniso_fft<T>;
  auto pyclass_name = std::string( "SFcalc_aniso_fft_" ) + dtype;
  nb::class_<Class, SFcalc_base<T>>( m, pyclass_name.c_str() )
      .def( nb::init<const ftype, const ftype, const ftype>(), nb::arg( "radius" ) = 2.5, nb::arg( "rate" ) = 1.5,
            nb::arg( "uadd" ) = 0.0 )
      .def( nb::init<HKL_data<F_phi<T>> &, const Atom_list &, const ftype, const ftype, const ftype>(),
            nb::arg( "fphidata_out" ), nb::arg( "atoms" ), nb::arg( "radius" ) = 2.5, nb::arg( "rate" ) = 1.5,
            nb::arg( "uadd" ) = 0.0 );
} // declare_sfcalc_aniso_fft

void init_sfcalc( nb::module_ &m ) {
  declare_sfcalc_base<ftype32>( m, "float" );
  declare_sfcalc_iso_sum<ftype32>( m, "float" );
  declare_sfcalc_aniso_sum<ftype32>( m, "float" );
  declare_sfcalc_iso_fft<ftype32>( m, "float" );
  declare_sfcalc_aniso_fft<ftype32>( m, "float" );

  declare_sfcalc_base<ftype64>( m, "double" );
  declare_sfcalc_iso_sum<ftype64>( m, "double" );
  declare_sfcalc_aniso_sum<ftype64>( m, "double" );
  declare_sfcalc_iso_fft<ftype64>( m, "double" );
  declare_sfcalc_aniso_fft<ftype64>( m, "double" );

} // init_sfcalc
