// Nanobind bindings for clipper edcalc
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <nanobind/trampoline.h>

using namespace clipper;

// Trampoline Class, from function_object_bases.h
template <class T> class PyEDcalc_base : public EDcalc_base<T> {
public:
  NB_TRAMPOLINE( EDcalc_base<T>, 2 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( Xmap<T> &xmap, const Atom_list &atoms ) const override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), xmap, atoms );
  }
  bool operator()( NXmap<T> &nxmap, const Atom_list &atoms ) const override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), nxmap, atoms );
  }
  //~PyEDcalc_base() { NB_OVERRIDE_PURE_NAME( "__del__", ~PyEDcalc_base() ); }
};

template <class T> void declare_edcalc_base( nb::module_ &m, const std::string &name ) {
  using Base = EDcalc_base<T>;
  nb::class_<EDcalc_base<T>, PyEDcalc_base<T>>( m, ( "_EDcalc_base_" + name ).c_str() )
      .def(
          "__call__", []( Base &self, Xmap<T> &xmap, const Atom_list &atoms ) { return self( xmap, atoms ); },
          nb::arg( "xmap" ), nb::arg( "atoms" ), "Calculate electron density. Function definition." )
      .def(
          "__call__", []( Base &self, NXmap<T> &nxmap, const Atom_list &atoms ) { return self( nxmap, atoms ); },
          nb::arg( "nxmap" ), nb::arg( "atoms" ), "Calculate electron density. Function definition." );
}

template <class T> void declare_edcalc_classes( nb::module_ &m, const std::string &name ) {
  using Mask = EDcalc_mask<T>;
  using Base = EDcalc_base<T>;
  nb::class_<Mask, Base>( m, ( "EDcalc_mask_" + name ).c_str() )
      .def( nb::init<const ftype>(), nb::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( Mask &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          nb::arg( "xmap" ), nb::arg( "atoms" ) )
      .def(
          "__call__", []( Mask &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          nb::arg( "nxmap" ), nb::arg( "atoms" ) )
      .doc() = "Atom mask calculation.\nAll points within the specified radius of an atom "
               "will be set to 1.0, all others will be set to 0.0.";

  using EDIso = EDcalc_iso<T>;
  nb::class_<EDIso, Base>( m, ( "EDcalc_iso_" + name ).c_str() )
      .def( nb::init<const ftype>(), nb::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( EDIso &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          nb::arg( "xmap" ), nb::arg( "atoms" ) )
      .def(
          "__call__", []( EDIso &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          nb::arg( "nxmap" ), nb::arg( "atoms" ) )
      .doc() = "Isotropic electron density calculation.";

  using EDAniso = EDcalc_aniso<T>;
  nb::class_<EDAniso, Base>( m, ( "EDcalc_aniso_" + name ).c_str() )
      .def( nb::init<const ftype>(), nb::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( EDAniso &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          nb::arg( "xmap" ), nb::arg( "atoms" ) )
      .def(
          "__call__", []( EDAniso &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          nb::arg( "nxmap" ), nb::arg( "atoms" ) )
      .doc() = "Anisotropic electron density calculation.";
}

void add_edcalc( nb::module_ &m ) {
  declare_edcalc_base<float>( m, "float" );
  declare_edcalc_classes<float>( m, "float" );

  declare_edcalc_base<double>( m, "double" );
  declare_edcalc_classes<double>( m, "double" );
}
