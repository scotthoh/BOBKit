
// PyBind11 Python bindings for Clipper edcalc
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/clipper-contrib.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

// Trampoline Class for EDcalc_base
template <class T> class PyEDcalc_base : public EDcalc_base<T> {
public:
  /* Inherit the constructor */
  using EDcalc_base<T>::EDcalc_base;
  /* Trampoline (need one for each virtual function) */
  bool operator()( Xmap<T> &xmap, const Atom_list &atoms ) const override {
    PYBIND11_OVERRIDE_PURE_NAME( bool, EDcalc_base<T>, "__call__", operator(), (xmap, atoms) );
  }
  bool operator()( NXmap<T> &nxmap, const Atom_list &atoms ) const override {
    PYBIND11_OVERRIDE_PURE_NAME( bool, EDcalc_base<T>, "__call__", operator(), (nxmap, atoms) );
  }
};

template <class T> void declare_edcalc_base( py::module &m, const std::string &name ) {
  using Base = EDcalc_base<T>;
  py::class_<Base, PyEDcalc_base<T>>( m, ( "_EDcalc_base_" + name ).c_str() )
      .def(
          "__call__", []( Base &self, Xmap<T> &xmap, const Atom_list &atoms ) { return self( xmap, atoms ); },
          py::arg( "xmap" ), py::arg( "atoms" ), "Calculate electron density. Function definition." )
      .def(
          "__call__", []( Base &self, NXmap<T> &nxmap, const Atom_list &atoms ) { return self( nxmap, atoms ); },
          py::arg( "nxmap" ), py::arg( "atoms" ), "Calculate electron density. Function definition." );
}

template <class T> void declare_edcalc_classes( py::module &m, const std::string &name ) {
  using Mask = EDcalc_mask<T>;
  using Base = EDcalc_base<T>;
  py::class_<Mask, Base>( m, ( "EDcalc_mask_" + name ).c_str() )
      .def( py::init<const ftype>(), py::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( Mask &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          py::arg( "xmap" ), py::arg( "atoms" ) )
      .def(
          "__call__", []( Mask &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          py::arg( "nxmap" ), py::arg( "atoms" ) )
      .doc() = "Atom mask calculation.\nAll points within the specified radius of an atom "
               "will be set to 1.0, all others will be set to 0.0.";

  using EDIso = EDcalc_iso<T>;
  py::class_<EDIso, Base>( m, ( "EDcalc_iso_" + name ).c_str() )
      .def( py::init<const ftype>(), py::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( EDIso &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          py::arg( "xmap" ), py::arg( "atoms" ) )
      .def(
          "__call__", []( EDIso &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          py::arg( "nxmap" ), py::arg( "atoms" ) )
      .doc() = "Isotropic electron density calculation.";

  using EDAniso = EDcalc_aniso<T>;
  py::class_<EDAniso, Base>( m, ( "EDcalc_aniso_" + name ).c_str() )
      .def( py::init<const ftype>(), py::arg( "radius" ) = 2.5, "Constructor" )
      .def(
          "__call__", []( EDAniso &self, Xmap<T> &xmap, Atom_list &atoms ) { return self( xmap, atoms ); },
          py::arg( "xmap" ), py::arg( "atoms" ) )
      .def(
          "__call__", []( EDAniso &self, NXmap<T> &nxmap, Atom_list &atoms ) { return self( nxmap, atoms ); },
          py::arg( "nxmap" ), py::arg( "atoms" ) )
      .doc() = "Anisotropic electron density calculation.";
}

void init_edcalc( py::module_ &m ) {
  declare_edcalc_base<float>( m, "float" );
  declare_edcalc_classes<float>( m, "float" );

  declare_edcalc_base<double>( m, "double" );
  declare_edcalc_classes<double>( m, "double" );
}