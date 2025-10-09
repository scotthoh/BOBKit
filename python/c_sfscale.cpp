// Nanobind bindings for clipper contrib_sfscale
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper sfscale by Tristan Croll

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <clipper/contrib/function_object_bases.h>
#include <nanobind/operators.h>
#include <nanobind/trampoline.h>

using namespace clipper;
using namespace clipper::datatypes;

// from function_object_bases.h
template <class T> class PySFscale_base : public SFscale_base<T> {
public:
  NB_TRAMPOLINE( SFscale_base<T>, 3 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( HKL_data<F_sigF<T> > &fo, const HKL_data<F_phi<T> > &fc ) override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), fo, fc );
  }
  bool operator()( HKL_data<F_phi<T> > &fc, const HKL_data<F_sigF<T> > &fo ) override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), fc, fo );
  }
  bool operator()( HKL_data<F_sigF<T> > &fo ) override { NB_OVERRIDE_PURE_NAME( "__call__", operator(), fo ); }
};

template <class T> void declare_sfscale_base( nb::module_ &m, const char *dtype ) {
  using Class = SFscale_base<T>;
  auto pyclass_name = std::string( "_SFscale_base_" ) + dtype;
  nb::class_<Class, PySFscale_base<T>>( m, pyclass_name.c_str(), "Base class for structure factor scaling methods" )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_sigF<T> > &, const HKL_data<F_phi<T> > & ) )&Class::operator(),
            nb::arg( "fo" ), nb::arg( "fc" ), "Scale Fo to Fc." )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_phi<T> > &, const HKL_data<F_sigF<T> > & ) )&Class::operator(),
            nb::arg( "fo" ), nb::arg( "fc" ), "Scale Fc to Fo." )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_sigF<T> > & ) )&Class::operator(), nb::arg( "fo" ),
            "Scale Fo to isotropic (approx.)." );

  //.def( "__call__", [](Class &self, HKL_data<F_sigF<T>> &fo, HKL_data<F_phi<T>> &fc) {
  //      return self(fo, fc);
  //    }, nb::arg("fo"), nb::arg("fc"), "Scale Fo to Fc." )
  //.def( "__call__", [](Class &self, HKL_data<F_phi<T>> &fc, HKL_data<F_sigF<T>> &fo) {
  //      return self(fc, fo);
  //    }, nb::arg("fc"), nb::arg("fo"), "Scale Fc to Fo." )
  //.def( "__call__", [](Class &self, HKL_data<F_sigF<T>> &fo) { return self(fo); },
  //    nb::arg("fo"), "Scale Fo to isotropic (approx.)." );
} // declare_sfscale_base

template <class Derived, class Base, class D, class T1, class T2, class S>
void add_scale_specialization( nb::class_<Derived, Base> &pyclass, const char *basis_type ) {
  auto fn_name = std::string( "scale_" ) + basis_type;
  pyclass.def( fn_name.c_str(), []( Derived &self, HKL_data<D> &fo, const ftype resfilter, const int npar_scl ) {
    return self.template scale<D, T1, T2, S>( fo, resfilter, npar_scl );
  } );
  // (bool (Derived::*)(HKL_data<D>&, const ftype, const int))
  // &Derived::scale<D, T1, T2, S>);
}

template <class T> void declare_sfscale_aniso( nb::module_ &m, const char *dtype ) {
  using Class = SFscale_aniso<T>;
  auto pyclass_name = std::string( "SFscale_aniso_" ) + dtype;
  nb::class_<Class, SFscale_base<T>> sfscale_aniso( m, pyclass_name.c_str() );
  sfscale_aniso.doc() = "Structure factor anisotropic scaling.\n"
                        "Perform structure factor anisotropic scaling, observed to calculated, "
                        "calculated to observed, or observed against itself.";
  nb::enum_<typename Class::TYPE>( sfscale_aniso, "TYPE", "U_aniso_orth return types." )
      .value( "F", Class::TYPE::F )
      .value( "I", Class::TYPE::I )
      .export_values();
  nb::enum_<typename Class::MODE>( sfscale_aniso, "MODE", "Scaling modes." )
      .value( "NORMAL", Class::MODE::NORMAL )
      .value( "SHARPEN", Class::MODE::SHARPEN )
      .value( "UNSHARPEN", Class::MODE::UNSHARPEN )
      .export_values();

  sfscale_aniso
      .def( nb::init<ftype, typename Class::MODE>(), nb::arg( "nsig" ) = 0.0, nb::arg( "mode" ) = Class::NORMAL,
            "Constructor: takes rejection criterion for F/sigF." )
      .def(
          "__call__",
          []( Class &self, HKL_data<F_sigF<T>> &fo, const ftype resfilter, const int npar_scl ) {
            return self( fo, resfilter, npar_scl );
          },
          nb::arg( "fo" ), nb::arg( "resfilter" ), nb::arg( "npar_scl" ), "Scale Fo to isotropic (approx.)." )
      .def(
          "__call__",
          []( Class &self, HKL_data<I_sigI<T>> &io, const ftype resfilter, const int npar_scl ) {
            return self( io, resfilter, npar_scl );
          },
          nb::arg( "io" ), nb::arg( "resfilter" ), nb::arg( "npar_scl" ), "Scale Io to isotropic (approx.)." )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_sigF<T> > &, const HKL_data<F_phi<T> > & ) )&Class::operator(),
            nb::arg( "fo" ), nb::arg( "fc" ), "Scale Fo to Fc." )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_phi<T> > &, const HKL_data<F_sigF<T> > & ) )&Class::operator(),
            nb::arg( "fo" ), nb::arg( "fc" ), "Scale Fc to Fo." )
      .def( "__call__", ( bool ( Class::* )( HKL_data<F_sigF<T> > & ) )&Class::operator(), nb::arg( "fo" ),
            "Scale Fo to isotropic (approx.)." )
      .def( "u_aniso_orth", ( const U_aniso_orth &( Class::* )( typename Class::TYPE ) const ) & Class::u_aniso_orth,
            nb::arg( "type" ), "Return aniso correction on F or I." );

  add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>, TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                           TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>, BasisFn_binner>( sfscale_aniso, "binned" );
  add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>, TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                           TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>, BasisFn_linear>( sfscale_aniso, "linear" );
  add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>, TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                           TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>, BasisFn_spline>( sfscale_aniso, "spline" );

  add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>, TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                           TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>, BasisFn_binner>( sfscale_aniso, "binned" );
  add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>, TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                           TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>, BasisFn_linear>( sfscale_aniso, "linear" );
  add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>, TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                           TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>, BasisFn_spline>( sfscale_aniso, "spline" );
}

void init_sfscale( nb::module_ &m ) {
  declare_sfscale_base<ftype32>( m, "float" );
  declare_sfscale_aniso<ftype32>( m, "float" );

  declare_sfscale_base<ftype64>( m, "double" );
  declare_sfscale_aniso<ftype64>( m, "double" );
} // init_sfscale