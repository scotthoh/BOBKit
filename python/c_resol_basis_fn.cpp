// Nanobind bindings for clipper resolfn, resol_basis, and resol_target functions methods
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <clipper/core/resol_fn.h>
#include <clipper/core/resol_basisfn.h>
#include <clipper/core/resol_targetfn.h>

using namespace clipper;

// from resol_fn.h

void declare_basisfn_base( nb::module_& m ) {
  nb::class_<BasisFn_base> basisfn_base(
      m, "_BasisFn_base" );

  nb::enum_<BasisFn_base::FNtype>( basisfn_base, "FNtype" )
      .value( "GENERAL", BasisFn_base::FNtype::GENERAL )
      .value( "LINEAR", BasisFn_base::FNtype::LINEAR );
  // do we need basisfn_base?
  //basisfn_base.def( "f", &BasisFn_base::f )
  //    .def_prop_ro( "num_params", &BasisFn_base::num_params );
}

void declare_basisfn_fderiv( nb::module_& m ) {
  using Class = BasisFn_base::Fderiv;
  nb::class_<Class>( m, "BasisFn_Fderiv" )
      .def( nb::init<>(), "Null constructor" )
      .def( nb::init<const int&>(), "Constructor with number of parameters", nb::arg("num_params"))
      .def_ro( "f", &Class::f) //[]( const Class& self ) { return self.f; } )
      .def_ro( "df", &Class::df) //[]( const Class& self ) { return self.df; } )
      .def_ro( "df2", &Class::df2); //[]( const Class& self ));
}

void declare_targetfn_base( nb::module_& m ) {
  nb::class_<TargetFn_base> targetfn_base( m, "_TargetFn_base" );

  nb::enum_<TargetFn_base::FNtype>( targetfn_base, "FNtype" )
      .value( "GENERAL", TargetFn_base::FNtype::GENERAL )
      .value( "QUADRATIC", TargetFn_base::FNtype::QUADRATIC );
}

void declare_resolution_fn( nb::module_& m ) {
  nb::class_<ResolutionFn>( m, "ResolutionFn" )
      /*
      As of 26/10/2021, binding of STL containers passed by const ref appears to be a no-no -
      the temporary vector produced by the implicit conversion machinery is out of scope by the
      time it makes it to the function. Replaced with a pass by value alternative.
      */
      // .def(nb::init<const HKL_info&, const BasisFn_base&, const TargetFn_base&,
      //         const std::vector<ftype>&, const ftype, const bool>(),
      //     nb::arg("hkl_info"), nb::arg("basis_fn"), nb::arg("target_fn"),
      //     nb::arg("params"), nb::arg("damp") = 0.0, nb::arg("debug") = false)
      .def( "__init__",  [](ResolutionFn *rfn, const HKL_info& hkl_info, const BasisFn_base& basis_fn,
                          const TargetFn_base& target_fn, std::vector<ftype> params, ftype damp,
                          bool debug ) {
              new(rfn) ResolutionFn( hkl_info, basis_fn, target_fn, params, damp, debug ) ;
            } ,
            nb::arg( "hkl_info" ), nb::arg( "basis_fn" ), nb::arg( "target_fn" ),
            nb::arg( "params" ), nb::arg( "damp" ) = 0.0, nb::arg( "debug" ) = false )
      .def( "f", &ResolutionFn::f )
      .def_prop_ro( "params", &ResolutionFn::params );
}

void declare_resolution_fn_nonlinear( nb::module_& m ) {
  nb::class_<ResolutionFn_nonlinear, ResolutionFn>( m, "ResolutionFn_nonlinear" )
      /*
      As of 26/10/2021, binding of STL containers passed by const ref appears to be a no-no -
      the temporary vector produced by the implicit conversion machinery is out of scope by the
      time it makes it to the function. Replaced with a pass by value alternative.
      */
      // .def(nb::init<const HKL_info&, const BasisFn_base&, const TargetFn_base&,
      //     const std::vector<ftype>&, const ftype, const bool >(),
      //     nb::arg("hkl_info"), nb::arg("basis_fn"), nb::arg("target_fn"),
      //     nb::arg("params"), nb::arg("damp") = 0.0, nb::arg("debug") = false)
      .def( "__init__", [](ResolutionFn_nonlinear* rfn, const HKL_info& hkl_info, const BasisFn_base& basis_fn,
                          const TargetFn_base& target_fn, std::vector<double> params, double damp,
                          bool debug ) {
              new (rfn) ResolutionFn_nonlinear(
                  hkl_info, basis_fn, target_fn, params, damp, debug );
            },
            nb::arg( "hkl_info" ), nb::arg( "basis_fn" ), nb::arg( "target_fn" ),
            nb::arg( "params" ), nb::arg( "damp" ) = 0.0, nb::arg( "debug" ) = false );
}

// from resol_basisfn.h
void declare_resolution_ordinal( nb::module_& m ) {
  nb::class_<Resolution_ordinal, Generic_ordinal>( m, "Resolution_ordinal" )
    .def( nb::init<>(), "Null constructor")
    .def( "__init__", [](Resolution_ordinal* ro, const HKL_info& hklinfo, const ftype& power ) {
      new(ro) Resolution_ordinal();
      ro->init( hklinfo, power );
    } ) 
    .def( "__init__", [](Resolution_ordinal* ro, const HKL_data_base& hkldata, const ftype& power ) {
      new(ro) Resolution_ordinal();
      ro->init( hkldata, power );
    } )
    .def( "__init__", [](Resolution_ordinal* ro, const HKL_data_base& hkldata, const Cell& cell, const ftype& power ) {
      new(ro) Resolution_ordinal();
      ro->init( hkldata, cell, power );
    } )
    .def( "init", ( void( Resolution_ordinal::* )( const HKL_info&, const ftype& ) ) &
                      Resolution_ordinal::init )
    .def( "init", ( void( Resolution_ordinal::* )( const HKL_data_base&, const ftype& ) ) &
                      Resolution_ordinal::init )
    .def( "init",
          ( void( Resolution_ordinal::* )( const HKL_data_base&, const Cell&, const ftype& ) ) &
              Resolution_ordinal::init );
}

template <class C, class B>
void add_basisfn_base_class_functions( nb::class_<C, B>& pyclass ) {
  pyclass.def( "f", &C::f )
      .def( "fderiv", &C::fderiv )
      .def_prop_ro( "num_params", &C::num_params )
      .def_prop_ro( "type", &C::type )
      .def_prop_ro( "num_diagonals", &C::num_diagonals );
}

void declare_basisfn_binner( nb::module_& m ) {
  nb::class_<BasisFn_binner, BasisFn_base> basis_binner( m, "BasisFn_binner" );
  basis_binner
      .def( nb::init<const HKL_info&, const int&, const ftype>(), nb::arg( "hkl_info" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( nb::init<const HKL_data_base&, const int&, const ftype>(), nb::arg( "hkl_data" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( "f_s", &BasisFn_binner::f_s )
      .def( "fderiv_s", &BasisFn_binner::fderiv_s );
  add_basisfn_base_class_functions<BasisFn_binner, BasisFn_base>( basis_binner );
}

void declare_basisfn_linear( nb::module_& m ) {
  nb::class_<BasisFn_linear, BasisFn_base> basis_linear( m, "BasisFn_linear" );
  basis_linear
      .def( nb::init<const HKL_info&, const int&, const ftype>(), nb::arg( "hkl_info" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( nb::init<const HKL_data_base&, const int&, const ftype>(), nb::arg( "hkl_data" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( "f_s", &BasisFn_linear::f_s )
      .def( "fderiv_s", &BasisFn_linear::fderiv_s );
  add_basisfn_base_class_functions<BasisFn_linear, BasisFn_base>( basis_linear );
}

void declare_basisfn_spline( nb::module_& m ) {
  nb::class_<BasisFn_spline, BasisFn_base> basis_spline( m, "BasisFn_spline" );
  basis_spline
      .def( nb::init<const HKL_info&, const int&, const ftype>(), nb::arg( "hkl_info" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( nb::init<const HKL_data_base&, const int&, const ftype>(), nb::arg( "hkl_data" ),
            nb::arg( "num_bins" ), nb::arg( "power" ) = 1.0 )
      .def( "f_s", &BasisFn_spline::f_s )
      .def( "fderiv_s", &BasisFn_spline::fderiv_s );
  add_basisfn_base_class_functions<BasisFn_spline, BasisFn_base>( basis_spline );
}

void declare_basisfn_gaussian( nb::module_& m ) {
  nb::class_<BasisFn_gaussian, BasisFn_base>( m, "BasisFn_gaussian" )
      .def( nb::init<>() )
      .def( "fderiv_s", &BasisFn_gaussian::fderiv_s )
      .def( "fderiv", &BasisFn_gaussian::fderiv )
      .def( "scale", &BasisFn_gaussian::scale )
      .def( "u_iso", &BasisFn_gaussian::u_iso );
}

void declare_basisfn_aniso_gaussian( nb::module_& m ) {
  nb::class_<BasisFn_aniso_gaussian, BasisFn_base>( m, "BasisFn_aniso_gaussian" )
      .def( nb::init<>() )
      .def( "fderiv_coord", &BasisFn_aniso_gaussian::fderiv_coord )
      .def( "fderiv", &BasisFn_aniso_gaussian::fderiv )
      .def( "scale", &BasisFn_aniso_gaussian::scale )
      .def( "u_aniso_orth", &BasisFn_aniso_gaussian::u_aniso_orth );
}

void declare_basisfn_log_gaussian( nb::module_& m ) {
  nb::class_<BasisFn_log_gaussian, BasisFn_base>( m, "BasisFn_log_gaussian" )
      .def( nb::init<>() )
      .def( "fderiv_s", &BasisFn_log_gaussian::fderiv_s )
      .def( "fderiv", &BasisFn_log_gaussian::fderiv )
      .def_prop_ro( "type", &BasisFn_log_gaussian::type )
      .def( "scale", &BasisFn_log_gaussian::scale )
      .def( "u_iso", &BasisFn_log_gaussian::u_iso );
}

void declare_basisfn_log_aniso_gaussian( nb::module_& m ) {
  nb::class_<BasisFn_log_aniso_gaussian, BasisFn_base>( m, "BasisFn_log_aniso_gaussian" )
      .def( nb::init<>() )
      .def( "fderiv_coord", &BasisFn_log_aniso_gaussian::fderiv_coord )
      .def( "fderiv", &BasisFn_log_aniso_gaussian::fderiv )
      .def_prop_ro( "type", &BasisFn_log_aniso_gaussian::type )
      .def( "scale", &BasisFn_log_aniso_gaussian::scale )
      .def( "u_aniso_orth", &BasisFn_log_aniso_gaussian::u_aniso_orth );
}

void declare_basisfn_expcubi( nb::module_& m ) {
  nb::class_<BasisFn_expcubic, BasisFn_base>( m, "BasisFn_expcubic" )
      .def( nb::init<>() , "Default constructor" )
      .def( "fderiv_s", &BasisFn_expcubic::fderiv_s)
      .def( "fderiv", &BasisFn_expcubic::fderiv );
}
// From resol_targetfn.h

template <class T>
void declare_targetfn_meanfnth( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_meanFnth_" ) + suffix;
  using Class = TargetFn_meanFnth<T>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T>&, const ftype&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T>
void declare_targetfn_meanenth( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_meanEnth_" ) + suffix;
  using Class = TargetFn_meanEnth<T>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T>&, const ftype&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T1, class T2>
void declare_targetfn_scalef1f2( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_scaleF1F2_" ) + suffix;
  using Class = TargetFn_scaleF1F2<T1, T2>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T1>&, const HKL_data<T2>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T1, class T2>
void declare_targetfn_scalelogf1f2( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_scaleLogF1F2_" ) + suffix;
  using Class = TargetFn_scaleLogF1F2<T1, T2>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T1>&, const HKL_data<T2>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T1, class T2>
void declare_targetfn_scalei1i2( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_scaleI1I2_" ) + suffix;
  using Class = TargetFn_scaleI1I2<T1, T2>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T1>&, const HKL_data<T2>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T1, class T2>
void declare_targetfn_scalelogi1i2( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_scaleLogI1I2_" ) + suffix;
  using Class = TargetFn_scaleLogI1I2<T1, T2>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T1>&, const HKL_data<T2>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T>
void declare_targetfn_scaleesq( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_scaleEsq_" ) + suffix;
  using Class = TargetFn_scaleEsq<T>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_prop_ro( "type", &Class::type );
}

template <class T>
void declare_targetfn_sigmaa_omegaa( nb::module_& m, const std::string& suffix ) {
  std::string pyname = std::string( "TargetFn_sigmaa_omegaa_" ) + suffix;
  using Class = TargetFn_sigmaa_omegaa<T>;
  nb::class_<Class, TargetFn_base>( m, pyname.c_str() )
      .def( nb::init<const HKL_data<T>&, const HKL_data<T>&>() )
      .def( "rderiv", &Class::rderiv )
      .def_static( "sigmaa", &Class::sigmaa );
}

void init_resol_fn( nb::module_& m ) {
  declare_basisfn_base( m );
  declare_basisfn_fderiv( m );
  declare_targetfn_base( m );
  declare_resolution_fn( m );
  declare_resolution_fn_nonlinear( m );
  declare_resolution_ordinal( m );
  declare_basisfn_binner( m );
  declare_basisfn_linear( m );
  declare_basisfn_spline( m );
  declare_basisfn_gaussian( m );
  declare_basisfn_aniso_gaussian( m );
  declare_basisfn_log_gaussian( m );
  declare_basisfn_log_aniso_gaussian( m );

  declare_targetfn_meanfnth<data32::F_sigF>( m, "F_sigF_float" );
  declare_targetfn_meanfnth<data32::F_sigF_ano>( m, "F_sigF_ano_float" );
  declare_targetfn_meanfnth<data32::F_phi>( m, "F_phi_float" );
  declare_targetfn_meanenth<data32::E_sigE>( m, "E_sigE_float" );

  declare_targetfn_scalef1f2<data32::F_sigF, data32::F_sigF>( m, "F_sigF_float" );
  declare_targetfn_scalef1f2<data32::F_phi, data32::F_phi>( m, "F_phi_float" );
  declare_targetfn_scalelogf1f2<data32::F_sigF, data32::F_sigF>( m, "F_sigF_float" );
  declare_targetfn_scalelogf1f2<data32::F_phi, data32::F_phi>( m, "F_phi_float" );
  declare_targetfn_scalelogf1f2<data32::F_phi, data32::F_sigF>( m, "F_phi_F_sigF_float" );
  declare_targetfn_scalei1i2<data32::I_sigI, data32::I_sigI>( m, "I_sigI_float" );
  declare_targetfn_scalelogi1i2<data32::I_sigI, data32::I_sigI>( m, "I_sigI_float" );

  declare_targetfn_scaleesq<data32::E_sigE>( m, "E_sigE_float" );

  declare_targetfn_sigmaa_omegaa<data32::E_sigE>( m, "E_sigE_float" );

  declare_targetfn_meanfnth<data64::F_sigF>( m, "F_sigF_double" );
  declare_targetfn_meanfnth<data64::F_sigF_ano>( m, "F_sigF_ano_double" );
  declare_targetfn_meanfnth<data64::F_phi>( m, "F_phi_double" );
  declare_targetfn_meanenth<data64::E_sigE>( m, "E_sigE_double" );

  declare_targetfn_scalef1f2<data64::F_sigF, data64::F_sigF>( m, "F_sigF_double" );
  declare_targetfn_scalelogf1f2<data64::F_sigF, data64::F_sigF>( m, "F_sigF_double" );
  declare_targetfn_scalef1f2<data64::F_phi, data64::F_phi>( m, "F_phi_double" );
  declare_targetfn_scalelogf1f2<data64::F_phi, data64::F_phi>( m, "F_phi_double" );

  declare_targetfn_scalei1i2<data64::I_sigI, data64::I_sigI>( m, "I_sigI_double" );
  declare_targetfn_scalelogi1i2<data64::I_sigI, data64::I_sigI>( m, "I_sigI_double" );

  declare_targetfn_scaleesq<data64::E_sigE>( m, "E_sigE_double" );

  declare_targetfn_sigmaa_omegaa<data64::E_sigE>( m, "E_sigE_double" );
}