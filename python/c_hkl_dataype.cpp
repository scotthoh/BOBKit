// Nanobind bindings for clipper hkl_datatypes
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/unique_ptr.h>

using namespace clipper;
using namespace clipper::datatypes;

#define GETSET( CLASS, FUNCNAME )                                                                                      \
  []( const CLASS &self ) { return self.FUNCNAME(); }, []( CLASS &self, ftype &val ) { self.FUNCNAME() = val; }

template <class C> void catch_null( const C &c ) {
  if ( c.is_null() )
    throw std::length_error( "Array is not initialised!" );
}

template <class HKLtype, class dtype> auto hkldatatype_to_array( const HKLtype &self ) {
  dtype *data = new dtype[self.data_size()];
  self.data_export( data );
  // Delete 'data' when the 'owner' capsule expires
  nb::capsule owner( data, []( void *p ) noexcept { delete[] ( dtype * )p; } );
  return nb::ndarray<nb::numpy, dtype, nb::ndim<1>>( data, { ( size_t )self.data_size() }, owner ).cast();
}

// template<class C, typename F, class dtype, class... Sources,  typename...
// Args> void safe_compute(C& target, const F& func, const Sources&... sources,
// const Args&... args)
// {
//     catch_null(target);
//     target.compute(sources, func(args));
// }

template <class Derived> void declare_base_methods( nb::class_<Derived> pyclass ) {
  pyclass.def( "set_null", &Derived::set_null, "Set data to null." )
      .def_prop_ro_static(
          "type", []( nb::handle /*unused*/ ) { return Derived::type(); },
          "Return name(or a list of names) for this data type." )
      .def( "friedel", &Derived::friedel, "Applies Friedel transformation." )
      .def( "shift_phase", &Derived::shift_phase, "Applies phase shift transformation." )
      .def( "missing", &Derived::missing, "Checks if data is present. Return true if data is missing." )
      .def_prop_ro_static(
          "data_size", []( nb::handle /*unused*/ ) { return Derived::data_size(); },
          "Return number of data elements in this type." )
      .def( "__len__", &Derived::data_size )
      .def_prop_ro_static(
          "data_names", []( nb::handle /*unused*/ ) { return Derived::data_names(); },
          "Names of data elements in this type." )
      // To/from numpy
      .def_prop_ro(
          "data",
          []( const Derived &self ) {
            return hkldatatype_to_array<Derived, xtype>( self );
            // return hkl_data_export_numpy<Derived, xtype>(self,
            //                                              self.data_size());
          },
          nb::rv_policy::reference_internal ) //,
      //    [](Derived &self, nb::array_t<xtype> vals) {
      //      hkl_data_import_numpy<Derived, xtype>(self, self.data_size(), vals);
      //    },
      //    "HKL data export/import to/from numpy.", )
      // Python methods common to all
      .def( "clone", []( const Derived &self ) { return new Derived( self ); }, "Return a copy." );
} // declare_base_methods

template <class T> void declare_i_sigi( nb::module_ &m, const char *dtype ) {
  using Class = I_sigI<T>;
  auto pyclass_name = std::string( "I_sigI_" ) + dtype;
  nb::class_<Class> isigi( m, pyclass_name.c_str() );
  isigi.def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "i" ), nb::arg( "sigi" ), "Constructor from I and sigI." )
      .def( "scale", &Class::scale, nb::arg( "s" ), "Apply magnitude scale factor." )
      .def_prop_rw( "i", GETSET( Class, I ), "Read/write access to I." )
      .def_prop_rw( "sigi", GETSET( Class, sigI ), "Read/write access to sigI." )
      .def_prop_ro( "i_pl", &Class::I_pl, "Read access to as anom." )
      .def_prop_ro( "sigi_pl", &Class::sigI_pl, "Read access as anom." )
      .def_prop_ro( "i_mi", &Class::I_mi, "Read access as anom." )
      .def_prop_ro( "sigi_mi", &Class::sigI_mi, "Read access as anom." )
      .def_prop_ro( "cov", &Class::cov, "Read access as anom." )
      .doc() = "Reflection data type: I + sigI.\n"
               "Note that I_sigI also has methods for returning I_pl(), "
               "sigI_pl(), I_mi, sigI_mi(), so you can use this type in any "
               "template type where you would use I_sigI_ano.";
  declare_base_methods<Class>( isigi );
} // declare_i_sigi

template <class T> void declare_i_sigi_ano( nb::module_ &m, const char *dtype ) {
  using Class = I_sigI_ano<T>;
  auto pyclass_name = std::string( "I_sigI_ano_" ) + dtype;
  nb::class_<Class> isigi_ano( m, pyclass_name.c_str() );
  isigi_ano.def( nb::init<>() )
      .def( "scale", &Class::scale, nb::arg( "s" ), "Apply magnitude scale factor." )
      .def_prop_rw( "i_pl", GETSET( Class, I_pl ), "Read/write accessor to I+." )
      .def_prop_rw( "sigi_pl", GETSET( Class, sigI_pl ) )
      .def_prop_rw( "i_mi", GETSET( Class, I_mi ), "Read/write accessor to I-." )
      .def_prop_rw( "sigi_mi", GETSET( Class, sigI_mi ), "Read/write accessor to sigI-." )
      .def_prop_rw( "cov", GETSET( Class, cov ), "Read/write accessor to covI+-." )
      .def_prop_ro( "i", &Class::I, "Read access as simple." )
      .def_prop_ro( "sigi", &Class::sigI, "Read access as simple." )
      .doc() = "Reflection data type: I(+) I(+) sigI(+) sigI(-) cov+- .\n"
               "Note that I_sigI_ano also has methods for returning I(), "
               "sigI(), so you can use this type in any template type "
               "where you would use I_sigI.";
  declare_base_methods<Class>( isigi_ano );
} // declare_i_sigi_ano

template <class T> void declare_f_sigf( nb::module_ &m, const char *dtype ) {
  using Class = F_sigF<T>;
  auto pyclass_name = std::string( "F_sigF_" ) + dtype;
  nb::class_<Class> fsigf( m, pyclass_name.c_str() );
  fsigf.def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "f" ), nb::arg( "sigf" ), "Constructor from F and sigF." )
      .def( "scale", &Class::scale, nb::arg( "s" ), "Apply magnitude scale factor." )
      .def_prop_rw( "f", GETSET( Class, f ), "Read/write access to F." )
      .def_prop_rw( "sigf", GETSET( Class, sigf ), "Read/write access to sigF." )
      .def_prop_ro( "f_pl", &Class::f_pl, "Read access as anom." )
      .def_prop_ro( "sigf_pl", &Class::sigf_pl, "Read access as anom." )
      .def_prop_ro( "f_mi", &Class::f_mi, "Read access as anom." )
      .def_prop_ro( "sigf_mi", &Class::sigf_mi, "Read access as anom." )
      .def_prop_ro( "cov", &Class::cov, "Read access as anom." )
      .doc() = "Reflection data type: F + sigF.\n"
               "Note that F_sigF also has methods for returning f_pl(), "
               "sigf_pl(), f_mi, sigf_mi(), so you can use this type in any "
               "template type where you would use F_sigF_ano.";
  declare_base_methods<Class>( fsigf );
} // declare_f_sigf

template <class T> void declare_f_sigf_ano( nb::module_ &m, const char *dtype ) {
  using Class = F_sigF_ano<T>;
  auto pyclass_name = std::string( "F_sigF_ano_" ) + dtype;
  nb::class_<Class> fsigf_ano( m, pyclass_name.c_str() );
  fsigf_ano.def( nb::init<>() )
      .def( "scale", &Class::scale, nb::arg( "s" ) )
      .def_prop_rw( "f_pl", GETSET( Class, f_pl ) )
      .def_prop_rw( "sigf_pl", GETSET( Class, sigf_pl ) )
      .def_prop_rw( "f_mi", GETSET( Class, f_mi ) )
      .def_prop_rw( "sigf_mi", GETSET( Class, sigf_mi ) )
      .def_prop_rw( "cov", GETSET( Class, cov ) )
      .def_prop_ro( "f", &Class::f )
      .def_prop_ro( "sigf", &Class::sigf )
      .doc() = "Reflection data type: F(+) F(+) sigF(+) sigF(-) cov+- .\n"
               "Note that F_sigF_ano also has methods for returning f(), "
               "sigf(), so you can use this type in any template type "
               "where you would use F_sigF. ";
  declare_base_methods<Class>( fsigf_ano );
} // declare_f_sigf_ano

template <class T> void declare_e_sige( nb::module_ &m, const char *dtype ) {
  using Class = E_sigE<T>;
  auto pyclass_name = std::string( "E_sigE_" ) + dtype;
  nb::class_<Class> esige( m, pyclass_name.c_str() );
  esige.def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "e" ), nb::arg( "sige" ), "Constructor from E, sigE." )
      .def_prop_rw( "e", GETSET( Class, E ), "Read/write accessor for E." )
      .def_prop_rw( "sige", GETSET( Class, sigE ), "Read/write accessor for sigE." )
      .def_prop_ro( "e_pl", &Class::E_pl, "Read access as anom." )
      .def_prop_ro( "sige_pl", &Class::sigE_pl, "Read access as anom." )
      .def_prop_ro( "e_mi", &Class::E_mi, "Read access as anom." )
      .def_prop_ro( "sige_mi", &Class::sigE_mi, "Read access as anom." )
      .def_prop_ro( "cov", &Class::cov, "Read access as anom." )
      .doc() = "Reflection data type: E + sigE.\n This is not strictly a type "
               "for storing E values, but rather a type for storing any "
               "structure factor magnitude-like quantity which has already had "
               "a symmetry enhancement factor (epsilon) removed from it. E's "
               "are most commonly stored in this form, wheras F's and U's are not.";
  declare_base_methods<Class>( esige );
} // declare_e_sige

// E_sigE_ano template is not currently instantiated in the Clipper core
// template<class T>
// void declare_e_sige_ano(nb::module_& m, const char* dtype)
// {
//     using Class=E_sigE_ano<T>;
//     auto pyclass_name=std::string("E_sigE_ano_") + dtype;
//     nb::class_<Class> esige_ano(m, pyclass_name.c_str());
//     esige_ano
//         .def(nb::init<>())
//         .def("scale", &Class::scale)
//         .def_prop_rw("e_pl", GETSET(Class, E_pl))
//         .def_prop_rw("sige_pl", GETSET(Class, sigE_pl))
//         .def_prop_rw("e_mi", GETSET(Class, E_mi))
//         .def_prop_rw("sige_mi", GETSET(Class, sigE_mi))
//         .def_prop_rw("cov", GETSET(Class, cov))
//         .def_prop_ro("e", &Class::E)
//         .def_prop_ro("sige", &Class::sigE)
//         ;
//     declare_base_methods<Class>(esige_ano);
// } // declare_e_sige_ano

template <class T> void declare_f_phi( nb::module_ &m, const char *dtype ) {
  using Class = F_phi<T>;
  auto pyclass_name = std::string( "F_phi_" ) + dtype;
  nb::class_<Class> fphi( m, pyclass_name.c_str() );
  fphi.def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "f" ), nb::arg( "phi" ), "Constructor from F,Phi." )
      .def( nb::init<const std::complex<T>>(), nb::arg( "c" ), "Convert from complex." )
      .def( "scale", &Class::scale, nb::arg( "s" ), "Apply magnitude scale factor." )
      .def_prop_rw( "f", GETSET( Class, f ), "Read/write accessor for F." )
      .def_prop_rw( "phi", GETSET( Class, phi ), "Read/write accessor for Phi." )
      .def_prop_ro( "a", &Class::a, "Read real part." )
      .def_prop_ro( "b", &Class::b, "Read imaginary part." )
      .def_prop_ro(
          "complex", []( const Class &self ) { return std::complex<T>( self ); }, "Convert to complex." )
      .def( "resolve", &Class::resolve, nb::arg( "phi" ), "Resolve along phase direction." )
      .def( "norm", &Class::norm, "Tidy up so that real part is positive and phase is 0...twopi." )
      // from hkl_operators.h
      .def( nb::self + nb::self )
      .def( nb::self - nb::self )
      .def( -nb::self )
      .doc() = "Reflection data type: F + phi model or map coeff "
               "(e.g. Fcalc, Fbest).";
  declare_base_methods<Class>( fphi );
} // declare_f_phi

template <class T> void declare_phi_fom( nb::module_ &m, const char *dtype ) {
  using Class = Phi_fom<T>;
  auto pyclass_name = std::string( "Phi_fom_" ) + dtype;
  nb::class_<Class> phifom( m, pyclass_name.c_str() );
  phifom.def( nb::init<>() )
      .def( nb::init<const T &, const T &>(), nb::arg( "phi" ), nb::arg( "fom" ), "Constructor from Phi,FOM." )
      .def_prop_rw( "phi", GETSET( Class, phi ), "Read/write accessor for Phi." )
      .def_prop_rw( "fom", GETSET( Class, fom ), "Read/write accessor for FOM." )
      .doc() = "Reflection data type: best phi + fom.";
  declare_base_methods<Class>( phifom );
} // declare_phi_fom

template <class T> void declare_abcd( nb::module_ &m, const char *dtype ) {
  using Class = ABCD<T>;
  auto pyclass_name = std::string( "ABCD_" ) + dtype;
  nb::class_<Class> abcd( m, pyclass_name.c_str() );
  abcd.def( nb::init<>() )
      .def( nb::init<const T &, const T &, const T &, const T &>(), nb::arg( "a" ), nb::arg( "b" ), nb::arg( "c" ),
            nb::arg( "d" ), "Constructor from ABCD." )
      .def_prop_rw( "a", GETSET( Class, a ), "Read/write accessor for A." )
      .def_prop_rw( "b", GETSET( Class, b ), "Read/write accessor for B." )
      .def_prop_rw( "c", GETSET( Class, c ), "Read/write accessor for C." )
      .def_prop_rw( "d", GETSET( Class, d ), "Read/write accessor for D." )
      .def( nb::self + nb::self ) // from hkl_operators.h
      .doc() = "Reflection data type: Hendrickson-Lattman coeff.";
  declare_base_methods<Class>( abcd );
} // declare_abcd

void declare_flag( nb::module_ &m ) {
  nb::class_<Flag> flag( m, "Flag", "Reflection data type: Free-R flag." );
  flag.def( nb::init<>() )
      .def( nb::init<const int &>(), nb::arg( "flag" ), "Constructor from flag (int)." )
      .def_prop_rw( "flag", GETSET( Flag, flag ), "Read/write accessor to flag." );
  declare_base_methods<Flag>( flag );
} // declare_flag

void declare_flag_bool( nb::module_ &m ) {
  nb::class_<Flag_bool> flag_bool( m, "Flag_bool", "Reflection data type: boolean (false = missing)." );
  flag_bool.def( nb::init<>() ).def_prop_rw( "flag", GETSET( Flag_bool, flag ), "Read/write acessor to flag bool." );
  declare_base_methods<Flag_bool>( flag_bool );
} // declare_flag_bool

template <class T> void declare_d_sigd( nb::module_ &m, const char *dtype ) {
  using Class = D_sigD<T>;
  auto pyclass_name = std::string( "D_sigD_" ) + dtype;
  nb::class_<Class> dsigd( m, pyclass_name.c_str(),
                           "Deprecated anomalous difference class for backward "
                           "compatibility only. Do not use." );
  dsigd.def( nb::init<>() )
      .def( nb::init<const T &, const T &>() )
      .def( "scale", &Class::scale )
      .def_prop_rw( "d", GETSET( Class, d ) )
      .def_prop_rw( "sigd", GETSET( Class, sigd ) );
  declare_base_methods<Class>( dsigd );
} // declare_flag_bool

// Actual wrapper initialisation
void init_hkl_datatypes( nb::module_ &m ) {
  // Non-floating-point datatypes go in the main module
  declare_flag( m );
  declare_flag_bool( m );
  {
    using namespace clipper::data32;
    // 32-bit floating point datatypes go in the data32 module
    declare_i_sigi<ftype32>( m, "float" );
    declare_i_sigi_ano<ftype32>( m, "float" );
    declare_f_sigf<ftype32>( m, "float" );
    declare_f_sigf_ano<ftype32>( m, "float" );
    declare_e_sige<ftype32>( m, "float" );
    declare_f_phi<ftype32>( m, "float" );
    declare_phi_fom<ftype32>( m, "float" );
    declare_abcd<ftype32>( m, "float" );
    declare_d_sigd<ftype32>( m, "float" );
  }
  {
    using namespace clipper::data64;
    // 64-bit floating point datatypes go in the data64 module
    declare_i_sigi<ftype64>( m, "double" );
    declare_i_sigi_ano<ftype64>( m, "double" );
    declare_f_sigf<ftype64>( m, "double" );
    declare_f_sigf_ano<ftype64>( m, "double" );
    declare_e_sige<ftype64>( m, "double" );
    declare_f_phi<ftype64>( m, "double" );
    declare_phi_fom<ftype64>( m, "double" );
    declare_abcd<ftype64>( m, "double" );
    declare_d_sigd<ftype64>( m, "double" );
  }
}