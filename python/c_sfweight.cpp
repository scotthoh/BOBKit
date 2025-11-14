// Nanobind bindings for clipper contrib sfweight
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper sfweight by Tristan Croll

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <nanobind/stl/vector.h>
#include <nanobind/trampoline.h>

using namespace clipper;
using namespace clipper::datatypes;

// Trampoline Class, from function_object_bases.h
template <class T> class PySFweight_base : public SFweight_base<T> {
public:
  NB_TRAMPOLINE( SFweight_base<T>, 1 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( HKL_data<datatypes::F_phi<T> > &fb, HKL_data<datatypes::F_phi<T> > &fd,
                   HKL_data<datatypes::Phi_fom<T> > &phiw, const HKL_data<datatypes::F_sigF<T> > &fo,
                   const HKL_data<datatypes::F_phi<T> > &fc, const HKL_data<datatypes::Flag> &usage ) override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), fb, fd, phiw, fo, fc, usage );
  }
}; // declare trampoline class for sfweight base

template <class T> void declare_sfweight_base( nb::module_ &m, const char *dtype ) {
  using Class = SFweight_base<T>;
  auto pyclass_name = std::string( "_SFweight_base_" ) + dtype;
  nb::class_<Class, PySFweight_base<T>> sfweightbase( m, pyclass_name.c_str() );
  nb::enum_<typename Class::TYPE>( sfweightbase, "TYPE", "Flag values for different reflection purposes" )
      .value( "NONE", Class::TYPE::NONE )
      .value( "SIGMAA", Class::TYPE::SIGMAA )
      .value( "SCALE", Class::TYPE::SCALE )
      .value( "BOTH", Class::TYPE::BOTH )
      .export_values();

  sfweightbase.def(
      "__call__",
      []( Class &self, HKL_data<F_phi<T>> &fb, HKL_data<F_phi<T>> &fd, HKL_data<Phi_fom<T>> &phiw,
          const HKL_data<F_sigF<T>> &fo, const HKL_data<F_phi<T>> &fc,
          const HKL_data<Flag> &usage ) { return self( fb, fd, phiw, fo, fc, usage ); },
      "Structure factor weighting (sigmaa) definition" );
} // declare_sfweight_base

template <class T> void declare_sfweight_spline( nb::module_ &m, const char *dtype ) {
  using Class = SFweight_spline<T>;
  auto pyclass_name = std::string( "SFweight_spline_" ) + dtype;
  nb::class_<Class, SFweight_base<T>> sfweight_spline( m, pyclass_name.c_str() );
  using ResultStruc = typename SFweight_spline<T>::TargetResult;
  nb::class_<ResultStruc>( sfweight_spline, "TargetResult", "Class containing target results" )
      .def( nb::init<>() )
      .def_rw( "r", &ResultStruc::r )
      .def_rw( "ds", &ResultStruc::ds )
      .def_rw( "dw", &ResultStruc::dw )
      .def_rw( "dss", &ResultStruc::dss )
      .def_rw( "dww", &ResultStruc::dww )
      .def_rw( "dsw", &ResultStruc::dsw );

  sfweight_spline
      .def( nb::init<const int, const int, const int>(), nb::arg( "n_reflns" ) = 1000, nb::arg( "n_params" ) = 20,
            nb::arg( "n_phases" ) = 24, "Constructor" )
      .def( nb::init<HKL_data<F_phi<T>> &, HKL_data<F_phi<T>> &, HKL_data<Phi_fom<T>> &, const HKL_data<F_sigF<T>> &,
                     const HKL_data<F_phi<T>> &, const HKL_data<Flag> &, const int, const int>(),
            nb::arg( "fb" ), nb::arg( "fd" ), nb::arg( "phiw" ), nb::arg( "fo" ), nb::arg( "fc" ), nb::arg( "usage" ),
            nb::arg( "n_reflns" ) = 1000, nb::arg( "n_params" ) = 20 )
      .def( "init", &Class::init, nb::arg( "n_reflns" ) = 1000, nb::arg( "n_params" ) = 20, nb::arg( "n_phases" ) = 24 )
      .def( "__call__", []( Class &self, HKL_data<datatypes::F_phi<T> > &fb, HKL_data<datatypes::F_phi<T> > &fd,
                            HKL_data<datatypes::Phi_fom<T> > &phiw, const HKL_data<datatypes::F_sigF<T> > &fo0,
                            const HKL_data<datatypes::F_phi<T> > &fc0,
                            const HKL_data<datatypes::Flag> &usage ) { return self( fb, fd, phiw, fo0, fc0, usage ); } )
      .def( "__call__", []( Class &self, HKL_data<F_phi<T>> &fb, HKL_data<F_phi<T>> &fd, HKL_data<Phi_fom<T>> &phiw,
                            HKL_data<ABCD<T>> &hl, const HKL_data<F_sigF<T>> &fo0, const HKL_data<ABCD<T>> &hl0,
                            const HKL_data<F_phi<T>> &fc0,
                            const HKL_data<Flag> &usage ) { return self( fb, fd, phiw, hl, fo0, hl0, fc0, usage ); } )
      .def_prop_ro( "params_scale", &Class::params_scale )
      .def_prop_ro( "params_error", &Class::params_error )
      .def_prop_ro( "log_likelihood_work", &Class::log_likelihood_work )
      .def_prop_ro( "log_likelihood_free", &Class::log_likelihood_free )
      .def( "targetfn",
            []( const Class &self, const HKL_class cls, const F_sigF<T> &fo0, const F_phi<T> &fc0, const ftype &s,
                const ftype &w )
            /*-> nb::tuple*/
            {
              return self.targetfn( cls, fo0, fc0, s, w );
              // auto r = self.targetfn(cls, fo0, fc0, s, w);
              // return nb::make_tuple(r.r, r.ds, r.dw, r.dss, r.dww, r.dsw);
            } )
      .def( "targethl",
            []( const Class &self, const HKL_class cls, const F_sigF<T> &fo0, const ABCD<T> &hl0, const F_phi<T> &fc0,
                const ftype &s, const ftype &w )
            /* -> nb::tuple*/
            {
              return self.targethl( cls, fo0, hl0, fc0, s, w );
              // auto r = self.targethl(cls, fo0, hl0, fc0, s, w);
              // return nb::make_tuple(r.r, r.ds, r.dw, r.dss, r.dww, r.dsw);
            } )
      .doc() = "Structure factor weighting by sigmaa-related spline method.\n"
               "Perform structure factor weighting to obtain likelihood "
               "weights for structure factors.\n"
               "This implementation uses a single list of reflections for both "
               "scaling and sigmaa, thus the only relevent usage flags are "
               "NONE/BOTH.\n"
               "The number of spline parameters or the number of reflections per "
               "parameter may be specified. If either is zero, the other takes "
               "priority. If both are non-zero, a compromise value is used.";
}

void init_sfweight( nb::module_ &m ) {
  declare_sfweight_base<ftype32>( m, "float" );
  declare_sfweight_spline<ftype32>( m, "float" );

  declare_sfweight_base<ftype64>( m, "double" );
  declare_sfweight_spline<ftype64>( m, "double" );
}