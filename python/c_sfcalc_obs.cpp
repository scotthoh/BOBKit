// Nanobind bindings for clipper contrib sfcalc obs
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper sfcalc obs by Tristan Croll

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <nanobind/trampoline.h>

using namespace clipper;
using namespace clipper::datatypes;

// Trampoline Class, from function_object_bases.h
template <class T> class PySFcalc_obs_base : public SFcalc_obs_base<T> {
public:
  NB_TRAMPOLINE( SFcalc_obs_base<T>, 1 );

  /* Trampoline (need one for each virtual function) */
  bool operator()( HKL_data<F_phi<T> > &fphidata, const HKL_data<F_sigF<T> > &fsig, const Atom_list &atoms ) override {
    NB_OVERRIDE_PURE_NAME( "__call__", operator(), fphidata, fsig, atoms );
  }
};

template <class T> void declare_sfcalc_obs_base( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_obs_base<T>;
  auto pyclass_name = std::string( "_SFcalc_obs_base_" ) + dtype;
  nb::class_<Class, PySFcalc_obs_base<T> >( m, pyclass_name.c_str() )
      .def(
          "__call__",
          []( Class &self, HKL_data<F_phi<T> > &fphi, const HKL_data<F_sigF<T> > &fsigf, const Atom_list &atoms ) {
            return self( fphi, fsigf, atoms );
          },
          nb::arg( "fphiout" ), nb::arg( "fsigf" ), nb::arg( "atoms" ),
          "Structure factor calculation function definition" )
      .doc() = "Base class for structure factor calculation vs observed methods";
} // declare_sfcalc_obs_base

template <class T> void declare_sfcalc_obs_bulk( nb::module_ &m, const char *dtype ) {
  using Class = SFcalc_obs_bulk<T>;
  auto pyclass_name = std::string( "SFcalc_obs_bulk_" ) + dtype;
  nb::class_<Class, SFcalc_obs_base<T>>( m, pyclass_name.c_str() )
      .def( nb::init<const int>(), nb::arg( "n_params" ) = 12, "Constructor" )
      .def( nb::init<HKL_data<F_phi<T>> &, const HKL_data<F_sigF<T>> &, const Atom_list &, const int>(),
            nb::arg( "fphiout" ), nb::arg( "fsigf" ), nb::arg( "atoms" ), nb::arg( "n_params" ) = 12,
            "Constructor: shorthand for constructor+operator" )
      .def(
          "__call__",
          []( Class &self, HKL_data<F_phi<T> > &fphi, const HKL_data<F_sigF<T> > &fsigf, const Atom_list &atoms ) {
            return self( fphi, fsigf, atoms );
          },
          nb::arg( "fphiout" ), nb::arg( "fsigf" ), nb::arg( "atoms" ), "Calculate structure factor" )
      .def_prop_ro( "bulk_frac", &Class::bulk_frac, "Return bulk fraction" )
      .def_prop_ro( "bulk_scale", &Class::bulk_scale, "Return bulk scale" )
      .doc() = "Structure factor calculation vs observed using bulk solvent.\n"
               "Perform structure factor calculation, adding an additional bulk "
               "solvent/missing correction to give best fit to observed data. "
               "A scaling is used internally, but not output.";
}

void init_sfcalc_obs( nb::module_ &m ) {
  declare_sfcalc_obs_base<ftype32>( m, "float" );
  declare_sfcalc_obs_bulk<ftype32>( m, "float" );

  declare_sfcalc_obs_base<ftype64>( m, "double" );
  declare_sfcalc_obs_bulk<ftype64>( m, "double" );

} // init_sfcalc_obs