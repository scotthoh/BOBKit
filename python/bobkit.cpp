// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

//#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "version.hpp" // for bobkit version
#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/stl/unique_ptr.h>

NB_MODULE( bobkit_ext, mbk_ ) {
  ( void ) mbk_;
  nb::module_ mbk = nb::module_::import_( "bobkit" );
  mbk.doc() =
      "Python bindings to Buccaneer and a subset of Clipper library for biomacromolecular "
      "model building kit (BOBKit).";
  mbk.attr("__version__") = BOBKIT_VERSION;
  nb::class_<clipper::Message_fatal>( mbk, "Message_fatal" )
    .def_prop_ro( "text", &clipper::Message_fatal::text )
    .def_prop_ro( "level", &clipper::Message_fatal::level );
  // auto package = pybind11::module::import("gemmi");
  // auto module = package.attr("Structure");
  // mbk.add_object("Structure", module);

  // try using handling custom exceptions
  // not sure how to solve this SystemError: nanobind::detail::nb_func_error_except(): exception could not be
  // translated!
  nb::register_exception_translator( [](const std::exception_ptr &p, void * ) {
    try {
      if ( p ) std::rethrow_exception( p );
    } catch ( const clipper::Message_fatal &e ) {
      PyErr_SetString( PyExc_RuntimeError, e.text().c_str() );
    } catch ( const std::system_error &e ) {
      const int errornum = e.code().value();
      PyErr_SetObject( PyExc_IOError, nb::make_tuple( errornum, e.what() ).ptr() );
    }
  } );

  mbk.def( "_long_running_func", []() {
    for (;;) {
      if ( PyErr_CheckSignals() != 0 )
        throw nb::python_error();
      // Long running iteration
    }
  } );

  auto mc = mbk.def_submodule( "clipper" );
  auto mb = mbk.def_submodule( "buccaneer" );
  auto mpdb = mbk.def_submodule( "protein_db" );
  auto mutil = mbk.def_submodule( "util" );
  auto mdata = mbk.def_submodule( "data", "Namespace for test data used in Test_core" );
  auto mm = mbk.def_submodule( "mm", "Dummy namespace to hold search modes." );
  add_clipper_util( mc );
  add_clipper_types( mc );
  add_thread_base( mc );
  add_coords( mc );
  add_atomlist( mc );
  add_cell( mc );
  add_clipper_memory( mc );
  add_ramachandran( mc );
  init_symop( mc );
  init_spacegroup( mc, mdata );
  init_clipper_stats( mc );
  init_resol_fn( mc );
  add_clipper_tests( mc, mdata );
  
  add_rotation( mc );
  add_xmap( mc );
  add_nxmap( mc );
  add_nxmap_operator( mc );
  add_mapfilters( mc );
  add_map_utils( mc );
  
  init_minimol(mc, mm);
  init_minimol_seq(mc);
  add_matomindex(mc);
  init_sfscale(mc);
  
  add_fffear(mc);
  add_edcalc(mc);



  add_hklinfo(mc);
  //init_containers(mc);
  init_hkl_datatypes(mc);
  init_hkl_data(mc);

  //init_gemmi_structure(mc);
  //init_map_io(mc);
  
  //init_clipper_util(mc);
  // init_ccp4_mtz_io(mc);

  add_simplex_lib(mb);
  add_map_simulate(mb);
  add_ca_build(mb);
  add_ca_join(mb);
  add_ca_filter(mb);
  add_buccaneer_lib(mb);
  add_buccaneer_util(mb);
  add_ca_prep(mb);
  add_buccaneer_prot(mb);
  add_protein_loop(mb);
  add_ca_merge(mb);
  add_ca_find(mb);
  add_ca_grow(mb);
  add_ca_link(mb);
  add_ca_sequence(mb);
  add_ca_correct(mb);
  add_ca_ncsbuild(mb);
  add_ca_prune(mb);
  add_knownstructure(mb);
  add_model_tidy(mb);

  add_proteindb(mpdb);

  add_utils( mutil );
}