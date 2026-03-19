// Nanobind bindings for clipper tests_core, test_contrib, test_minimol
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper.h>
#include <clipper/contrib/test_contrib.h>
#include <clipper/core/test_core.h>
#include <clipper/core/test_data.h>
// #include <clipper/minimol/test_minimol_gemmi.h>

using namespace clipper;

// Self test for clipper_core and clipper_contrib,
// useful to check if fftw3 compiled correctly
void add_clipper_tests( nb::module_ &m, nb::module_ &mdata ) {
  nb::class_<Test_core>( m, "Test_core", "Class test clipper core methods." )
      .def( nb::init<>() )
      .def( "__call__", &Test_core::operator() )
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_core::set_stream)
      .def( "__repr__", []( const Test_core &self ) { return "<clipper.Test_core class.>"; } );

  nb::class_<Test_contrib>( m, "Test_contrib", "Class test clipper contrib methods." )
      .def( nb::init<>() )
      .def( "__call__", &Test_contrib::operator() )
      // if needed, have to find a way for the type conversion std::ostream
      //.def("set_stream", &Test_contrib::set_stream)
      .def( "__repr__", []( const Test_contrib &self ) { return "<clipper.Test_contrib class.>"; } );

  nb::class_<data::Test_data>( mdata, "Test_data", "Class to return test data." )
      .def( nb::init<>(), "Null constructor, fills the arrays" )
      .def_prop_ro( "hkl_data_f_sigf", &data::Test_data::hkl_data_f_sigf, "Return HKL_data class for F and SIGF." )
      .def_prop_ro( "hkl_data_abcd", &data::Test_data::hkl_data_abcd, "Return HKL_data class for ABCD." )
      .def_prop_ro( "atom_list", &data::Test_data::atom_list, "Return atom list." )
      .def( "__repr__", []( const data::Test_data &self ) { return "<clipper.Test_data class.>"; } );
  // nb::class_<Test_minimol_gemmi>(m, "Test_minimol_gemmi", "Class test clipper minimol-gemmi methods.")
  //     .def(nb::init<>())
  //     .def("run", &Test_minimol_gemmi::run, nb::arg("input_file"))
  //     .def("__repr__", [](const Test_minimol_gemmi &self) {
  //       return "<clipper.Test_minimol_gemmi class.>";
  //     });
}