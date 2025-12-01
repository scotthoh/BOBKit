// Nanobind bindings for buccaneer-build
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-build.h"
#include "commons.h"
#include <nanobind/operators.h>

void add_ca_build( nb::module_ &m ) {
  nb::class_<Ca_build>( m, "Ca_build" )
      .def( nb::init<clipper::String, bool>(), nb::arg( "newrestype" ) = "ALA", nb::arg( "flexible" ) = false,
            "Constructor for Ca_build." )
      //.def_static( "build", &Ca_build::build, nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "newrestype" ) = "ALA",
      //    nb::arg( "flexible" ) = false, "Build Ca chains using density." )
      .def_static(
          "build",
          []( clipper::MiniMol &mol, const clipper::Xmap<float> &xmap, clipper::String newrestype,
              bool flexible, const nb::object &pystream ) -> bool {
            bool success = Ca_build::build( mol, xmap, newrestype, flexible );
            std::string m = " C-alphas after rebuilding : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
            if ( pystream.is_valid() ) to_pystream(m, pystream);
            else std::cout << m;
            return success;
          },
          /*&Ca_build::build,*/ nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "newrestype" ) = "ALA",
          nb::arg( "flexible" ) = false, nb::arg( "stdout" ) = nullptr,
          "Build Ca chains using density, static function with an option to print summary." )
      .def( "__call__", &Ca_build::operator(), nb::arg( "mol" ), nb::arg( "xmap" ), "Build Ca chains using density." )
      .def( "__repr__", []( const Ca_build &self ) { return "<buccaneer.Ca_build class>"; } )
      .doc() = "Class for building Ca chains using density.";
}