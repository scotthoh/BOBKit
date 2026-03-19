// Nanobind bindings for buccaneer-ncsbuild
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-ncsbuild.h"
#include "commons.h"
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

void add_ca_ncsbuild( nb::module_ &m ) {
  nb::class_<Ca_ncsbuild>( m, "Ca_ncsbuild" )
      .def( nb::init<double, double, int>(), nb::arg( "reliability" ) = 0.5, nb::arg( "rmsd" ) = 1.0,
            nb::arg( "nmin" ) = 12, "Constructor for Ca_ncsbuild class." )
      .def( "__call__", &Ca_ncsbuild::operator(), nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
            nb::arg( "seq" ), "Build NCS related chains." )
      .def_static(
          "ncsbuild",
          []( clipper::MiniMol &mol, const clipper::Xmap<float> &xmap, const std::vector<LLK_map_target> &llkcls,
              const clipper::MMoleculeSequence &seq, double &reliability, double &rmsd, int &nmin, const nb::object &pystream ) -> bool {
            Ca_ncsbuild cancsbuild( reliability, rmsd, nmin );
            std::string m = " C-alphas after NCS build  : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
            bool success = cancsbuild( mol, xmap, llkcls, seq );
            if ( pystream.is_valid() ) to_pystream(m, pystream);
            std::cout << m;
            return success;
          },
          nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktarget" ), nb::arg( "seq" ), nb::arg( "reliability" ) = 0.5,
          nb::arg( "rmsd" ) = 1.0, nb::arg( "nmin" ) = 12, nb::arg( "stdout" ) = nullptr,
          "Static function to build NCS related chains, with option to print summary." )
      .def( "__repr__", []( const Ca_ncsbuild &self ) { return "<buccaneer.Ca_ncsbuild class.>"; } )
      .doc() = "Class for build NCS related chains using density.";
}