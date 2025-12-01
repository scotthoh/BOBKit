// Nanobind bindings for buccaneer-lib
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-link.h"
#include "commons.h"
#include <nanobind/operators.h>

using namespace clipper;

void add_ca_link( nb::module_ &m ) {
  nb::class_<Ca_link>( m, "Ca_link" )
      .def( nb::init<ftype, int>(), nb::arg( "rad_link" ) = 5.0, nb::arg( "torsion_sampling" ) = 24,
            "Constructor for Ca_link class." )
      .def( "__call__", &Ca_link::operator(), nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
            "Merge overlapped Ca chains." )
      .def_prop_ro( "num_linked", &Ca_link::num_linked, "Get number of Ca alphas linked." )
      .def_static( "link",
        []( clipper::MiniMol &mol, const clipper::Xmap<float> &xmap, const LLK_map_target &llktarget,
            clipper::ftype radlink, int torsamp, const nb::object &pystream ) -> bool {
          Ca_link calink( radlink, torsamp );
          bool success = calink( mol, xmap, llktarget );
          std::string m = " C-alphas linked           : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
          if ( pystream.is_valid() ) to_pystream(m, pystream);
          else std::cout << m;
          return success;
      }, nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktarget" ), nb::arg("radlink")=5.0, nb::arg("torsion_sampling")=24, nb::arg("stdout")=nullptr )
      .def( "__repr__",
            []( const Ca_link &self ) {
              return "<buccaneer.Ca_link with " + String( self.num_linked() ) + " C-alphas linked.>";
            } )
      .doc() = "Class for merging overlapped Ca chains and grouping by symmetry.";
}