// Nanobind bindings for buccaneer-correct
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-correct.h"
#include "commons.h"
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

void add_ca_correct( nb::module_ &m ) {
  nb::class_<Ca_correct>( m, "Ca_correct" )
      .def( nb::init<int>(), nb::arg( "torsion_sampling" ) = 12, "Constructor for Ca correct class." )
      .def( "__call__", &Ca_correct::operator(), nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktargets" ),
            nb::arg( "seq" ), "Rebuild chain to fix insertion/deletions." )
      .def_prop_ro( "num_corrected", &Ca_correct::num_corrected, "Get number of C-alphas corrected." )
      .def_static( "correct", []( clipper::MiniMol &mol, const clipper::Xmap<float> &xmap, const std::vector<LLK_map_target> &llkcls,
        const clipper::MMoleculeSequence &seq, int &torsamp, const nb::object &pystream ) -> bool {
        Ca_correct cacor(torsamp);
        bool success = cacor( mol, xmap, llkcls, seq );
        std::string m = " C-alphas corrected        : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
        if ( pystream.is_valid() ) to_pystream(m, pystream);  
        else std::cout << m;
        return success;
      }, nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "llktargets" ), nb::arg( "seq" ), nb::arg("torsion_sampling")=12,
      nb::arg("stdout")=nullptr, "Rebuild chain to fix insertion/deletions. Static function with option to print summary.")
      .def( "__repr__",
            []( const Ca_correct &self ) {
              return "<buccaneer.Ca_correct with " + String( self.num_corrected() ) + " C-alphas corrected.>";
            } )
      .doc() = "Class for correcting Ca chains using density.";
}