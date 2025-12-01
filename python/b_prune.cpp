// Nanobind bindings for buccaneer-prune
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prune.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>
#include <iostream>


void add_ca_prune( nb::module_ &m ) {
  nb::class_<Ca_prune>( m, "Ca_prune", "Class for pruning clashing Ca chains using density." )
      .def( nb::init<double>(), nb::arg( "rad" ) = 3.0, "Constructor with radius." )
      //.def_static("prune", &Ca_prune::prune, nb::arg("mol"), nb::arg("xmap"),
      //            nb::arg("rad") = 3.0,
      //            "Static function to prune clashing Ca chains using density "
      //            "and specified radius.")
      .def_static(
          "prune",
          []( clipper::MiniMol &mol, const clipper::Xmap<float> &xmap, double &rad, const nb::object &pystream ) -> bool {
            bool success = Ca_prune::prune( mol, xmap, rad );
            std::string m = " C-alphas after pruning    : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
            if ( pystream.is_valid() ) to_pystream(m, pystream);
            else std::cout << m;
            return success;
          },
          nb::arg( "mol" ), nb::arg( "xmap" ), nb::arg( "rad" ) = 3.0, nb::arg("stdout")=nullptr,
          "Static function to prune clashing Ca chains using density and specified radius, with option to print summary." )
      .def( "__call__", &Ca_prune::operator(), nb::arg( "mol" ), nb::arg( "xmap" ),
            "Prune clashing Ca chains using density." )
      .def( "__repr__", []( const Ca_prune &self ) { return "<buccaneer.Ca_prune class.>"; } );
}