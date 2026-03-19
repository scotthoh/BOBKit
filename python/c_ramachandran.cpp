// Nanobind bindings for clipper Ramachandran
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "arrays.h"
#include "commons.h"
#include <clipper/core/ramachandran.h>
// #include <nanobind/stl/vector.h>
#include <nanobind/stl/array.h>

using namespace clipper;

void add_ramachandran( nb::module_ &m ) {
  nb::class_<Prob_phi_2d>( m, "Prob_phi_2d" )
      .def( nb::init<>() )
      .def( "init", &Prob_phi_2d::init, nb::arg( "size" ), "Initialise: with sampling" )
      .def(
          "accumulate", []( Prob_phi_2d &self, const cpu_c_array<float> &table ) { self.accumulate( table.data() ); },
          nb::arg( "table" ), "Accumulate new table (C-contiguous numpy array) of samples to probability" )
      .def( "accumulate", ( void ( Prob_phi_2d::* )( const double &, const double &, double ) )&Prob_phi_2d::accumulate,
            nb::arg( "phi1" ), nb::arg( "phi2" ), nb::arg( "wgt" ) = 1.0, "Accumulate new sample to probability." )
      .def( "normalise", &Prob_phi_2d::normalise, "Normalise to integrate to 1/(2pi)^2" )
      .def( "probability", &Prob_phi_2d::probability, nb::arg( "phi1" ), nb::arg( "phi2" ),
            "Get probability for a particular pair of angles." )
      .def( "format", &Prob_phi_2d::format,
            "Return formatted string representation (as C++ code), values rounded to integer." )
      .def( "get_data", ( const double &( Prob_phi_2d::* )( const int &, const int & ) const ) & Prob_phi_2d::data,
            nb::arg( "i" ), nb::arg( "j" ), "2d read access." )
      .def( "set_data", ( double &( Prob_phi_2d::* )( const int &, const int & )) & Prob_phi_2d::data, nb::arg( "i" ),
            nb::arg( "j" ), "2d write access" )
      .def( "__repr__", []( const Prob_phi_2d &self ) { return "<clipper.Prob_phi_2d class>"; } )
      .doc() = "2-d angular probability distibution class.\n"
               "Base for Ramachandran class (and other similar classes, such as "
               "a pseudo-ramachandran plot or the JPD of two phases).";

  nb::class_<Ramachandran> ramachandran( m, "Ramachandran" );

  nb::enum_<Ramachandran::TYPE>( ramachandran, "TYPE", "Enumeration of built-in Ramachandran tables." )
      .value( "Gly", Ramachandran::TYPE::Gly )
      .value( "Pro", Ramachandran::TYPE::Pro )
      .value( "NonGlyPro", Ramachandran::TYPE::NonGlyPro )
      .value( "NonGly", Ramachandran::TYPE::NonGly )
      .value( "All", Ramachandran::TYPE::All )
      .value( "Gly5", Ramachandran::TYPE::Gly5 )
      .value( "Pro5", Ramachandran::TYPE::Pro5 )
      .value( "NonGlyPro5", Ramachandran::TYPE::NonGlyPro5 )
      .value( "NonGly5", Ramachandran::TYPE::NonGly5 )
      .value( "All5", Ramachandran::TYPE::All5 )
      .value( "All2", Ramachandran::TYPE::All2 )
      .value( "Gly2", Ramachandran::TYPE::Gly2 )
      .value( "Pro2", Ramachandran::TYPE::Pro2 )
      .value( "PrePro2", Ramachandran::TYPE::PrePro2 )
      .value( "IleVal2", Ramachandran::TYPE::IleVal2 )
      .value( "NoGPIVpreP2", Ramachandran::TYPE::NoGPIVpreP2 )
      .export_values();

  ramachandran.def( nb::init<>(), "Null constructor." )
      .def( nb::init<Ramachandran::TYPE>(), nb::arg( "type" ), "Constructor from standard plot." )
      .def( "init", &Ramachandran::init, nb::arg( "type" ), "Initialise from standard plot." )
      .def( "set_thresholds", &Ramachandran::set_thresholds, nb::arg( "prob_favoured" ) = 0.01,
            nb::arg( "prob_allowed" ) = 0.0005, "Change thresholds to difference values." )
      .def( "probability", &Ramachandran::probability, nb::arg( "phi" ), nb::arg( "psi" ),
            "Get probability for a particular pair of angles." )
      .def( "favoured", &Ramachandran::favored, nb::arg( "phi" ), nb::arg( "psi" ),
            "Test if a pair of angles are in the favoured region." )
      .def( "favored", &Ramachandran::favored, nb::arg( "phi" ), nb::arg( "psi" ),
            "Test if a pair of angles are in the favoured region." )
      .def( "allowed", &Ramachandran::allowed, nb::arg( "phi" ), nb::arg( "psi" ),
            "Test if a pair of angles are in an allowed region." )
      .def( "__repr__", []( const Ramachandran &self ) { return "<clipper.Ramachandran plot class.>"; } )
      .doc() = "Ramachandran plot class.\nThis class provides a reference Ramachandran plot for "
               "Gly, Pro, other, and combinations of those types of residues. The source data comes from the "
               "best residues from the 'top500' best-determined structures list of D. C. and J. S. Richardson, "
               "http://kinemage.biochem.duke.edu/index.html \nThe Ramachandran plot is normalised in inverse "
               "radians squared, so the mean value of a probability is 1/(2 pi)<sup>2</sup>.";
}
