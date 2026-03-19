// Nanobind bindings for clipper atomsf
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

void add_atomsf( nb::module_ &m ) {
  nb::enum_<ScatteringFactorsType>( m, "ScatteringFactorsType" )
      .value( "SF_WAASMAIER_KIRFEL", ScatteringFactorsType::SF_WAASMAIER_KIRFEL )
      .value( "SF_ELECTRON", ScatteringFactorsType::SF_ELECTRON );

  nb::class_<ScatteringFactors>( m, "ScatteringFactors" )
      .def_static( "select_SFType", &ScatteringFactors::selectScattteringFactorsType );

  nb::class_<AtomShapeFn> atomsf( m, "AtomShapeFn" );
  nb::enum_<AtomShapeFn::TYPE>( atomsf, "TYPE", "Atomic parameters." )
      .value( "X", AtomShapeFn::TYPE::X )
      .value( "Y", AtomShapeFn::TYPE::Y )
      .value( "Z", AtomShapeFn::TYPE::Z )
      .value( "Uiso", AtomShapeFn::TYPE::Uiso )
      .value( "Occ", AtomShapeFn::TYPE::Occ )
      .value( "U11", AtomShapeFn::TYPE::U11 )
      .value( "U22", AtomShapeFn::TYPE::U22 )
      .value( "U33", AtomShapeFn::TYPE::U33 )
      .value( "U12", AtomShapeFn::TYPE::U12 )
      .value( "U13", AtomShapeFn::TYPE::U13 )
      .value( "U23", AtomShapeFn::TYPE::U23 )
      .export_values();

  atomsf.def( nb::init<>(), "Null constructor." )
      .def( nb::init<const Atom &>(), nb::arg( "atom" ), "Constructor: from atom object." )
      .def( nb::init<const Coord_orth &, const String &, const ftype, const ftype>(), nb::arg( "pos" ),
            nb::arg( "element" ), nb::arg( "u_iso" ) = 0.0, nb::arg( "occ" ) = 1.0,
            "Constructor: from coord, element, isotropic U, occupancy." )
      .def( nb::init<const Coord_orth &, const String &, const U_aniso_orth &, const ftype>(), nb::arg( "pos" ),
            nb::arg( "element" ), nb::arg( "u_aniso" ), nb::arg( "occ" ) = 1.0,
            "Constructor: from coord, element, anisotropic U, occupancy." )
      .def( "init", ( void ( AtomShapeFn::* )( const Atom & ) )&AtomShapeFn::init, nb::arg( "atom" ),
            "Initialiser: from atom object." )
      .def(
          "init",
          ( void ( AtomShapeFn::* )( const Coord_orth &, const String &, const ftype, const ftype ) )&AtomShapeFn::init,
          nb::arg( "pos" ), nb::arg( "element" ), nb::arg( "u_iso" ) = 0.0, nb::arg( "occ" ) = 1.0,
          "Initialiser: from coord, element, isotropic U, occupancy." )
      .def( "init",
            ( void ( AtomShapeFn::* )( const Coord_orth &, const String &, const U_aniso_orth &,
                                       const ftype ) )&AtomShapeFn::init,
            nb::arg( "pos" ), nb::arg( "element" ), nb::arg( "u_aniso" ), nb::arg( "occ" ) = 1.0,
            "Initialiser: from coord, element, anisotropic U, occupancy." )
      .def( "f", ( ftype ( AtomShapeFn::* )( const Coord_reci_orth & ) const ) & AtomShapeFn::f,
            nb::arg( "coord_reci_orth" ), "Return scattering factor as a function of reflection position." )
      .def( "rho", ( ftype ( AtomShapeFn::* )( const Coord_orth & ) const ) & AtomShapeFn::rho, nb::arg( "pos" ),
            "Return electron density as a function of coordinate." )
      .def( "f_iso", ( ftype ( AtomShapeFn::* )( const ftype & ) const ) & AtomShapeFn::f, nb::arg( "invresolsq" ),
            "Return isotropic scattering factor as a function of resolution." )
      .def( "rho_iso", ( ftype ( AtomShapeFn::* )( const ftype & ) const ) & AtomShapeFn::rho, nb::arg( "rad" ),
            "Return isotropic electron density as a function of radius." )
      // need to use lambda and return tuple of results
      .def(
          "rho_grad",
          []( const AtomShapeFn &self, const Coord_orth &pos ) {
            ftype rho;
            std::vector<ftype> grad;
            self.rho_grad( pos, rho, grad );
            return nb::make_tuple( rho, grad );
          },
          nb::arg( "pos" ), "Return a tuple of Agarwal density and gradients as a function of coordinate." )
      .def(
          "rho_curv",
          []( const AtomShapeFn &self, const Coord_orth &pos, Matrix<ftype> &curv ) {
            ftype rho;
            std::vector<ftype> grad;
            self.rho_curv( pos, rho, grad, curv );
            return nb::make_tuple( rho, grad );
          },
          nb::arg( "pos" ), nb::arg( "curv" ),
          "Updates given curvature Matrix and returns a tuple of Agarwal density and gradient as a function of "
          "coordinate." )
      .def_prop_rw(
          "agarwal_params",
          []( AtomShapeFn &self ) -> std::vector<AtomShapeFn::TYPE> & { return self.agarwal_params(); },
          []( AtomShapeFn &self, const std::vector<AtomShapeFn::TYPE> &params ) { self.agarwal_params() = params; },
          nb::for_getter( "Get parameters for Agarwal gradient/curvature calculations." ),
          nb::for_setter( "Set parameters for Agarwal gradient/curvature calculations." ) )
      .doc() = "Atomic shape function object\n"
               "The atomic scattering factor object is instantiated for each "
               "atom in turn, giving the atom parameters: position, element, "
               "occupancy and the isotropic or anisotropic U-value. (See "
               "clipper::Util for conversion from B-factors.). The methods of the "
               "class may then be called to return the scattering in reciprocal "
               "space or density in real space using either isotropic or "
               "anistropic models as required.\n"
               "If the atom only has an isotropic U, the faster isotropic methods"
               "will be used where available.\n"
               "This implementation uses the coefficients from Waasmaier & Kirfel"
               "(1995), Acta Cryst. A51, 416-431. The source data can be found at:"
               "ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat"
               "\nTo use electron scattering factors, call :\n"
               "ScatteringFactors.select_SFType(ScatteringFactorsType.SF_ELECTRON)\n"
               "then the next instantiation of AtomShapeFn will use electron scattering factors.";

  // atomsf
  //   .def(nb::init<>(), "Null constructor")
  //   .def(nb::init<const String &, const ftype &, const ftype &, const ftype &>(),
  //        nb::arg("label"), nb::arg("f0"), nb::arg("f_prime"),
  //        nb::arg("f_double_prime"),
  //        "Constructor: takes label, f0, f_prime, f_double_prime.")
  //   .def("__init__", [](Atom_scat* atom_scat, const std::string &label,
  //                      const double &f0, const double &f_prime,
  //                      const double &f_double_prime) {
  //         new(atom_scat) Atom_scat(label, f0, f_prime, f_double_prime);
  //       },
  //       nb::arg("label"), nb::arg("f0"), nb::arg("f_prime"),
  //       nb::arg("f_double_prime"),
  //       "Constructor: takes label, f0, f_prime, f_double_prime.")
  //   .def("label", &Atom_scat::label, "Get the label")
  //   .def("f0", &Atom_scat::f0, "Get the f0")
  //   .def("f_prime", &Atom_scat::f_prime, "Get the f_prime")
  //   .def("f_double_prime", &Atom_scat::f_double_prime, "Get the f_double_prime")
  //   .def("__repr__",
  //        [](const Atom_scat &self) { return "<clipper.Atom_scat class.>"; });
}