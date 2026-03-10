// Nanobind bindings for clipper coord.h
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "arrays.h"
#include "commons.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/core/coords.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
//#include <gemmi/math.hpp>

using namespace clipper;

void add_coords( nb::module_ &m ) {
  //auto tr = nb::module_::import_("gemmi").attr("Transform");

  nb::class_<Resolution>( m, "Resolution" )
      .def( nb::init<>() )
      .def( nb::init<const double &>(), nb::arg( "resolution" ), "Constructor." )
      .def( "init", &Resolution::init, nb::arg( "resolution" ), "Initialiser." )
      .def( "limit", &Resolution::limit, "Get resolution limit." )
      .def( "invresolsq_limit", &Resolution::invresolsq_limit, "Get inverse resolution square limit." )
      .def( "is_null", &Resolution::is_null, "Test if value has been initialised." )
      .def( "__str__", []( const Resolution &self ) { return String( self.limit(), 6, 4 ); } )
      .def( "__repr__",
            []( const Resolution &self ) { return "<clipper.Resolution " + String( self.limit(), 6, 4 ) + " Å.>"; } )
      .doc() = "Resolution in angstroms\n"
               "This object represents a resolution limit which will be used for "
               "all aspects of a calculation. This is a base for a donor type.";

  nb::class_<HKL_class>( m, "HKL_class" )
      .def( nb::init<>() )
      .def( nb::init<const Spacegroup &, const HKL &>(), nb::arg( "spacegroup" ), nb::arg( "hkl" ),
            "Constructor from spacegroup and HKL" )
      .def(
          "__init__",
          []( HKL_class *p, const Spacegroup &sg, const std::array<int, 3> &hkl ) {
            new ( p ) HKL_class( sg, HKL( hkl[0], hkl[1], hkl[2] ) );
          },
          "Constructor from spacegroup and HKL(list or numpy array)" )
      .def_prop_ro( "epsilon", &HKL_class::epsilon, "Get epsilon." )
      .def( "epsilonc", &HKL_class::epsilonc, "Get epsilon for acentric, 2x epsilon for centric." )
      .def( "allowed", &HKL_class::allowed, "Get allowed phase." )
      .def( "is_centric", &HKL_class::centric, "Is it centric?" )
      .def( "is_sys_abs", &HKL_class::sys_abs, "Is there systematic absences?" )
      .def( "__repr__",
            []( const HKL_class &self ) {
              return "<clipper.HKL_class: Describes the type of reflection in a "
                     "given spacegroup.>";
            } )
      .doc() = "Reflection class\n"
               "This describes the type of a reflection in a given spacegroup, "
               "including centricity, systematic absence, phase restriction, and "
               "multiplicity.";

  nb::class_<RTop_orth, RTop<>> rtop_orth( m, "RTop_orth" );
  rtop_orth.def( nb::init<>() )
      .def( nb::init<const RTop<> &>(), nb::arg( "rtop" ), "Constructor: copy/convert." )
      .def( nb::init<const Mat33<> &>(), nb::arg( "rot" ), "Constructor from rotation." )
      .def( nb::init<const Mat33<> &, const Vec3<> &>(), nb::arg( "rot" ), nb::arg( "trn" ),
            "Constructor from rotation and translation" )
      .def(
          "__init__",
          []( RTop_orth *rt, const cpu_c_2darray<double, 3, 3> &m, cpu_c_array<double> &v ) {
            clipper::Mat33<> rotation( m( 0, 0 ), m( 0, 1 ), m( 0, 2 ), m( 1, 0 ), m( 1, 1 ), m( 1, 2 ), m( 2, 0 ),
                                       m( 2, 1 ), m( 2, 2 ) );
            clipper::Vec3<> translation( v( 0 ), v( 1 ), v( 2 ) );
            new ( rt ) RTop_orth( rotation, translation );
          },
          "Constructor from rotation and translation Ndarrays" )
      .def( nb::init<const std::vector<Coord_orth> &, const std::vector<Coord_orth> &>(), nb::arg( "src" ),
            nb::arg( "tgt" ), "Constructor from two lists of Coord_orth." )
      .def( nb::init<const std::vector<Coord_orth> &, const std::vector<Coord_orth> &, const std::vector<double> &>(),
            nb::arg( "src" ), nb::arg( "tgt" ), nb::arg( "weight" ),
            "Constructor from two lists of Coord_orth and weight." )
      .def( nb::init<const Atom_list &, const Atom_list &>(), nb::arg( "src" ), nb::arg( "tgt" ),
            "Constructor from two Atom_list type objects." )
      .def(
          "__init__",
          []( RTop_orth *rt, std::array<std::array<double, 3>, 3> &rot, std::array<double, 3> &trn ) {
            auto rotation = array_to_mat33( rot );
            auto translation = array_to_vec3( trn );
            new ( rt ) RTop_orth( *rotation, *translation );
          },
          nb::arg( "rot" ), nb::arg( "trn" ), "Constructor from rotation and translation(list/arrays)." )
      .def_static(
          "from_gemmi_transform", []( const gemmi::Transform &rtop ) { return GEMMI::transform( rtop ); },
          nb::arg( "rtop" ), //nb::sig( "def from_gemmi_transform(self, rtop: gemmi.Transform, /) -> bobkit.clipper.RTop_orth" ),
          "Convert gemmi::Transform to RTop_orth." )
      .def_static(
          "to_gemmi_transform", []( const RTop_orth &rtop ) { return GEMMI::transform( rtop ); }, nb::arg( "rtop" ),
          "Convert RTop_orth to gemmi::Transform." )
      .def( "rtop_frac", &RTop_orth::rtop_frac, nb::arg( "cell" ), "Orthogonal-fractional conversion." )
      .def( "inverse", &RTop_orth::inverse, "Inverse operator." )
      .def( "axis_coordinate_near", &RTop_orth::axis_coordinate_near, nb::arg( "centre" ),
            "Return point on axis near the specified coordinate." )
      .def( "screw_translation", &RTop_orth::screw_translation, "Return screw translation." )
      .def_static( "identity", &RTop_orth::identity, "Return identity operator." )
      .def_static( "null", &RTop_orth::null, "Return null (uninitialised) operator." )
      .def( "__repr__", []( const RTop_orth &self ) { return "<clipper.RTop_orth class.>"; } )
      .doc() = "Orthogonal operator class.\n"
               "This class is used for any RT-operator which operates on orthogonal "
               "coordinates. For a full list of methods, see clipper::RTop";

  nb::class_<HKL, Vec3<int>> hkl( m, "HKL", "Reflection 'Miller index." );
  hkl.def( nb::init<>() )
      .def( nb::init<const Vec3<int> &>(), "Constructor copy/convert." )
      .def( nb::init<const int &, const int &, const int &>(), nb::arg( "h" ), nb::arg( "k" ), nb::arg( "l" ),
            "Constructor from H,K,L." )
      .def(
          "__init__",
          []( HKL *p, const std::array<int, 3> &hkl ) {
            new ( p ) HKL( hkl[0], hkl[1], hkl[2] );
            // std::unique_ptr<HKL>(
            //     new HKL(hkl[0], hkl[1], hkl[2]));
          },
          nb::arg( "hkl" ), "Constructor from list/array." )
      .def( "__init__", ( []( HKL *p, const gemmi::Miller &hkl ) { new ( p ) HKL( hkl[0], hkl[1], hkl[2] ); } ),
            "Constructor from gemmi::Miller." )
      .def_static(
          "from_gemmi_Miller", []( const gemmi::Miller &hkl ) { return GEMMI::Hkl( hkl ); }, nb::arg( "hkl" ),
          "Convert gemmi::Miller to HKL." )
      .def_static(
          "to_gemmi_Miller", []( const HKL &hkl ) { return GEMMI::Hkl( hkl ); }, nb::arg( "hkl" ),
          "Convert HKL to gemmi::Miller." )
      .def_prop_rw(
          "h", []( const HKL &self ) -> const int & { return self.h(); },
          []( HKL &self, const int &val ) { self.h() = val; }, "Get/set h." )
      .def_prop_rw(
          "k", []( const HKL &self ) -> const int & { return self.k(); },
          []( HKL &self, const int &val ) { self.k() = val; }, "Get/set k." )
      .def_prop_rw(
          "l", []( const HKL &self ) -> const int & { return self.l(); },
          []( HKL &self, const int &val ) { self.l() = val; }, "Get/set l." )
      .def_prop_rw(
          "hkl", []( const HKL &self ) -> std::array<int, 3> { return { self.h(), self.k(), self.l() }; },
          []( HKL &self, std::array<int, 3> hkl ) {
            self.h() = hkl[0];
            self.k() = hkl[1];
            self.l() = hkl[2];
          },
          "Get/set hkl." )
      .def( "invresolsq", &HKL::invresolsq, nb::arg( "cell" ),
            "Return inverse resolution squared for this reflection in given "
            "cell." )
      .def( "coord_reci_frac", &HKL::coord_reci_frac, "Return fractional reciprocal coordinate (i.e. non-integer HKL)" )
      .def( "coord_reci_orth", &HKL::coord_reci_orth, nb::arg( "cell" ),
            "Orthogonal-fractional reciprocal space coordinate conversion" )
      .def( "transform", ( HKL ( HKL::* )( const Symop & ) const ) & HKL::transform, nb::arg( "symop" ),
            "Return transformed hkl." )
      .def( "transform", ( HKL ( HKL::* )( const Isymop & ) const ) & HKL::transform, nb::arg( "isymop" ),
            "Return transformed hkl." )
      .def( "sym_phase_shift", &HKL::sym_phase_shift, nb::arg( "op" ),
            "Return symmerty phase shift for this HKL under op." )
      .def( "format", &HKL::format, "Return formatted String representation." )
      .def( "__repr__", []( const HKL &self ) { return "<clipper." + self.format() + ">"; } )
      .def( -nb::self )
      .def( nb::self - nb::self )
      .def( nb::self + nb::self )
      .def( nb::self += nb::self, nb::rv_policy::none )
      .def( operator-=( nb::self, nb::self ), nb::rv_policy::none )
      .def(
          "__mul__", []( const HKL &self, const int &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const HKL &self, const int &s ) { return s * self; }, nb::is_operator() )
      .def( "__rmul__", []( const HKL &self, const Isymop &op ) { return op * self; } )
      .def_prop_ro( "array", []( HKL &self) {
        return nb::ndarray<int, nb::numpy>(&self[0], {3}, nb::handle());
      });

  nb::class_<Coord_reci_orth, Vec3<>> coord_reci_orth(
      m, "Coord_reci_orth", "Orthogonal reciprocal coordinate (length of which is invresolsq)" );
  coord_reci_orth.def( nb::init<>() )
      .def( nb::init<const Vec3<> &>(), nb::arg( "v" ), "Constructor copy/convert." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "xs" ), nb::arg( "ys" ),
            nb::arg( "zs" ), "Constructor from x*, y*, z*." )
      .def(
          "__init__",
          []( Coord_reci_orth *cro, const std::array<double, 3> &a ) {
            new ( cro ) Coord_reci_orth( a[0], a[1], a[2] );
          },
          "Constructor from list/array of x*,y*,z*." )
      .def_prop_ro( "xs", &Coord_reci_orth::xs, "Get x*" )
      .def_prop_ro( "ys", &Coord_reci_orth::ys, "Get y*" )
      .def_prop_ro( "zs", &Coord_reci_orth::zs, "Get z*" )
      .def( "invresolsq", &Coord_reci_orth::invresolsq, "Return inverse resolution squared for this coordinate." )
      .def( "coord_reci_frac", &Coord_reci_orth::coord_reci_frac, nb::arg( "cell" ),
            "Orthogonal-fractional reciprocal space coordinate conversion" )
      .def( "transform", &Coord_reci_orth::transform, nb::arg( "op" ), "Return transformed coordinate." )
      .def( "format", &Coord_reci_orth::format, "Return formatted string representation." )
      .def( "__str__", &Coord_reci_orth::format );

  nb::class_<Coord_reci_frac, Vec3<>> coord_reci_frac( m, "Coord_reci_frac",
                                                       "Fractional reciprocal coordinate (i.e. non-integer hkl)" );
  coord_reci_frac.def( nb::init<>() )
      .def( nb::init<const Vec3<> &>(), nb::arg( "v" ), "Constructor copy/convert." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "us" ), nb::arg( "vs" ),
            nb::arg( "ws" ), "Constructor from u*,v*,w*." )
      .def(
          "__init__",
          []( Coord_reci_frac *crf, const std::array<double, 3> &a ) {
            new ( crf ) Coord_reci_frac( a[0], a[1], a[2] );
          },
          "Constructor from list/array of u*,v*,w*." )
      .def( nb::init<const HKL &>(), nb::arg( "hkl" ), "Constructor from HKL." )
      .def( "hkl", &Coord_reci_frac::hkl, "Round to HKL." )
      .def( "invresolsq", &Coord_reci_frac::invresolsq, nb::arg( "cell" ),
            "Return inverse resolution squared for this reflection in given "
            "cell." )
      .def_prop_ro( "us", &Coord_reci_frac::us, "Get u*" )
      .def_prop_ro( "vs", &Coord_reci_frac::vs, "Get v*" )
      .def_prop_ro( "ws", &Coord_reci_frac::ws, "Get w*" )
      .def( "coord_reci_orth", &Coord_reci_frac::coord_reci_orth, nb::arg( "cell" ),
            "Fractional-orthogonal reciprocal space coordinate conversion" )
      .def( "transform", &Coord_reci_frac::transform, nb::arg( "op" ), "Returned transformed coordinate." )
      .def( "__str__", &Coord_reci_frac::format )
      .def( "format", &Coord_reci_frac::format, "Return formatted string representation." );

  nb::class_<Coord_grid, Vec3<int>> coord_grid( m, "Coord_grid", "Grid coordinate." );
  coord_grid.def( nb::init<>() )
      .def( nb::init<const Vec3<int> &>(), "Constructor copy/convert." )
      .def( nb::init<const int &, const int &, const int &>(), nb::arg( "u" ), nb::arg( "v" ), nb::arg( "w" ),
            "Constructor from u,v,w" )
      .def( nb::init<const Grid &, const int &>(), nb::arg( "grid" ), nb::arg( "index" ),
            "Constructor from a grid and an index in that grid." )
      .def_prop_rw(
          "u", []( const Coord_grid &self ) -> const int & { return self.u(); },
          []( Coord_grid &self, const int &val ) { self.u() = val; }, "Get/set u." )
      .def_prop_rw(
          "v", []( const Coord_grid &self ) -> const int & { return self.v(); },
          []( Coord_grid &self, const int &val ) { self.v() = val; }, "Get/set v." )
      .def_prop_rw(
          "w", []( const Coord_grid &self ) -> const int & { return self.w(); },
          []( Coord_grid &self, const int &val ) { self.w() = val; }, "Get/set w." )
      .def( "coord_map", &Coord_grid::coord_map, "Convert to Coord_map." )
      .def( "coord_frac", &Coord_grid::coord_frac, nb::arg( "grid_sampling" ),
            "Convert to Coord_fract using given Grid_sampling." )
      .def( "transform", &Coord_grid::transform, nb::arg( "op" ), "Returned transformed coordinate." )
      .def( "unit", &Coord_grid::unit, nb::arg( "grid_sampling" ), "Reduced to unit box: (0..nu-1, 0..nv-1, 0..nw-1)" )
      .def( "next", ( const Coord_grid &( Coord_grid::* )( const Grid & )) & Coord_grid::next, nb::arg( "grid" ),
            "Increment in storage order (see index())." )
      .def( "next", ( const Coord_grid &( Coord_grid::* )( const Grid_range & )) & Coord_grid::next,
            nb::arg( "grid_range" ), "Increment in storage order (see index())." )
      .def( "last", ( bool ( Coord_grid::* )( const Grid & ) const ) & Coord_grid::last, nb::arg( "grid" ),
            "Test if done in storage order (see index())." )
      .def( "last", ( bool ( Coord_grid::* )( const Grid_range & ) const ) & Coord_grid::last, nb::arg( "grid_range" ),
            "Test if done in storage order (see index())." )
      .def( "index", &Coord_grid::index, nb::arg( "grid" ), "Grid indexing operator." )
      .def( "deindex", [](Coord_grid &self, const Grid &g, const itype64& ind) {
        return self.deindex(g, ind);
      },
        //&Coord_grid::deindex, nb::arg( "grid" ), nb::arg( "index" ),
        "Grid deindexing operator." )
      .def( "__str__", &Coord_grid::format )
      .def( "format", &Coord_grid::format, "Return formatted string representation." )
      .def( -nb::self )
      .def( nb::self + nb::self )
      .def( nb::self += nb::self, nb::rv_policy::none )
      .def( nb::self - nb::self )
      .def( operator-=( nb::self, nb::self ), nb::rv_policy::none )
      .def( nb::self == nb::self )
      .def( nb::self != nb::self )
      .def(
          "__mul__", []( const Coord_grid &self, const int &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Coord_grid &self, const int &s ) { return s * self; }, nb::is_operator() )
      //.def(
      //    "__eq__", []( const Coord_grid &self, const Coord_grid &other ) { return self == other; }, nb::is_operator() )
      //.def(
      //    "__ne__", []( const Coord_grid &self, const Coord_grid &other ) { return self != other; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Coord_grid &self, const Isymop &op ) { return op * self; }, nb::is_operator() )
      .def(
          "__mul__", []( const Coord_grid &self, const Coord_grid &other ) { return self * other; },
          nb::is_operator() );

  nb::class_<Coord_orth, Vec3<>> coord_orth( m, "Coord_orth", "Orthogonal (Angstrom) coordinates" );
  coord_orth.def( nb::init<>() )
      .def( nb::init<const Vec3<> &>(), nb::arg( "v" ), "Constructor copy/convert." )
      .def( nb::init<const ftype &, const ftype &, const ftype &>(), nb::arg( "x" ), nb::arg( "y" ), nb::arg( "z" ),
            "Constructor from x,y,z." )
      .def( nb::init<const Coord_orth &, const Coord_orth &, const Coord_orth &, const double &, const double &,
                     const double &>(),
            nb::arg( "x1" ), nb::arg( "x2" ), nb::arg( "x3" ), nb::arg( "length" ), nb::arg( "angle" ),
            nb::arg( "torsion" ), "Constructor: from 3 coords and bond length, angle, torsion" )
      .def(
          "__init__",
          []( Coord_orth *co, const std::array<double, 3> &a ) { new ( co ) Coord_orth( a[0], a[1], a[2] ); },
          "Constructor from list/array of x,y,z." )
      .def_prop_ro( "x", &Coord_orth::x, "Get x." )
      .def_prop_ro( "y", &Coord_orth::y, "Get y." )
      .def_prop_ro( "z", &Coord_orth::z, "Get z." )
      .def( "lengthsq", &Coord_orth::lengthsq, "Return square of length of vector in Angstroms." )
      .def( "coord_frac", &Coord_orth::coord_frac, nb::arg( "cell" ), "Orthogonal-fraction coordinate conversion." )
      .def( "transform", &Coord_orth::transform, nb::arg( "op" ), "Return transformed coordinate." )
      .def( "format", &Coord_orth::format, "Return formatted string representation." )
      .def( "__str__", &Coord_orth::format )
      .def( "__getstate__", [](const Coord_orth &self) {
        return nb::make_tuple(self.x(), self.y(), self.z());
      })
      .def( "__setstate__", [](Coord_orth &self, const std::tuple<ftype, ftype, ftype> &t) {
        new (&self) Coord_orth(std::get<0>(t), std::get<1>(t), std::get<2>(t));
      })
      .def_static( "length", &Coord_orth::length, nb::arg( "x1" ), nb::arg( "x2" ),
                   "Return length of vector between two coordinates." )
      .def_static( "angle", &Coord_orth::angle, nb::arg( "x1" ), nb::arg( "x2" ), nb::arg( "x3" ),
                   "Return angle between three orthogonal coordinates." )
      .def_static( "torsion", &Coord_orth::torsion, nb::arg( "x1" ), nb::arg( "x2" ), nb::arg( "x3" ), nb::arg( "x4" ),
                   "Return torsion betwee four orthogonal coordinates." )
      .def( -nb::self )
      .def( nb::self - nb::self )
      .def( nb::self + nb::self )
      .def( operator-=( nb::self, nb::self ), nb::rv_policy::none )
      .def( nb::self += nb::self, nb::rv_policy::none )
      .def(
          "__mul__", []( const Coord_orth &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Coord_orth &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def( "__rmul__", []( const Coord_orth &self, const RTop_orth &op ) { return op * self; }, nb::is_operator() );
  // inherited from Vec3: __iter__, as_array, from_array, buffer_protocol

  nb::class_<Coord_frac, Vec3<>> coord_frac( m, "Coord_frac", "Fractional (cell) coordinates" );
  coord_frac.def( nb::init<>() )
      .def( nb::init<const Vec3<> &>(), "Constructor copy/convert." )
      .def( nb::init<const ftype &, const ftype &, const ftype &>(), "Constructor from u,v,w." )
      .def(
          "__init__",
          []( Coord_frac *cf, const std::array<double, 3> &a ) { new ( cf ) Coord_frac( a[0], a[1], a[2] ); },
          "Constructor from list/array of u,v,w." )
      .def_prop_ro( "u", &Coord_frac::u, "Get u." )
      .def_prop_ro( "v", &Coord_frac::v, "Get v." )
      .def_prop_ro( "w", &Coord_frac::w, "Get w." )
      .def( "lengthsq", &Coord_frac::lengthsq, nb::arg( "cell" ), "Return square of length of vector in Angstrom." )
      .def( "coord_orth", &Coord_frac::coord_orth, nb::arg( "cell" ), "Fractional-orthogonal coordinate conversion." )
      .def( "coord_map", &Coord_frac::coord_map, nb::arg( "grid" ), "Fractional-map coordinate conversion." )
      .def( "coord_grid", &Coord_frac::coord_grid, nb::arg( "grid" ), "Fractional-grid coordinate conversion." )
      .def( "transform", &Coord_frac::transform, nb::arg( "op" ), "Return transformed coordinate" )
      .def( "lattice_copy_zero", &Coord_frac::lattice_copy_zero, "Return lattice copy nearest origin." )
      .def( "lattice_copy_unit", &Coord_frac::lattice_copy_unit,
            "Return lattice copy in unit box (0...1,0...1,0...1)." )
      .def( "lattice_copy_near", &Coord_frac::lattice_copy_near, "Return lattice copy near the specified coordinate." )
      .def( "symmetry_copy_near", &Coord_frac::symmetry_copy_near, nb::arg( "spacegroup" ), nb::arg( "cell" ),
            nb::arg( "cf" ), "Return symmetry copy near the specified coordinate." )
      .def( "__str__", &Coord_frac::format )
      .def( "format", &Coord_frac::format, "Return formatted string representation." )
      .def( "__getstate__", [](const Coord_frac &self) {
        return nb::make_tuple(self.u(), self.v(), self.w());
      })
      .def( "__setstate__", [](Coord_frac &self, const std::tuple<ftype, ftype, ftype> &t) {
        new (&self) Coord_frac(std::get<0>(t), std::get<1>(t), std::get<2>(t));
      })
      .def( -nb::self )
      .def( nb::self + nb::self )
      .def( nb::self += nb::self, nb::rv_policy::none )
      .def( nb::self - nb::self )
      .def( operator-=( nb::self, nb::self ), nb::rv_policy::none )
      .def(
          "__mul__", []( const Coord_frac &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Coord_frac &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Coord_frac &self, const RTop_frac &op ) { return op * self; },
          nb::is_operator() ); // RTop_frac is from core/symop.h

  nb::class_<Coord_map, Vec3<>> coord_map( m, "Coord_map",
                                           "Map coordinate.\nThis is like Coord_grid, but non-integer." );
  coord_map.def( nb::init<>() )
      .def( nb::init<const Vec3<> &>(), "Constructor copy/convert." )
      .def( nb::init<const Coord_grid &>(), nb::arg( "c" ), "Constructor from Coord_grid." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "u" ), nb::arg( "v" ), nb::arg( "w" ),
            "Constructor from u,v,w." )
      .def( "coord_frac", &Coord_map::coord_frac, nb::arg( "grid" ), "Grid-fractional coordinate conversion." )
      .def( "coord_grid", &Coord_map::coord_grid, "Return interger Coord_grid nearest to this coordinate." )
      .def( "floor", &Coord_map::floor, "Return integer Coord_grid below this coordinate." )
      .def( "ceil", &Coord_map::ceil, "Return integer Coord_grid above this coordinate." )
      .def_prop_ro( "u", &Coord_map::u, "Get u." )
      .def_prop_ro( "v", &Coord_map::v, "Get v." )
      .def_prop_ro( "w", &Coord_map::w, "Get w." )
      .def( "__str__", &Coord_map::format )
      .def( "format", &Coord_map::format, "Return formatted string representation." )
      .def( -nb::self )
      .def( nb::self - nb::self )
      .def( operator-=( nb::self, nb::self ), nb::rv_policy::none )
      .def( nb::self + nb::self )
      .def( nb::self += nb::self, nb::rv_policy::none )
      .def(
          "__mul__", []( const Coord_map &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def( "__rmul__", []( const Coord_map &self, const ftype &s ) { return s * self; }, nb::is_operator() );

  nb::class_<U_aniso_orth, Mat33sym<>> uaniso_orth( m, "U_aniso_orth" );
  uaniso_orth.def( nb::init<>() )
      .def( nb::init<const Mat33sym<> &>(), nb::arg( "m" ), "Constructor from Mat33sym." )
      .def( nb::init<const double &>(), nb::arg( "u" ), "Constructor from isotropic U." )
      .def( nb::init<const double &, const double &, const double &, const double &, const double &, const double &>(),
            nb::arg( "u11" ), nb::arg( "u22" ), nb::arg( "u33" ), nb::arg( "u12" ), nb::arg( "u13" ), nb::arg( "u23" ),
            "Constructor from Uij." )
      .def( "u_iso", &U_aniso_orth::u_iso, "Return nearest isotropic U." )
      .def( "u_aniso_frac", &U_aniso_orth::u_aniso_frac, nb::arg( "cell" ), "Orthogonal-fractional conversion." )
      .def( "transform", &U_aniso_orth::transform, nb::arg( "op" ), "Return transformed U_aniso." )
      .def( -nb::self )
      .def( nb::self + nb::self )
      .def(
          "__mul__", []( const U_aniso_orth &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const U_aniso_orth &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def( "__getstate__", [](const U_aniso_orth &u) {
        return nb::make_tuple(u.mat00(), u.mat11(), u.mat22(), u.mat01(), u.mat02(), u.mat12());
      })
      .def( "__setstate__", [](U_aniso_orth &self, const nb::tuple &t) {
        if (t.size() != 6) throw std::runtime_error("Invalid size!");
        
        new (&self) U_aniso_orth(nb::cast<ftype>(t[0]), nb::cast<ftype>(t[1]), nb::cast<ftype>(t[2]),
                                nb::cast<ftype>(t[3]), nb::cast<ftype>(t[4]), nb::cast<ftype>(t[5]));
      })
      .doc() = "Anisotropic orthogonal atomic displacement parameters.\n"
               "These are defined on orthogonal atomic coordinates in "
               "A<sup>-2</sup>, i.e. they are anisotropic U values.";

  nb::class_<U_aniso_frac, Mat33sym<>> uaniso_frac( m, "U_aniso_frac" );
  uaniso_frac.def( nb::init<>() )
      .def( nb::init<const Mat33sym<> &>(), nb::arg( "m" ), "Constructor from Mat33sym." )
      .def( nb::init<const double &, const double &, const double &, const double &, const double &, const double &>(),
            nb::arg( "u11" ), nb::arg( "u22" ), nb::arg( "u33" ), nb::arg( "u12" ), nb::arg( "u13" ), nb::arg( "u23" ),
            "Constructor from Uij." )
      .def( "u_aniso_orth", &U_aniso_frac::u_aniso_orth, nb::arg( "cell" ), "Fractional-orthogonal conversion." )
      .def( "transform", &U_aniso_frac::transform, nb::arg( "op" ), "Return transformed U_aniso." )
      .def( nb::self + nb::self )
      .def( -nb::self )
      .def(
          "__mul__", []( const U_aniso_frac &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__rmul__", []( const U_aniso_frac &self, const ftype &s ) { return s * self; }, nb::is_operator() )
      .doc() = "Anisotropic fractional atomic displacement parameters.\n"
               "These are defined on fractional atomic coordinates in "
               "A<sup>-2</sup>, i.e. they are anisotropic U values.";

  nb::class_<Grid, Vec3<int>> grid( m, "Grid" );
  grid.def( nb::init<>() )
      .def( nb::init<const int &, const int &, const int &>(), nb::arg( "nu" ), nb::arg( "nv" ), nb::arg( "nw" ),
            "Constructor from nu,nv,nw." )
      .def_prop_ro( "nu", &Grid::nu, "Get nu." )
      .def_prop_ro( "nv", &Grid::nv, "Get nv." )
      .def_prop_ro( "nw", &Grid::nw, "Get nw." )
      .def( "size", &Grid::size, "Return size of grid array." )
      .def( "in_grid", &Grid::in_grid, nb::arg( "grid" ), "Determin if a point is in the grid." )
      .def( "index", &Grid::index, nb::arg( "grid" ), "Grid indexing operator." )
      .def( "deindex", [](Grid &self, const itype64& ind) {
        return self.deindex(ind);
      },
         //&Grid::deindex, nb::arg( "index" ),
         "Grid deindexing operator." )
      .def( "__str__", &Grid::format )
      .def( "format", &Grid::format, "Return formatted string representation." )
      .def( "debug", &Grid::debug, "Output debug details." )
      .doc() = "Generic grid.\nThis holds the dimensions of a 3D array, "
               "indexed from 0 along each dimension.";

  nb::class_<Grid_sampling, Grid> grid_sampling( m, "Grid_sampling" );
  grid_sampling.def( nb::init<>() )
      .def( nb::init<const int &, const int &, const int &>(), nb::arg( "nu" ), nb::arg( "nv" ), nb::arg( "nw" ),
            "Constructor from nu,nv,nw." )
      .def( nb::init<const Spacegroup &, const Cell &, const Resolution &, const ftype>(), nb::arg( "spacegroup" ),
            nb::arg( "cell" ), nb::arg( "resolution" ), nb::arg( "rate" ) = 1.5,
            "Constructor from spacegroup, cell, resolution, Shannon rate." )
      .def( "init", &Grid_sampling::init, nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "resolution" ),
            nb::arg( "rate" ) = 1.5, "Initialiser from spacegroup, cell, resolution, Shannon rate." )
      .def( "matrix_grid_frac", &Grid_sampling::matrix_grid_frac,
            "Return matrix which converts grid to fractional coordinates." )
      .def( "matrix_frac_grid", &Grid_sampling::matrix_frac_grid,
            "Return matrix which converts fractional to grid coordinates." )
      .def( "is_null", &Grid_sampling::is_null, "Test if object has been initialised." )
      .def( "__setstate__", [](Grid_sampling &self, const std::tuple<int, int, int> &t) {
        new (&self) Grid_sampling(std::get<0>(t), std::get<1>(t), std::get<2>(t));
      })
      .doc() = "Grid sampling of a unit cell.\nThis class represents the grid "
               "sampling of a unit cell. It is otherwise identical to its parent, "
               "clipper::Grid_cell, but has an additional constructor which takes a "
               "spacegroup, cell and resolution and produces an appropriate grid "
               "obeying all of the symmetry constraints, and using efficient "
               "factors for the calculation of FFTs.";

  nb::class_<HKL_sampling> hkl_sampling( m, "HKL_sampling" );
  hkl_sampling.def( nb::init<>() )
      .def( nb::init<const Cell &, const Resolution &>(), nb::arg( "cell" ), nb::arg( "resolution" ),
            "Constructor from cell (normal or inverse) and resolution." )
      .def( "hkl_limit", &HKL_sampling::hkl_limit, "Return limiting values of H, K, L." )
      .def( "resolution", &HKL_sampling::resolution, nb::arg( "cell" ), "Return approximate resolution given cell." )
      .def( "in_resolution", &HKL_sampling::in_resolution, nb::arg( "hkl" ),
            "Test if a reflection is within the resolution limit." )
      .def( "is_null", &HKL_sampling::is_null, "Test if object has been initialised. " )
      .def( "__str__", &HKL_sampling::format )
      .def( "format", &HKL_sampling::format, "Return formatted string representation." )
      .def( nb::self == nb::self )
      //.def(
      //    "__eq__", []( const HKL_sampling &self, const HKL_sampling &other ) { return self == other; },
      //    nb::is_operator() )
      .doc() = "HKL sampling of reciprocal space.\nThe HKL_sampling class "
               "uniquely describes a P0 reflection list bounded by some "
               "resolution limit in reciprocal space. It is described in "
               "terms of large integers, and so immune from rounding errors "
               "once the object is constructed.";

  nb::class_<Grid_range, Grid> grid_range( m, "Grid_range" );
  grid_range.def( nb::init<>() )
      .def( nb::init<const Coord_grid &, const Coord_grid &>(), nb::arg( "min" ), nb::arg( "max" ),
            "Constructor from grid limits (Coord_grid)" )
      .def( nb::init<const Grid &, const Coord_frac &, const Coord_frac &>(), nb::arg( "grid" ), nb::arg( "min" ),
            nb::arg( "max" ), "Constructor from cell grid and fractional limits." )
      .def( nb::init<const Cell &, const Grid &, const ftype &>(), nb::arg( "cell" ), nb::arg( "grid" ),
            nb::arg( "radius" ), "Constructor: make grid to hold a sphere from cell, grid and radius." )
      .def_prop_ro( "min", &Grid_range::min, "Access grid limits." )
      .def_prop_ro( "max", &Grid_range::max, "Access grid limits." )
      .def( "add_border", &Grid_range::add_border, nb::arg( "b" ), "Border: increase grid to include given border." )
      .def( "in_grid", &Grid_range::in_grid, nb::arg( "grid" ), "Determine if a point is in the grid." )
      .def( "index", &Grid_range::index, nb::arg( "c" ), "Grid indexing operator." )
      .def( "deindex", [](Grid_range &self, const itype64& ind) {
        return self.deindex(ind);
      },
        //&Grid_range::deindex, nb::arg( "index" ),
        "Grid deindexing operator." )
      .doc() = "Grid range class: defines array limits for a grid.\nThis class "
               "is used for describing 3D grids covering an arbitrary "
               "part of the 3D space, i.e. which do not start from (0,0,0).";

  // probably need something more;
}