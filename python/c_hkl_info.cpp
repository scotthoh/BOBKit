// Nanobind bindings for clipper hkl_info
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include "arrays.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <nanobind/operators.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/tuple.h>
// #include <gemmi/symmetry.hpp>
// #include <gemmi/unitcell.hpp>

using namespace clipper;

void add_hklinfo( nb::module_ &m ) {
  nb::class_<HKL_info> hklinfo( m, "HKL_info" );
  hklinfo.def( nb::init<>(), "Null constructor" )
      .def( nb::init<const Spacegroup &, const Cell &, const Resolution &, const bool &>(), nb::arg( "spacegroup" ),
            nb::arg( "cell" ), nb::arg( "resolution" ), nb::arg( "generate" ) = false,
            "Constructor: takes Spacegroup, Cell and Resolution types." )
      .def(
          "__init__",
          []( HKL_info *hklinfo, const Spacegroup &sg, const Cell &c, const double &res, const bool &gen ) {
            clipper::Resolution reso( res );
            new ( hklinfo ) HKL_info( sg, c, reso, gen );
          },
          nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "resolution" ), nb::arg( "generate" ) = false,
          "Constructor: takes Spacegroup, Cell and resolution(double)." )
      .def(
          "init",
          ( void ( HKL_info::* )( const Spacegroup &, const Cell &, const Resolution &, const bool & ) )&HKL_info::init,
          nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "resolution" ), nb::arg( "generate" ) = false,
          "Initialiser: takes Spacegroup, Cell, and Resolution types." )
      .def( "init",
            ( void ( HKL_info::* )( const Spacegroup &, const Cell &, const HKL_sampling &,
                                    const bool & ) )&HKL_info::init,
            nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "hkl_sampling" ), nb::arg( "generate" ) = true,
            "Initialiser: takes Spacegroup, Cell, and HKL_sampling" )
      // from gemmi types
      .def(
          "init",
          []( HKL_info &self, const gemmi::SpaceGroup &sg, const gemmi::UnitCell &uc, const double dmin,
              const double tol, const bool &generate ) {
            self.init( GEMMI::spacegroup( sg ), GEMMI::cell( uc ), Resolution( dmin - tol ), generate );
          },
          nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "dmin" ), nb::arg( "tol" ) = 1.e-8,
          nb::arg( "generate" ) = false )
      .def( "is_null", &HKL_info::is_null, "Test if object has been initialised" )
      .def_prop_ro( "cell", &HKL_info::cell, "Get the cell" )
      .def_prop_ro( "spacegroup", &HKL_info::spacegroup, "Get the spacegroup" )
      .def_prop_ro( "hkl_sampling", &HKL_info::hkl_sampling, "Get HKL_sampling" )
      .def_prop_ro( "resolution", &HKL_info::resolution, "Get the resolution" )

      // from gemmi::Mtz
      .def_static(
          "from_gemmi_mtz",
          []( const gemmi::Mtz &mtzobj, const double tol, const bool &generate ) {
            return GEMMI::as_HKL_info( mtzobj, tol, generate );
          },
          nb::arg( "mtz" ), nb::arg( "tol" ) = 1.e-8, nb::arg( "generate" ) = false )
      .def( "generate_hkl_list", &HKL_info::generate_hkl_list, "Synthesize hkl list" )
      //.def("add_hkl_list", &HKL_info::add_hkl_list, nb::arg("hkls"))

      //.def("add_hkl_list",
      //     (void(HKL_info::*)(const std::vector<HKL> &add)) &
      //         HKL_info::add_hkl_list,
      //     nb::arg("hkls"))
      .def(
          "add_hkl_list", []( HKL_info &self, std::vector<HKL> &hkl ) { self.add_hkl_list( hkl ); },
          "Add new reflections to the list" )
      .def_prop_ro( "hkl_array", [](const HKL_info &self) {
        int* arr = new int[self.num_reflections()*3];
        std::initializer_list<int64_t> strides={3,1};
        for ( int i = 0; i < self.num_reflections(); i++ ) {
            auto hkl = self.hkl_of(i);
            arr[i*3] = hkl.h();
            arr[i*3+1] = hkl.k();
            arr[i*3+2] = hkl.l();
        }
        nb::capsule owner(arr, [](void *p) noexcept {delete[] static_cast<int*>(p);});
        return nb::ndarray<int, nb::numpy, nb::ndim<2>, nb::device::cpu, nb::c_contig>(arr, {(size_t)self.num_reflections(), 3}, owner, strides);

      },nb::rv_policy::automatic)
      .def( "num_reflections", &HKL_info::num_reflections, "Get number of reflections in the object" )
      .def( "hkl_of", &HKL_info::hkl_of, nb::arg( "index" ), "Return the corresponding HKL to the index given." )
      .def( "index_of", &HKL_info::index_of, nb::arg( "hkl" ),
            "Reflection index from hkl. This does not check symmetry equivalence." )
      .def( "invresolsq", &HKL_info::invresolsq, nb::arg( "index" ), "Get reflection resolution using lookup" )
      .def_prop_ro( "invresolsq_range", &HKL_info::invresolsq_range, "Get resolution limits of the list" )
      .def( "hkl_class", &HKL_info::hkl_class, nb::arg( "index" ), "Get reflection class using lookup" )
      .def( "find_sym", &HKL_info::find_sym, nb::arg( "hkl" ), nb::arg( "sym" ), nb::arg( "friedel" ),
            "Find symop number and friedel to bring an HKL into ASU" )
      .def( "first", &HKL_info::first, "Return HKL_reference_index pointing to first reflection" )
      .def( "debug", &HKL_info::debug, "Debug: print number of reflections" )
      .def( "__repr__",
            []( const HKL_info &self ) {
              return "<clipper.HKL_info with spacegroup " + self.spacegroup().symbol_hm() + ", " +
                     clipper::String( self.num_reflections() ) + " reflections>";
            } )
      .doc() = "HKL list container and tree root.\n"
               "This object contains contains a reflection list, and all the "
               "properties on which such a list depends, i.e. spacegroup, cell, "
               "resolution. It also keeps a fast reflection lookup list and lookup "
               "lists for resolutions and reflection classes.";

  using HKLB = HKL_info::HKL_reference_base;
  nb::class_<HKLB>( hklinfo, "HKL_reference_base" )
      .def_prop_ro( "base_hkl_info", &HKLB::base_hkl_info, "Return the parent HKL_info" )
      .def_prop_ro( "index", &HKLB::index, "Return the current index(-1 if invalid)" )
      .def( "invresolsq", ( ftype ( HKLB::* )( const HKL_data_base & ) const ) & HKLB::invresolsq, nb::arg( "hkldata" ),
            "Return the inverse resolution squared for the reflections (assumes index valid)" )
      .def( "invresolsq", ( ftype ( HKLB::* )() const ) & HKLB::invresolsq,
            "Return the inverse resolution squared for the reflection of current index (assumes index valid)" )
      .def( "last", &HKLB::last, "Test if index has gone past last reflection" )
      .doc() = "HKL reference base class\n"
               "This is a reference to an HKL. It forms a base class for "
               "index-like and coordinate-like HKL references. If you write a "
               "method which will work with either, then specify this instead of "
               "either of the derived classed. ";

  using HKLI = HKL_info::HKL_reference_index;
  nb::class_<HKLI, HKLB>( hklinfo, "HKL_reference_index" )
      .def( nb::init<>(), "Null constructor" )
      .def( nb::init<const HKL_info &, const int &>(), nb::arg( "hklinfo" ), nb::arg( "index" ),
            "Constructor: takes parent HKL_info and initial index" )
      .def( "hkl", &HKLI::hkl, "Return the current HKL" )
      .def( "hkl_class", &HKLI::hkl_class, "Return the reflection class for the reflection" )
      // note from Tristan: avoid creating new Python objects when incrementing
      .def(
          "next", []( HKLI &self ) { self.next(); }, "Increment to the next reflection" )
      .doc() = "HKL reference with index-like behaviour\n"
               "This is a reference to an HKL. It behaves like a simple index "
               "into the reflection list, but can be easily converted into an "
               "HKL as and when required. It also implements methods for "
               "iterating through a reflection list. "
               "NOTE: The following methods are inherited from "
               "HKL_reference_base but are documented here for convenience: "
               "base_hkl_info(), index(), invresolsq(), last(). ";

  using HKLC = HKL_info::HKL_reference_coord;
  nb::class_<HKLC, HKLB>( hklinfo, "HKL_reference_coord" )
      .def( nb::init<>(), "Null constructor" )
      .def( nb::init<const HKL_info &, const HKL &>(), nb::arg( "hklinfo" ), nb::arg( "hkl" ),
            "Constructor: takes parent HKL_info and initial HKL" )
      .def_prop_rw(
          "hkl", &HKLC::hkl, []( HKLC &self, const HKL &hkl ) { self.set_hkl( hkl ); }, "Return/Set current HKL" )
      .def_prop_ro( "sym", &HKLC::sym, "Get current symop number" )
      .def_prop_ro( "friedel", &HKLC::friedel, "Get current friedel flag" )
      // note from Tristan: avoid creating new Python objects when incrementing
      .def(
          "next", []( HKLC &self ) { self.next(); }, "Increment to next reflection" )
      .def(
          "next_h", []( HKLC &self ) { self.next_h(); }, "Increment to next h" )
      .def(
          "next_k", []( HKLC &self ) { self.next_k(); }, "Increment to next k" )
      .def(
          "next_l", []( HKLC &self ) { self.next_l(); }, "Increment to next l" )
      .def(
          "prev_h", []( HKLC &self ) { self.prev_h(); }, "Decrement to previous h" )
      .def(
          "prev_k", []( HKLC &self ) { self.prev_k(); }, "Decrement to previous h" )
      .def(
          "prev_l", []( HKLC &self ) { self.prev_l(); }, "Decrement to previous h" )
      .doc() = "HKL reference with coord-like behaviour.\nThis is a reference to an HKL. "
               "It behaves like an HKL, but also stores the index of the corresponding reflection in the "
               "reflection list, if it exists, and the symmetry and friedel operators required to get there. "
               "NOTE: The following methods are inherited from HKL_reference_base but are documented here "
               "for convenience: base_hkl_info(), index(), invresolsq(), last(). ";
}