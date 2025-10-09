// Nanobind bindings for clipper cell.h
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper nxmap by Tristan Croll

#include "commons.h"
#include "arrays.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
// #include <nanobind/ndarray.h>
#include <nanobind/operators.h>
//#include <nanobind/stl/array.h>

using namespace clipper;

template <class T> std::string display_repr( const NXmap<T> nxmap ) {
  if ( strcmp( typeid( T ).name(), "f" ) == 0 )
    return "<clipper.NXmap_float: Grid " + nxmap.grid().format() + ">";
  else if ( strcmp( typeid( T ).name(), "d" ) == 0 )
    return "<clipper.NXmap_double: Grid " + nxmap.grid().format() + ">";
  else if ( strcmp( typeid( T ).name(), "i" ) == 0 )
    return "<clipper.NXmap_int: Grid " + nxmap.grid().format() + ">";
  // return "<clipper.Xmap_int>";
  else
    return "<clipper.NXmap>";
}

template<typename T>
auto nxmap_to_array(NXmap<T> &xm) {
  auto grid = xm.grid();
  return nb::ndarray<nb::numpy, T>(&xm[xm.first()],
                                   {(size_t)grid.nu(), (size_t)grid.nv(), (size_t)grid.nw()},
                                   nb::handle(),
                                   {grid.nv()*grid.nw(), grid.nw(), 1});
                                   //{1, grid.nu(), grid.nu()*grid.nv()});
}

// numpy array to xmap
template <typename T>
void numpy_to_nxmap(const nb::ndarray<nb::numpy, T, nb::ndim<3>> & data, const Cell &cell, NXmap<T> &nxmap,
                   const std::array<int, 3> &axis_pos) {
  auto d = data.view();
  if ( d.ndim() != 3 ) throw std::runtime_error("Array must be 3-dimension!");
  
  // set map params
  Grid_sampling grid_sam(int(d.shape(0)), int(d.shape(1)), int(d.shape(2)));
  Grid_range grid_extent(Coord_grid(0, 0, 0), Coord_grid(int(d.shape(0)-1), int(d.shape(1)-1), int(d.shape(2)-1)));
  nxmap.init(cell, grid_sam, grid_extent);
  int g[3], gfms[3], orderxyz[3];
  for ( int i = 0; i < 3; ++i ) {
    gfms[i] = d.shape(i) - 1;
  }

  // copy map data
  for (g[0] = 0; g[0] <= gfms[0]; g[0]++)
    for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
      for (g[2] = 0; g[2] <= gfms[2]; g[2]++) {
        nxmap.set_data(Coord_grid(g[0], g[1], g[2]), T(d(g[0], g[1], g[2])));
        //ix.set_coord(clipper::Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
        //xmap[ix] = T(d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
        //xmap[ix] = d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]);
        //std::cout << ix.coord().format() << " = " << String(xmap[ix]) << " ; " << String(T(d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]))) << std::endl;
      }
  //std::cout << "TEST" << std::endl;
  //for (g[0] = 0; g[0] <= gfms[0]; g[0]++)
  //  for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
  //    for (g[2] = 0; g[2] <= gfms[2]; g[2]++) {
  //      ix.set_coord(clipper::Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
  //      //std::cout << ix.coord().format() << " = " << String(xmap[ix]) << std::endl;
  //    }
  //std::cout << "TEST" << std::endl;
  // copy map data
  //Xmap_base::Map_reference_coord x(xmap);
  //for (g[2] = 0; g[2] <= gfms[2]; g[2]++)
  //  for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
  //    for (g[0] = 0; g[0] <= gfms[0]; g[0]++) {
  //      x.set_coord(clipper::Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
  //      xmap[x] = T(d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
  //      //std::cout << x.coord().format() << " = " << String(xmap[x]) << " ; " << String(T(d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]))) << std::endl;
  //    }
}

void declare_map_base_methods( nb::module_ &m, nb::class_<NXmap_base> &nxmapbase ) {
  // nb::class_<NXmap_base> nxmapbase(m, "NXmap_base");
  nxmapbase.def( "is_null", &NXmap_base::is_null, "Test if object has been initialised." )
      .def_prop_ro( "grid", &NXmap_base::grid, "Return the grid dimensions for this map." )
      .def_prop_ro( "operator_orth_grid", &NXmap_base::operator_orth_grid,
                    "Return the orthogonal-to-grid coordinate operator." )
      .def_prop_ro( "operator_grid_orth", &NXmap_base::operator_grid_orth,
                    "Return the grid-to-orthogonal coordinate operator." )
      .def( "coord_orth", &NXmap_base::coord_orth, nb::arg( "coord_map" ), "Convert map coordinate to orthogonal." )
      .def( "coord_map", &NXmap_base::coord_map, nb::arg( "coord_orth" ), "Convert orthogonal coordinate to map." )
      .def( "in_map", ( bool ( NXmap_base::* )( const clipper::Coord_grid & ) const ) & NXmap_base::in_map,
            nb::arg( "coord_grid" ), "Check if given coord available in the map?" )
      // in_map using different interpolators
      .def( "in_map_interp_nearest",
            ( bool ( NXmap_base::* )( const clipper::Coord_map & ) const ) & NXmap_base::in_map<Interp_nearest>,
            nb::arg( "coord_map" ), "Check if given coord available in the map?" )
      .def( "in_map_interp_linear",
            ( bool ( NXmap_base::* )( const clipper::Coord_map & ) const ) & NXmap_base::in_map<Interp_linear>,
            nb::arg( "coord_map" ), "Check if given coord available in the map?" )
      .def( "in_map_interp_cubic",
            ( bool ( NXmap_base::* )( const clipper::Coord_map & ) const ) & NXmap_base::in_map<Interp_cubic>,
            nb::arg( "coord_map" ), "Check if given coord available in the map?" )
      // TODO: check if the inmap methods are working
      .def( "multiplicity", &NXmap_base::multiplicity, nb::arg( "coord_grid" ) = nullptr,
            "Get multiplicity of a map grid point (always 1 for NXmap)." )
      .def( "first", &NXmap_base::first, "Return a basic Map_reference_index for this map." )
      .def( "first_coord", &NXmap_base::first_coord, "Return a coord Map_reference_index for this map" )
      .doc() = "NXmap_base: base for non-crystallographic map class.\n"
               "The non-crystallographic map class stores a map of arbitrary "
               "data type. Unlike an Xmap it is finite in extent and has no "
               "symmetry. An RT operator provides mapping onto an arbitrary "
               "orthogonal coordinate frame. Iterators provide efficient "
               "access to data. This base contains everything except the data, "
               "which is templated in the derived type clipper::NXmap<T>.";
}

void declare_map_reference_base( nb::module_ &m, nb::class_<NXmap_base> &nxmapbase ) {
  using MRB = NXmap_base::Map_reference_base;
  nb::class_<MRB>( nxmapbase, "NXmap_reference_base" )
      .def_prop_ro( "base_nxmap", &MRB::base_nxmap, "Return the parent NXmap." )
      .def_prop_ro( "index", &MRB::index, "Get the index into the map data array." )
      .def( "last", &MRB::last, "Check for end of map." )
      .doc() = "Map reference base class.\n"
               "This is a reference to an Map. It forms a base class for "
               "index-like and coordinate-like Map references.";
}

void declare_map_reference_index( nb::module_ &m, nb::class_<NXmap_base> &nxmapbase ) {
  using MRI = NXmap_base::Map_reference_index;
  nb::class_<MRI, NXmap_base::Map_reference_base>( nxmapbase, "NXmap_reference_index" )
      .def_prop_rw(
          "coord", &MRI::coord, []( MRI &self, const Coord_grid &pos ) -> void { self.set_coord( pos ); },
          nb::for_getter( "Get current grid coordinate." ), nb::for_setter( "Set current grid coordinate." ) )
      .def( "coord_orth", &MRI::coord_orth )
      .def( "next", []( MRI &self ) -> void { self.next(); } )
      .def( "index_offset", &MRI::index_offset, nb::arg( "du" ), nb::arg( "dv" ), nb::arg( "dw" ) )
      .doc() = "Map reference with index-like behaviour.\n This is a reference "
               "to a map coordinate. It behaves like a simple index into the "
               "map, but can be easily converted into a coordinate as and when "
               "required. It also implements methods for iterating through the "
               "map. It is very compact, but coord() involves some overhead.";
}

void declare_map_reference_coord( nb::module_ &m, nb::class_<NXmap_base> &nxmapbase ) {
  using MRC = NXmap_base::Map_reference_coord;
  nb::class_<MRC, NXmap_base::Map_reference_base>( nxmapbase, "NXmap_reference_coordinate" )
      .def_prop_rw(
          "coord", &MRC::coord, []( MRC &self, const Coord_grid &pos ) -> void { self.set_coord( pos ); },
          nb::for_getter( "Get current grid coordinate." ), nb::for_setter( "Set current grid coordinate." ) )
      .def( "coord_orth", &MRC::coord_orth )
      .def( "next", []( MRC &self ) -> void { self.next(); } )
      .def( "next_u", []( MRC &self ) -> void { self.next_u(); } )
      .def( "next_v", []( MRC &self ) -> void { self.next_v(); } )
      .def( "next_w", []( MRC &self ) -> void { self.next_w(); } )
      .def( "prev_u", []( MRC &self ) -> void { self.prev_u(); } )
      .def( "prev_v", []( MRC &self ) -> void { self.prev_v(); } )
      .def( "prev_w", []( MRC &self ) -> void { self.prev_w(); } );
}

// void declare_map_reference_base(nb::module_ &m) {}
//
// void declare_map_reference_index(nb::module_ &m) {}
//
// void declare_map_reference_coord(nb::module_ &m) {}

template <class Derived, class T> void apply_nxmap_interpolation_methods( nb::class_<Derived, NXmap_base> &pyclass ) {
  // need to check if they are working
  pyclass
      .def(
          "interp_linear",
          []( const Derived &self, const Coord_map &pos ) { return self.template interp<Interp_linear>( pos ); },
          nb::arg( "pos" ), "Get map value for map coordinate using linear interpolator." )
      .def(
          "interp_cubic",
          []( const Derived &self, const Coord_map &pos ) { return self.template interp<Interp_cubic>( pos ); },
          nb::arg( "pos" ), "Get map value for map coordinate using cubic interpolator." )
      .def(
          "interp_grad_cubic",
          []( const Derived &self, const Coord_map &pos ) {
            T val;
            Grad_map<T> grad;
            if ( self.template in_map<Interp_cubic>( pos ) ) {
              self.template interp_grad<Interp_cubic>( pos, val, grad );
            } else {
              Message::message( Message_fatal( "Interpolation point outside map - " + pos.format() ) );
            }
            return nb::make_tuple( val, grad );
          },
          nb::arg( "pos" ), "Get map value and gradient for map coordinate using cubic interpolator." )
      .def(
          "interp_curv_cubic",
          []( const Derived &self, const Coord_map &pos ) {
            T val;
            Grad_map<T> grad;
            Curv_map<T> curv;
            if ( self.template in_map<Interp_cubic>( pos ) ) {
              self.template interp_curv<Interp_cubic>( pos, val, grad, curv );
            } else {
              Message::message( Message_fatal( "Interpolation point outside map - " + pos.format() ) );
            }
            return nb::make_tuple( val, grad, curv );
          },
          nb::arg( "pos" ), "Get map value, gradient and curvatures for map coordinate using cubic interpolator." );
}

template <class T> void declare_nxmap( nb::module_ &m, const std::string &name ) {
  using MRI = NXmap_base::Map_reference_index;
  using MRC = NXmap_base::Map_reference_coord;
  using NXMClass = NXmap<T>;
  std::string PyClass_name = std::string( "NXmap_" ) + name;
  nb::class_<NXMClass, NXmap_base> nxmap( m, PyClass_name.c_str() );
  // nb::buffer_protocol());
  nxmap.def( nb::init<>() )
      .def( nb::init<const Grid &, const RTop<> &>() )
      .def( nb::init<const Cell &, const Grid_sampling &, const Grid_range &>() )
      .def( "init", ( void ( NXMClass::* )( const Grid &, const RTop<> & ) )&NXMClass::init, nb::arg( "grid" ),
            nb::arg( "rtop" ) )
      .def( "init", ( void ( NXMClass::* )( const Cell &, const Grid_sampling &, const Grid_range & ) )&NXMClass::init,
            nb::arg( "cell" ), nb::arg( "grid_sampling" ), nb::arg( "grid_range" ) )
      .def( "__init__", []( NXMClass *nxmap, const nb::ndarray<nb::numpy, T, nb::ndim<3>> &data, 
                            const Cell &cell) {//}, const std::array<int, 3> &axis_pos) {
        new ( nxmap ) NXMClass();
        std::array<int, 3> axis_pos = {0,1,2};
        numpy_to_nxmap<T>(data, cell, *nxmap, axis_pos);
      }, nb::arg("array").noconvert(), nb::arg("cell")) //, nb::arg("axis_pos"))
      .def_prop_ro( "array", [](NXMClass& xm) {
        return nxmap_to_array<T>(xm); }, nb::rv_policy::reference_internal)
      //.def_buffer([=](NXMClass &self) -> nb::buffer_info {
      //  return nb::buffer_info(
      //      &self[self.first()],
      //      {self.grid().nu(), self.grid().nv(),
      //       self.grid().nw()}, // dimensions
      //      {sizeof(T),
      //       sizeof(T) * self.grid().nu(), // strides
      //       sizeof(T) * self.grid().nu() * self.grid().nv()});
      //})
      .def( "__repr__", []( const NXMClass &self ) { return display_repr<T>( self ); } )
      .def( "__str__", []( const NXMClass &self ) { return display_repr<T>( self ); } )
      .def( "__getitem__", []( const NXMClass &self, const MRI &ix ) { return self[ix]; } )
      .def( "__setitem__", []( NXMClass &self, const MRI &ix, const T &val ) { self[ix] = val; } )
      .def( "__getitem__", []( const NXMClass &self, const MRC &ix ) { return self[ix]; } )
      .def( "__setitem__", []( NXMClass &self, const MRC &ix, const T &val ) { self[ix] = val; } )
      .def( "get_data", []( const NXMClass &self, const Coord_grid &pos ) { return self.get_data( pos ); } )
      .def( "set_data", []( NXMClass &self, const Coord_grid &pos, const T &val ) { self.set_data( pos, val ); } )
      .def( "fill_map_with", []( NXMClass &self, const T &val ) { self = val; } )
      .def( nb::self += nb::self )
      .def( nb::self -= nb::self )
      // 
      .def(
        "import_from_gemmi",
        [](NXMClass &self, const gemmi::Ccp4<float> &mapobj) {
          GEMMI::import_nxmap(self, mapobj);
        },
        nb::arg("target"), "Import NXmap from gemmi.Ccp4Map.")
        .def(
          "export_to_gemmi",
          [](const NXMClass &self, gemmi::Ccp4<float> &mapobj,
            const Cell &unitcell) {
              GEMMI::export_nxmap(self, mapobj, unitcell);
            },
            nb::arg("target"), nb::arg("cell"),
            "Export NXmap from gemmi.Ccp4Map.");
  apply_nxmap_interpolation_methods<NXMClass, T>( nxmap );
}

void add_nxmap( nb::module_ &m ) {
  nb::class_<NXmap_base> nxmapbase( m, "NXmap_base" );
  declare_map_base_methods( m, nxmapbase );
  declare_map_reference_base( m, nxmapbase );
  declare_map_reference_index( m, nxmapbase );
  declare_map_reference_coord( m, nxmapbase );
  declare_nxmap<int>( m, "int" );
  declare_nxmap<ftype32>( m, "float" );
  declare_nxmap<ftype64>( m, "double" );
}