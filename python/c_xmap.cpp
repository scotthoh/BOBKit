// Nanobind bindings for clipper xmap
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York
// Adapted and rewritten from pybind11 bindings for clipper xmap by Tristan Croll

#include "commons.h"
#include "arrays.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <nanobind/operators.h>
// #include <nanobind/stl/stl_bind.h>

using namespace clipper;

// orthogonal cut out maps
// rotational function on cutout box and interpolation

// template <class T> std::string display_repr(const Xmap<T> xmapclass) {
template <class T> std::string display_repr( const Xmap<T> xmap ) {
  if ( !xmap.is_null() ) {
    if ( strcmp( typeid( T ).name(), "f" ) == 0 )
      return "<clipper.Xmap_float: ASU grid " + xmap.grid_asu().format() + ">";
    else if ( strcmp( typeid( T ).name(), "d" ) == 0 )
      return "<clipper.Xmap_double: ASU grid " + xmap.grid_asu().format() + ">";
    else if ( strcmp( typeid( T ).name(), "i" ) == 0 )
      return "<clipper.Xmap_int: ASU grid " + xmap.grid_asu().format() + ">";
    else
      return "<clipper.Xmap>";
  } else
    return "<clipper.Xmap uninitialised.>";
}

// xmap to numpy array
template<typename T>
auto xmap_to_array(const Xmap<T> &xm, const bool &asu=false) { //}, nb::ndarray<T, nb::numpy, nb::ndim<3>> &arr) {

  //auto grid_asu = xm.grid_asu();
  std::array<int,3> orderxyz = { 0, 1, 2 };
  // bounds
  std::array<int, 3> gfms0, gfms1, g, gs;
  if (asu) {
    for ( int i = 0; i < 3; i++ ) {
      gfms0[orderxyz[i]] = xm.grid_asu().min()[i];
      gfms1[orderxyz[i]] = xm.grid_asu().max()[i];
      gs[i] = gfms1[orderxyz[i]] - gfms0[orderxyz[i]] + 1;
    }
  }
  else {
    for ( int i = 0; i < 3; i++ ) {
      gfms0[orderxyz[i]] = 0;
      gfms1[orderxyz[i]] = xm.grid_sampling()[i] - 1;
      gs[i] = gfms1[orderxyz[i]] - gfms0[orderxyz[i]] + 1;
    }
  } 
  // size_t array_size = grid_asu.nu() * grid_asu.nv() * grid_asu.nw();
  T* arr = new T[gs[0]*gs[1]*gs[2]];
  std::initializer_list<int64_t> strides={gs[1]*gs[2], gs[2], 1};
  Xmap_base::Map_reference_coord ix( xm );
  for ( g[0]=gfms0[0]; g[0] <= gfms1[0]; g[0]++ )
    for ( g[1]=gfms0[1]; g[1] <= gfms1[1]; g[1]++ )
      for ( g[2]=gfms0[2]; g[2] <= gfms1[2]; g[2]++ )
      {
        ix.set_coord(Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
        //const int ind = xm.index_of(Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
        // xm.index_of and ix.index does not work
        const int ind = ((g[orderxyz[0]]*gs[1]+g[orderxyz[1]]) * gs[2]) + g[2];
        arr[ind] = xm[ix];
      }
  nb::capsule owner(arr, [](void *p) noexcept {delete[] static_cast<T*>(p);});
  return nb::ndarray<T, nb::numpy, nb::ndim<3>, nb::device::cpu, nb::c_contig>(arr, {(size_t)gs[0], (size_t)gs[1], (size_t)gs[2]}, owner, strides);
}

// numpy array to xmap
template <typename T>
void numpy_to_xmap(const nb::ndarray<nb::numpy, T, nb::ndim<3>> & data, /*np_cpu_c_3darray<T> &data,*/ const Spacegroup &sg, const Cell &cell, Xmap<T> &xmap,
                   const std::array<int, 3> &axis_pos) {
  auto d = data.view();
  if ( d.ndim() != 3 ) throw std::runtime_error("Array must be 3-Dimension!");
  
  // set map params
  Grid_sampling grid_sam(int(d.shape(0)), int(d.shape(1)), int(d.shape(2)));
  xmap.init(sg, cell, grid_sam);
  const auto &asu_grid = xmap.grid_asu();
  auto grid_size = asu_grid.max()-asu_grid.min();
  int g[3], gfms[3], orderxyz[3];
  for ( int i = 0; i < 3; ++i ) {
    gfms[i] = grid_size[i];
    orderxyz[axis_pos[i]] = i; 
  }
  Xmap_base::Map_reference_coord ix(xmap);
  // copy map data
  for (g[0] = 0; g[0] <= gfms[0]; g[0]++)
    for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
      for (g[2] = 0; g[2] <= gfms[2]; g[2]++) {
        ix.set_coord(clipper::Coord_grid(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
        xmap[ix] = T(d(g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]]));
      }
}

nb::class_<Xmap_base> declare_xmap_base( nb::module_ &m ) {
  nb::class_<Xmap_base> xmap_base( m, "Xmap_base" );
  nb::enum_<Xmap_base::FFTtype>( xmap_base, "FFTtype", "FFT backend selection." )
      .value( "Default", Xmap_base::FFTtype::Default )
      .value( "Normal", Xmap_base::FFTtype::Normal )
      .value( "Sparse", Xmap_base::FFTtype::Sparse )
      .export_values();

  xmap_base.def( "is_null", &Xmap_base::is_null )
      .def_prop_ro( "cell", &Xmap_base::cell, "Get the cell." )
      .def_prop_ro( "spacegroup", &Xmap_base::spacegroup, "Get the spacegroup." )
      .def_prop_ro( "grid_sampling", &Xmap_base::grid_sampling, "Get cell grid." )
      .def_prop_ro( "grid_asu", &Xmap_base::grid_asu, "Get ASU grid." )
      .def( "coord_of", &Xmap_base::coord_of, nb::arg( "index" ), "Map coordinate from index." )
      .def( "index_of", &Xmap_base::index_of, nb::arg( "coord_grid" ), "Map index from coordinate." )
      .def( "to_map_unit", &Xmap_base::to_map_unit, nb::arg( "coord_grid" ),
            "Function to pick right cell repeat for any grid coord." )
      .def_prop_ro( "operator_orth_grid", &Xmap_base::operator_orth_grid,
                    "Return the orthogonal-to-grid coordinate "
                    "operator (translation is zero)." )
      .def_prop_ro( "operator_grid_orth", &Xmap_base::operator_grid_orth,
                    "Return the grid-to-orthogonal coordinate "
                    "operator (translation is zero)." )
      .def( "coord_orth", &Xmap_base::coord_orth, nb::arg( "coord_map" ), "Convert map coordinate to orthogonal." )
      .def( "coord_map", &Xmap_base::coord_map, nb::arg( "coord_orth" ), "Convert orthogonal coordinate to map." )
      .def( "in_map", ( bool ( Xmap_base::* )( const clipper::Coord_grid & ) const ) & Xmap_base::in_map,
            nb::arg( "Coord_grid" ) = nullptr,
            "(This method is for compatibility with NXmap - it always returns "
            "true)." )
      //.def("in_map_interp_nearest",
      //     (bool(NXmap_base::*)(const clipper::Coord_map &) const) &
      //         NXmap_base::in_map<Interp_nearest>,
      //     nb::arg("coord_map"), "Check if given coord available in the map?")
      //.def("in_map_interp_linear",
      //     (bool(NXmap_base::*)(const clipper::Coord_map &) const) &
      //         NXmap_base::in_map<Interp_linear>,
      //     nb::arg("coord_map"), "Check if given coord available in the map?")
      //.def("in_map_interp_cubic",
      //     (bool(NXmap_base::*)(const clipper::Coord_map &) const) &
      //         NXmap_base::in_map<Interp_cubic>,
      //     nb::arg("coord_map"), "Check if given coord available in the map?")
      .def( "multiplicity", &Xmap_base::multiplicity, nb::arg( "coord_grid" ), "Get multiplicity of map grid point." )
      .def( "first", &Xmap_base::first, "Return Map_reference_index for this map." )
      .def( "first_coord", &Xmap_base::first_coord, "Return a Map_reference_coord for this map." )
      // Tristan in pybind11 version of clipper xmap
      // - let xmap take control of creation of Map_reference... types
      .def(
          "map_reference_index",
          []( const Xmap_base &self, const Coord_grid &pos ) { return Xmap_base::Map_reference_index( self, pos ); },
          nb::arg( "pos" ), "Return Map_reference_index for this position in map." )
      .def(
          "map_reference_coord",
          []( const Xmap_base &self, const Coord_grid &pos ) { return Xmap_base::Map_reference_coord( self, pos ); },
          nb::arg( "pos" ), "Return Map_reference_coord for this position in map." )
      .doc() = "Xmap_base: base for crystallographic map class.\n"
               "The crystallographic map class stores a map of arbitrary data "
               "type. Its main difference from a 3-d array is that the data extent "
               "appears to be infinite, and yet internally only a unique ASU is "
               "stored. Iterators provide efficient access to data. "
               "This base contains everything except the data, which is templated "
               "in the derived type Xmap<T>";
  return xmap_base;
}

void declare_xmap_reference_base( nb::module_ &m, nb::class_<Xmap_base> &xmap_base ) {
  using MRB = Xmap_base::Map_reference_base;
  nb::class_<MRB>( xmap_base, "Xmap_reference_base" )
      .def_prop_ro( "base_xmap", &MRB::base_xmap, "Return parent Xmap." )
      .def_prop_ro( "index", &MRB::index, "Get the index into the map data array." )
      .def( "last", &MRB::last, "Check for end of map." )
      .doc() = "Map reference base class.\n"
               "This is a reference to an Map. It forms a base class for "
               "index-like and coordinate-like Map references. If you write a "
               "method which will work with either, then specify this instead of "
               "either of the derived classed.";
}

void declare_xmap_reference_index( nb::module_ &m, nb::class_<Xmap_base> &xmap_base ) {
  using MRI = Xmap_base::Map_reference_index;
  nb::class_<MRI, Xmap_base::Map_reference_base>( xmap_base, "Xmap_reference_index" )
      // Let xmap take control of creation ...
      //.def(nb::init<>(), "Null constructor.")
      //.def(nb::init<const Xmap_base &>())
      //.def(nb::init<const Xmap_base &, const Coord_grid &>())
      .def_prop_rw(
          "coord", &MRI::coord, []( MRI &self, const Coord_grid &pos ) -> void { self.set_coord( pos ); },
          nb::for_getter( "Get current grid coordinate." ),
          nb::for_setter( "Set current grid coordinate - optimised for nearby coords" ) )
      .def( "coord_orth", &MRI::coord_orth, "Get current value of orthogonal coordinate." )
      .def(
          "next", []( MRI &self ) -> void { self.next(); }, "Simpe increment." )
      .def( "index_offset", &MRI::index_offset, nb::arg( "du" ), nb::arg( "dv" ), nb::arg( "dw" ),
            "Index of neighbouring point." )
      .doc() = "Xmap reference with index-like behaviour.\n"
               "This is a reference to a map coordinate. It behaves like a "
               "simple index into the map, but can be easily converted into a "
               "coordinate as and when required. It also implements methods for "
               "iterating through the unique portion of a map. It is very "
               "compact, but coord() involves some overhead and loses any "
               "information concerning symmetry equivelents.";
}

void declare_xmap_reference_coord( nb::module_ &m, nb::class_<Xmap_base> &xmap_base ) {
  using MRC = Xmap_base::Map_reference_coord;
  nb::class_<MRC, Xmap_base::Map_reference_base>( xmap_base, "Xmap_reference_coord" )
      // Let xmap take control of creation ...
      //.def(nb::init<>(), "Null constructor.")
      //.def(nb::init<const Xmap_base &>())
      //.def(nb::init<const Xmap_base &, const Coord_grid &>())
      .def_prop_rw(
          "coord", &MRC::coord, []( MRC &self, const Coord_grid &pos ) -> void { self.set_coord( pos ); },
          nb::for_getter( "Get current value of coordinate." ), nb::for_setter( "Set current value of coordinate." ) )
      .def( "coord_orth", &MRC::coord_orth )
      .def_prop_ro( "sym", &MRC::sym, "Get current symmetry operator." )
      .def(
          "next", []( MRC &self ) -> void { self.next(); },
          "Simple increment. Use of this function resets the stored "
          "coordinat and sym." )
      .def(
          "next_u", []( MRC &self ) -> void { self.next_u(); }, "Increment u." )
      .def(
          "next_v", []( MRC &self ) -> void { self.next_v(); }, "Increment v." )
      .def(
          "next_w", []( MRC &self ) -> void { self.next_w(); }, "Increment w." )
      .def(
          "prev_u", []( MRC &self ) -> void { self.prev_u(); }, "Decrement u." )
      .def(
          "prev_v", []( MRC &self ) -> void { self.prev_v(); }, "Decrement v." )
      .def(
          "prev_w", []( MRC &self ) -> void { self.prev_w(); }, "Decrement w." )
      .doc() = "Xmap reference with coord-like behaviour.\n"
               "This is a reference to an HKL. It behaves like an HKL, but "
               "also stores the index of the corresponding reflection in the "
               "reflection list, if it exists, and the symmetry and friedel "
               "operators required to get there. ";
}

template <class Derived, class T, class H> void apply_xmap_fft_methods( nb::class_<Derived, Xmap_base> &pyclass ) {
  pyclass
      .def(
          "fft_from",
          []( Derived &self, const H &phidata, const Xmap_base::FFTtype type ) { self.fft_from( phidata ); },
          nb::arg( "fphidata" ), nb::arg( "ffttype" ) = Xmap_base::FFTtype::Default,
          "FFT from reflection list to map." )
      .def(
          "fft_to", []( const Derived &self, H &phidata, const Xmap_base::FFTtype type ) { self.fft_to( phidata ); },
          nb::arg( "fphidata" ), nb::arg( "ffttype" ) = Xmap_base::FFTtype::Default,
          "FFT from map to reflection list." );
}

template <class Derived, class T> void apply_xmap_interpolation_methods( nb::class_<Derived, Xmap_base> &pyclass ) {
  pyclass
      .def( "interp_linear_frac",
            []( const Derived &self, const Coord_frac &pos ) { return self.template interp<Interp_linear>( pos ); } )
      .def( "interp_cubic_frac",
            []( const Derived &self, const Coord_frac &pos ) { return self.template interp<Interp_cubic>( pos ); } )
      .def( "interp_linear_orth",
            []( const Derived &self, const Coord_orth &xyz ) {
              return self.template interp<Interp_linear>( xyz.coord_frac( self.cell() ) );
            } )
      .def( "interp_cubic_orth",
            []( const Derived &self, const Coord_orth &xyz ) {
              return self.template interp<Interp_cubic>( xyz.coord_frac( self.cell() ) );
            } )
      .def( "interp_cubic_grad_frac",
            []( const Derived &self, const Coord_frac &pos ) -> nb::tuple {
              T val;
              Grad_frac<T> grad;
              self.template interp_grad<Interp_cubic>( pos, val, grad );
              return nb::make_tuple( val, grad );
            } )
      .def( "interp_cubic_curv_frac", []( const Derived &self, const Coord_frac &pos ) -> nb::tuple {
        T val;
        Grad_frac<T> grad;
        Curv_frac<T> curv;
        self.template interp_curv<Interp_cubic>( pos, val, grad, curv );
        return nb::make_tuple( val, grad, curv );
      } );
} // apply_xmap_interpolation_methods

template <class T> void declare_xmap( nb::module_ &m, const std::string &name ) {
  using MRI = Xmap_base::Map_reference_index;
  using MRC = Xmap_base::Map_reference_coord;
  using XMClass = Xmap<T>;
  std::string PyClass = std::string( "Xmap_" ) + name;
  nb::class_<XMClass, Xmap_base> xmap( m, PyClass.c_str() );
  // nb::buffer_protocol());
  xmap.def( nb::init<>(), "Null constructor" )
      .def( nb::init<const Spacegroup &, const Cell &, const Grid_sampling &>(), nb::arg( "spacegroup" ),
            nb::arg( "cell" ), nb::arg( "grid_sampling" ), "Constructor: from spacegroup, cell, and grid." )
      .def(
          "__init__", []( XMClass *self, const XMClass &xmap ) { new ( self ) XMClass( xmap ); }, nb::arg( "xmap" ),
          "Copy constructor." )
      //.def( "__init__", []( XMClass *xmap, const nb::ndarray<nb::numpy, T, nb::ndim<3>> &data, const Cell &cell,
      .def( "__init__", []( XMClass *xmap, const nb::ndarray<nb::numpy, T, nb::ndim<3>> &data /*const np_cpu_c_3darray<T> &data*/, const Spacegroup &sg, 
                            const Cell &cell) {//}, const std::array<int, 3> &axis_pos) {
        new ( xmap ) XMClass();
        std::array<int, 3> axis_pos = {0,1,2};
        numpy_to_xmap<T>(data, sg, cell, *xmap, axis_pos);
      }, nb::arg("array").noconvert(), nb::arg("spacegroup"), nb::arg("cell")) //, nb::arg("axis_pos"))
      .def( "init", &XMClass::init, nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "grid_sampling" ),
            "Initialise from spacegroup, cell, and grid." )
      //.def_prop_ro( "asu_array", [](XMClass &xm) {
      //  
      //})
      //.def( "find_sym", [](XMClass &self, const Coord_grid &cg) {
      //  int ind, sym;
      //  self.find_sym_1(cg, ind, sym);
      //  return nb::make_tuple(ind, sym);
      //})
      .def_prop_ro( "asu_array", []( XMClass &xm) { return xmap_to_array(xm, true); }, nb::rv_policy::automatic)
      .def_prop_ro( "array", []( XMClass &xm) { return xmap_to_array(xm); }, nb::rv_policy::automatic)
      .def( "__repr__", []( const XMClass &self ) { return display_repr<T>( self ); } )
      .def( "__str__", []( const XMClass &self ) { return display_repr<T>( self ); } )
      .def( "__getitem__", []( const XMClass &self, const MRI &ix ) { return self[ix]; } )
      .def( "__setitem__", []( XMClass &self, const MRI &ix, const T &val ) { self[ix] = val; } )
      .def( "__getitem__", []( const XMClass &self, const MRC &ix ) { return self[ix]; } )
      .def( "__setitem__", []( XMClass &self, const MRC &ix, const T &val ) { self[ix] = val; } )
      .def(
          "get_data", []( const XMClass &self, const Coord_grid &pos ) { return self.get_data( pos ); },
          nb::arg( "pos" ), "Get data." )
      .def(
          "set_data", []( XMClass &self, const Coord_grid &pos, const T &val ) { self.set_data( pos, val ); },
          nb::arg( "pos" ), nb::arg( "val" ), "Set data." )
      .def(
          "fill_map_with", []( XMClass &self, const T &val ) { self = val; }, nb::arg( "val" ),
          "Fill map with given value." )
      .def( nb::self += nb::self )
      .def( nb::self -= nb::self )
      // because python's gemmi map object's data is float type.
      .def(
          "import_from_gemmi",
          []( XMClass& self, const gemmi::Ccp4<float>& mapobj ) {
            GEMMI::import_xmap(self, mapobj);
          },
          nb::arg( "target" ), "Import Xmap from gemmi.Ccp4Map." )
      .def(
          "export_to_gemmi",
          []( const XMClass& self, gemmi::Ccp4<float>& mapobj ) {
            GEMMI::export_xmap(self, mapobj);
          },
          nb::arg( "target" ), "Export Xmap to gemmi.Ccp4Map." );
  // interpolators
  apply_xmap_interpolation_methods<XMClass, T>( xmap );
  // ffts
  apply_xmap_fft_methods<XMClass, T, HKL_data<clipper::data32::F_phi>>( xmap );
  apply_xmap_fft_methods<XMClass, T, HKL_data<clipper::data64::F_phi>>( xmap );
}

void add_xmap( nb::module_ &m ) {
  auto xmapbase = declare_xmap_base( m );
  declare_xmap_reference_base( m, xmapbase );
  declare_xmap_reference_index( m, xmapbase );
  declare_xmap_reference_coord( m, xmapbase );
  declare_xmap<int>( m, "int" );
  declare_xmap<ftype32>( m, "float" );
  declare_xmap<ftype64>( m, "double" );
}