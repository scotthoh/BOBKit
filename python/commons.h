#pragma once

#include <clipper/core/clipper_types.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <gemmi/math.hpp>
#include <memory>

namespace nb = nanobind;
constexpr auto rv_ri = nb::rv_policy::reference_internal;

// declare bindings
void add_clipper_util( nb::module_ &m );  //1
void add_clipper_types( nb::module_ &m );  //1
void add_thread_base( nb::module_ &m );  //1
void add_coords( nb::module_ &m );  //1
void add_atomlist( nb::module_ &m );  //1
void add_cell( nb::module_ &m );  //1
void add_clipper_memory( nb::module_ &m );  //1
void add_ramachandran( nb::module_ &m );  //1
void add_hklinfo( nb::module_ &m );  //1
void init_hkl_data( nb::module_ &m );  //1
void init_hkl_datatypes( nb::module_ &m );  //1
void init_symop( nb::module_ &m );  //1
void init_spacegroup( nb::module_ &m, nb::module_ &mdata );  //1
void init_clipper_stats( nb::module_ &m );  //1
void init_resol_fn( nb::module_ &m );  //1
void add_atomsf( nb::module_ &m );  //1
void init_sfcalc( nb::module_ &m );
void init_sfcalc_obs( nb::module_ &m );
void init_sfweight( nb::module_ &m );
void add_clipper_tests( nb::module_ &m, nb::module_ &mdata );  //1
void add_cifdata_io( nb::module_ &m );

void add_rotation( nb::module_ &m );  //1
void add_xmap( nb::module_ &m );  //1
void add_nxmap( nb::module_ &m );  //1
void add_nxmap_operator( nb::module_ &m );  //1
void add_mapfilters( nb::module_ &m );  //1
void add_map_utils( nb::module_ &m );  //1

void init_minimol( nb::module_ &m, nb::module_ &mm );  //1
void init_minimol_seq( nb::module_ &m );  //1
void add_matomindex( nb::module_ &m );  //1
void init_sfscale( nb::module_ &m );  //1

void add_fffear( nb::module_ &m );  //1
void add_edcalc( nb::module_ &m );  //1

void init_gemmi_structure( nb::module_ &m );
void add_minimol_io_gemmi( nb::module_ &m );  

void add_ca_build( nb::module_ &m );
void add_ca_correct( nb::module_ &m );
void add_ca_filter( nb::module_ &m );
void add_ca_find( nb::module_ &m );
void add_ca_grow( nb::module_ &m );
void add_ca_join( nb::module_ &m );
void add_knownstructure( nb::module_ &m );
void add_buccaneer_lib( nb::module_ &m );
void add_ca_link( nb::module_ &m );
void add_ca_merge( nb::module_ &m );
void add_ca_ncsbuild( nb::module_ &m );
void add_ca_prep( nb::module_ &m );
void add_buccaneer_prot( nb::module_ &m );
void add_proteindb( nb::module_ &m );
void add_protein_loop( nb::module_ &m );
void add_ca_prune( nb::module_ &m );
void add_ca_sequence( nb::module_ &m );
void add_simplex_lib( nb::module_ &m );
void add_map_simulate( nb::module_ &m );
void add_model_tidy( nb::module_ &m );
void add_buccaneer_util( nb::module_ &m );
void add_utils( nb::module_ &m );

namespace clipper {
  class MiniMol;
  class MModel;
}
//void add_convert_structure2mmol(nb::module_ &m, nb::class_<clipper::MiniMol, clipper::MModel>& mmol);
//void declare_gemmi_atom( nb::module_ &m );

namespace nanobind {
namespace detail {
template <> struct type_caster<clipper::String> {
  NB_TYPE_CASTER( clipper::String, const_name( "String" ) );
  // Python -> C++
  bool from_python( handle src, uint8_t, cleanup_list * ) noexcept {
    PyObject *source = src.ptr();
    PyObject *tmp = PyUnicode_AsUTF8String( source );
    if ( !tmp ) {
      PyErr_Clear();
      return false;
    }
    char *buffer;
    ssize_t length;
    if ( PyBytes_AsStringAndSize( tmp, &buffer, &length ) )
      return false;
    value = clipper::String( std::string( buffer, ( size_t )length ) );
    Py_DECREF( tmp );
    return true;
  }
  // C++ -> Python
  static handle from_cpp( const clipper::String &src, rv_policy policy, cleanup_list * ) noexcept {
    return PyUnicode_FromString( src.c_str() );
  }
};

} // namespace detail
} // namespace nanobind

// return Vec3 from array
template <typename T> std::unique_ptr< clipper::Vec3<T> > array_to_vec3( std::array<T, 3> vals ) {
  return std::unique_ptr<clipper::Vec3<T>>( new clipper::Vec3<T>( T( vals[0] ), T( vals[1] ), T( vals[2] ) ) );
}

// return Mat33 from array
template <typename T> std::unique_ptr< clipper::Mat33<T> > array_to_mat33( std::array<std::array<T, 3>, 3> &val ) {

  return std::unique_ptr< clipper::Mat33<T> >(
      new clipper::Mat33<T>( T( val[0][0] ), T( val[0][1] ), T( val[0][2] ), T( val[1][0] ), T( val[1][1] ),
                             T( val[1][2] ), T( val[2][0] ), T( val[2][1] ), T( val[2][2] ) ) );
}

// check index from python and return C++ index
inline int normalise_index( int index, const size_t container_size ) {
  if ( index < 0 )
    index += ( int )container_size;
  if ( ( size_t )index >= container_size )
    throw nb::index_error( "Index out of bounds" );
  return index;
}

template <typename Item>
void delitem_at_index(std::vector<Item> &items, ssize_t index) {
  items.erase(items.begin() + index);
}

template <typename T> void delitem_at_index(T &container, ssize_t index) {
  container.erase(container.begin() + index);
}

template <typename T>
void delitem_range(T &container, ssize_t start, ssize_t end) {
  container.erase(container.begin() + start, container.begin() + end);
}

template <typename T, typename C>
C &add_item(T &container, C child, int pos)
{
  if ((size_t)pos > container.size()) // true for negative pos
    pos = (int)container.size();
  return *container.insert(container.begin() + pos, std::move(child));
}

inline void to_pystream( std::string &msg, const nb::object &pystream ) {
  if ( nb::hasattr( pystream, "write" ) && nb::hasattr( pystream, "flush" ) ) {
    pystream.attr( "write" )( msg.c_str() );
    pystream.attr( "flush" )();
  }
}
