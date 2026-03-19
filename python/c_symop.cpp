// Nanobind bindings for clipper symop
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>

// #include <clipper/core/clipper_types.h>
// #include <clipper/core/symop.h>
// #include "c_symop.h"

using namespace clipper;

template <typename T, size_t N> using mat44 = nb::ndarray<T, nb::shape<-1, N>, nb::device::cpu, nb::c_contig>;

void declare_rtopfrac( nb::module_ &m ) {
  nb::class_<RTop_frac, RTop<>>( m, "RTop_frac" )
      .def( nb::init<>() )
      .def( nb::init<const RTop<> &>(), nb::arg( "rtop" ), "Constructor: copy/convert" )
      .def( nb::init<const Mat33<> &>(), nb::arg( "rot" ), "Constructor from rotation. " )
      .def( nb::init<const String &>(), nb::arg( "strop" ), "Constructor from string description." )
      .def( nb::init<const Mat33<> &, const Vec3<> &>(), nb::arg( "rot" ), nb::arg( "trn" ),
            "Constructor from rotation and translation." )
      .def(
          "__init__",
          []( RTop_frac *rtf, std::array<std::array<double, 3>, 3> &rot, std::array<double, 3> &trn ) {
            auto rotation = array_to_mat33( rot );
            auto translation = array_to_vec3( trn );
            new ( rtf ) RTop_frac( *rotation, *translation );
          },
          nb::arg( "rot" ), nb::arg( "trn" ), "Constructor from list/array of rotation and translation." )
      .def( "rtop_orth", &RTop_frac::rtop_orth, nb::arg( "cell" ), "Fractional-orthogonal conversion." )
      .def( "inverse", &RTop_frac::inverse, "Inverse operator." )
      .def_static( "identity", &RTop_frac::identity, "Return identity operator." )
      .def_static( "null", &RTop_frac::null, "Return null (uninitialised) operator." )
      .def( "__repr__", []( const RTop_frac &self ) { return "<clipper.RTop_frac: Fractional operator class.>"; } )
      .doc() = "Fractional operator class.\nThis class is used "
               "for any RT-operator which operates on fractional "
               "coordinates. For a full list of methods, see clipper::RTop.";
}

void declare_symop( nb::module_ &m ) {
  nb::class_<Symop, RTop_frac>( m, "Symop" )
      .def( nb::init<>() )
      .def( nb::init<const RTop<> &>(), nb::arg( "rtop" ), "Constructor from RTop." )
      // 4x4 array
      .def( "__init__",
            []( Symop *self, const mat44<double, 4> &values ) {
              auto v = values.view();
              if ( v.shape( 0 ) != 4 )
                throw std::length_error( "Array must be 4x4" );
              // is there a better way than to use reinterpret_cast?
              new ( self ) Symop( *reinterpret_cast<double ( * )[4][4]>( v.data() ) );
            } )
      //.def(nb::init<const std::array<std::array<double,4>,4>>(),  nb::arg("mat"))
      //.def("__init__", [](const std:array<std::array<double,4>,4>& a) {
      //       check_array_shape(a, {4, 4}, true);
      //       auto buf = a.request();
      //       return std::unique_ptr<Symop>(
      //           new Symop(*reinterpret_cast<ftype(*)[4][4]>(buf.ptr)));
      //     },
      //     nb::arg("rtop"), "Constructor from 4x4 matrix (list/array).")
      .def( "format", &Symop::format, "Return formatted string representation." )
      .def( "__str__", &Symop::format )
      .def( "__repr__",
            []( const Symop &self ) { return "<clipper.Symop: Crystallographic symmetry operator class.>"; } )
      .doc() = "Crystallographic symmetry operator.\nThis is identical to "
               "a fractional RTop, but has its own class since not all "
               "fractional  RTops are symops. For a full list of methods, "
               "see clipper::RTop and clipper::RTop_frac.";
}

void declare_isymop( nb::module_ &m ) {
  nb::class_<Isymop, RTop<int>>( m, "Isymop" )
      .def( nb::init<>() )
      .def( nb::init<const RTop<int> &>(), nb::arg( "rtop" ), "Constructor from RTop." )
      .def( nb::init<const Symop &, const Grid &>(), nb::arg( "symop" ), nb::arg( "grid" ),
            "Constructor from symop and grid." )
      .def( "__repr__", []( const Isymop &self ) { return "<clipper.Isymop: Integerised symmetry matrix class.>"; } )
      .doc() = "Integerised symmetry matrix.\nThis is used for "
               "optimised calculations in real and reciprocal space.";
}

void declare_symop_code( nb::module_ &m ) {
  nb::class_<Symop_code>( m, "Symop_code" )
      .def( nb::init<>() )
      .def( nb::init<const int &>(), nb::arg( "code" ), "Constructor from int." )
      .def( nb::init<const Symop &>(), nb::arg( "symop" ), "Constructor from Symop." )
      .def( nb::init<const Isymop &>(), nb::arg( "isymop" ), "Contructor from Isymop." )
      .def( "init", &Symop_code::init, nb::arg( "isymop" ), "Initialiser from Isymop." )
      .def( "code_rot", &Symop_code::code_rot, "Return code for rotational part." )
      .def( "code_trn", &Symop_code::code_trn, "Return code for translational part." )
      .def( "symop", &Symop_code::symop, "Convert to Symop." )
      .def( "isymop", &Symop_code::isymop, "Convert to integerised symop." )
      .def_static( "identity", &Symop_code::identity, "Identity code." )
      .def(
          "__int__", []( const Symop_code &self ) { return int( self ); }, "Convert to integer." )
      .def( "__repr__",
            []( const Symop_code &self ) {
              return "<clipper.Symop_code: Compressed encoded symmetry operator "
                     "class.>";
            } )
      .doc() = "Compressed encoded symmetry operator.\n"
               "This is a compresses representation of a crystallographic "
               "symmetry operator, stored as a single 32-bit integer. It may be "
               "converted to or from a symop or an int and compared, sorted, etc. "
               "The following guarantees are made concerning the code: "
               "- The identity operator has a code of zero.\n"
               "- Operators with non-identity rotations will have higher codes "
               "than operators with identity rotations, for the same translation.\n"
               "- Operators with non-zero translations will have higher codes "
               "than operators with zero translations, for the same rotation.";
}

void init_symop( nb::module_ &m ) {
  declare_rtopfrac( m );
  declare_symop( m );
  declare_isymop( m );
  declare_symop_code( m );
}
