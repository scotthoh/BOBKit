// Nanobind bindings for clipper types
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "arrays.h"
#include "commons.h"
#include <clipper/core/clipper_types.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/string.h>

// #include <clipper/core/clipper_types.h>
// change printVec to printVec with template<class T>
// change format() for mat33, rtop to deal with int

using namespace clipper;

//! Return formatted string representation for Vec3<>
template <class T> std::string printVec( Vec3<T> v ) {
  return "(" + String( double( v[0] ), 10, 4 ) + "," + String( double( v[1] ), 10, 4 ) + "," +
         String( double( v[2] ), 10, 4 ) + ")";
}

//! Return formatted string representation for Mat33<>
template <class T> std::string printMat( Mat33<T> m ) {
  return "|" + String( double( m( 0, 0 ) ), 10, 4 ) + "," + String( double( m( 0, 1 ) ), 10, 4 ) + "," +
         String( double( m( 0, 2 ) ), 10, 4 ) + "|\n|" + String( double( m( 1, 0 ) ), 10, 4 ) + "," +
         String( double( m( 1, 1 ) ), 10, 4 ) + "," + String( double( m( 1, 2 ) ), 10, 4 ) + "|\n|" +
         String( double( m( 2, 0 ) ), 10, 4 ) + "," + String( double( m( 2, 1 ) ), 10, 4 ) + "," +
         String( double( m( 2, 2 ) ), 10, 4 ) + "|";
}

//! Return formatted string representation for Mat33sym<>
template <class T> std::string printSMat( Mat33sym<T> m ) {
  return "|" + String( double( m.mat00() ), 10, 4 ) + "," + String( double( m.mat01() ), 10, 4 ) + "," +
         String( double( m.mat02() ), 10, 4 ) + "|\n|" + String( double( m.mat01() ), 10, 4 ) + "," +
         String( double( m.mat11() ), 10, 4 ) + "," + String( double( m.mat12() ), 10, 4 ) + "|\n|" +
         String( double( m.mat02() ), 10, 4 ) + "," + String( double( m.mat12() ), 10, 4 ) + "," +
         String( double( m.mat22() ), 10, 4 ) + "|";
}

template <class T> auto mat33_to_array( Mat33<T> &self ) {
  return nb::ndarray<nb::numpy, T, nb::shape<3, 3>, nb::c_contig>( &self( 0, 0 ), { 3, 3 }, nb::handle() );
}

void check_index( const int i, const int j ) {
  if ( i < 0 || i > 2 )
    throw nb::index_error( "Index out of bounds" );
  if ( j < 0 || j > 2 )
    throw nb::index_error( "Index out of bounds" );
}

template <class T> void declare_vec3( nb::module_ &m, const std::string &name ) {
  using V3 = Vec3<T>;
  nb::class_<V3> vec3( m, ( "Vec3_" + name ).c_str(), "3-vector Class" );
  vec3.def( nb::init<>(), "Null constructor." )
      .def( nb::init<const T &, const T &, const T &>(), nb::arg( "x" ), nb::arg( "y" ), nb::arg( "z" ),
            "Constructor from individual values." )
      // init from array or list
      .def(
          "__init__", []( V3 *v, std::array<T, 3> &a ) { new ( v ) Vec3<T>( a[0], a[1], a[2] ); }, nb::arg( "a" ),
          "Constructor from list/array." )
      .def( nb::init<const Vec3<float> &>(), nb::arg( "v" ), "Constructor: copy/convert." )
      .def( nb::init<const Vec3<double> &>(), nb::arg( "v" ), "Constructor: copy/convert." )
      .def( nb::init<const Vec3<int> &>(), nb::arg( "v" ), "Constructor: copy/convert." )
      //.def_buffer([](V3 &self) -> py::buffer_info {
      //  return py::buffer_info(&self[0], sizeof(T),
      //                         py::format_descriptor<T>::format(), 1, {3},
      //                         {sizeof(T)});
      //})
      .def( "equals", &V3::equals, nb::arg( "vec" ), nb::arg( "tol" ), "Test equality." )
      .def( "__getitem__",
            []( const V3 &self, const int i ) {
              check_index(i, 0);
              return self[i];
            } )
      .def( "__setitem__",
            []( V3 &self, const int &i, const T &val ) {
              check_index(i, 0);
              self[i] = val;
            } )
      //.def(
      //    "as_array",
      //    [](const V3 &self) {
      //      return make_array_1d<V3, T>(self, 3);
      //    },
      //    "Return as array.")
      //.def(
      //    "from_array",
      //    [](V3 &self, py::array_t<T> values) {
      //      fill_array_1d<V3, T>(self, 3, values);
      //    },
      //    "Import values from array.")
      .def( "unit", &V3::unit, "Return unit vector with same direction as this vector." )
      .def_static( "zero", &V3::zero, "Return zero vector." )
      .def_static( "null", &V3::null, "Return null vector (only valid for floating point types)." )
      // ugly hack to bypass error due to non-ftype argument
      // passed to clipper::Util::is_nan(ftype)
      // or we can edit clipper::Util to add is_nan(int)???
      .def(
          "is_null",
          []( const V3 &self ) {
            Vec3<double> vtmp( self );
            return vtmp.is_null();
          },
          "Test for null vector (only valid for floating point types)." )
      .def_static( "dot", &V3::dot, nb::arg( "v1" ), nb::arg( "v2" ), "Vector dot product. (Equivalent to v1 * v2)" )
      .def_static( "cross", &V3::cross, nb::arg( "v1" ), nb::arg( "v2" ), "Vector cross product." )
      .def(
          "__iter__",
          []( const V3 &self ) {
            return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<V3>(), "iterator", &self[0],
                                                                         &self[3] );
          },
          nb::keep_alive<0, 1>() )
      // don't understand why copy constructor requires len() in python side
      // thus this is defined to avoid error when calling constructor
      // Vec3_T(Vec3_T)
      .def( "__len__", []( const V3 &self ) { return 3; } )
      .def( "__str__", []( const V3 &self ) { return printVec( self ); } )
      .def(
          "format", []( const V3 &self ) { return printVec( self ); }, "Return formatted string representation." )
      .def( "__repr__", [=]( const V3 &self ) { return "<clipper.Vec3_" + name + " " + printVec( self ) + ">"; } )
      .def(
          "__iadd__", []( V3 &self, const V3 &other ) { return self += other; }, nb::is_operator(),
          nb::rv_policy::none )
      .def(
          "__isub__", []( V3 &self, const V3 &other ) { return self -= other; }, nb::is_operator(),
          nb::rv_policy::none )
      .def(
          "__eq__", []( const V3 &self, const V3 &other ) { return self == other; }, nb::is_operator() )
      .def(
          "__ne__", []( const V3 &self, const V3 &other ) { return self != other; }, nb::is_operator() )
      .def(
          "__neg__", []( const V3 &self ) { return -self; }, nb::is_operator() )
      .def(
          "__add__", []( const V3 &self, const V3 &other ) { return self + other; }, nb::is_operator() )
      .def(
          "__sub__", []( const V3 &self, const V3 &v ) { return self - v; }, nb::is_operator() )
      .def(
          "__mul__", []( const V3 &self, const T &s ) { return self * s; }, nb::is_operator() )
      .def(
          "__rmul__", []( const V3 &self, const T &s ) { return s * self; }, nb::is_operator() )
      .def(
          "__mul__", []( const V3 &self, const V3 &other ) { return Vec3<T>::dot( self, other ); }, nb::is_operator() );
}

template <class T> void declare_mat33( nb::module_ m, const std::string &name ) {
  using M33 = Mat33<T>;
  nb::class_<M33> mat33( m, ( "Mat33_" + name ).c_str(), "3x3 matrix class." );
  mat33.def( nb::init<>() )
      .def(
          nb::init<const T &, const T &, const T &, const T &, const T &, const T &, const T &, const T &, const T &>(),
          "Constructor from individual values." )
      .def( nb::init<const Mat33<float> &>(), "Constructor: copy/convert." )
      .def( nb::init<const Mat33<double> &>(), "Constructor: copy/convert." )
      .def( nb::init<const Mat33sym<float> &>(), "Constructor: copy/convert from symmetric matrix." )
      .def( nb::init<const Mat33sym<double> &>(), "Constructor: copy/convert from symmetric matrix." )
      .def(
          "__init__",
          []( M33 *m, const cpu_c_2darray<T, 3, 3> &a ) {
            new ( m ) Mat33<T>( a( 0, 0 ), a( 0, 1 ), a( 0, 2 ), a( 1, 0 ), a( 1, 1 ), a( 1, 2 ), a( 2, 0 ), a( 2, 1 ),
                                a( 2, 2 ) );
          },
          nb::arg( "array" ) )
      .def_prop_ro( "array", &mat33_to_array<T>, nb::rv_policy::reference_internal )
      //.def("__array__", [](nb:handle_t<Mat33>& h, nb::handle dtype, nb::handle copy) {
      //  return hand
      //})
      //.def(nb::init([](py::array_t<T> vals) { return numpy_to_mat33(vals); }),
      //     nb::arg("a"), "Constructor from 3x3 numpy array.")
      //.def_buffer([](M33 &self) -> py::buffer_info {
      //  return py::buffer_info(&self(0, 0), {3, 3}, {sizeof(T) * 3, sizeof(T)});
      //})
      .def( "det", &M33::det, "Determinant." )
      .def( "inverse", &M33::inverse, "Inverse." )
      .def( "transpose", &M33::transpose, "Transpose." )
      .def( "equals", &M33::equals, nb::arg( "m" ), nb::arg( "tol" ), "Test equality." )
      .def(
          "get",
          []( M33 &self, const int i, const int j ) {
            check_index( i, j );
            return self( i, j );
          },
          "Get element." )
      .def(
          "set",
          []( M33 &self, const int i, const int j, const T val ) {
            check_index( i, j );
            self( i, j ) = val;
          },
          "Set element." )
      //.def("format", &M33::format)
      .def(
          "format", []( const M33 &self ) { return printMat( self ); }, "Return formatted string representation." )
      .def( "__str__", []( const M33 &self ) { return printMat( self ); } )
      .def( "__repr__", [=]( const M33 &self ) { return "<clipper.Matt33_" + name + " class.>"; } )
      .def_static( "identity", &M33::identity, "Return identity matrix." )
      .def_static( "null", &M33::null, "Return null matrix (only valid for floating point types)." )
      // ugly hack to bypass error due to non-ftype argument
      // passed to clipper::Util::is_nan(ftype)
      // or we can edit clipper::Util to add is_nan(int)???
      .def(
          "is_null",
          []( const M33 &self ) {
            Mat33<double> mtmp( self );
            return mtmp.is_null();
          },
          "Test for null vector (only valid for floating point types)." )
      // operator Matrix-vector, assumes column vector
      .def(
          "__mul__", []( const M33 &self, const Vec3<T> &v ) { return self * v; }, nb::is_operator() )
      .def(
          "__rmul__", []( const Vec3<T> &v, const M33 &self ) { return v * self; }, nb::is_operator() )
      .def(
          "__mul__", []( const M33 &self, const M33 &other ) { return self * other; }, nb::is_operator() )
      .def(
          "__add__", []( const M33 &self, const M33 &other ) { return self + other; }, nb::is_operator() )
      .def(
          "__sub__", []( const M33 &self, const M33 &other ) { return self + ( -other ); }, nb::is_operator() )
      .def( "__neg__", []( const M33 &self ) { return -self; }, nb::is_operator() );
}

template <class T> void declare_mat33sym( nb::module_ m, const std::string &name ) {
  using SMat = Mat33sym<T>;
  nb::class_<SMat> smat( m, ( "Mat33sym_" + name ).c_str(), "Compressed form for 3x3 symmetric matrix class." );
  smat.def( nb::init<>() )
      .def( nb::init<const Mat33<float> &>(), "Constructor from Mat33 (does not check for symmetry)." )
      .def( nb::init<const Mat33<double> &>(), "Constructor from Mat33 (does not check for symmetry)." )
      .def( nb::init<const Mat33sym<float> &>(), "Constructor from Mat33sym." )
      .def( nb::init<const Mat33sym<double> &>(), "Constructor from Mat33sym." )
      .def( nb::init<const T &, const T &, const T &, const T &, const T &, const T &>(), nb::arg( "c00" ),
            nb::arg( "c11" ), nb::arg( "c22" ), nb::arg( "c01" ), nb::arg( "c02" ), nb::arg( "c12" ),
            "Constructor from coefficients." )
      .def(
          "format", []( const SMat &self ) { return printSMat( self ); }, "Return formatted string representation." )
      .def( "__str__", []( const SMat &self ) { return printSMat( self ); } )
      .def( "__repr__", [=]( const SMat &self ) { return "<clipper.Mat33sym_" + name + " class.>"; } )
      .def_static( "identity", &SMat::identity, "Return identity matrix." )
      .def_static( "null", &SMat::null, "Return null matrix (only valid for floating point types)." )
      // ugly hack to bypass error due to non-ftype argument
      // passed to clipper::Util::is_nan(ftype)
      .def(
          "is_null",
          []( const SMat &self ) {
            Mat33sym<double> mtmp( self );
            return mtmp.is_null();
          },
          "Test for null vector (only valid for floating point types)." )
      .def( "quad_form", &SMat::quad_form, nb::arg( "v" ), "Return quadratic form with vector." )
      .def( "det", &SMat::det, "Determinant." )
      .def( "sqrt", &SMat::sqrt, "Square root." )
      .def( "inverse", &SMat::inverse, "Inverse." )
      .def_prop_ro( "m00", &SMat::mat00, "Return element (0,0)." )
      .def_prop_ro( "m11", &SMat::mat11, "Return element (1,1)." )
      .def_prop_ro( "m22", &SMat::mat22, "Return element (2,2)." )
      .def_prop_ro( "m01", &SMat::mat01, "Return element (0,1)." )
      .def_prop_ro( "m02", &SMat::mat02, "Return element (0,2)." )
      .def_prop_ro( "m12", &SMat::mat12, "Return element (1,2)." )
      .def(
          "__mul__", []( const SMat &self, const Vec3<T> &v ) { return self * v; }, nb::is_operator() )
      .def(
          "__add__", []( const SMat &self, const SMat &other ) { return self + other; }, nb::is_operator() )
      .def(
          "__sub__", []( const SMat &self, const SMat &other ) { return self + ( -other ); }, nb::is_operator() )
      .def( "__neg__", []( const SMat &self ) { return -self; }, nb::is_operator() );
}

template <class T> void declare_rtop( nb::module_ m, const std::string &name ) {
  using RT = RTop<T>;
  nb::class_<RT> rtop( m, ( "RTop_" + name ).c_str(), "Rotation-translation operator." );
  rtop.def( nb::init<>() )
      .def( nb::init<const Mat33<T> &>(), nb::arg( "rot" ), "Constructor from rotation (Mat33)." )
      .def( nb::init<const Mat33<T> &, const Vec3<T> &>(), nb::arg( "rot" ), nb::arg( "trn" ),
            "Constructor from rotation (Mat33) and translation (Vec3)." )
      .def(
          "__init__",
          []( RT *rt, std::array<std::array<T, 3>, 3> &rot, std::array<T, 3> &trn ) {
            auto rotation = array_to_mat33( rot );
            auto translation = array_to_vec3( trn );
            new ( rt ) RT( *rotation, *translation );
          },
          nb::arg( "rot" ), nb::arg( "trn" ), "Constructor from rotation and translation(list/arrays)." )
      .def( "inverse", &RT::inverse, "Inverse." )
      .def( "equals", &RT::equals, nb::arg( "m" ), nb::arg( "tol" ), "Test equality with some tolerance." )
      .def_prop_rw(
          "rot", []( const RT &self ) { return self.rot(); }, []( RT &self, const Mat33<T> &rot ) { self.rot() = rot; },
          nb::for_getter( "Get rotation. Use .rot().as_array() to return as numpy array." ),
          nb::for_setter( "Set rotation." ) )
      .def_prop_rw(
          "trn", []( const RT &self ) { return self.trn(); }, []( RT &self, const Vec3<T> &trn ) { self.trn() = trn; },
          nb::for_getter( "Get translation. Use .trn().as_array() to return as numpy array." ),
          nb::for_setter( "Set translation." ) )
      //.def(
      //    "trn", [](RT &self, const Vec3<T> &trn) { self.trn() = trn; },
      //    "Set translation.")
      // bind set rot and trn with array
      .def_static( "identity", &RT::identity, "Return identity operator." )
      .def_static( "null", &RT::null, "Return null operator." )
      // ugly hack to bypass error due to non-ftype argument
      .def(
          "is_null",
          []( const RT &self ) {
            // Mat33<ftype> m(self.rot());
            // Vec3<ftype> v(self.trn());
            return Util::is_nan( ftype( self.rot()( 0, 0 ) ) ) && Util::is_nan( ftype( self.trn()[0] ) );
          },
          "Test for null operator (only valid for floating point types)." )
      .def(
          "format", []( const RT &self ) { return printMat( self.rot() ) + "\n" + printVec( self.trn() ); },
          "Return formatted string representation." )
      .def( "__str__", []( const RT &self ) { return printMat( self.rot() ) + "\n" + printVec( self.trn() ); } )
      .def( "__repr__", [=]( const RT &self ) { return "<clipper.RTop_" + name + " class.>"; } )
      .def(
          "__mul__", []( const RT &self, const Vec3<T> &v ) { return self * v; }, nb::is_operator() )
      .def( "__mul__", []( const RT &self, const RT &other ) { return self * other; }, nb::is_operator() );
}

template <class T> void declare_array2d( nb::module_ m, const std::string &name ) {
  using Arr2d = Array2d<T>;
  nb::class_<Arr2d> array2d( m, ( "Array2d_" + name ).c_str(), "Simple 2-d array class." );
  array2d.def( nb::init<>() )
      .def( nb::init<const int &, const int &>(), nb::arg( "rows" ), nb::arg( "cols" ),
            "Constructor with number of rows and columns." )
      .def( nb::init<const int &, const int &, T>(), nb::arg( "rows" ), nb::arg( "cols" ), nb::arg( "val" ),
            "Constructor with number of rows, columns and fill value." )
      .def( "resize", ( void ( Arr2d::* )( const int &, const int & ) )&Arr2d::resize, nb::arg( "rows" ),
            nb::arg( "cols" ), "Resize array." )
      .def( "resize", ( void ( Arr2d::* )( const int &, const int &, const T & ) )&Arr2d::resize, nb::arg( "rows" ),
            nb::arg( "cols" ), nb::arg( "val" ), "Resize array and fill with value." )
      .def_prop_ro( "size", &Arr2d::size, "Size of array." )
      .def_prop_ro( "rows", &Arr2d::rows, "Number of rows." )
      .def_prop_ro( "cols", &Arr2d::cols, "Number of columns." )
      .def(
          "shape",
          []( const Arr2d &self ) {
            // return std::array<int, 2>(self.rows(), self.cols());
            //  return [ self.rows(), self.cols() ];
            return nb::make_tuple( self.rows(), self.cols() );
          },
          "Return shape of array." )
      .def(
          "get", []( const Arr2d &self, const int &i, const int &j ) { return self( i, j ); }, nb::arg( "row" ),
          nb::arg( "col" ), "Read accessor." )
      .def(
          "set", []( Arr2d &self, const int &i, const int &j, const T &val ) { self( i, j ) = val; }, nb::arg( "row" ),
          nb::arg( "col" ), nb::arg( "val" ), "Write accessor." )
      //.def(
      //    "from_numpy",
      //    [](Arr2d &self, const py::array_t<T> vals) {
      //      numpy_to_array2d(vals, self);
      //    },
      //    "Import values from numpy array.")
      //.def("as_numpy")
      .def( "__repr__", [=]( const Arr2d &self ) {
        return "<clipper.Array2d_" + name + " class with shape (" + String( self.rows() ) + ", " +
               String( self.cols() ) + ").>";
      } );
}

template <class T> void declare_matrix( nb::module_ m, const std::string &name ) {
  using Mat = Matrix<T>;
  nb::class_<Mat, Array2d<T>> matrix( m, ( "Matrix_" + name ).c_str(),
                                      "General matrix class: like Array2d but with numerical methods" );
  matrix.def( nb::init<>() )
      .def( nb::init<const int &, const int &>(), nb::arg( "rows" ), nb::arg( "cols" ),
            "Constructor with number of rows and columns." )
      .def( nb::init<const int &, const int &, T>(), nb::arg( "rows" ), nb::arg( "cols" ), nb::arg( "val" ),
            "Constructor with number of rows, columns and fill value." )
      .def( "solve", &Mat::solve, nb::arg( "b" ), "Equation solver (square matrices only)." )
      .def( "eigen", &Mat::eigen, nb::arg( "sort" ) = true, "Eigenvalue calculation (square symmetric matrices only)." )
      .def( "__mul__", []( const Mat &self, const std::vector<T> &v ) { return self * v; } )
      .def( "__repr__", [=]( const Mat &self ) {
        return "<clipper.Matrix_" + name + " class with shape (" + String( self.rows() ) + ", " +
               String( self.cols() ) + ").>";
      } );
}

void add_clipper_types( nb::module_ &m ) {
  declare_vec3<int>( m, "int" );
  declare_mat33<int>( m, "int" );
  declare_rtop<int>( m, "int" );
  declare_array2d<int>( m, "int" );

  declare_vec3<double>( m, "double" );
  declare_mat33<double>( m, "double" );
  declare_mat33sym<double>( m, "double" );
  declare_rtop<double>( m, "double" );
  declare_array2d<double>( m, "double" );
  declare_matrix<double>( m, "double" );

  declare_vec3<float>( m, "float" );
  declare_mat33<float>( m, "float" );
  declare_mat33sym<float>( m, "float" );
  declare_rtop<float>( m, "float" );
  declare_array2d<float>( m, "float" );
  declare_matrix<float>( m, "float" );
}