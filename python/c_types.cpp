// Wrapper for clipper types
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/clipper_types.h>
#include <clipper/core/clipper_util.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <fstream>
#include <utility>

// #include <clipper/core/clipper_types.h>
// change printVec to printVec with template<class T>
// change format() for mat33, rtop to deal with int

namespace py = pybind11;
using namespace clipper;

//! Return formatted string representation for Vec3<>
template <class T> std::string printVec(Vec3<T> v) {
  return "(" + String(ftype64(v[0]), 10, 4) + "," +
         String(ftype64(v[1]), 10, 4) + "," + String(ftype64(v[2]), 10, 4) +
         ")";
}

//! Return formatted string representation for Mat33<>
template <class T> std::string printMat(Mat33<T> m) {
  return "|" + String(ftype64(m(0, 0)), 10, 4) + "," +
         String(ftype64(m(0, 1)), 10, 4) + "," +
         String(ftype64(m(0, 2)), 10, 4) + "|\n|" +
         String(ftype64(m(1, 0)), 10, 4) + "," +
         String(ftype64(m(1, 1)), 10, 4) + "," +
         String(ftype64(m(1, 2)), 10, 4) + "|\n|" +
         String(ftype64(m(2, 0)), 10, 4) + "," +
         String(ftype64(m(2, 1)), 10, 4) + "," +
         String(ftype64(m(2, 2)), 10, 4) + "|";
}

//! Return formatted string representation for Mat33sym<>
template <class T> std::string printSMat(Mat33sym<T> m) {
  return "|" + String(ftype64(m.mat00()), 10, 4) + "," +
         String(ftype64(m.mat01()), 10, 4) + "," +
         String(ftype64(m.mat02()), 10, 4) + "|\n|" +
         String(ftype64(m.mat01()), 10, 4) + "," +
         String(ftype64(m.mat11()), 10, 4) + "," +
         String(ftype64(m.mat12()), 10, 4) + "|\n|" +
         String(ftype64(m.mat02()), 10, 4) + "," +
         String(ftype64(m.mat12()), 10, 4) + "," +
         String(ftype64(m.mat22()), 10, 4) + "|";
}

template <class T> void declare_vec3(py::module &m, const std::string &name) {
  using VecClass = Vec3<T>;
  std::string PyClass = std::string("Vec3_") + name;
  py::class_<VecClass> vec3(m, PyClass.c_str(), py::buffer_protocol(),
                            "3-vector Class");
  vec3.def(py::init<>(), "Null constructor.")
      .def(py::init<const T &, const T &, const T &>(), py::arg("x"),
           py::arg("y"), py::arg("z"), "Constructor from individual values.")
      // init from array or list
      .def(py::init([](std::array<T, 3> &a) {
             return std::unique_ptr<Vec3<T>>(new Vec3<T>(a[0], a[1], a[2]));
           }),
           py::arg("a"), "Constructor from list/array.")
      .def(py::init<const Vec3<ftype32> &>(), py::arg("v"),
           "Constructor: copy/convert.")
      .def(py::init<const Vec3<ftype64> &>(), py::arg("v"),
           "Constructor: copy/convert.")
      .def(py::init<const Vec3<int> &>(), py::arg("v"),
           "Constructor: copy/convert.")
      .def_buffer([](VecClass &self) -> py::buffer_info {
        return py::buffer_info(&self[0], sizeof(T),
                               py::format_descriptor<T>::format(), 1, {3},
                               {sizeof(T)});
      })
      .def("equals", &VecClass::equals, py::arg("vec"), py::arg("tol"),
           "Test equality.")
      .def("__getitem__",
           [](const VecClass &self, const int i) {
             return self[normalise_index(i, 3)];
           })
      .def("__setitem__",
           [](VecClass &self, const int &i, const T &val) {
             self[normalise_index(i, 3)] = val;
           })
      .def(
          "as_array",
          [](const VecClass &self) {
            return make_array_1d<VecClass, T>(self, 3);
          },
          "Return as array.")
      .def(
          "from_array",
          [](VecClass &self, py::array_t<T> values) {
            fill_array_1d<VecClass, T>(self, 3, values);
          },
          "Import values from array.")
      .def("unit", &VecClass::unit,
           "Return unit vector with same direction as this vector.")
      .def_static("zero", &VecClass::zero, "Return zero vector.")
      .def_static("null", &VecClass::null,
                  "Return null vector (only valid for floating point types).")
      // ugly hack to bypass error due to non-ftype argument
      // passed to clipper::Util::is_nan(ftype)
      // or we can edit clipper::Util to add is_nan(int)???
      .def(
          "is_null",
          [](const VecClass &self) {
            Vec3<ftype64> vtmp(self);
            return vtmp.is_null();
          },
          "Test for null vector (only valid for floating point types).")
      .def_static("dot", &VecClass::dot, py::arg("v1"), py::arg("v2"),
                  "Vector dot product. (Equivalent to v1 * v2)")
      .def_static("cross", &VecClass::cross, py::arg("v1"), py::arg("v2"),
                  "Vector cross product.")
      .def(
          "__iter__",
          [](VecClass &self) { return py::make_iterator(&self[0], &self[3]); },
          py::return_value_policy::reference_internal)
      // don't understand why copy constructor requires len() in python side
      // thus this is defined to avoid error when calling constructor
      // Vec3_T(Vec3_T)
      .def("__len__", [](const VecClass &self) { return 3; })
      // also hack for error because pybind11 don't like
      // the non-existent function String(int,w,p) called in format()
      .def("__str__",
           [](const VecClass &self) {
             return printVec(self);
             // return printVec(ftype64(self[0]), ftype64(self[1]),
             //                 ftype64(self[2]));
             //"(" + String(ftype64(self[0]), 10, 4) + "," +
             //        String(ftype64(self[1]), 10, 4) + "," +
             //        String(ftype64(self[2]), 10, 4) + ")";
           })
      .def(
          "format",
          [](const VecClass &self) {
            return printVec(self);
            // return printVec(ftype64(self[0]), ftype64(self[1]),
            //                 ftype64(self[2]));
            //  return "(" + String(ftype64(self[0]), 10, 4) + "," +
            //         String(ftype64(self[1]), 10, 4) + "," +
            //         String(ftype64(self[2]), 10, 4) + ")";
          },
          "Return formatted string representation.")
      .def("__repr__",
           [=](const VecClass &self) {
             return "<clipper.Vec3_" + name + " " + printVec(self) + ">";
             // return "<clipper.Vec3 (" + String(ftype64(self[0]), 10, 4) + ","
             // +
             //        String(ftype64(self[1]), 10, 4) + "," +
             //        String(ftype64(self[2]), 10, 4) + ")>";
           })
      //.def("__iadd__", &VecClass::operator+=)
      //.def("__isub__", &VecClass::operator-=)
      .def("__iadd__",
           [](VecClass &self, const VecClass &other) { return self += other; })
      .def("__isub__",
           [](VecClass &self, const VecClass &other) { return self -= other; })
      //.def(py::self += py::self)
      //.def(py::self -= py::self)
      .def(
          "__eq__",
          [](const VecClass &self, const VecClass &other) {
            return self == other;
          },
          py::is_operator())
      .def(
          "__ne__",
          [](const VecClass &self, const VecClass &other) {
            return self != other;
          },
          py::is_operator())
      .def(
          "__neg__", [](const VecClass &self) { return -self; },
          py::is_operator())
      .def(
          "__add__",
          [](const VecClass &self, const VecClass &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const VecClass &self, const VecClass &v) { return self - v; },
          py::is_operator())
      .def(
          "__mul__", [](const VecClass &self, const T &s) { return self * s; },
          py::is_operator())
      .def(
          "__rmul__", [](const VecClass &self, const T &s) { return s * self; },
          py::is_operator())
      .def(
          "__mul__",
          [](const VecClass &self, const VecClass &other) {
            return Vec3<T>::dot(self, other);
          },
          py::is_operator());
}

template <class T> void declare_mat33(py::module m, const std::string &name) {
  using Mat33Class = Mat33<T>;
  std::string PyClass = std::string("Mat33_") + name;
  py::class_<Mat33Class> mat33(m, PyClass.c_str(), py::buffer_protocol(),
                               "3x3 matrix class.");
  mat33.def(py::init<>())
      .def(py::init<const T &, const T &, const T &, const T &, const T &,
                    const T &, const T &, const T &, const T &>(),
           "Constructor from individual values.")
      .def(py::init<const Mat33<ftype32> &>(), "Constructor: copy/convert.")
      .def(py::init<const Mat33<ftype64> &>(), "Constructor: copy/convert.")
      .def(py::init<const Mat33sym<ftype32> &>(),
           "Constructor: copy/convert from symmetric matrix.")
      .def(py::init<const Mat33sym<ftype64> &>(),
           "Constructor: copy/convert from symmetric matrix.")
      .def(py::init([](py::array_t<T> vals) { return numpy_to_mat33(vals); }),
           py::arg("a"), "Constructor from 3x3 numpy array.")
      .def_buffer([](Mat33Class &self) -> py::buffer_info {
        return py::buffer_info(&self(0, 0), {3, 3}, {sizeof(T) * 3, sizeof(T)});
      })
      .def("det", &Mat33Class::det, "Determinant.")
      .def("inverse", &Mat33Class::inverse, "Inverse.")
      .def("transpose", &Mat33Class::transpose, "Transpose.")
      .def("equals", &Mat33Class::equals, py::arg("m"), py::arg("tol"),
           "Test equality.")
      .def(
          "get",
          [](Mat33Class &self, const int i, const int j) { return self(i, j); },
          "Get element.")
      .def(
          "set",
          [](Mat33Class &self, const int i, const int j, const T val) {
            self(i, j) = val;
          },
          "Set element.")
      //.def("format", &Mat33Class::format)
      .def(
          "format", [](const Mat33Class &self) { return printMat(self); },
          "Return formatted string representation.")
      .def("__str__", [](const Mat33Class &self) { return printMat(self); })
      .def("__repr__",
           [=](const Mat33Class &self) {
             return "<clipper.Matt33_" + name + " class.>";
           })
      //.def("__str__", &Mat33Class::format)
      //.def("__repr__", [](const Mat33Class &self)
      //     { return "<clipper.Mat33 |" + String(self(0, 0), 10, 4) + "," +
      //     String(self(0, 1), 10, 4) + "," + String(self(0, 2), 10, 4) + "|\n"
      //     +
      //              "               |" + String(self(1, 0), 10, 4) + "," +
      //              String(self(1, 1), 10, 4) + "," + String(self(1, 2), 10,
      //              4) + "|\n" + "               |" + String(self(2, 0), 10,
      //              4) + "," + String(self(2, 1), 10, 4) + "," +
      //              String(self(2, 2), 10, 4) + "|>"; })
      .def_static("identity", &Mat33Class::identity, "Return identity matrix.")
      .def_static("null", &Mat33Class::null,
                  "Return null matrix (only valid for floating point types).")
      // ugly hack to bypass error due to non-ftype argument
      // passed to clipper::Util::is_nan(ftype)
      // or we can edit clipper::Util to add is_nan(int)???
      .def(
          "is_null",
          [](const Mat33Class &self) {
            Mat33<ftype64> mtmp(self);
            return mtmp.is_null();
          },
          "Test for null vector (only valid for floating point types).")
      // operator Matrix-vector, assumes column vector
      .def(
          "__mul__",
          [](const Mat33Class &self, const Vec3<T> &v) { return self * v; },
          py::is_operator())
      .def(
          "__rmul__",
          [](const Vec3<T> &v, const Mat33Class &self) { return v * self; },
          py::is_operator())
      .def(
          "__mul__",
          [](const Mat33Class &self, const Mat33Class &other) {
            return self * other;
          },
          py::is_operator())
      .def(
          "__add__",
          [](const Mat33Class &self, const Mat33Class &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const Mat33Class &self, const Mat33Class &other) {
            return self + (-other);
          },
          py::is_operator())
      .def(
          "__neg__", [](const Mat33Class &self) { return -self; },
          py::is_operator());
}

template <class T>
void declare_mat33sym(py::module m, const std::string &name) {
  using SMatClass = Mat33sym<T>;
  std::string PyClass = std::string("Mat33sym_") + name;
  py::class_<SMatClass> smat(m, PyClass.c_str(),
                             "Compressed form for 3x3 symmetric matrix class.");
  smat.def(py::init<>())
      .def(py::init<const Mat33<ftype32> &>(),
           "Constructor from Mat33 (does not check for symmetry).")
      .def(py::init<const Mat33<ftype64> &>(),
           "Constructor from Mat33 (does not check for symmetry).")
      .def(py::init<const Mat33sym<ftype32> &>(), "Constructor from Mat33sym.")
      .def(py::init<const Mat33sym<ftype64> &>(), "Constructor from Mat33sym.")
      .def(py::init<const T &, const T &, const T &, const T &, const T &,
                    const T &>(),
           py::arg("c00"), py::arg("c11"), py::arg("c22"), py::arg("c01"),
           py::arg("c02"), py::arg("c12"), "Constructor from coefficients.")
      //.def("__str__", &SMatClass::format)
      //.def("format", &SMatClass::format)
      .def(
          "format", [](const SMatClass &self) { return printSMat(self); },
          "Return formatted string representation.")
      .def("__str__", [](const SMatClass &self) { return printSMat(self); })
      .def("__repr__",
           [=](const SMatClass &self) {
             return "<clipper.Mat33sym_" + name + " class.>";
           })
      .def_static("identity", &SMatClass::identity, "Return identity matrix.")
      .def_static("null", &SMatClass::null,
                  "Return null matrix (only valid for floating point types).")
      //.def("is_null", &SMatClass::is_null)
      .def("quad_form", &SMatClass::quad_form, py::arg("v"),
           "Return quadratic form with vector.")
      .def("det", &SMatClass::det, "Determinant.")
      .def("sqrt", &SMatClass::sqrt, "Square root.")
      .def("inverse", &SMatClass::inverse, "Inverse.")
      .def_property_readonly("m00", &SMatClass::mat00, "Return element (0,0).")
      .def_property_readonly("m11", &SMatClass::mat11, "Return element (1,1).")
      .def_property_readonly("m22", &SMatClass::mat22, "Return element (2,2).")
      .def_property_readonly("m01", &SMatClass::mat01, "Return element (0,1).")
      .def_property_readonly("m02", &SMatClass::mat02, "Return element (0,2).")
      .def_property_readonly("m12", &SMatClass::mat12, "Return element (1,2).")
      .def(
          "get",
          [](const SMatClass &self, const int &i, const int &j) {
            return self(i, j);
          },
          "Get element (inefficient).")
      .def(
          "__mul__",
          [](const SMatClass &self, const Vec3<T> &v) { return self * v; },
          py::is_operator())
      .def(
          "__add__",
          [](const SMatClass &self, const SMatClass &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const SMatClass &self, const SMatClass &other) {
            return self + (-other);
          },
          py::is_operator())
      .def(
          "__neg__", [](const SMatClass &self) { return -self; },
          py::is_operator());
}

template <class T> void declare_rtop(py::module m, const std::string &name) {
  using RTopClass = RTop<T>;
  std::string PyClass = std::string("RTop_") + name;
  py::class_<RTopClass> rtop(m, PyClass.c_str(),
                             "Rotation-translation operator.");
  rtop.def(py::init<>())
      .def(py::init<const Mat33<T> &>(), py::arg("rot"),
           "Constructor from rotation (Mat33).")
      .def(py::init<const Mat33<T> &, const Vec3<T> &>(), py::arg("rot"),
           py::arg("trn"),
           "Constructor from rotation (Mat33) and translation (Vec3).")
      .def(py::init([](py::array_t<T> rot, py::array_t<T> trn) {
             auto rotation = numpy_to_mat33(rot);
             auto translation = numpy_to_vec3(trn);
             return std::unique_ptr<RTopClass>(
                 new RTopClass(*rotation, *translation));
           }),
           py::arg("rot"), py::arg("trn"),
           "Constructor from rotation and translation(list/arrays).")
      .def("inverse", &RTopClass::inverse, "Inverse.")
      .def("equals", &RTopClass::equals, py::arg("m"), py::arg("tol"),
           "Test equality with some tolerance.")
      .def(
          "rot", [](const RTopClass &self) { return self.rot(); },
          "Get rotation. Use .rot().as_array() to return as numpy array.")
      .def(
          "rot", [](RTopClass &self, const Mat33<T> &rot) { self.rot() = rot; },
          "Set rotation.")
      .def(
          "trn", [](const RTopClass &self) { return self.trn(); },
          "Get translation. Use .trn().as_array() to return as numpy array.")
      .def(
          "trn", [](RTopClass &self, const Vec3<T> &trn) { self.trn() = trn; },
          "Set translation.")
      // bind set rot and trn with array
      .def_static("identity", &RTopClass::identity, "Return identity operator.")
      .def_static("null", &RTopClass::null, "Return null operator.")
      //.def("is_null", &RTopClass::is_null)
      //.def("format", &RTopClass::format)
      .def(
          "format",
          [](const RTopClass &self) {
            return printMat(self.rot()) + "\n" + printVec(self.trn());
          },
          "Return formatted string representation.")
      .def("__str__",
           [](const RTopClass &self) {
             return printMat(self.rot()) + "\n" + printVec(self.trn());
           })
      .def("__repr__",
           [=](const RTopClass &self) {
             return "<clipper.RTop_" + name + " class.>";
           })
      .def(
          "__mul__",
          [](const RTopClass &self, const Vec3<T> &v) { return self * v; },
          py::is_operator())
      .def(
          "__mul__",
          [](const RTopClass &self, const RTopClass &other) {
            return self * other;
          },
          py::is_operator());
}

template <class T> void declare_array2d(py::module m, const std::string &name) {
  using Array2dClass = Array2d<T>;
  std::string PyClass = std::string("Array2d_") + name;
  py::class_<Array2dClass> array2d(m, PyClass.c_str(),
                                   "Simple 2-d array class.");
  array2d.def(py::init<>())
      .def(py::init<const int &, const int &>(), py::arg("rows"),
           py::arg("cols"), "Constructor with number of rows and columns.")
      .def(py::init<const int &, const int &, T>(), py::arg("rows"),
           py::arg("cols"), py::arg("val"),
           "Constructor with number of rows, columns and fill value.")
      .def("resize",
           (void(Array2dClass::*)(const int &, const int &)) &
               Array2dClass::resize,
           py::arg("rows"), py::arg("cols"), "Resize array.")
      .def("resize",
           (void(Array2dClass::*)(const int &, const int &, const T &)) &
               Array2dClass::resize,
           py::arg("rows"), py::arg("cols"), py::arg("val"),
           "Resize array and fille with value.")
      .def("size", &Array2dClass::size, "Size of array.")
      .def("rows", &Array2dClass::rows, "Number of rows.")
      .def("cols", &Array2dClass::cols, "Number of columns.")
      .def(
          "shape",
          [](const Array2dClass &self) {
            // return std::array<int, 2>(self.rows(), self.cols());
            //  return [ self.rows(), self.cols() ];
            return py::make_tuple(self.rows(), self.cols());
          },
          "Return shape of array.")
      .def(
          "get",
          [](const Array2dClass &self, const int &i, const int &j) {
            return self(i, j);
          },
          py::arg("row"), py::arg("col"), "Read accessor.")
      .def(
          "set",
          [](Array2dClass &self, const int &i, const int &j, const T &val) {
            self(i, j) = val;
          },
          py::arg("row"), py::arg("col"), py::arg("val"), "Write accessor.")
      .def(
          "from_numpy",
          [](Array2dClass &self, const py::array_t<T> vals) {
            numpy_to_array2d(vals, self);
          },
          "Import values from numpy array.")
      //.def("as_numpy")
      .def("__repr__", [=](const Array2dClass &self) {
        std::stringstream stream;
        stream << "<clipper.Array2d_" << name << " class with shape (";
        stream << self.rows() << ", " << self.cols() << ").>";
        return stream.str();
      });
}

template <class T> void declare_matrix(py::module m, const std::string &name) {
  using MatrixClass = Matrix<T>;
  std::string PyClass = std::string("Matrix_") + name;
  py::class_<MatrixClass, Array2d<T>> matrix(
      m, PyClass.c_str(),
      "General matrix class: like Array2d but with numerical methods");
  matrix.def(py::init<>())
      .def(py::init<const int &, const int &>(), py::arg("rows"),
           py::arg("cols"), "Constructor with number of rows and columns.")
      .def(py::init<const int &, const int &, T>(), py::arg("rows"),
           py::arg("cols"), py::arg("val"),
           "Constructor with number of rows, columns and fill value.")
      .def("solve", &MatrixClass::solve, py::arg("b"),
           "Equation solver (square matrices only).")
      .def("eigen", &MatrixClass::eigen, py::arg("sort") = true,
           "Eigenvalue calculation (square symmetric matrices only).")
      .def("__mul__", [](const MatrixClass &self,
                         const std::vector<T> &v) { return self * v; })
      .def("__repr__", [=](const MatrixClass &self) {
        std::stringstream stream;
        stream << "<clipper.Matrix_" << name << " class with shape (";
        stream << self.rows() << ", " << self.cols() << ").>";
        return stream.str();
      });
}

void init_clipper_types(py::module &m) {
  declare_vec3<int>(m, "int");
  declare_mat33<int>(m, "int");
  declare_rtop<int>(m, "int");
  declare_array2d<int>(m, "int");

  declare_vec3<ftype64>(m, "double");
  declare_mat33<ftype64>(m, "double");
  declare_mat33sym<ftype64>(m, "double");
  declare_rtop<ftype64>(m, "double");
  declare_array2d<ftype64>(m, "double");
  declare_matrix<ftype64>(m, "double");

  declare_vec3<ftype32>(m, "float");
  declare_mat33<ftype32>(m, "float");
  declare_mat33sym<ftype32>(m, "float");
  declare_rtop<ftype32>(m, "float");
  declare_array2d<ftype32>(m, "float");
  declare_matrix<ftype32>(m, "float");
}