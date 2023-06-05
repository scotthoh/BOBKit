#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "type_conversions.h"
#include <clipper/clipper.h>
// #include <clipper/core/clipper_types.h>
#include "helper_functions.h"

namespace py = pybind11;
using namespace clipper;

template <class T>
void declare_vec3(py::module &m, const std::string &name)
{
  using VecClass = Vec3<T>;
  std::string PyClass = std::string("Vec3_") + name;
  py::class_<VecClass> vec3(m, PyClass.c_str());
  vec3
      .def(py::init<>([]()
                      {
        Vec3<T> *v = new Vec3<T>();
        //Util::set_null(v[0]);
        return std::unique_ptr<Vec3<T>>(v); }))
      .def(py::init<const T &, const T &, const T &>())
      // init from array or list
      .def(py::init([](std::array<T, 3> &a)
                    {
                      return std::unique_ptr<Vec3<T>>(new Vec3<T>(a[0], a[1], a[2]));

                      //(a[0], a[1], a[2]);
                      // return v;
                    }))
      .def(py::init<const Vec3<ftype32> &>())
      .def(py::init<const Vec3<ftype64> &>())
      .def(py::init<const Vec3<int> &>())
      .def("equals", &VecClass::equals, py::arg("v"), py::arg("tol"))
      .def("__getitem__", [](const VecClass &self, const int i)
           { return self[normalise_index(i, 3)]; })
      .def("__setitem__", [](VecClass &self, const int &i, const T &val)
           { self[normalise_index(i, 3)] = val; })
      .def("as_array", [](const VecClass &self)
           { return make_array_1d<VecClass, T>(self, 3); })
      .def("from_array", [](VecClass &self, py::array_t<T> values)
           { fill_array_1d<VecClass, T>(self, 3, values); })
      .def("unit", &VecClass::unit)
      .def_static("zero", &VecClass::zero)
      .def_static("null", &VecClass::null)
      .def("is_null", &VecClass::is_null)
      .def_static("dot", &VecClass::dot)
      .def_static("cross", &VecClass::cross)
      .def("__str__", [](const VecClass &self)
           { return self.format(); })
      .def("__repr__", [](const VecClass &self)
           { return "<clipper.Vec3 " + self.format() + ">"; })
      .def("format", &VecClass::format)

      //.def("__iadd__", &VecClass::operator+=)
      //.def("__isub__", &VecClass::operator-=)
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(
          "__eq__", [](const VecClass &self, const VecClass &v)
          { return self == v; },
          py::is_operator())
      .def(
          "__ne__", [](const VecClass &self, const VecClass &v)
          { return self != v; },
          py::is_operator())
      .def(
          "__neg__", [](const VecClass &self)
          { return -self; },
          py::is_operator())
      .def(
          "__add__", [](const VecClass &self, const VecClass &v)
          { return self + v; },
          py::is_operator())
      .def(
          "__sub__", [](const VecClass &self, const VecClass &v)
          { return self - v; },
          py::is_operator())
      .def(
          "__mul__", [](const VecClass &self, const T &s)
          { return self * s; },
          py::is_operator())
      .def(
          "__rmul__", [](const VecClass &self, const T &s)
          { return s * self; },
          py::is_operator())
      .def(
          "__mul__", [](const VecClass &self, const VecClass &v)
          { return Vec3<T>::dot(self, v); },
          py::is_operator());
}

template <class T>
void declare_mat33(py::module m, const std::string &name)
{
  using Mat33Class = Mat33<T>;
  std::string PyClass = std::string("Mat33_") + name;
  py::class_<Mat33Class> mat33(m, PyClass.c_str(), py::buffer_protocol());
  mat33
      .def(py::init<>())
      .def(py::init<const T &, const T &, const T &, const T &, const T &, const T &, const T &, const T &, const T &>())
      .def(py::init<const Mat33<ftype32> &>())
      .def(py::init<const Mat33<ftype64> &>())
      .def(py::init<const Mat33sym<ftype32> &>())
      .def(py::init<const Mat33sym<ftype64> &>())
      .def_buffer([](Mat33Class &self) -> py::buffer_info
                  { return py::buffer_info(
                        &self(0, 0),
                        {3, 3},
                        {sizeof(T) * 3, sizeof(T)}); })
      .def("det", &Mat33Class::det)
      .def("inverse", &Mat33Class::inverse)
      .def("transpose", &Mat33Class::transpose)
      .def("equals", &Mat33Class::equals, py::arg("m"), py::arg("tol"))
      .def("get", [](Mat33Class &self, const int i, const int j)
           { return self(i, j); })
      .def("set", [](Mat33Class &self, const int i, const int j, const T val)
           { self(i, j) = val; })
      .def("format", &Mat33Class::format)
      .def("__str__", &Mat33Class::format)
      .def("__repr__", [](const Mat33Class &self)
           { return "<clipper.Mat33 |" + String(self(0, 0), 10, 4) + "," + String(self(0, 1), 10, 4) + "," + String(self(0, 2), 10, 4) + "|\n" +
                    "               |" + String(self(1, 0), 10, 4) + "," + String(self(1, 1), 10, 4) + "," + String(self(1, 2), 10, 4) + "|\n" +
                    "               |" + String(self(2, 0), 10, 4) + "," + String(self(2, 1), 10, 4) + "," + String(self(2, 2), 10, 4) + "|>"; })
      .def_static("identity", &Mat33Class::identity)
      .def_static("null", &Mat33Class::null)
      // something weird with the null constructor that is_null is not returning true
      .def("is_null", &Mat33Class::is_null)
      // operator Matrix-vector, assumes column vector
      .def(
          "__mul__", [](const Mat33Class &self, const Vec3<T> &v)
          { return self * v; },
          py::is_operator())
      .def(
          "__rmul__", [](const Vec3<T> &v, const Mat33Class &self)
          { return v * self; },
          py::is_operator())
      .def(
          "__mul__", [](const Mat33Class &self, const Mat33Class &m)
          { return self * m; },
          py::is_operator())
      .def(
          "__add__", [](const Mat33Class &self, const Mat33Class &m)
          { return self + m; },
          py::is_operator())
      .def(
          "__sub__", [](const Mat33Class &self, const Mat33Class &m)
          { return self + (-m); },
          py::is_operator())
      .def(
          "__neg__", [](const Mat33Class &self)
          { return -self; },
          py::is_operator());
}

template <class T>
void declare_mat33sym(py::module m, const std::string &name)
{
  using SMatClass = Mat33sym<T>;
  std::string PyClass = std::string("Mat33sym_") + name;
  py::class_<SMatClass> smat(m, PyClass.c_str());
  smat
      .def(py::init<>())
      .def(py::init<const Mat33<ftype32> &>())
      .def(py::init<const Mat33<ftype64> &>())
      .def(py::init<const Mat33sym<ftype32> &>())
      .def(py::init<const Mat33sym<ftype64> &>())
      .def(py::init<const T &, const T &, const T &, const T &, const T &, const T &>())
      .def("__str__", &SMatClass::format)
      .def("format", &SMatClass::format)
      .def_static("identity", &SMatClass::identity)
      .def_static("null", &SMatClass::null)
      .def("is_null", &SMatClass::is_null)
      .def("quad_form", &SMatClass::quad_form)
      .def("det", &SMatClass::det)
      .def("sqrt", &SMatClass::sqrt)
      .def("inverse", &SMatClass::inverse)
      .def("m00", &SMatClass::mat00)
      .def("m11", &SMatClass::mat11)
      .def("m22", &SMatClass::mat22)
      .def("m01", &SMatClass::mat01)
      .def("m02", &SMatClass::mat02)
      .def("m12", &SMatClass::mat12)
      .def("get", [](const SMatClass &self, const int &i, const int &j)
           { return self(i, j); })
      .def(
          "__mul__", [](const SMatClass &self, const Vec3<T> &v)
          { return self * v; },
          py::is_operator())
      .def(
          "__add__", [](const SMatClass &self, const SMatClass &m)
          { return self + m; },
          py::is_operator())
      .def(
          "__sub__", [](const SMatClass &self, const SMatClass &m)
          { return self + (-m); },
          py::is_operator())
      .def(
          "__neg__", [](const SMatClass &self)
          { return -self; },
          py::is_operator());
}

template <class T>
void declare_rtop(py::module m, const std::string &name)
{
  using RTopClass = RTop<T>;
  std::string PyClass = std::string("RTop_") + name;
  py::class_<RTopClass> rtop(m, PyClass.c_str());
  rtop
      .def(py::init<>())
      .def(py::init<const Mat33<T> &>(), py::arg("rot"))
      .def(py::init<const Mat33<T> &, const Vec3<T> &>(), py::arg("rot"), py::arg("trn"))
      .def("inverse", &RTopClass::inverse)
      .def("equals", &RTopClass::equals, py::arg("m"), py::arg("tol"))
      .def(
          "rot", [](const RTopClass &self)
          { return self.rot(); })
      .def(
          "rot", [](RTopClass &self, const Mat33<T> &rot)
          { self.rot() = rot; })
      .def(
          "trn", [](const RTopClass &self)
          { return self.trn(); })
      .def("trn", [](RTopClass &self, const Vec3<T> &trn)
           { self.trn() = trn; })
      .def_static("identity", &RTopClass::identity)
      .def_static("null", &RTopClass::null)
      .def("is_null", &RTopClass::is_null)
      .def("format", &RTopClass::format)
      .def("__str__", &RTopClass::format)
      .def(
          "__mul__", [](const RTopClass &self, const Vec3<T> &v)
          { return self * v; },
          py::is_operator())
      .def(
          "__mul__", [](const RTopClass &self, const RTopClass &rt2)
          { return self * rt2; },
          py::is_operator());
}

template <class T>
void declare_array2d(py::module m, const std::string &name)
{
  using Array2dClass = Array2d<T>;
  std::string PyClass = std::string("Array2d_") + name;
  py::class_<Array2dClass> array2d(m, PyClass.c_str());
  array2d
      .def(py::init<>())
      .def(py::init<const int &, const int &>(), py::arg("d1"), py::arg("d2"))
      .def(py::init<const int &, const int &, T>(), py::arg("d1"), py::arg("d2"), py::arg("val"))
      .def("resize", (void(Array2dClass::*)(const int &, const int &)) & Array2dClass::resize)
      .def("resize", (void(Array2dClass::*)(const int &, const int &, const T &)) & Array2dClass::resize)
      .def("size", &Array2dClass::size)
      .def("rows", &Array2dClass::rows)
      .def("cols", &Array2dClass::cols)
      .def("get", [](const Array2dClass &self, const int &i, const int &j)
           { return self(i, j); })
      .def("set", [](Array2dClass &self, const int &i, const int &j, const T &val)
           { self(i, j) = val; });
}

template <class T>
void declare_matrix(py::module m, const std::string &name)
{
  using MatrixClass = Matrix<T>;
  std::string PyClass = std::string("Matrix_") + name;
  py::class_<MatrixClass, Array2d<T>> matrix(m, PyClass.c_str());
  matrix
      .def(py::init<>())
      .def(py::init<const int &, const int &>(), py::arg("d1"), py::arg("d2"))
      .def(py::init<const int &, const int &, T>(), py::arg("d1"), py::arg("d2"), py::arg("val"))
      .def("solve", &MatrixClass::solve)
      .def("eigen", &MatrixClass::eigen)
      .def("__mul__", [](const MatrixClass &self, const std::vector<T> &v)
           { return self * v; });
}

void init_clipper_types(py::module &m)
{
  declare_vec3<int>(m, "int");
  declare_mat33<int>(m, "int");
  declare_rtop<int>(m, "int");
  declare_array2d<int>(m, "int");

  declare_vec3<ftype64>(m, "float64");
  declare_mat33<ftype64>(m, "float64");
  declare_mat33sym<ftype64>(m, "float64");
  declare_rtop<ftype64>(m, "float64");
  declare_array2d<ftype64>(m, "float64");
  declare_matrix<ftype64>(m, "float64");

  declare_vec3<ftype32>(m, "float32");
  declare_mat33<ftype32>(m, "float32");
  declare_mat33sym<ftype32>(m, "float32");
  declare_rtop<ftype32>(m, "float32");
  declare_array2d<ftype32>(m, "float32");
  declare_matrix<ftype32>(m, "float32");

  // declare_rtop<ftype32>(m, "float")
  //  declare_vec3<
}