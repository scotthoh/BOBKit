// Wrapper for clipper symop
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/symop.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void declare_rtopfrac(py::module &m) {
  py::class_<RTop_frac, RTop<>>(m, "RTop_frac")
      .def(py::init<>())
      .def(py::init<const RTop<> &>(), py::arg("rtop"),
           "Constructor: copy/convert")
      .def(py::init<const Mat33<> &>(), py::arg("rot"),
           "Constructor from rotation. ")
      .def(py::init<const String &>(), py::arg("strop"),
           "Constructor from string description.")
      .def(py::init<const Mat33<> &, const Vec3<> &>(), py::arg("rot"),
           py::arg("trn"), "Constructor from rotation and translation.")
      .def(py::init([](py::array_t<ftype> rot, py::array_t<ftype> trn) {
             auto rotation = numpy_to_mat33(rot);
             auto translation = numpy_to_vec3(trn);
             return std::unique_ptr<RTop_frac>(
                 new RTop_frac(*rotation, *translation));
           }),
           py::arg("rot"), py::arg("trn"),
           "Constructor from list/array of rotation and translation.")
      .def("rtop_orth", &RTop_frac::rtop_orth, py::arg("cell"),
           "Fractional-orthogonal conversion.")
      .def("inverse", &RTop_frac::inverse, "Inverse operator.")
      .def_static("identity", &RTop_frac::identity, "Return identity operator.")
      .def_static("null", &RTop_frac::null,
                  "Return null (uninitialised) operator.")
      .def("__repr__",
           [](const RTop_frac &self) {
             return "<clipper.RTop_frac: Fractional operator class.>";
           })
      .doc() = "Fractional operator class.\nThis class is used "
               "for any RT-operator which operates on fractional "
               "coordinates. For a full list of methods, see clipper::RTop.";
}

void declare_symop(py::module &m) {
  py::class_<Symop, RTop_frac>(m, "Symop")
      .def(py::init<>())
      .def(py::init<const RTop<> &>(), py::arg("rtop"),
           "Constructor from RTop.")
      // 4x4 array
      .def(py::init([](const py::array_t<ftype> a) {
             check_array_shape(a, {4, 4}, true);
             auto buf = a.request();
             return std::unique_ptr<Symop>(
                 new Symop(*reinterpret_cast<ftype(*)[4][4]>(buf.ptr)));
           }),
           py::arg("rtop"), "Constructor from 4x4 matrix (list/array).")
      .def("format", &Symop::format, "Return formatted string representation.")
      .def("__str__", &Symop::format)
      .def(
          "__repr__",
          [](const Symop &self) {
            return "<clipper.Symop: Crystallographic symmetry operator class.>";
          })
      .doc() = "Crystallographic symmetry operator.\nThis is identical to "
               "a fractional RTop, but has its own class since not all "
               "fractional  RTops are symops. For a full list of methods, "
               "see clipper::RTop and clipper::RTop_frac.";
}

void declare_isymop(py::module &m) {
  py::class_<Isymop, RTop<int>>(m, "Isymop")
      .def(py::init<>())
      .def(py::init<const RTop<int> &>(), py::arg("rtop"),
           "Constructor from RTop.")
      .def(py::init<const Symop &, const Grid &>(), py::arg("symop"),
           py::arg("grid"), "Constructor from symop and grid.")
      .def("__repr__",
           [](const Isymop &self) {
             return "<clipper.Isymop: Integerised symmetry matrix class.>";
           })
      .doc() = "Integerised symmetry matrix.\nThis is used for "
               "optimised calculations in real and reciprocal space.";
}

void declare_symop_code(py::module &m) {
  py::class_<Symop_code>(m, "Symop_code")
      .def(py::init<>())
      .def(py::init<const int &>(), py::arg("code"), "Constructor from int.")
      .def(py::init<const Symop &>(), py::arg("symop"),
           "Constructor from Symop.")
      .def(py::init<const Isymop &>(), py::arg("isymop"),
           "Contructor from Isymop.")
      .def("init", &Symop_code::init, py::arg("isymop"),
           "Initialiser from Isymop.")
      .def("code_rot", &Symop_code::code_rot,
           "Return code for rotational part.")
      .def("code_trn", &Symop_code::code_trn,
           "Return code for translational part.")
      .def("symop", &Symop_code::symop, "Convert to Symop.")
      .def("isymop", &Symop_code::isymop, "Convert to integerised symop.")
      .def_static("identity", &Symop_code::identity, "Identity code.")
      .def(
          "__int__", [](const Symop_code &self) { return int(self); },
          "Convert to integer.")
      .def("__repr__",
           [](const Symop_code &self) {
             return "<clipper.Symop_code: Compressed encoded symmetry operator "
                    "class.>";
           })
      .doc() =
      "Compressed encoded symmetry operator.\n"
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

void init_symop(py::module &m) {
  declare_rtopfrac(m);
  declare_symop(m);
  declare_isymop(m);
  declare_symop_code(m);
}
