#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void declare_rtopfrac(py::module &m) {
  py::class_<RTop_frac, RTop<>>(m, "RTop_frac")
      .def(py::init<>())
      .def(py::init<const RTop<> &>(), py::arg("rtop"))
      .def(py::init<const Mat33<> &>(), py::arg("rot"))
      .def(py::init<const String &>(), py::arg("strop"))
      .def(py::init<const Mat33<> &, const Vec3<> &>(), py::arg("rot"),
           py::arg("trn"))
      .def(py::init([](py::array_t<ftype> rot, py::array_t<ftype> trn) {
             auto rotation = new_mat33_from_numpy(rot);
             auto translation = new_vec3_from_numpy(trn);
             return std::unique_ptr<RTop_frac>(
                 new RTop_frac(*rotation, *translation));
           }),
           py::arg("rot"), py::arg("trn"))
      .def("rtop_orth", &RTop_frac::rtop_orth, py::arg("cell"))
      .def("inverse", &RTop_frac::inverse)
      .def_static("identity", &RTop_frac::identity)
      .def_static("null", &RTop_frac::null)
      .def("__repr__", [](const RTop_frac &self) {
        return "<clipper.RTop_frac: Fractional operator class.>";
      });
}

void declare_symop(py::module &m) {
  py::class_<Symop, RTop_frac>(m, "Symop")
      .def(py::init<>())
      .def(py::init<const RTop<> &>(), py::arg("rtop"))
      // 4x4 array
      .def(py::init([](const py::array_t<ftype> a) {
             check_array_shape(a, {4, 4}, true);
             auto buf = a.request();
             return std::unique_ptr<Symop>(
                 new Symop(*reinterpret_cast<ftype(*)[4][4]>(buf.ptr)));
           }),
           py::arg("a"))
      .def("format", &Symop::format)
      .def("__str__", &Symop::format)
      .def("__repr__", [](const Symop &self) {
        return "<clipper.Symop: Crystallographic symmetry operator class.>";
      });
}

void declare_isymop(py::module &m) {
  py::class_<Isymop, RTop<int>>(m, "Isymop")
      .def(py::init<>())
      .def(py::init<const RTop<int> &>(), py::arg("rtop"))
      .def(py::init<const Symop &, const Grid &>(), py::arg("symop"),
           py::arg("grid"))
      .def("__repr__", [](const Isymop &self) {
        return "<clipper.Isymop: Integerised symmetry matrix class.>";
      });
}

void declare_symop_code(py::module &m) {
  py::class_<Symop_code>(m, "Symop_code")
      .def(py::init<>())
      .def(py::init<const int &>(), py::arg("code"))
      .def(py::init<const Symop &>(), py::arg("symop"))
      .def(py::init<const Isymop &>(), py::arg("isymop"))
      .def("init", &Symop_code::init, py::arg("isymop"))
      .def("code_rot", &Symop_code::code_rot)
      .def("code_trn", &Symop_code::code_trn)
      .def("symop", &Symop_code::symop)
      .def("isymop", &Symop_code::isymop)
      .def_static("identity", &Symop_code::identity)
      .def("__int__", [](const Symop_code &self) { return int(self); })
      .def("__repr__", [](const Symop_code &self) {
        return "<clipper.Symop_code: Compressed encoded symmetry operator "
               "class.>";
      });
}

void init_symop(py::module &m) {
  declare_rtopfrac(m);
  declare_symop(m);
  declare_isymop(m);
  declare_symop_code(m);
}
