// Wrapper for clipper coord
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <gemmi/math.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void init_coord_orth(py::module &m) {
  py::class_<Resolution>(m, "Resolution")
      .def(py::init<>())
      .def(py::init<const ftype &>(), py::arg("resolution"),
           "Constructor from ftype.")
      .def("init", &Resolution::init, py::arg("resolution"),
           "Initialiser from ftype.")
      .def("limit", &Resolution::limit, "Get resolution limit.")
      .def("invresolsq_limit", &Resolution::invresolsq_limit,
           "Get inverse resolution square limit.")
      .def("is_null", &Resolution::is_null,
           "Test if value has been initialised.")
      .def("__str__",
           [](const Resolution &self) { return String(self.limit(), 6, 4); })
      .def("__repr__",
           [](const Resolution &self) {
             return "<clipper.Resolution " + String(self.limit(), 6, 4) +
                    " Å.>";
           })
      .doc() =
      "Resolution in angstroms\n"
      "This object represents a resolution limit which will be used for "
      "all aspects of a calculation. This is a base for a donor type.";

  py::class_<HKL_class>(m, "HKL_class")
      .def(py::init<>())
      .def(py::init<const Spacegroup &, const HKL &>(), py::arg("spacegroup"),
           py::arg("hkl"), "Constructor from spacegroup and HKL")
      .def(py::init([](const Spacegroup &sg, const py::array_t<int> hkl) {
             check_array_shape(hkl, {3}, true);
             return std::unique_ptr<HKL_class>(
                 new HKL_class(sg, HKL(hkl.at(0), hkl.at(1), hkl.at(2))));
           }),
           "Constructor from spacegroup and HKL(list or numpy array)")
      .def("epsilon", &HKL_class::epsilon, "Get epsilon.")
      .def("epsilonc", &HKL_class::epsilonc,
           "Get epsilon for acentric, 2x epsilon for centric.")
      .def("allowed", &HKL_class::allowed, "Get allowed phase.")
      .def("centric", &HKL_class::centric, "Is centric?")
      .def("sys_abs", &HKL_class::sys_abs, "Is sys abs?")
      .def("__repr__",
           [](const HKL_class &self) {
             return "<clipper.HKL_class: Describes the type of reflection in a "
                    "given spacegroup.>";
           })
      .doc() =
      "Reflection class\n"
      "This describes the type of a reflection in a given spacegroup, "
      "including centricity, systematic absence, phase restriction, and "
      "multiplicity.";

  py::class_<RTop_orth, RTop<>> rtop_orth(m, "RTop_orth");
  rtop_orth.def(py::init<>())
      .def(py::init<const RTop<> &>(), py::arg("rtop"),
           "Constructor: copy/convert.")
      .def(py::init<const Mat33<> &>(), py::arg("rot"),
           "Constructor from rotation.")
      .def(py::init<const Mat33<> &, const Vec3<> &>(), py::arg("rot"),
           py::arg("trn"), "Constructor from rotation and translation")
      .def(py::init<const std::vector<Coord_orth> &,
                    const std::vector<Coord_orth> &>(),
           py::arg("src"), py::arg("tgt"),
           "Constructor from two lists of Coord_orth.")
      .def(py::init<const std::vector<Coord_orth> &,
                    const std::vector<Coord_orth> &,
                    const std::vector<ftype> &>(),
           py::arg("src"), py::arg("tgt"), py::arg("weight"),
           "Constructor from two lists of Coord_orth and weight.")
      .def(py::init<const Atom_list &, const Atom_list &>(), py::arg("src"),
           py::arg("tgt"), "Constructor from two Atom_list type objects.")
      .def(py::init([](py::array_t<ftype> rot, py::array_t<ftype> trn) {
             auto rotation = numpy_to_mat33(rot);
             auto translation = numpy_to_vec3(trn);
             return std::unique_ptr<RTop_orth>(
                 new RTop_orth(*rotation, *translation));
           }),
           py::arg("rot"), py::arg("trn"),
           "Constructor from rotation and translation(list/arrays).")
      .def_static(
          "from_gemmi_transform",
          [](const gemmi::Transform &rtop) { return GEMMI::transform(rtop); },
          py::arg("rtop"), "gemmi::Transform-RTop_orth conversion.")
      .def_static(
          "to_gemmi_transform",
          [](const RTop_orth &rtop) { return GEMMI::transform(rtop); },
          py::arg("rtop"), "RTop_orth-gemmi::Transform conversion.")
      .def("rtop_frac", &RTop_orth::rtop_frac, py::arg("cell"),
           "Orthogonal-fractional conversion.")
      .def("inverse", &RTop_orth::inverse, "Inverse operator.")
      .def("axis_coordinate_near", &RTop_orth::axis_coordinate_near,
           py::arg("centre"),
           "Return point on axis near the specified coordinate.")
      .def("screw_translation", &RTop_orth::screw_translation,
           "Return screw translation.")
      .def_static("identity", &RTop_orth::identity, "Return identity operator.")
      .def_static("null", &RTop_orth::null,
                  "Return null (uninitialised) operator.")
      .def("__repr__",
           [](const RTop_orth &self) { return "<clipper.RTop_orth class.>"; })
      .doc() =
      "Orthogonal operator class.\n"
      "This class is used for any RT-operator which operates on orthogonal "
      "coordinates. For a full list of methods, see clipper::RTop";

  py::class_<HKL, Vec3<int>> hkl(m, "HKL", "Reflection 'Miller index.");
  hkl.def(py::init<>())
      .def(py::init<const Vec3<int> &>(), "Constructor copy/convert.")
      .def(py::init<const int &, const int &, const int &>(), py::arg("h"),
           py::arg("k"), py::arg("l"), "Constructor from H,K,L.")
      .def(py::init([](const py::array_t<int> &hkl) {
             check_array_shape(hkl, {3}, true);
             return std::unique_ptr<HKL>(
                 new HKL(hkl.at(0), hkl.at(1), hkl.at(2)));
           }),
           py::arg("hkl"), "Constructor from list/array.")
      .def(py::init([](HKL &self, const gemmi::Miller &hkl) {
             return std::unique_ptr<HKL>(
                 new HKL(hkl.at(0), hkl.at(1), hkl.at(2)));
           }),
           "Constructor from gemmi::Miller.")
      .def_static(
          "from_gemmi_Miller",
          [](const gemmi::Miller &hkl) { return GEMMI::Hkl(hkl); },
          py::arg("hkl"), "gemmi::Miller-HKL conversion.")
      .def_static(
          "to_gemmi_Miller", [](const HKL &hkl) { return GEMMI::Hkl(hkl); },
          py::arg("hkl"), "HKL-gemmi::Miller conversion.")
      .def_property(
          "h", [](const HKL &self) -> const int & { return self.h(); },
          [](HKL &self, const int &val) { self.h() = val; }, "Get/set h.")
      .def_property(
          "k", [](const HKL &self) -> const int & { return self.k(); },
          [](HKL &self, const int &val) { self.k() = val; }, "Get/set k.")
      .def_property(
          "l", [](const HKL &self) -> const int & { return self.l(); },
          [](HKL &self, const int &val) { self.l() = val; }, "Get/set l.")
      .def_property(
          "hkl",
          [](const HKL &self) { return make_array_1d<HKL, int>(self, 3); },
          [](HKL &self, py::array_t<int> hkl) {
            fill_array_1d<HKL, int>(self, 3, hkl);
          },
          "Get/set hkl.")
      .def("invresolsq", &HKL::invresolsq, py::arg("cell"),
           "Return inverse resolution squared for this reflection in given "
           "cell.")
      .def("coord_reci_frac", &HKL::coord_reci_frac,
           "Return fractional reciprocal coordinate (i.e. non-integer HKL)")
      .def("coord_reci_orth", &HKL::coord_reci_orth, py::arg("cell"),
           "Orthogonal-fractional reciprocal space coordinate conversion")
      .def("transform", (HKL(HKL::*)(const Symop &) const) & HKL::transform,
           py::arg("symop"), "Return transformed hkl.")
      .def("transform", (HKL(HKL::*)(const Isymop &) const) & HKL::transform,
           py::arg("isymop"), "Return transformed hkl.")
      .def("sym_phase_shift", &HKL::sym_phase_shift, py::arg("op"),
           "Return symmerty phase shift for this HKL under op.")
      .def("format", &HKL::format, "Return formatted String representation.")
      .def("__repr__",
           [](const HKL &self) { return "<clipper." + self.format() + ">"; })
      .def("__neg__", [](const HKL &self) { return -self; })
      .def(
          "__add__",
          [](const HKL &self, const HKL &other) { return self + other; },
          py::is_operator())
      .def(
          "__sub__",
          [](const HKL &self, const HKL &other) { return self - other; },
          py::is_operator())
      .def(
          "__mul__", [](const HKL &self, const int &s) { return s * self; },
          py::is_operator())
      .def(
          "__rmul__", [](const HKL &self, const int &s) { return s * self; },
          py::is_operator())
      .def("__rmul__",
           [](const HKL &self, const Isymop &op) { return op * self; });

  py::class_<Coord_reci_orth, Vec3<>> coord_reci_orth(
      m, "Coord_reci_orth",
      "Orthogonal reciprocal coordinate (length of which is invresolsq)");
  coord_reci_orth.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"),
           "Constructor copy/convert.")
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("xs"), py::arg("ys"), py::arg("zs"),
           "Constructor from x*, y*, z*.")
      .def(py::init([](const py::array_t<ftype> &c) {
             check_array_shape(c, {3}, true);
             return std::unique_ptr<Coord_reci_orth>(
                 new Coord_reci_orth(c.at(0), c.at(1), c.at(2)));
           }),
           "Constructor from list/array of x*,y*,z*.")
      .def_property_readonly("xs", &Coord_reci_orth::xs, "Get x*")
      .def_property_readonly("ys", &Coord_reci_orth::ys, "Get y*")
      .def_property_readonly("zs", &Coord_reci_orth::zs, "Get z*")
      .def("invresolsq", &Coord_reci_orth::invresolsq,
           "Return inverse resolution squared for this coordinate.")
      .def("coord_reci_frac", &Coord_reci_orth::coord_reci_frac,
           py::arg("cell"),
           "Orthogonal-fractional reciprocal space coordinate conversion")
      .def("transform", &Coord_reci_orth::transform, py::arg("op"),
           "Return transformed coordinate.")
      .def("format", &Coord_reci_orth::format,
           "Return formatted string representation.")
      .def("__str__", &Coord_reci_orth::format);

  py::class_<Coord_reci_frac, Vec3<>> coord_reci_frac(
      m, "Coord_reci_frac",
      "Fractional reciprocal coordinate (i.e. non-integer hkl)");
  coord_reci_frac.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"),
           "Constructor copy/convert.")
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("us"), py::arg("vs"), py::arg("ws"),
           "Constructor from u*,v*,w*.")
      .def(py::init([](const py::array_t<ftype> &c) {
             check_array_shape(c, {3}, true);
             return std::unique_ptr<Coord_reci_frac>(
                 new Coord_reci_frac(c.at(0), c.at(1), c.at(2)));
           }),
           "Constructor from list/array of u*,v*,w*.")
      .def(py::init<const HKL &>(), py::arg("hkl"), "Constructor from HKL.")
      .def("hkl", &Coord_reci_frac::hkl, "Round to HKL.")
      .def("invresolsq", &Coord_reci_frac::invresolsq, py::arg("cell"),
           "Return inverse resolution squared for this reflection in given "
           "cell.")
      .def_property_readonly("us", &Coord_reci_frac::us, "Get u*")
      .def_property_readonly("vs", &Coord_reci_frac::vs, "Get v*")
      .def_property_readonly("ws", &Coord_reci_frac::ws, "Get w*")
      .def("coord_reci_orth", &Coord_reci_frac::coord_reci_orth,
           py::arg("cell"),
           "Fractional-orthogonal reciprocal space coordinate conversion")
      .def("transform", &Coord_reci_frac::transform, py::arg("op"),
           "Returned transformed coordinate.")
      .def("__str__", &Coord_reci_frac::format)
      .def("format", &Coord_reci_frac::format,
           "Return formatted string representation.");

  py::class_<Coord_grid, Vec3<int>> coord_grid(m, "Coord_grid",
                                               "Grid coordinate.");
  coord_grid.def(py::init<>())
      .def(py::init<const Vec3<int> &>(), "Constructor copy/convert.")
      .def(py::init<const int &, const int &, const int &>(), py::arg("u"),
           py::arg("v"), py::arg("w"), "Constructor from u,v,w")
      .def(py::init<const Grid &, const int &>(), py::arg("grid"),
           py::arg("index"),
           "Constructor from a grid and an index in that grid.")
      .def_property(
          "u", [](const Coord_grid &self) -> const int & { return self.u(); },
          [](Coord_grid &self, const int &val) { self.u() = val; },
          "Get/set u.")
      .def_property(
          "v", [](const Coord_grid &self) -> const int & { return self.v(); },
          [](Coord_grid &self, const int &val) { self.v() = val; },
          "Get/set v.")
      .def_property(
          "w", [](const Coord_grid &self) -> const int & { return self.w(); },
          [](Coord_grid &self, const int &val) { self.w() = val; },
          "Get/set w.")
      .def("coord_map", &Coord_grid::coord_map, "Convert to Coord_map.")
      .def("coord_frac", &Coord_grid::coord_frac, py::arg("grid_sampling"),
           "Convert to Coord_fract using given Grid_sampling.")
      .def("transform", &Coord_grid::transform, py::arg("op"),
           "Returned transformed coordinate.")
      .def("unit", &Coord_grid::unit, py::arg("grid_sampling"),
           "Reduced to unit box: (0..nu-1, 0..nv-1, 0..nw-1)")
      .def("next",
           (const Coord_grid &(Coord_grid::*)(const Grid &)) & Coord_grid::next,
           py::arg("grid"), "Increment in storage order (see index()).")
      .def("next",
           (const Coord_grid &(Coord_grid::*)(const Grid_range &)) &
               Coord_grid::next,
           py::arg("grid_range"), "Increment in storage order (see index()).")
      .def("last", (bool(Coord_grid::*)(const Grid &) const) & Coord_grid::last,
           py::arg("grid"), "Test if done in storage order (see index()).")
      .def("last",
           (bool(Coord_grid::*)(const Grid_range &) const) & Coord_grid::last,
           py::arg("grid_range"),
           "Test if done in storage order (see index()).")
      .def("index", &Coord_grid::index, py::arg("grid"),
           "Grid indexing operator.")
      .def("deindex", &Coord_grid::deindex, py::arg("grid"), py::arg("index"),
           "Grid deindexing operator.")
      .def("__str__", &Coord_grid::format)
      .def("format", &Coord_grid::format,
           "Return formatted string representation.")
      .def(
          "__neg__", [](const Coord_grid &self) { return -self; },
          py::is_operator())
      .def(
          "__add__",
          [](const Coord_grid &self, const Coord_grid &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const Coord_grid &self, const Coord_grid &other) {
            return self - other;
          },
          py::is_operator())
      .def(
          "__mul__",
          [](const Coord_grid &self, const int &s) { return s * self; },
          py::is_operator())
      .def(
          "__rmul__",
          [](const Coord_grid &self, const int &s) { return s * self; },
          py::is_operator())
      .def(
          "__eq__",
          [](const Coord_grid &self, const Coord_grid &other) {
            return self == other;
          },
          py::is_operator())
      .def(
          "__ne__",
          [](const Coord_grid &self, const Coord_grid &other) {
            return self != other;
          },
          py::is_operator())
      .def(
          "__rmul__",
          [](const Coord_grid &self, const Isymop &op) { return op * self; },
          py::is_operator());

  py::class_<Coord_orth, Vec3<>> coord_orth(
      m, "Coord_orth", "Orthogonal (Angstrom) coordinates");
  coord_orth.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"),
           "Constructor copy/convert.")
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("x"), py::arg("y"), py::arg("z"), "Constructor from x,y,z.")
      .def(py::init<const Coord_orth &, const Coord_orth &, const Coord_orth &,
                    const ftype &, const ftype &, const ftype &>(),
           py::arg("x1"), py::arg("x2"), py::arg("x3"), py::arg("length"),
           py::arg("angle"), py::arg("torsion"),
           "Constructor: from 3 coords and bond length, angle, torsion")
      .def(py::init([](const py::array_t<ftype> &a) {
             check_array_shape(a, {3}, true);
             return std::unique_ptr<Coord_orth>(
                 new Coord_orth(a.at(0), a.at(1), a.at(2)));
           }),
           "Constructor from list/array of x,y,z.")
      .def_property_readonly("x", &Coord_orth::x, "Get x.")
      .def_property_readonly("y", &Coord_orth::y, "Get y.")
      .def_property_readonly("z", &Coord_orth::z, "Get z.")
      .def("lengthsq", &Coord_orth::lengthsq,
           "Return square of length of vector in Angstroms.")
      .def("coord_frac", &Coord_orth::coord_frac, py::arg("cell"),
           "Orthogonal-fraction coordinate conversion.")
      .def("transform", &Coord_orth::transform, py::arg("op"),
           "Return transformed coordinate.")
      .def("format", &Coord_orth::format,
           "Return formatted string representation.")
      .def("__str__", &Coord_orth::format)
      .def_static("length", &Coord_orth::length, py::arg("x1"), py::arg("x2"),
                  "Return length of vector between two coordinates.")
      .def_static("angle", &Coord_orth::angle, py::arg("x1"), py::arg("x2"),
                  py::arg("x3"),
                  "Return angle between three orthogonal coordinates.")
      .def_static("torsion", &Coord_orth::torsion, py::arg("x1"), py::arg("x2"),
                  py::arg("x3"), py::arg("x4"),
                  "Return torsion betwee four orthogonal coordinates.")
      .def(
          "__neg__", [](const Coord_orth &self) { return -self; },
          py::is_operator())
      .def(
          "__add__",
          [](const Coord_orth &self, const Coord_orth &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const Coord_orth &self, const Coord_orth &other) {
            return self - other;
          },
          py::is_operator())
      .def(
          "__mul__", [](const Coord_orth &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_orth &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__",
          [](const Coord_orth &self, const RTop_orth &op) { return op * self; },
          py::is_operator());
  // inherited from Vec3: __iter__, as_array, from_array, buffer_protocol

  py::class_<Coord_frac, Vec3<>> coord_frac(m, "Coord_frac",
                                            "Fractional (cell) coordinates");
  coord_frac.def(py::init<>())
      .def(py::init<const Vec3<> &>(), "Constructor copy/convert.")
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           "Constructor from u,v,w.")
      .def(py::init([](const py::array_t<ftype> &a) {
             check_array_shape(a, {3}, true);
             return std::unique_ptr<Coord_frac>(
                 new Coord_frac(a.at(0), a.at(1), a.at(2)));
           }),
           "Constructor from list/array of u,v,w.")
      .def_property_readonly("u", &Coord_frac::u, "Get u.")
      .def_property_readonly("v", &Coord_frac::v, "Get v.")
      .def_property_readonly("w", &Coord_frac::w, "Get w.")
      .def("lengthsq", &Coord_frac::lengthsq, py::arg("cell"),
           "Return square of length of vector in Angstrom.")
      .def("coord_orth", &Coord_frac::coord_orth, py::arg("cell"),
           "Fractional-orthogonal coordinate conversion.")
      .def("coord_map", &Coord_frac::coord_map, py::arg("grid"),
           "Fractional-map coordinate conversion.")
      .def("coord_grid", &Coord_frac::coord_grid, py::arg("grid"),
           "Fractional-grid coordinate conversion.")
      .def("transform", &Coord_frac::transform, py::arg("op"),
           "Return transformed coordinate")
      .def("lattice_copy_zero", &Coord_frac::lattice_copy_zero,
           "Return lattice copy nearest origin.")
      .def("lattice_copy_unit", &Coord_frac::lattice_copy_unit,
           "Return lattice copy in unit box (0...1,0...1,0...1).")
      .def("lattice_copy_near", &Coord_frac::lattice_copy_near,
           "Return lattice copy near the specified coordinate.")
      .def("symmetry_copy_near", &Coord_frac::symmetry_copy_near,
           py::arg("spacegroup"), py::arg("cell"), py::arg("cf"),
           "Return symmetry copy `near the specified coordinate.")
      .def("__str__", &Coord_frac::format)
      .def("format", &Coord_frac::format,
           "Return formatted string representation.")
      .def(
          "__neg__", [](const Coord_frac &self) { return -self; },
          py::is_operator())
      .def("__add__", [](const Coord_frac &self,
                         const Coord_frac &other) { return self + other; })
      .def("__sub__", [](const Coord_frac &self,
                         const Coord_frac &other) { return self - other; })
      .def(
          "__mul__", [](const Coord_frac &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_frac &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__",
          [](const Coord_frac &self, const RTop_frac &op) { return op * self; },
          py::is_operator()); // RTop_frac is from core/symop.h

  py::class_<Coord_map, Vec3<>> coord_map(
      m, "Coord_map",
      "Map coordinate.\nThis is like Coord_grid, but non-integer.");
  coord_map.def(py::init<>())
      .def(py::init<const Vec3<> &>(), "Constructor copy/convert.")
      .def(py::init<const Coord_grid &>(), py::arg("c"),
           "Constructor from Coord_grid.")
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("u"), py::arg("v"), py::arg("w"), "Constructor from u,v,w.")
      .def("coord_frac", &Coord_map::coord_frac, py::arg("grid"),
           "Grid-fractional coordinate conversion.")
      .def("coord_grid", &Coord_map::coord_grid,
           "Return interger Coord_grid nearest to this coordinate.")
      .def("floor", &Coord_map::floor,
           "Return integer Coord_grid below this coordinate.")
      .def("ceil", &Coord_map::ceil,
           "Return integer Coord_grid above this coordinate.")
      .def_property_readonly("u", &Coord_map::u, "Get u.")
      .def_property_readonly("v", &Coord_map::v, "Get v.")
      .def_property_readonly("w", &Coord_map::w, "Get w.")
      .def("__str__", &Coord_map::format)
      .def("format", &Coord_map::format,
           "Return formatted string representation.")
      .def(
          "__neg__", [](const Coord_map &self) { return -self; },
          py::is_operator())
      .def("__add__", [](const Coord_map &self,
                         const Coord_map &other) { return self + other; })
      .def("__sub__", [](const Coord_map &self,
                         const Coord_map &other) { return self - other; })
      .def(
          "__mul__", [](const Coord_map &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_map &self, const ftype &s) { s *self; },
          py::is_operator());

  py::class_<U_aniso_orth, Mat33sym<>> uaniso_orth(m, "U_aniso_orth");
  uaniso_orth.def(py::init<>())
      .def(py::init<const Mat33sym<> &>(), py::arg("m"),
           "Constructor from Mat33sym.")
      .def(py::init<const ftype &>(), py::arg("u"),
           "Constructor from isotropic U.")
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("u11"), py::arg("u22"), py::arg("u33"), py::arg("u12"),
           py::arg("u13"), py::arg("u23"), "Constructor from Uij.")
      .def("u_iso", &U_aniso_orth::u_iso, "Return nearest isotropic U.")
      .def("u_aniso_frac", &U_aniso_orth::u_aniso_frac, py::arg("cell"),
           "Orthogonal-fractional conversion.")
      .def("transform", &U_aniso_orth::transform, py::arg("op"),
           "Return transformed U_aniso.")
      .def(
          "__neg__", [](const U_aniso_orth &self) { return -self; },
          py::is_operator())
      .def(
          "__add__",
          [](const U_aniso_orth &self, const U_aniso_orth &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__mul__", [](const U_aniso_orth &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const U_aniso_orth &self, const ftype &s) { s *self; },
          py::is_operator())
      .doc() = "Anisotropic orthogonal atomic displacement parameters.\n"
               "These are defined on orthogonal atomic coordinates in "
               "A<sup>-2</sup>, i.e. they are anisotropic U values.";

  py::class_<U_aniso_frac, Mat33sym<>> uaniso_frac(m, "U_aniso_frac");
  uaniso_frac.def(py::init<>())
      .def(py::init<const Mat33sym<> &>(), py::arg("m"),
           "Constructor from Mat33sym.")
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("u11"), py::arg("u22"), py::arg("u33"), py::arg("u12"),
           py::arg("u13"), py::arg("u23"), "Constructor from Uij.")
      .def("u_aniso_orth", &U_aniso_frac::u_aniso_orth, py::arg("cell"),
           "Fractional-orthogonal conversion.")
      .def("transform", &U_aniso_frac::transform, py::arg("op"),
           "Return transformed U_aniso.")
      //.def(py::self + py::self)
      .def(
          "__neg__", [](const U_aniso_frac &self) { return -self; },
          py::is_operator())
      .def(
          "__add__",
          [](const U_aniso_frac &self, const U_aniso_frac &other) {
            return self + other;
          },
          py::is_operator())
      .def(
          "__mul__", [](const U_aniso_frac &self, const ftype &s) { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const U_aniso_frac &self, const ftype &s) { s *self; },
          py::is_operator())
      .doc() = "Anisotropic fractional atomic displacement parameters.\n"
               "These are defined on fractional atomic coordinates in "
               "A<sup>-2</sup>, i.e. they are anisotropic U values.";

  py::class_<Grid, Vec3<int>> grid(m, "Grid");
  grid.def(py::init<>())
      .def(py::init<const int &, const int &, const int &>(), py::arg("nu"),
           py::arg("nv"), py::arg("nw"), "Constructor from nu,nv,nw.")
      .def_property_readonly("nu", &Grid::nu, "Get nu.")
      .def_property_readonly("nv", &Grid::nv, "Get nv.")
      .def_property_readonly("nw", &Grid::nw, "Get nw.")
      .def("size", &Grid::size, "Return size of grid array.")
      .def("in_grid", &Grid::in_grid, py::arg("grid"),
           "Determin if a point is in the grid.")
      .def("index", &Grid::index, py::arg("grid"), "Grid indexing operator.")
      .def("deindex", &Grid::deindex, py::arg("index"),
           "Grid deindexing operator.")
      .def("__str__", &Grid::format)
      .def("format", &Grid::format, "Return formatted string representation.")
      .def("debug", &Grid::debug, "Output debug details.")
      .doc() = "Generic grid.\nThis holds the dimensions of a 3D array, "
               "indexed from 0 along each dimension.";

  py::class_<Grid_sampling, Grid> grid_sampling(m, "Grid_sampling");
  grid_sampling.def(py::init<>())
      .def(py::init<const int &, const int &, const int &>(), py::arg("nu"),
           py::arg("nv"), py::arg("nw"), "Constructor from nu,nv,nw.")
      .def(py::init<const Spacegroup &, const Cell &, const Resolution &,
                    const ftype>(),
           py::arg("spacegroup"), py::arg("cell"), py::arg("resolution"),
           py::arg("rate") = 1.5,
           "Constructor from spacegroup, cell, resolution, Shannon rate.")
      .def("init", &Grid_sampling::init, py::arg("spacegroup"), py::arg("cell"),
           py::arg("resolution"), py::arg("rate") = 1.5,
           "Initialiser from spacegroup, cell, resolution, Shannon rate.")
      .def("matrix_grid_frac", &Grid_sampling::matrix_grid_frac,
           "Return matrix which converts grid to fractional coordinates.")
      .def("matrix_frac_grid", &Grid_sampling::matrix_frac_grid,
           "Return matrix which converts fractional to grid coordinates.")
      .def("is_null", &Grid_sampling::is_null,
           "Test if object has been initialised.")
      .doc() =
      "Grid sampling of a unit cell.\nThis class represents the grid "
      "sampling of a unit cell. It is otherwise identical to its parent, "
      "clipper::Grid_cell, but has an additional constructor which takes a "
      "spacegroup, cell and resolution and produces an appropriate grid "
      "obeying all of the symmetry constraints, and using efficient "
      "factors for the calculation of FFTs.";

  py::class_<HKL_sampling> hkl_sampling(m, "HKL_sampling");
  hkl_sampling.def(py::init<>())
      .def(py::init<const Cell &, const Resolution &>(), py::arg("cell"),
           py::arg("resolution"),
           "Constructor from cell (normal or inverse) and resolution.")
      .def("hkl_limit", &HKL_sampling::hkl_limit,
           "Return limiting values of H, K, L.")
      .def("resolution", &HKL_sampling::resolution, py::arg("cell"),
           "Return approximate resolution given cell.")
      .def("in_resolution", &HKL_sampling::in_resolution, py::arg("hkl"),
           "Test if a reflection is within the resolution limit.")
      .def("is_null", &HKL_sampling::is_null,
           "Test if object has been initialised. ")
      .def("__str__", &HKL_sampling::format)
      .def("format", &HKL_sampling::format,
           "Return formatted string representation.")
      .def(
          "__eq__",
          [](const HKL_sampling &self, const HKL_sampling &other) {
            return self == other;
          },
          py::is_operator())
      .doc() = "HKL sampling of reciprocal space.\nThe HKL_sampling class "
               "uniquely describes a P0 reflection list bounded by some "
               "resolution limit in reciprocal space. It is described in "
               "terms of large integers, and so immune from rounding errors "
               "once the object is constructed.";

  py::class_<Grid_range, Grid> grid_range(m, "Grid_range");
  grid_range.def(py::init<>())
      .def(py::init<const Coord_grid &, const Coord_grid &>(), py::arg("min"),
           py::arg("max"), "Constructor from grid limits (Coord_grid)")
      .def(py::init<const Grid &, const Coord_frac &, const Coord_frac &>(),
           py::arg("grid"), py::arg("min"), py::arg("max"),
           "Constructor from cell grid and fractional limits.")
      .def(
          py::init<const Cell &, const Grid &, const ftype &>(),
          py::arg("cell"), py::arg("grid"), py::arg("radius"),
          "Constructor: make grid to hold a sphere from cell, grid and radius.")
      .def("min", &Grid_range::min, "Access grid limits.")
      .def("max", &Grid_range::max, "Access grid limits.")
      .def("add_border", &Grid_range::add_border, py::arg("b"),
           "Border: increase grid to include given border.")
      .def("in_grid", &Grid_range::in_grid, py::arg("grid"),
           "Determine if a point is in the grid.")
      .def("index", &Grid_range::index, py::arg("c"), "Grid indexing operator.")
      .def("deindex", &Grid_range::deindex, py::arg("index"),
           "Grid deindexing operator.")
      .doc() = "Grid range class: defines array limits for a grid.\nThis class "
               "is used for describing 3D grids covering an arbitrary "
               "part of the 3D space, i.e. which do not start from (0,0,0).";

  // probably need something more;
}