#include <clipper/clipper-gemmi.h>
#include <clipper/core/coords.h>
#include <gemmi/math.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// #include "coord.h"
#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
namespace py = pybind11;
using namespace clipper;

void declare_atom(py::module &m) {
  py::class_<Atom> atom(m, "Atom");
  atom.def(py::init<>())
      .def(py::init<const Atom &>())
      .def_property("element", &Atom::element, &Atom::set_element)
      .def_property("pos", &Atom::coord_orth, &Atom::set_coord_orth)
      .def_property("occupancy", &Atom::occupancy, &Atom::set_occupancy)
      .def_property("u_iso", &Atom::u_iso, &Atom::set_u_iso)
      .def_property("u_aniso_orth", &Atom::u_aniso_orth,
                    [](Atom &self, py::array_t<ftype> u) {
                      std::vector<ftype> uval(6);
                      fill_array_1d<std::vector<ftype>, ftype>(uval, 6, u);
                      auto uaniso = U_aniso_orth(uval[0], uval[1], uval[2],
                                                 uval[3], uval[4], uval[5]);
                      self.set_u_aniso_orth(uaniso);
                    })
      .def_property(
          "b_iso",
          [](const Atom &self) { return clipper::Util::u2b(self.u_iso()); },
          [](Atom &self, const ftype &value) {
            self.set_u_iso(clipper::Util::b2u(value));
          })
      .def("x", [](const Atom &self) { return self.coord_orth().x(); })
      .def("y", [](const Atom &self) { return self.coord_orth().y(); })
      .def("z", [](const Atom &self) { return self.coord_orth().z(); })
      .def("transform", &Atom::transform)
      .def("is_null", &Atom::is_null)
      .def_static("null", &Atom::null)
      .def("__repr__", [](const Atom &self) {
        return "<clipper.Atom " + self.element().trim() + " " +
               self.coord_orth().format() + ">";
      });

  py::class_<Atom_list> atomlist(m, "Atom_list");
  atomlist.def(py::init<>())
      .def(py::init<const std::vector<Atom> &>())
      //.def(py::init<const T &>())
      .def("clear", [](Atom_list &self) { self.clear(); })
      .def("pop_back", [](Atom_list &self) { self.pop_back(); })
      .def("push_back",
           [](Atom_list &self, Atom &atom) { self.push_back(atom); })
      .def("__repr__",
           [](const Atom_list &self) {
             std::stringstream stream;
             stream << "<clipper.Atom_list containing " << self.size()
                    << " atom(s).>";
             return stream.str();
           })
      .def("__len__", [](const Atom_list &self) { return self.size(); })
      .def("size", [](const Atom_list &self) { return self.size(); })
      .def(
          "__getitem__",
          [](const Atom_list &self, const int &i) -> Atom {
            return self.at(i);
            // return self.at(normalise_index(i, self.size()));
          },
          "Returns a copy of Atom at given index.")
      .def(
          "__setitem__",
          [](Atom_list &self, const int &i, const Atom &atom) {
            self.at(normalise_index(i, self.size())) = atom;
          },
          "Replaces Atom at given index.")
      .def(
          "delete_atom",
          [](Atom_list &self, const int &pos) { delitem_at_index(self, pos); },
          py::arg("index"))
      .def(
          "add_atom",
          [](Atom_list &self, const int &pos, const Atom &atom) {
            add_item(self, atom, pos);
          },
          py::arg("index"), py::arg("Atom"))
      .def(
          "__iter__",
          [](Atom_list &self) {
            return py::make_iterator(&self[0], &self[self.size()]);
          },
          py::keep_alive<0, 1>());
}

void declare_coord_orth(py::module &m) {
  py::class_<Resolution>(m, "Resolution")
      .def(py::init<>())
      .def(py::init<const ftype &>(), py::arg("resolution"))
      .def("init", &Resolution::init, py::arg("resolution"))
      .def("limit", &Resolution::limit)
      .def("invresolsq_limit", &Resolution::invresolsq_limit)
      .def_static("test_precision64",
                  []() { return ftype64(1.8501448145505965); })
      .def_static("test_precision32",
                  []() { return ftype32(1.8501448145505965); })
      .def("is_null", &Resolution::is_null)
      .def("__str__",
           [](const Resolution &self) { return String(self.limit(), 6, 4); })
      .def("__repr__", [](const Resolution &self) {
        return "<clipper.Resolution " + String(self.limit(), 6, 4) + " Å.>";
      });

  py::class_<HKL_class>(m, "HKL_class")
      .def(py::init<>())
      .def(py::init<const Spacegroup &, const HKL &>(), py::arg("spacegroup"),
           py::arg("hkl"))
      .def(py::init([](const Spacegroup &sg, const py::array_t<int> hkl) {
        check_array_shape(hkl, {3}, true);
        return std::unique_ptr<HKL_class>(
            new HKL_class(sg, HKL(hkl.at(0), hkl.at(1), hkl.at(2))));
      }))
      .def("epsilon", &HKL_class::epsilon)
      .def("epsilonc", &HKL_class::epsilonc)
      .def("allowed", &HKL_class::allowed)
      .def("centric", &HKL_class::centric)
      .def("sys_abs", &HKL_class::sys_abs)
      .def("__repr__", [](const HKL_class &self) {
        return "<clipper.HKL_class: Describes the type of reflection in a "
               "given spacegroup.>";
      });

  py::class_<RTop_orth, RTop<>> rtop_orth(m, "RTop_orth");
  rtop_orth.def(py::init<>())
      .def(py::init<const RTop<> &>())
      .def(py::init<const Mat33<> &>())
      .def(py::init<const Mat33<> &, const Vec3<> &>(), py::arg("rot"),
           py::arg("trn"))
      .def(py::init<const std::vector<Coord_orth> &,
                    const std::vector<Coord_orth> &>(),
           py::arg("src"), py::arg("tgt"))
      .def(py::init<const std::vector<Coord_orth> &,
                    const std::vector<Coord_orth> &,
                    const std::vector<ftype> &>(),
           py::arg("src"), py::arg("tgt"), py::arg("weight"))
      .def(py::init<const Atom_list &, const Atom_list &>(), py::arg("src"),
           py::arg("tgt"))
      .def(py::init([](py::array_t<ftype> rot, py::array_t<ftype> trn) {
             auto rotation = new_mat33_from_numpy(rot);
             auto translation = new_vec3_from_numpy(trn);
             return std::unique_ptr<RTop_orth>(
                 new RTop_orth(*rotation, *translation));
           }),
           py::arg("rot"), py::arg("trn"))
      .def_static(
          "from_gemmi_transform",
          [](const gemmi::Transform &rtop) { return GEMMI::transform(rtop); })
      .def_static("to_gemmi_transform",
                  [](const RTop_orth &rtop) { return GEMMI::transform(rtop); })
      .def("rtop_frac", &RTop_orth::rtop_frac, py::arg("cell"))
      .def("inverse", &RTop_orth::inverse)
      .def("axis_coordinate_near", &RTop_orth::axis_coordinate_near,
           py::arg("centre"))
      .def("screw_translation", &RTop_orth::screw_translation)
      .def_static("identity", &RTop_orth::identity)
      .def_static("null", &RTop_orth::null)
      .def("__repr__",
           [](const RTop_orth &self) { return "<clipper.RTop_orth class.>"; });

  py::class_<HKL, Vec3<int>> hkl(m, "HKL");
  hkl.def(py::init<>())
      .def(py::init<const Vec3<int> &>())
      .def(py::init<const int &, const int &, const int &>())
      .def(py::init([](const py::array_t<int> &hkl) {
        check_array_shape(hkl, {3}, true);
        return std::unique_ptr<HKL>(new HKL(hkl.at(0), hkl.at(1), hkl.at(2)));
      }))
      .def("init",
           [](HKL &self, const gemmi::Miller &hkl) {
             self.h() = hkl.at(0);
             self.k() = hkl.at(1);
             self.l() = hkl.at(2);
           })
      .def_static("from_gemmi_Miller",
                  [](const gemmi::Miller &hkl) { return GEMMI::Hkl(hkl); })
      .def_static("to_gemmi_Miller",
                  [](const HKL &hkl) { return GEMMI::Hkl(hkl); })
      .def_property(
          "h", [](const HKL &self) -> const int & { return self.h(); },
          [](HKL &self, const int &val) { self.h() = val; })
      .def_property(
          "k", [](const HKL &self) -> const int & { return self.k(); },
          [](HKL &self, const int &val) { self.k() = val; })
      .def_property(
          "l", [](const HKL &self) -> const int & { return self.l(); },
          [](HKL &self, const int &val) { self.l() = val; })
      .def_property(
          "hkl",
          [](const HKL &self) { return make_array_1d<HKL, int>(self, 3); },
          [](HKL &self, py::array_t<int> hkl) {
            fill_array_1d<HKL, int>(self, 3, hkl);
          })
      .def("invresolsq", &HKL::invresolsq, py::arg("cell"))
      .def("coord_reci_frac", &HKL::coord_reci_frac)
      .def("coord_reci_orth", &HKL::coord_reci_orth, py::arg("cell"))
      .def("transform", (HKL(HKL::*)(const Symop &) const) & HKL::transform,
           py::arg("op"))
      .def("transform", (HKL(HKL::*)(const Isymop &) const) & HKL::transform,
           py::arg("op"))
      .def("sym_phase_shift", &HKL::sym_phase_shift, py::arg("op"))
      .def("format", &HKL::format)
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

  py::class_<Coord_reci_orth, Vec3<>> coord_reci_orth(m, "Coord_reci_orth");
  coord_reci_orth.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"))
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("xs"), py::arg("ys"), py::arg("zs"))
      .def(py::init([](const py::array_t<ftype> &c) {
        check_array_shape(c, {3}, true);
        return std::unique_ptr<Coord_reci_orth>(
            new Coord_reci_orth(c.at(0), c.at(1), c.at(2)));
      }))
      .def_property_readonly("xs", &Coord_reci_orth::xs)
      .def_property_readonly("ys", &Coord_reci_orth::ys)
      .def_property_readonly("zs", &Coord_reci_orth::zs)
      .def("invresolsq", &Coord_reci_orth::invresolsq)
      .def("coord_reci_frac", &Coord_reci_orth::coord_reci_frac,
           py::arg("cell"))
      .def("transform", &Coord_reci_orth::transform, py::arg("op"))
      .def("format", &Coord_reci_orth::format)
      .def("__str__", &Coord_reci_orth::format);

  py::class_<Coord_reci_frac, Vec3<>> coord_reci_frac(m, "Coord_reci_frac");
  coord_reci_frac.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"))
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("us"), py::arg("vs"), py::arg("ws"))
      .def(py::init([](const py::array_t<ftype> &c) {
        check_array_shape(c, {3}, true);
        return std::unique_ptr<Coord_reci_frac>(
            new Coord_reci_frac(c.at(0), c.at(1), c.at(2)));
      }))
      .def(py::init<const HKL &>(), py::arg("hkl"))
      .def("hkl", &Coord_reci_frac::hkl)
      .def("invresolsq", &Coord_reci_frac::invresolsq, py::arg("cell"))
      .def_property_readonly("us", &Coord_reci_frac::us)
      .def_property_readonly("vs", &Coord_reci_frac::vs)
      .def_property_readonly("ws", &Coord_reci_frac::ws)
      .def("coord_reci_orth", &Coord_reci_frac::coord_reci_orth,
           py::arg("cell"))
      .def("transform", &Coord_reci_frac::transform, py::arg("op"))
      .def("__str__", &Coord_reci_frac::format)
      .def("format", &Coord_reci_frac::format);

  py::class_<Coord_grid, Vec3<int>> coord_grid(m, "Coord_grid");
  coord_grid.def(py::init<>())
      .def(py::init<const Vec3<int> &>())
      .def(py::init<const int &, const int &, const int &>(), py::arg("u"),
           py::arg("v"), py::arg("w"))
      .def(py::init<const Grid &, const int &>(), py::arg("grid"),
           py::arg("index"))
      .def_property(
          "u", [](const Coord_grid &self) -> const int & { return self.u(); },
          [](Coord_grid &self, const int &val) { self.u() = val; })
      .def_property(
          "v", [](const Coord_grid &self) -> const int & { return self.v(); },
          [](Coord_grid &self, const int &val) { self.v() = val; })
      .def_property(
          "w", [](const Coord_grid &self) -> const int & { return self.w(); },
          [](Coord_grid &self, const int &val) { self.w() = val; })
      .def("coord_map", &Coord_grid::coord_map)
      .def("coord_frac", &Coord_grid::coord_frac, py::arg("grid_sampling"))
      .def("transform", &Coord_grid::transform, py::arg("op"))
      .def("unit", &Coord_grid::unit, py::arg("grid_sampling"))
      .def("next",
           (const Coord_grid &(Coord_grid::*)(const Grid &)) & Coord_grid::next,
           py::arg("grid"))
      .def("next",
           (const Coord_grid &(Coord_grid::*)(const Grid_range &)) &
               Coord_grid::next,
           py::arg("grid_range"))
      .def("last", (bool(Coord_grid::*)(const Grid &) const) & Coord_grid::last,
           py::arg("grid"))
      .def("last",
           (bool(Coord_grid::*)(const Grid_range &) const) & Coord_grid::last,
           py::arg("grid_range"))
      .def("index", &Coord_grid::index, py::arg("grid"))
      .def("deindex", &Coord_grid::deindex, py::arg("grid"), py::arg("index"))
      .def("__str__", &Coord_grid::format)
      .def("format", &Coord_grid::format)
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

  py::class_<Coord_orth, Vec3<>> coord_orth(m, "Coord_orth");
  coord_orth.def(py::init<>())
      .def(py::init<const Vec3<> &>(), py::arg("v"))
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("x"), py::arg("y"), py::arg("z"))
      .def(py::init<const Coord_orth &, const Coord_orth &, const Coord_orth &,
                    const ftype &, const ftype &, const ftype &>(),
           py::arg("x1"), py::arg("x2"), py::arg("x3"), py::arg("length"),
           py::arg("angle"), py::arg("torsion"))
      .def(py::init([](const py::array_t<ftype> &a) {
        check_array_shape(a, {3}, true);
        return std::unique_ptr<Coord_orth>(
            new Coord_orth(a.at(0), a.at(1), a.at(2)));
      }))
      .def_property_readonly("x", &Coord_orth::x)
      .def_property_readonly("y", &Coord_orth::y)
      .def_property_readonly("z", &Coord_orth::z)
      .def("lengthsq", &Coord_orth::lengthsq)
      .def("coord_frac", &Coord_orth::coord_frac, py::arg("cell"))
      .def("transform", &Coord_orth::transform, py::arg("op"))
      .def("format", &Coord_orth::format)
      .def("__str__", &Coord_orth::format)
      .def_static("length", &Coord_orth::length, py::arg("x1"), py::arg("x2"))
      .def_static("angle", &Coord_orth::angle, py::arg("x1"), py::arg("x2"),
                  py::arg("x3"))
      .def_static("torsion", &Coord_orth::torsion, py::arg("x1"), py::arg("x2"),
                  py::arg("x3"), py::arg("x4"))
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
          py::is_operator())
      .def("__iter__",
           [](Coord_orth &self) {
             return py::make_iterator(&self[0], &self[2]);
           })
      .def(
          "to_numpy",
          [](const Coord_orth &self) {
            return to_numpy_1d<Coord_orth, ftype>(self, 3);
          },
          py::return_value_policy::reference_internal);
  // py::bind_vector<std::vector<Coord_orth>>(m, "CoordList");
  py::class_<Coord_frac, Vec3<>> coord_frac(m, "Coord_frac");
  coord_frac.def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>())
      .def(py::init([](const py::array_t<ftype> &a) {
        check_array_shape(a, {3}, true);
        return std::unique_ptr<Coord_frac>(
            new Coord_frac(a.at(0), a.at(1), a.at(2)));
      }))
      .def_property_readonly("u", &Coord_frac::u)
      .def_property_readonly("v", &Coord_frac::v)
      .def_property_readonly("w", &Coord_frac::w)
      .def("lengthsq", &Coord_frac::lengthsq, py::arg("cell"))
      .def("coord_orth", &Coord_frac::coord_orth, py::arg("cell"))
      .def("coord_map", &Coord_frac::coord_map, py::arg("grid"))
      .def("coord_grd", &Coord_frac::coord_grid, py::arg("grid"))
      .def("transform", &Coord_frac::transform, py::arg("op"))
      .def("lattice_copy_zero", &Coord_frac::lattice_copy_zero)
      .def("lattice_copy_unit", &Coord_frac::lattice_copy_unit)
      .def("lattice_copy_near", &Coord_frac::lattice_copy_near)
      .def("symmetry_copy_near", &Coord_frac::symmetry_copy_near,
           py::arg("spacegroup"), py::arg("cell"), py::arg("cf"))
      .def("__str__", &Coord_frac::format)
      .def("format", &Coord_frac::format)
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

  py::class_<Coord_map, Vec3<>> coord_map(m, "Coord_map");
  coord_map.def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const Coord_grid &>(), py::arg("c"))
      .def(py::init<const ftype &, const ftype &, const ftype &>(),
           py::arg("u"), py::arg("v"), py::arg("w"))
      .def("coord_frac", &Coord_map::coord_frac, py::arg("grid"))
      .def("coord_grid", &Coord_map::coord_grid)
      .def("floor", &Coord_map::floor)
      .def("ceil", &Coord_map::ceil)
      .def_property_readonly("u", &Coord_map::u)
      .def_property_readonly("v", &Coord_map::v)
      .def_property_readonly("w", &Coord_map::w)
      .def("__str__", &Coord_map::format)
      .def("format", &Coord_map::format)
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
      .def(py::init<const Mat33sym<> &>(), py::arg("m"))
      .def(py::init<const ftype &>(), py::arg("u"))
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("u11"), py::arg("u22"), py::arg("u33"), py::arg("u12"),
           py::arg("u13"), py::arg("u23"))
      .def("u_iso", &U_aniso_orth::u_iso)
      .def("u_aniso_frac", &U_aniso_orth::u_aniso_frac, py::arg("cell"))
      .def("transform", &U_aniso_orth::transform, py::arg("op"))
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
          py::is_operator());

  py::class_<U_aniso_frac, Mat33sym<>> uaniso_frac(m, "U_aniso_frac");
  uaniso_frac.def(py::init<>())
      .def(py::init<const Mat33sym<> &>(), py::arg("m"))
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &,
                    const ftype &, const ftype &>(),
           py::arg("u11"), py::arg("u22"), py::arg("u33"), py::arg("u12"),
           py::arg("u13"), py::arg("u23"))
      .def("u_aniso_orth", &U_aniso_frac::u_aniso_orth, py::arg("cell"))
      .def("transform", &U_aniso_frac::transform, py::arg("op"))
      .def(py::self + py::self)
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
          py::is_operator());

  py::class_<Grid, Vec3<int>> grid(m, "Grid");
  grid.def(py::init<>())
      .def(py::init<const int &, const int &, const int &>(), py::arg("nu"),
           py::arg("nv"), py::arg("nw"))
      .def_property_readonly("nu", &Grid::nu)
      .def_property_readonly("nv", &Grid::nv)
      .def_property_readonly("nw", &Grid::nw)
      .def("size", &Grid::size)
      .def("in_grid", &Grid::in_grid, py::arg("grid"))
      .def("index", &Grid::index, py::arg("grid"))
      .def("deindex", &Grid::deindex, py::arg("index"))
      .def("__str__", &Grid::format)
      .def("format", &Grid::format)
      .def("debug", &Grid::debug);

  py::class_<Grid_sampling, Grid> grid_sampling(m, "Grid_sampling");
  grid_sampling.def(py::init<>())
      .def(py::init<const int &, const int &, const int &>(), py::arg("nu"),
           py::arg("nv"), py::arg("nw"))
      .def(py::init<const Spacegroup &, const Cell &, const Resolution &,
                    const ftype>(),
           py::arg("spacegroup"), py::arg("cell"), py::arg("resolution"),
           py::arg("rate") = 1.5)
      .def("init", &Grid_sampling::init, py::arg("spacegroup"), py::arg("cell"),
           py::arg("resolution"), py::arg("rate") = 1.5)
      .def("matrix_grid_frac", &Grid_sampling::matrix_grid_frac)
      .def("matrix_frac_grid", &Grid_sampling::matrix_frac_grid)
      .def("is_null", &Grid_sampling::is_null);

  py::class_<HKL_sampling> hkl_sampling(m, "HKL_sampling");
  hkl_sampling.def(py::init<>())
      .def(py::init<const Cell &, const Resolution &>(), py::arg("cell"),
           py::arg("resolution"))
      .def("hkl_limit", &HKL_sampling::hkl_limit)
      .def("resolution", &HKL_sampling::resolution, py::arg("cell"))
      .def("in_resolution", &HKL_sampling::in_resolution, py::arg("hkl"))
      .def("is_null", &HKL_sampling::is_null)
      .def("__str__", &HKL_sampling::format)
      .def("format", &HKL_sampling::format)
      .def(
          "__eq__",
          [](const HKL_sampling &self, const HKL_sampling &other) {
            return self == other;
          },
          py::is_operator());

  py::class_<Grid_range, Grid> grid_range(m, "Grid_range");
  grid_range.def(py::init<>())
      .def(py::init<const Coord_grid &, const Coord_grid &>(), py::arg("min"),
           py::arg("max"))
      .def(py::init<const Grid &, const Coord_frac &, const Coord_frac &>(),
           py::arg("grid"), py::arg("min"), py::arg("max"))
      .def(py::init<const Cell &, const Grid &, const ftype &>(),
           py::arg("cell"), py::arg("grid"), py::arg("radius"))
      .def("min", &Grid_range::min)
      .def("max", &Grid_range::max)
      .def("add_border", &Grid_range::add_border, py::arg("b"))
      .def("in_grid", &Grid_range::in_grid, py::arg("grid"))
      .def("index", &Grid_range::index, py::arg("c"))
      .def("deindex", &Grid_range::deindex, py::arg("index"));

  // probably need something more;
}

void init_coord_orth(py::module &m) {
  declare_atom(m);
  declare_coord_orth(m);
}