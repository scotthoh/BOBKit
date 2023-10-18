#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <clipper/clipper.h>
#include "type_conversions.h"
#include "helper_functions.h"

namespace py = pybind11;
using namespace clipper;
void init_coord_orth(py::module &m)
{
  py::class_<Resolution>(m, "Resolution")
      .def(py::init<>())
      .def(py::init<const ftype &>(), py::arg("resolution"))
      .def("init", &Resolution::init, py::arg("resolution"))
      .def("limit", &Resolution::limit)
      .def("invresolsq_limit", &Resolution::invresolsq_limit)
      .def("is_null", &Resolution::is_null)
      .def("__str__", [](const Resolution &self)
           { return String(self.limit(), 6, 4); })
      .def("__repr__", [](const Resolution &self)
           { return "<clipper.Resolution " + String(self.limit(), 6, 4) + "Å>"; });

  py::class_<HKL_class>(m, "HKL_class")
      .def(py::init<>())
      .def(py::init<const Spacegroup &, const HKL &>())
      .def("epsilon", &HKL_class::epsilon)
      .def("epsilonc", &HKL_class::epsilonc)
      .def("allowed", &HKL_class::allowed)
      .def("centric", &HKL_class::centric)
      .def("sys_abs", &HKL_class::sys_abs);

  py::class_<RTop_orth, RTop<>> rtop_orth(m, "RTop_orth");
  rtop_orth
      .def(py::init<>())
      .def(py::init<const RTop<> &>())
      .def(py::init<const Mat33<> &>())
      .def(py::init<const Mat33<> &, const Vec3<> &>(), py::arg("rot"), py::arg("trn"))
      .def(py::init<const std::vector<Coord_orth> &, const std::vector<Coord_orth> &>(),
           py::arg("src"), py::arg("tgt"))
      .def(py::init<const std::vector<Coord_orth> &, const std::vector<Coord_orth> &, const std::vector<ftype> &>(),
           py::arg("src"), py::arg("tgt"), py::arg("weight"))
      // need to bind atomlist for the atomlist type constructor
      .def(py::init<const Atom_list &, const Atom_list &>(), py::arg("src"), py::arg("tgt"))
      .def(py::init([](py::array_t<ftype> rot, py::array_t<ftype> trn)
                    {
        auto rotation = new_mat33_from_numpy(rot);
        auto translation = new_vec3_from_numpy(trn);
        return std::unique_ptr<RTop_orth>(new RTop_orth(*rotation, *translation)); }))
      .def("rtop_frac", &RTop_orth::rtop_frac)
      .def("inverse", &RTop_orth::inverse)
      .def("axis_coordinate_near", &RTop_orth::axis_coordinate_near)
      .def("screw_translation", &RTop_orth::screw_translation)
      .def_static("identity", &RTop_orth::identity)
      .def_static("null", &RTop_orth::null);

  py::class_<HKL, Vec3<int>> hkl(m, "HKL");
  hkl
      .def(py::init<>())
      .def(py::init<const Vec3<int> &>())
      .def(py::init<const int &, const int &, const int &>())
      .def_property(
          "h", [](const HKL &self) -> const int &
          { return self.h(); },
          [](HKL &self, const int &val)
          { self.h() = val; })
      .def_property(
          "k", [](const HKL &self) -> const int &
          { return self.k(); },
          [](HKL &self, const int &val)
          { self.k() = val; })
      .def_property(
          "l", [](const HKL &self) -> const int &
          { return self.l(); },
          [](HKL &self, const int &val)
          { self.l() = val; })
      .def_property(
          "hkl", [](const HKL &self)
          { return make_array_1d<HKL, int>(self, 3); },
          [](HKL &self, py::array_t<int> hkl)
          { fill_array_1d<HKL, int>(self, 3, hkl); })
      .def("invresolsq", &HKL::invresolsq)
      .def("coord_reci_frac", &HKL::coord_reci_frac)
      .def("coord_reci_orth", &HKL::coord_reci_orth)
      .def("transform", (HKL(HKL::*)(const Symop &) const) & HKL::transform)
      .def("transform", (HKL(HKL::*)(const Isymop &) const) & HKL::transform)
      .def("sym_phase_shift", &HKL::sym_phase_shift)
      .def("format", &HKL::format)
      .def("__repr__", [](const HKL &self)
           { return "<clipper." + self.format() + ">"; })
      .def("__neg__", [](const HKL &self)
           { return -self; })
      .def(
          "__add__", [](const HKL &self, const HKL &other)
          { return self + other; },
          py::is_operator())
      .def(
          "__sub__", [](const HKL &self, const HKL &other)
          { return self - other; },
          py::is_operator())
      .def(
          "__mul__", [](const HKL &self, const int &s)
          { return s * self; },
          py::is_operator())
      .def(
          "__rmul__", [](const HKL &self, const int &s)
          { return s * self; },
          py::is_operator())
      .def("__rmul__", [](const HKL &self, const Isymop &op)
           { return op * self; });

  py::class_<Coord_reci_orth, Vec3<>> coord_reci_orth(m, "Coord_reci_orth");
  coord_reci_orth
      .def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>())
      .def_property_readonly("xs", &Coord_reci_orth::xs)
      .def_property_readonly("ys", &Coord_reci_orth::ys)
      .def_property_readonly("zs", &Coord_reci_orth::zs)
      .def("invresolsq", &Coord_reci_orth::invresolsq)
      .def("coord_reci_frac", &Coord_reci_orth::coord_reci_frac)
      .def("transform", &Coord_reci_orth::transform)
      .def("format", &Coord_reci_orth::format)
      .def("__str__", &Coord_reci_orth::format);

  py::class_<Coord_reci_frac, Vec3<>> coord_reci_frac(m, "Coord_reci_frac");
  coord_reci_frac
      .def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>())
      .def(py::init<const HKL &>())
      .def("hkl", &Coord_reci_frac::hkl)
      .def("invresolsq", &Coord_reci_frac::invresolsq)
      .def_property_readonly("us", &Coord_reci_frac::us)
      .def_property_readonly("vs", &Coord_reci_frac::vs)
      .def_property_readonly("ws", &Coord_reci_frac::ws)
      .def("coord_reci_orth", &Coord_reci_frac::coord_reci_orth)
      .def("transform", &Coord_reci_frac::transform)
      .def("__str__", &Coord_reci_frac::format)
      .def("format", &Coord_reci_frac::format);

  py::class_<Coord_grid, Vec3<int>> coord_grid(m, "Coord_grid");
  coord_grid
      .def(py::init<>())
      .def(py::init<const Vec3<int> &>())
      .def(py::init<const int &, const int &, const int &>())
      .def(py::init<const Grid &, const int &>())
      .def_property(
          "u", [](const Coord_grid &self) -> const int &
          { return self.u(); },
          [](Coord_grid &self, const int &val)
          { self.u() = val; })
      .def_property(
          "v", [](const Coord_grid &self) -> const int &
          { return self.v(); },
          [](Coord_grid &self, const int &val)
          { self.v() = val; })
      .def_property(
          "w", [](const Coord_grid &self) -> const int &
          { return self.w(); },
          [](Coord_grid &self, const int &val)
          { self.w() = val; })
      .def("coord_map", &Coord_grid::coord_map)
      .def("coord_frac", &Coord_grid::coord_frac)
      .def("transform", &Coord_grid::transform)
      .def("unit", &Coord_grid::unit)
      .def("next", (const Coord_grid &(Coord_grid::*)(const Grid &)) & Coord_grid::next)
      .def("next", (const Coord_grid &(Coord_grid::*)(const Grid_range &)) & Coord_grid::next)
      .def("last", (bool(Coord_grid::*)(const Grid &) const) & Coord_grid::last)
      .def("last", (bool(Coord_grid::*)(const Grid_range &) const) & Coord_grid::last)
      .def("index", &Coord_grid::index)
      .def("deindex", &Coord_grid::deindex)
      .def("__str__", &Coord_grid::format)
      .def("format", &Coord_grid::format)
      .def(
          "__neg__", [](const Coord_grid &self)
          { return -self; },
          py::is_operator())
      .def(
          "__add__", [](const Coord_grid &self, const Coord_grid &other)
          { return self + other; },
          py::is_operator())
      .def(
          "__sub__", [](const Coord_grid &self, const Coord_grid &other)
          { return self - other; },
          py::is_operator())
      .def(
          "__mul__", [](const Coord_grid &self, const int &s)
          { return s * self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_grid &self, const int &s)
          { return s * self; },
          py::is_operator())
      .def(
          "__eq__", [](const Coord_grid &self, const Coord_grid &other)
          { return self == other; },
          py::is_operator())
      .def(
          "__ne__", [](const Coord_grid &self, const Coord_grid &other)
          { return self != other; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_grid &self, const Isymop &op)
          { return op * self; },
          py::is_operator());

  py::class_<Coord_orth, Vec3<>> coord_orth(m, "Coord_orth");
  coord_orth
      .def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>(), py::arg("x"), py::arg("y"), py::arg("z"))
      .def(py::init<const Coord_orth &, const Coord_orth &, const Coord_orth &, const ftype &, const ftype &, const ftype &>())
      .def(py::init([](const std::array<ftype, 3> &a)
                    { return std::unique_ptr<Coord_orth>(new Coord_orth(a[0], a[1], a[2])); }))
      .def_property_readonly("x", &Coord_orth::x)
      .def_property_readonly("y", &Coord_orth::y)
      .def_property_readonly("z", &Coord_orth::z)
      .def("lengthsq", &Coord_orth::lengthsq)
      .def("coord_frac", &Coord_orth::coord_frac)
      .def("transform", &Coord_orth::transform)
      .def("format", &Coord_orth::format)
      .def("__str__", &Coord_orth::format)
      .def("length", &Coord_orth::length)
      .def("angle", &Coord_orth::angle)
      .def("torsion", &Coord_orth::torsion)
      .def(
          "__neg__", [](const Coord_orth &self)
          { return -self; },
          py::is_operator())
      .def("__add__", [](const Coord_orth &self, const Coord_orth &other)
           { return self + other; })
      .def("__sub__", [](
                          const Coord_orth &self, const Coord_orth &other)
           { return self - other; })
      .def(
          "__mul__", [](const Coord_orth &self, const ftype &s)
          { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_orth &self, const ftype &s)
          { s *self; },
          py::is_operator())
      .def(
          "__rmul__", [](const Coord_orth &self, const RTop_orth &op)
          { return op * self; },
          py::is_operator())
      .def("__iter__", [](Coord_orth &self)
           { return py::make_iterator(&self[0], &self[2]); })
      .def(
          "to_numpy", [](const Coord_orth &self)
          { return to_numpy_1d<Coord_orth, ftype>(self, 3); },
          py::return_value_policy::reference_internal);

  py::class_<Coord_frac, Vec3<>> coord_frac(m, "Coord_frac");
  coord_frac
      .def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>())
      .def_property_readonly("u", &Coord_frac::u)
      .def_property_readonly("v", &Coord_frac::v)
      .def_property_readonly("w", &Coord_frac::w)
      .def("lengthsq", &Coord_frac::lengthsq)
      .def("coord_orth", &Coord_frac::coord_orth)
      .def("coord_map", &Coord_frac::coord_map)
      .def("coord_grd", &Coord_frac::coord_grid)
      .def("transform", &Coord_frac::transform)
      .def("lattice_copy_zero", &Coord_frac::lattice_copy_zero)
      .def("lattice_copy_unit", &Coord_frac::lattice_copy_unit)
      .def("lattice_copy_near", &Coord_frac::lattice_copy_near)
      .def("symmetry_copy_near", &Coord_frac::symmetry_copy_near)
      .def("__str__", &Coord_frac::format)
      .def("format", &Coord_frac::format)
      .def(
          "__neg__", [](const Coord_frac &self)
          { return -self; },
          py::is_operator())
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(ftype() * py::self)
      .def(RTop_frac() * py::self); // RTop_frac is from core/symop.h

  py::class_<Coord_map, Vec3<>> coord_map(m, "Coord_map");
  coord_map
      .def(py::init<>())
      .def(py::init<const Vec3<> &>())
      .def(py::init<const Coord_grid &>())
      .def(py::init<const ftype &, const ftype &, const ftype &>())
      .def("coord_frac", &Coord_map::coord_frac)
      .def("coord_grid", &Coord_map::coord_grid)
      .def("floor", &Coord_map::floor)
      .def("ceil", &Coord_map::ceil)
      .def_property_readonly("u", &Coord_map::u)
      .def_property_readonly("v", &Coord_map::v)
      .def_property_readonly("w", &Coord_map::w)
      .def("__str__", &Coord_map::format)
      .def("format", &Coord_map::format)
      .def(
          "__neg__", [](const Coord_map &self)
          { return -self; },
          py::is_operator())
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(ftype() * py::self);

  py::class_<U_aniso_orth, Mat33sym<>> uaniso_orth(m, "U_aniso_orth");
  uaniso_orth
      .def(py::init<>())
      .def(py::init<const Mat33sym<> &>())
      .def(py::init<const ftype &>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &, const ftype &, const ftype &>())
      .def("u_iso", &U_aniso_orth::u_iso)
      .def("u_aniso_frac", &U_aniso_orth::u_aniso_frac)
      .def("transform", &U_aniso_orth::transform)
      .def(py::self + py::self)
      .def(
          "__neg__", [](const U_aniso_orth &self)
          { return -self; },
          py::is_operator())
      .def(ftype() * py::self);

  py::class_<U_aniso_frac, Mat33sym<>> uaniso_frac(m, "U_aniso_frac");
  uaniso_frac
      .def(py::init<>())
      .def(py::init<const Mat33sym<> &>())
      .def(py::init<const ftype &, const ftype &, const ftype &, const ftype &, const ftype &, const ftype &>())
      .def("u_aniso_orth", &U_aniso_frac::u_aniso_orth)
      .def("transofrorm", &U_aniso_frac::transform)
      .def(py::self + py::self)
      .def(
          "__neg__", [](const U_aniso_frac &self)
          { return -self; },
          py::is_operator())
      .def(ftype() * py::self);

  py::class_<Grid, Vec3<int>> grid(m, "Grid");
  grid
      .def(py::init<>())
      .def(py::init<const int &, const int &, const int &>())
      .def_property_readonly("nu", &Grid::nu)
      .def_property_readonly("nv", &Grid::nv)
      .def_property_readonly("nw", &Grid::nw)
      .def("size", &Grid::size)
      .def("in_grid", &Grid::in_grid)
      .def("index", &Grid::index)
      .def("deindex", &Grid::deindex)
      .def("__str__", &Grid::format)
      .def("format", &Grid::format)
      .def("debug", &Grid::debug);

  py::class_<Grid_sampling, Grid> grid_sampling(m, "Grid_sampling");
  grid_sampling
      .def(py::init<>())
      .def(py::init<const int &, const int &, const int &>())
      .def(py::init<const Spacegroup &, const Cell &, const Resolution &, const ftype>())
      .def("init", &Grid_sampling::init)
      .def("matrix_grid_frac", &Grid_sampling::matrix_grid_frac)
      .def("matrix_frac_grid", &Grid_sampling::matrix_frac_grid)
      .def("is_null", &Grid_sampling::is_null);

  py::class_<HKL_sampling> hkl_sampling(m, "HKL_sampling");
  hkl_sampling
      .def(py::init<>())
      .def(py::init<const Cell &, const Resolution &>())
      .def("hkl_limit", &HKL_sampling::hkl_limit)
      .def("resolution", &HKL_sampling::resolution)
      .def("in_resolution", &HKL_sampling::in_resolution)
      .def("is_null", &HKL_sampling::is_null)
      .def("__str__", &HKL_sampling::format)
      .def("format", &HKL_sampling::format)
      .def(py::self == py::self);

  py::class_<Grid_range, Grid> grid_range(m, "Grid_range");
  grid_range
      .def(py::init<>())
      .def(py::init<const Coord_grid &, const Coord_grid &>())
      .def(py::init<const Grid &, const Coord_frac &, const Coord_frac &>())
      .def(py::init<const Cell &, const Grid &, const ftype &>())
      .def("min", &Grid_range::min)
      .def("max", &Grid_range::max)
      .def("add_border", &Grid_range::add_border)
      .def("in_grid", &Grid_range::in_grid)
      .def("index", &Grid_range::index)
      .def("deindex", &Grid_range::deindex);

  py::class_<Atom> atom(m, "Atom");
  atom
      .def(py::init<>())
      .def(py::init<const Atom &>())
      .def_property("element", &Atom::element, &Atom::set_element)
      .def_property("pos", &Atom::coord_orth, &Atom::set_coord_orth)
      .def_property("occupancy", &Atom::occupancy, &Atom::set_occupancy)
      .def_property("u_iso", &Atom::u_iso, &Atom::set_u_iso)
      .def_property(
          "u_aniso_orth",
          &Atom::u_aniso_orth,
          [](Atom &self, py::array_t<ftype> u)
          {
    std::vector<ftype> uval(6);
    fill_array_1d<std::vector<ftype>, ftype>(uval, 6, u);
    auto uaniso = U_aniso_orth(uval[0], uval[1], uval[2], uval[3], uval[4], uval[5]);
    self.set_u_aniso_orth(uaniso); })
      .def_property(
          "b_iso",
          [](const Atom &self)
          { return clipper::Util::u2b(self.u_iso()); },
          [](Atom &self, const ftype &value)
          { self.set_u_iso(clipper::Util::b2u(value)); })
      .def("x", [](const Atom &self)
           { return self.coord_orth().x(); })
      .def("y", [](const Atom &self)
           { return self.coord_orth().y(); })
      .def("z", [](const Atom &self)
           { return self.coord_orth().z(); })
      .def("transform", &Atom::transform)
      .def("is_null", &Atom::is_null)
      .def_static("null", &Atom::null)
      .def("__repr__", [](const Atom &self)
           { return "<clipper.Atom " + self.element().trim() + " " + self.coord_orth().format() + ">"; });

  py::class_<Atom_list> atomlist(m, "Atom_list", "test");
  atomlist
      .def(py::init<>())
      .def(py::init<const std::vector<Atom> &>())
      //.def(py::init<const T &>())
      .def("__repr__", [](const Atom_list &self)
           {  std::stringstream stream;
           stream << "<clipper.Atom_list containin " << self.size() << " atom(s).>";
           return stream.str(); })
      .def("__len__", [](const Atom_list &self)
           { return self.size(); })
      .def("size", [](const Atom_list &self)
           { return self.size(); })
      .def("__getitem__", [](const Atom_list &self, const int &i) -> Atom
           { return self.at(normalise_index(i, self.size())); })
      .def("__setitem__", [](Atom_list &self, const int &i, const Atom &atom)
           { self.at(normalise_index(i, self.size())) = atom; })
      .def(
          "delete_atom", [](Atom_list &self, const int &pos)
          { delitem_at_index(self, pos); },
          py::arg("index"))
      .def(
          "add_atom", [](Atom_list &self, const int &pos, const Atom &atom)
          { add_item(self, atom, pos); },
          py::arg("index"), py::arg("Atom"))
      .def(
          "__iter__", [](Atom_list &self)
          { return py::make_iterator(&self[0], &self[self.size()]); },
          py::keep_alive<0, 1>());

  // probably need something more;
}
