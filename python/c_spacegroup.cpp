#include <clipper/clipper-gemmi.h>
#include <clipper/core/spacegroup.h>
#include <gemmi/symmetry.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_spgr_descr(py::module &m) {
  py::class_<Spgr_descr> spgr_descr(m, "Spgr_descr");

  py::enum_<Spgr_descr::TYPE>(spgr_descr, "TYPE")
      .value("Hall", Spgr_descr::TYPE::Hall)
      .value("HM", Spgr_descr::TYPE::HM)
      .value("XHM", Spgr_descr::TYPE::XHM)
      .value("Symops", Spgr_descr::TYPE::Symops)
      .value("Number", Spgr_descr::TYPE::Number)
      .value("Unknown", Spgr_descr::TYPE::Unknown);

  spgr_descr.def(py::init<>())
      .def(py::init<const String &, Spgr_descr::TYPE>(), py::arg("symbol"),
           py::arg("type") = Spgr_descr::TYPE::Unknown)
      .def(py::init<const int &>(), py::arg("num"))
      .def(py::init<const Spgr_descr::Symop_codes &>(),
           py::arg("symop_codes_list"))
      .def("spacegroup_number", &Spgr_descr::spacegroup_number)
      .def("symbol_hall", &Spgr_descr::symbol_hall)
      .def("symbol_hm", &Spgr_descr::symbol_hm)
      .def("symbol_xhm", &Spgr_descr::symbol_xhm)
      .def("symbol_hm_ext", &Spgr_descr::symbol_hm_ext)
      .def_static("set_preferred", &Spgr_descr::set_preferred)
      .def_property_readonly("generator_ops", &Spgr_descr::generator_ops)
      .def("hash", &Spgr_descr::hash)
      .def("__hash__", &Spgr_descr::hash)
      .def("__repr__", [](const Spgr_descr &self) {
        return "<clipper.Spgr_descr " + self.symbol_hm() + " >";
      });

  using SymopC = Spgr_descr::Symop_codes;
  py::class_<SymopC>(spgr_descr, "Symop_codes")
      .def(py::init<>())
      .def("init_hall", &SymopC::init_hall, py::arg("symbol"))
      .def("init_symops", &SymopC::init_symops, py::arg("symops"))
      .def("expand", &SymopC::expand)
      .def("primitive_noninversion_ops", &SymopC::primitive_noninversion_ops)
      .def("inversion_ops", &SymopC::inversion_ops)
      .def("primitive_ops", &SymopC::primitive_ops)
      .def("centering_ops", &SymopC::centering_ops)
      .def("laue_ops", &SymopC::laue_ops)
      .def("pgrp_ops", &SymopC::pgrp_ops)
      .def("patterson_ops", &SymopC::patterson_ops)
      .def("generator_ops", &SymopC::generator_ops)
      .def("product", &SymopC::product, py::arg("symop_codes_list"))
      .def("hash", &SymopC::hash)
      .def("__hash__", &SymopC::hash)
      .def("__repr__",
           [](const SymopC &self) {
             return "<clipper.Symop_codes containing " +
                    clipper::String(int(self.size())) +
                    " compressed encoded symmetry operator(s).>";
           })
      // need getitem
      //  iterators
      .def("__getitem__",
           [](const SymopC &self, const int &i) { return self[i]; })
      .def(
          "__iter__",
          [](SymopC &self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def(
          "__reversed__",
          [](SymopC &self) {
            return py::make_iterator(self.rbegin(), self.rend());
          },
          py::keep_alive<0, 1>())
      // modifiers, some std::vector methods
      .def(
          "append",
          [](SymopC &self, const clipper::Symop_code &c) { self.push_back(c); },
          "Append item to the end.")
      .def(
          "clear", [](SymopC &self) { self.clear(); }, "Clear list.")
      .def(
          "insert",
          [](SymopC &self, const int &pos, const clipper::Symop_code &c) {
            add_item(self, c, pos);
          },
          "Insert item at given index.")
      .def("pop",
           [](SymopC &self) {
             auto ret = self.back();
             self.pop_back();
             return ret;
           })
      .def(
          "remove",
          [](SymopC &self, const int &pos) { delitem_at_index(self, pos); },
          "Delete item at given index.")
      .def(
          "remove",
          [](SymopC &self, const int &start, const int &end) {
            delitem_range(self, start, end);
          },
          "Delete items at range of indices given [start, end).");
}

// void declare_symop_codes(py::module &m) {}

void declare_spacegroup(py::module &m) {
  py::class_<Spacegroup, Spgr_descr> spacegroup(m, "Spacegroup");

  py::enum_<Spacegroup::TYPE>(spacegroup, "SGTYPE")
      .value("Null", Spacegroup::TYPE::Null)
      .value("P1", Spacegroup::TYPE::P1);

  py::enum_<Spacegroup::AXIS>(spacegroup, "AXIS")
      .value("A", Spacegroup::AXIS::A)
      .value("B", Spacegroup::AXIS::B)
      .value("C", Spacegroup::AXIS::C);

  spacegroup.def(py::init<>())
      .def(py::init<Spacegroup::TYPE>())
      .def(py::init<const Spgr_descr &>())
      .def(
          "init", [](Spacegroup &self, const Spgr_descr &sd) { self.init(sd); },
          py::arg("spgr_descr"))
      .def(
          "init",
          [](Spacegroup &self, const gemmi::SpaceGroup &sg) {
            self.init(GEMMI::spacegroup(sg).descr());
          },
          py::arg("spacegroup"))
      .def(
          "init",
          [](Spacegroup &self, const Spacegroup &sg) { self.init(sg.descr()); },
          py::arg("spacegroup"))
      .def("is_null", &Spacegroup::is_null)
      .def("descr", &Spacegroup::descr)
      .def("num_symops", &Spacegroup::num_symops)
      .def("num_primops", &Spacegroup::num_primops)
      .def("num_primitive_symops", &Spacegroup::num_primitive_symops)
      .def("num_centering_symops", &Spacegroup::num_centering_symops)
      .def("num_inversion_symops", &Spacegroup::num_inversion_symops)
      .def("num_primitive_noninversion_symops",
           &Spacegroup::num_primitive_noninversion_symops)
      .def("symop", &Spacegroup::symop, py::arg("sym_no"))
      .def("primitive_symop", &Spacegroup::primitive_symop, py::arg("sym_no"))
      .def("inversion_symop", &Spacegroup::inversion_symop, py::arg("sym_no"))
      .def("centering_symop", &Spacegroup::centering_symop, py::arg("sym_no"))
      .def("order_of_symmetry_about_axis",
           &Spacegroup::order_of_symmetry_about_axis, py::arg("axis"))
      .def("hkl_class", &Spacegroup::hkl_class, py::arg("hkl"))
      .def("recip_asu", &Spacegroup::recip_asu, py::arg("hkl"))
      .def("product_op", &Spacegroup::product_op, py::arg("s1"), py::arg("s2"))
      .def("inverse_op", &Spacegroup::inverse_op, py::arg("s"))
      .def_property_readonly("asu_max", &Spacegroup::asu_max)
      .def_property_readonly("asu_min", &Spacegroup::asu_min)
      .def_property_readonly("invariant_under_change_of_hand",
                             &Spacegroup::invariant_under_change_of_hand)
      .def_property_readonly("symbol_laue", &Spacegroup::symbol_laue)
      .def_static("p1", &Spacegroup::p1)
      .def_static("null", &Spacegroup::null)
      // from clipper-gemmi
      .def_static(
          "from_gemmi_spacegroup",
          [](const gemmi::SpaceGroup &sg) { return GEMMI::spacegroup(sg); })
      .def_static("to_gemmi_spacegroup",
                  [](const Spacegroup &sg) { return GEMMI::spacegroup(sg); })
      .def("__repr__",
           [](const Spacegroup &self) {
             return "<clipper.Spacegroup " + self.symbol_hm() + " >";
           })
      .def("__str__", &Spacegroup::symbol_hm);
}

void init_spacegroup(py::module &m) {
  declare_spgr_descr(m);
  declare_spacegroup(m);
}