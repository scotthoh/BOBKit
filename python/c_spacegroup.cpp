// Wrapper for clipper spacegroup
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/core/spacegroup.h>
#include <gemmi/symmetry.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace clipper;

void declare_spgr_descr(py::module &m) {
  py::class_<Spgr_descr> spgr_descr(m, "Spgr_descr");

  py::enum_<Spgr_descr::TYPE>(spgr_descr, "TYPE", "Spacegroup symbol type.")
      .value("Hall", Spgr_descr::TYPE::Hall)
      .value("HM", Spgr_descr::TYPE::HM)
      .value("XHM", Spgr_descr::TYPE::XHM)
      .value("Symops", Spgr_descr::TYPE::Symops)
      .value("Number", Spgr_descr::TYPE::Number)
      .value("Unknown", Spgr_descr::TYPE::Unknown);

  spgr_descr.def(py::init<>())
      .def(py::init<const String &, Spgr_descr::TYPE>(), py::arg("symbol"),
           py::arg("type") = Spgr_descr::TYPE::Unknown,
           "Constructor from symbol or operators.")
      .def(py::init<const int &>(), py::arg("num"), "Constructor from number.")
      .def(py::init<const Spgr_descr::Symop_codes &>(),
           py::arg("symop_codes_list"), "Constructor from symop list.")
      .def("spacegroup_number", &Spgr_descr::spacegroup_number,
           "Return spacegroup number.")
      .def("symbol_hall", &Spgr_descr::symbol_hall, "Return Hall symbol.")
      .def("symbol_hm", &Spgr_descr::symbol_hm, "Return H-M symbol.")
      .def("symbol_xhm", &Spgr_descr::symbol_xhm, "Return extended H-M symbol.")
      .def("symbol_hm_ext", &Spgr_descr::symbol_hm_ext,
           "Return extension H-M symbol.")
      .def_static("set_preferred", &Spgr_descr::set_preferred,
                  "Set preferred default spacegroup choice.")
      .def_property_readonly("generator_ops", &Spgr_descr::generator_ops,
                             "Return the generators for the spacegroup.")
      .def("hash", &Spgr_descr::hash,
           "Return the hash code for the spacegroup.")
      .def("__hash__", &Spgr_descr::hash)
      .def("__repr__",
           [](const Spgr_descr &self) {
             return "<clipper.Spgr_descr " + self.symbol_hm() + " >";
           })
      .doc() = "Spacegroup description.\nThe spacegroup description "
               "is a compact description of a spacegroup. It may be "
               "initialised from Hall or H-M symbols, a string of symops "
               "or a number. Internally a hash code is used to refer to "
               "the spacegroup, so this object is only 32 bits in size.";

  using SymopC = Spgr_descr::Symop_codes;
  py::class_<SymopC>(spgr_descr, "Symop_codes",
                     "Vector of symop codes and associated methods.")
      .def(py::init<>())
      .def("init_hall", &SymopC::init_hall, py::arg("symbol"),
           "Initialise from Hall symbol.")
      .def("init_symops", &SymopC::init_symops, py::arg("symops"),
           "Initialise from symops.")
      .def("expand", &SymopC::expand, "Expand (incomplete) list of symops.")
      .def("primitive_noninversion_ops", &SymopC::primitive_noninversion_ops,
           "Return primitive non-inversion ops (by computation).")
      .def("inversion_ops", &SymopC::inversion_ops,
           "Return inversion ops (by computation).")
      .def("primitive_ops", &SymopC::primitive_ops,
           "Return primitive incl. inversion ops (by computation).")
      .def("centering_ops", &SymopC::centering_ops,
           "Return lattice centering ops (by computation).")
      .def("laue_ops", &SymopC::laue_ops, "Return Laue ops.")
      .def("pgrp_ops", &SymopC::pgrp_ops, "Return point groups ops.")
      .def("patterson_ops", &SymopC::patterson_ops, "Return Patterson ops.")
      .def("generator_ops", &SymopC::generator_ops,
           "Return minimal list of generator ops.")
      .def("product", &SymopC::product, py::arg("symop_codes_list"),
           "Return product of this (expanded) list by another (expanded) list.")
      .def("hash", &SymopC::hash, "Return hash code of symop list.")
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

  // py::enum_<Spacegroup::TYPE>(
  //     spacegroup, "SGTYPE",
  //     "enumeration for fast construction of Null or P1 spacegroup.")
  //     .value("Null", Spacegroup::TYPE::Null)
  //     .value("P1", Spacegroup::TYPE::P1)
  //     .export_values();

  py::enum_<Spacegroup::AXIS>(spacegroup, "AXIS", "enumeration for cell axes.")
      .value("A", Spacegroup::AXIS::A)
      .value("B", Spacegroup::AXIS::B)
      .value("C", Spacegroup::AXIS::C)
      .export_values();

  spacegroup
      .def(py::init<>())
      //.def(py::init<Spacegroup::TYPE>(), "Fast constructor for Null or P1
      // spacegroup.")
      .def(py::init<const Spgr_descr &>(),
           "Constructor from spacegroup description.")
      .def_static("p1", &Spacegroup::p1, "Return P1 spacegroup.")
      .def_static("null", &Spacegroup::null, "Return Null spacegroup.")
      .def(
          "init", [](Spacegroup &self, const Spgr_descr &sd) { self.init(sd); },
          py::arg("spgr_descr"), "Initialiser from spacegroup")
      .def(
          "init",
          [](Spacegroup &self, const gemmi::SpaceGroup &sg) {
            self.init(GEMMI::spacegroup(sg).descr());
          },
          py::arg("spacegroup"), "Initialiser from gemmi spacegroup.")
      .def(
          "init",
          [](Spacegroup &self, const Spacegroup &sg) { self.init(sg.descr()); },
          py::arg("spacegroup"), "Initialiser from spacegroup.")
      .def("is_null", &Spacegroup::is_null,
           "Test if object has been initialised.")
      .def("descr", &Spacegroup::descr, "Return spacegroup description.")
      .def("num_symops", &Spacegroup::num_symops, "Get number of symops.")
      .def("num_primops", &Spacegroup::num_primops,
           "Get number of primitive symops (identical to "
           "num_primitive_symops()).")
      .def("num_primitive_symops", &Spacegroup::num_primitive_symops,
           "Get number of primitive symops(inc. identity and inversion).")
      .def("num_centering_symops", &Spacegroup::num_centering_symops,
           "Get number of centering symops(inc. identity).")
      .def("num_inversion_symops", &Spacegroup::num_inversion_symops,
           "Get number of inversion symops(inc. identity).")
      .def("num_primitive_noninversion_symops",
           &Spacegroup::num_primitive_noninversion_symops,
           "Get number of primitive non-inversion symops(inc. identity).")
      .def("symop", &Spacegroup::symop, py::arg("sym_no"), "Get n'th symop.")
      .def("primitive_symop", &Spacegroup::primitive_symop, py::arg("sym_no"),
           "Get n'th primitive symop (identival to symop(sym_no)).")
      .def("inversion_symop", &Spacegroup::inversion_symop, py::arg("sym_no"),
           "Get n'th inversion symop (0...1 max)")
      .def("centering_symop", &Spacegroup::centering_symop, py::arg("sym_no"),
           "Get n'th centering symop (0...3 max)")
      .def("order_of_symmetry_about_axis",
           &Spacegroup::order_of_symmetry_about_axis, py::arg("axis"),
           "Get the order of rotational symmetry about a given axis.")
      .def("hkl_class", &Spacegroup::hkl_class, py::arg("hkl"),
           "Get \'class\' of reflection: multiplicity, allowed phase, absence.")
      .def("recip_asu", &Spacegroup::recip_asu, py::arg("hkl"),
           "Test if hkl is in default reciprocal ASU.")
      .def("product_op", &Spacegroup::product_op, py::arg("s1"), py::arg("s2"),
           "Get symop number corresponding to the product of two symop.")
      .def("inverse_op", &Spacegroup::inverse_op, py::arg("s"),
           "Get symop number corresponding to the inverse of a symop.")
      .def_property_readonly("asu_max", &Spacegroup::asu_max,
                             "Get map ASU, upper bound.")
      .def_property_readonly("asu_min", &Spacegroup::asu_min,
                             "Get map ASU, lower bound.")
      .def_property_readonly("invariant_under_change_of_hand",
                             &Spacegroup::invariant_under_change_of_hand,
                             "Test if change of hand preserves spacegroup.")
      .def_property_readonly("symbol_laue", &Spacegroup::symbol_laue,
                             "Return Laue group symbol.")

      // from clipper-gemmi
      .def_static(
          "from_gemmi_spacegroup",
          [](const gemmi::SpaceGroup &sg) { return GEMMI::spacegroup(sg); },
          "Convert gemmi spacegroup to clipper object.")
      .def_static(
          "to_gemmi_spacegroup",
          [](const Spacegroup &sg) { return GEMMI::spacegroup(sg); },
          "Convert clipper spacegroup to gemmi object.")
      .def("__repr__",
           [](const Spacegroup &self) {
             return "<clipper.Spacegroup " + self.symbol_hm() + " >";
           })
      .def("__str__", &Spacegroup::symbol_hm)
      .def("debug", &Spacegroup::debug, "Output debug details.")
      .doc() = "Spacegroup object.\nThe spacegroup object is a full "
               "description of a spacegroup, including all the most "
               "regularly used information in an efficient form. It may "
               "be initialised from a clipper::Spgr_descr. This object.";
}

void init_spacegroup(py::module &m) {
  declare_spgr_descr(m);
  declare_spacegroup(m);
}