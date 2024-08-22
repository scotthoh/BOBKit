// bind maps stuff
#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <clipper/core/xmap.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// template <class T> std::string display_repr(const Xmap<T> xmapclass) {
template <class T> std::string display_repr(const Xmap<T> xmap) {
  if (!xmap.is_null()) {
    if (strcmp(typeid(T).name(), "f") == 0)
      return "<clipper.Xmap_float: ASU grid " + xmap.grid_asu().format() + ">";
    else if (strcmp(typeid(T).name(), "d") == 0)
      return "<clipper.Xmap_double: ASU grid " + xmap.grid_asu().format() + ">";
    else if (strcmp(typeid(T).name(), "i") == 0)
      return "<clipper.Xmap_int: ASU grid " + xmap.grid_asu().format() + ">";
    else
      return "<clipper.Xmap>";
  } else
    return "<clipper.Xmap uninitialised.>";
}

void declare_xmap_base(py::module &m) {
  py::class_<Xmap_base> xmap_base(m, "Xmap_base");
  py::enum_<Xmap_base::FFTtype>(xmap_base, "FFTtype", "FFT backend selection.")
      .value("Default", Xmap_base::FFTtype::Default)
      .value("Normal", Xmap_base::FFTtype::Normal)
      .value("Sparse", Xmap_base::FFTtype::Sparse)
      .export_values();

  xmap_base.def("is_null", &Xmap_base::is_null)
      .def_property_readonly("cell", &Xmap_base::cell, "Get the cell.")
      .def_property_readonly("spacegroup", &Xmap_base::spacegroup,
                             "Get the spacegroup.")
      .def_property_readonly("grid_sampling", &Xmap_base::grid_sampling,
                             "Get cell grid.")
      .def_property_readonly("grid_asu", &Xmap_base::grid_asu, "Get ASU grid.")
      .def("coord_of", &Xmap_base::coord_of, py::arg("index"),
           "Map coordinate from index.")
      .def("index_of", &Xmap_base::index_of, py::arg("c"),
           "Map index from coordinate.")
      .def("to_map_unit", &Xmap_base::to_map_unit, py::arg("pos"),
           "Function to pick right cell repeat for any grid coord.")
      .def_property_readonly("operator_orth_grid",
                             &Xmap_base::operator_orth_grid,
                             "Return the orthogonal-to-grid coordinate "
                             "operator (translation is zero).")
      .def_property_readonly("operator_grid_orth",
                             &Xmap_base::operator_grid_orth,
                             "Return the grid-to-orthogonal coordinate "
                             "operator (translation is zero).")
      .def("coord_orth", &Xmap_base::coord_orth, py::arg("cm"),
           "Convert map coordinate to orthogonal.")
      .def("coord_map", &Xmap_base::coord_map, py::arg("co"),
           "Convert orthogonal coordinate to map.")
      .def("in_map",
           (bool(Xmap_base::*)(const clipper::Coord_grid &) const) &
               Xmap_base::in_map,
           py::arg("Coord_grid") = nullptr,
           "(This method is for compatibility with NXmap - it always returns "
           "true).")
      .def("multiplicity", &Xmap_base::multiplicity, py::arg("pos"),
           "Get multiplicity of map grid point.")
      .def("first", &Xmap_base::first,
           "Return Map_reference_index for this map.")
      .def("first_coord", &Xmap_base::first_coord,
           "Return a Map_reference_index for this map.")
      .doc() =
      "Xmap_base: base for crystallographic map class.\n"
      "The crystallographic map class stores a map of arbitrary data "
      "type. Its main difference from a 3-d array is that the data extent "
      "appears to be infinite, and yet internally only a unique ASU is "
      "stored. Iterators provide efficient access to data. "
      "This base contains everything except the data, which is templated "
      "in the derived type Xmap<T>";

  using MRB = Xmap_base::Map_reference_base;
  py::class_<MRB>(xmap_base, "Xmap_reference_base")
      .def_property_readonly("base_xmap", &MRB::base_xmap,
                             "Return parent Xmap.")
      .def_property_readonly("index", &MRB::index,
                             "Get the index into the map data array.")
      .def("last", &MRB::last, "Check for end of map.");

  using MRI = Xmap_base::Map_reference_index;
  py::class_<MRI, Xmap_base::Map_reference_base>(xmap_base,
                                                 "Xmap_reference_index")
      .def_property(
          "coord", &MRI::coord,
          [](MRI &self, const Coord_grid &pos) -> void { self.set_coord(pos); },
          "Get set current grid coordinate.")
      .def("coord_orth", &MRI::coord_orth,
           "Get current value of orthogonal coordinate.")
      .def(
          "next", [](MRI &self) -> void { self.next(); }, "Simpe increment.")
      .def("index_offset", &MRI::index_offset, py::arg("du"), py::arg("dv"),
           py::arg("dw"), "Index of neighbouring point.");

  using MRC = Xmap_base::Map_reference_coord;
  py::class_<MRC, Xmap_base::Map_reference_base>(xmap_base,
                                                 "Xmap_reference_coord")
      .def_property(
          "coord", &MRC::coord,
          [](MRC &self, const Coord_grid &pos) -> void { self.set_coord(pos); },
          "Get/set current value of coordinate.")
      .def("coord_orth", &MRC::coord_orth)
      .def_property_readonly("sym", &MRC::sym, "Get current symmetry operator.")
      .def(
          "next", [](MRC &self) -> void { self.next(); },
          "Simple increment. Use of this function resets the stored coordinate "
          "and sym.")
      .def(
          "next_u", [](MRC &self) -> void { self.next_u(); }, "Increment u.")
      .def(
          "next_v", [](MRC &self) -> void { self.next_v(); }, "Increment v.")
      .def(
          "next_w", [](MRC &self) -> void { self.next_w(); }, "Increment w.")
      .def(
          "prev_u", [](MRC &self) -> void { self.prev_u(); }, "Decrement u.")
      .def(
          "prev_v", [](MRC &self) -> void { self.prev_v(); }, "Decrement v.")
      .def(
          "prev_w", [](MRC &self) -> void { self.prev_w(); }, "Decrement w.");
}
// void declare_xmap_reference_base(py::module &m) {}
// void declare_xmap_reference_index(py::module &m) {}
//
// void declare_xmap_reference_coord(py::module &m) {}

template <class Derived, class T, class H>
void apply_xmap_fft_methods(py::class_<Derived, Xmap_base> &pyclass) {
  pyclass
      .def(
          "fft_from",
          [](Derived &self, const H &phidata, const Xmap_base::FFTtype type) {
            self.fft_from(phidata);
          },
          py::arg("fphidata"), py::arg("type") = Xmap_base::FFTtype::Default,
          "FFT from reflection list to map.")
      .def(
          "fft_to",
          [](const Derived &self, H &phidata, const Xmap_base::FFTtype type) {
            self.fft_to(phidata);
          },
          py::arg("fphidata"), py::arg("type") = Xmap_base::FFTtype::Default,
          "FFT from map to reflection list.");
}

template <class Derived, class T>
void apply_xmap_interpolation_methods(py::class_<Derived, Xmap_base> &pyclass) {
  pyclass
      .def("interp_linear_frac",
           [](const Derived &self, const Coord_frac &pos) {
             return self.template interp<Interp_linear>(pos);
           })
      .def("interp_cubic_frac",
           [](const Derived &self, const Coord_frac &pos) {
             return self.template interp<Interp_cubic>(pos);
           })
      .def("interp_linear_orth",
           [](const Derived &self, const Coord_orth &xyz) {
             return self.template interp<Interp_linear>(
                 xyz.coord_frac(self.cell()));
           })
      .def("interp_cubic_orth",
           [](const Derived &self, const Coord_orth &xyz) {
             return self.template interp<Interp_cubic>(
                 xyz.coord_frac(self.cell()));
           })
      //.def("interp_cubic_grad_frac", [](const Derived& self, const Coord_frac&
      // pos) -> py::tuple
      //{
      //    T val;
      //    Grad_frac<T> grad;
      //    self.template interp_grad<Interp_cubic>(pos, val, grad);
      //    return py::make_tuple(val, grad);
      //})
      //.def("interp_cubic_curv_frac", [](const Derived& self, const Coord_frac&
      // pos) -> py::tuple
      //{
      //    T val;
      //    Grad_frac<T> grad;
      //    Curv_frac<T> curv;
      //    self.template interp_curv<Interp_cubic>(pos, val, grad, curv);
      //    return py::make_tuple(val, grad, curv);
      //})
      ;
} // apply_xmap_interpolation_methods

template <class T> void declare_xmap(py::module &m, const std::string &name) {
  using MRI = Xmap_base::Map_reference_index;
  using MRC = Xmap_base::Map_reference_coord;
  using XMClass = Xmap<T>;
  std::string PyClass = std::string("Xmap_") + name;
  py::class_<XMClass, Xmap_base> xmap(m, PyClass.c_str(),
                                      py::buffer_protocol());
  xmap.def(py::init<>())
      .def(py::init<const Spacegroup &, const Cell &, const Grid_sampling &>())
      .def("init", &XMClass::init)
      .def_buffer([=](XMClass &self) -> py::buffer_info {
        return py::buffer_info(
            &self[self.first()],
            {self.grid_sampling().nu(), self.grid_sampling().nv(),
             self.grid_sampling().nw()},                       // dimensions
            {sizeof(T), sizeof(T) * self.grid_sampling().nu(), // strides
             sizeof(T) * self.grid_sampling().nu() *
                 self.grid_sampling().nv()});
      })
      .def("__repr__",
           [](const XMClass &self) { return display_repr<T>(self); })
      .def("__str__", [](const XMClass &self) { return display_repr<T>(self); })
      .def("__getitem__",
           [](const XMClass &self, const MRI &ix) { return self[ix]; })
      .def("__setitem__",
           [](XMClass &self, const MRI &ix, const T &val) { self[ix] = val; })
      .def("__getitem__",
           [](const XMClass &self, const MRC &ix) { return self[ix]; })
      .def("__setitem__",
           [](XMClass &self, const MRC &ix, const T &val) { self[ix] = val; })
      .def(
          "get_data",
          [](const XMClass &self, const Coord_grid &pos) {
            return self.get_data(pos);
          },
          py::arg("pos"), "Get data.")
      .def(
          "set_data",
          [](XMClass &self, const Coord_grid &pos, const T &val) {
            self.set_data(pos, val);
          },
          py::arg("pos"), py::arg("val"), "Set data.")
      .def(
          "fill_map_with", [](XMClass &self, const T &val) { self = val; },
          py::arg("val"), "Fill map with given value.")
      .def(py::self += py::self)
      .def(py::self -= py::self)
      // because python's gemmi map object's data is float type.
      .def(
          "import_from_gemmi",
          [](XMClass &self, const gemmi::Ccp4<float> &mapobj) {
            GEMMI::import_xmap(self, mapobj);
          },
          py::arg("target"), "Import Xmap from gemmi.Ccp4Map.")
      .def(
          "export_to_gemmi",
          [](const XMClass &self, gemmi::Ccp4<float> &mapobj) {
            GEMMI::export_xmap(self, mapobj);
          },
          py::arg("target"), "Export Xmap from gemmi.Ccp4Map.");
  // interpolators
  apply_xmap_interpolation_methods<XMClass, T>(xmap);
  // ffts
  apply_xmap_fft_methods<XMClass, T, HKL_data<clipper::data32::F_phi>>(xmap);
  apply_xmap_fft_methods<XMClass, T, HKL_data<clipper::data64::F_phi>>(xmap);
}

void init_maps(py::module &m) {
  declare_xmap_base(m);
  // declare_xmap_reference_base(m);
  // declare_xmap_reference_index(m);
  // declare_xmap_reference_coord(m);
  declare_xmap<int>(m, "int");
  declare_xmap<ftype32>(m, "float");
  declare_xmap<ftype64>(m, "double");
}