// Wrapper for clipper NXmap
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

template <class T> std::string display_repr(const NXmap<T> nxmap) {
  if (strcmp(typeid(T).name(), "f") == 0)
    return "<clipper.NXmap_float: Grid " + nxmap.grid().format() + ">";
  else if (strcmp(typeid(T).name(), "d") == 0)
    return "<clipper.NXmap_double: Grid " + nxmap.grid().format() + ">";
  else if (strcmp(typeid(T).name(), "i") == 0)
    return "<clipper.NXmap_int: Grid " + nxmap.grid().format() + ">";
  // return "<clipper.Xmap_int>";
  else
    return "<clipper.NXmap>";
}

void declare_map_base(py::module &m) {
  py::class_<NXmap_base> nxmapbase(m, "NXmap_base");
  nxmapbase
      .def("is_null", &NXmap_base::is_null,
           "Test if object has been initialised.")
      .def_property_readonly("grid", &NXmap_base::grid,
                             "Return the grid dimensions for this map.")
      .def_property_readonly(
          "operator_orth_grid", &NXmap_base::operator_orth_grid,
          "Return the orthogonal-to-grid coordinate operator.")
      .def_property_readonly(
          "operator_grid_orth", &NXmap_base::operator_grid_orth,
          "Return the grid-to-orthogonal coordinate operator.")
      .def("coord_orth", &NXmap_base::coord_orth, py::arg("coord_map"),
           "Convert map coordinate to orthogonal.")
      .def("coord_map", &NXmap_base::coord_map, py::arg("coord_orth"),
           "Convert orthogonal coordinate to map.")
      .def("in_map",
           (bool(NXmap_base::*)(const clipper::Coord_grid &) const) &
               NXmap_base::in_map,
           py::arg("pos"), "Is the given coord available in the map?")
      .def("multiplicity", &NXmap_base::multiplicity,
           py::arg("coord_grid") = nullptr,
           "Get multiplicity of a map grid point (always 1 for NXmap).")
      .def("first", &NXmap_base::first,
           "Return a basic Map_reference_index for this map.")
      .def("first_coord", &NXmap_base::first_coord,
           "Return a coord Map_reference_index for this map")
      .doc() = "NXmap_base: base for non-crystallographic map class.\n"
               "The non-crystallographic map class stores a map of arbitrary "
               "data type. Unlike an Xmap it is finite in extent and has no "
               "symmetry. An RT operator provides mapping onto an arbitrary "
               "orthogonal coordinate frame. Iterators provide efficient "
               "access to data. This base contains everything except the data, "
               "which is templated in the derived type clipper::NXmap<T>.";

  using MRB = NXmap_base::Map_reference_base;
  py::class_<MRB>(nxmapbase, "NXmap_reference_base")
      .def_property_readonly("base_nxmap", &MRB::base_nxmap,
                             "Return the parent NXmap.")
      .def_property_readonly("index", &MRB::index,
                             "Get the index into the map data array.")
      .def("last", &MRB::last, "Check for end of map.")
      .doc() = "Map reference base class.\n"
               "This is a reference to an Map. It forms a base class for "
               "index-like and coordinate-like Map references.";

  using MRI = NXmap_base::Map_reference_index;
  py::class_<MRI, MRB>(nxmapbase, "NXmap_reference_index")
      .def_property(
          "coord", &MRI::coord,
          [](MRI &self, const Coord_grid &pos) -> void { self.set_coord(pos); },
          "Get/set current grid coordinate.")
      .def("coord_orth", &MRI::coord_orth)
      .def("next", [](MRI &self) -> void { self.next(); })
      .def("index_offset", &MRI::index_offset, py::arg("du"), py::arg("dv"),
           py::arg("dw"))
      .doc() = "Map reference with index-like behaviour.\n This is a reference "
               "to a map coordinate. It behaves like a simple index into the "
               "map, but can be easily converted into a coordinate as and when "
               "required. It also implements methods for iterating through the "
               "map. It is very compact, but coord() involves some overhead.";

  using MRC = NXmap_base::Map_reference_coord;
  py::class_<MRC, MRB>(nxmapbase, "NXmap_reference_coordinate")
      .def_property(
          "coord", &MRC::coord,
          [](MRC &self, const Coord_grid &pos) -> void { self.set_coord(pos); })
      .def("coord_orth", &MRC::coord_orth)
      .def("next", [](MRC &self) -> void { self.next(); })
      .def("next_u", [](MRC &self) -> void { self.next_u(); })
      .def("next_v", [](MRC &self) -> void { self.next_v(); })
      .def("next_w", [](MRC &self) -> void { self.next_w(); })
      .def("prev_u", [](MRC &self) -> void { self.prev_u(); })
      .def("prev_v", [](MRC &self) -> void { self.prev_v(); })
      .def("prev_w", [](MRC &self) -> void { self.prev_w(); });
}

// void declare_map_reference_base(py::module &m) {}
//
// void declare_map_reference_index(py::module &m) {}
//
// void declare_map_reference_coord(py::module &m) {}

template <class T> void declare_nxmap(py::module &m, const std::string &name) {
  using MRI = NXmap_base::Map_reference_index;
  using MRC = NXmap_base::Map_reference_coord;
  using NXMClass = NXmap<T>;
  std::string PyClass_name = std::string("NXmap_") + name;
  py::class_<NXMClass, NXmap_base> nxmap(m, PyClass_name.c_str(),
                                         py::buffer_protocol());
  nxmap.def(py::init<>())
      .def(py::init<const Grid &, const RTop<> &>())
      .def(py::init<const Cell &, const Grid_sampling &, const Grid_range &>())
      .def("init",
           (void(NXMClass::*)(const Grid &, const RTop<> &)) & NXMClass::init,
           py::arg("grid"), py::arg("rtop"))
      .def("init",
           (void(NXMClass::*)(const Cell &, const Grid_sampling &,
                              const Grid_range &)) &
               NXMClass::init,
           py::arg("cell"), py::arg("grid_sampling"), py::arg("grid_range"))
      .def_buffer([=](NXMClass &self) -> py::buffer_info {
        return py::buffer_info(
            &self[self.first()],
            {self.grid().nu(), self.grid().nv(),
             self.grid().nw()}, // dimensions
            {sizeof(T),
             sizeof(T) * self.grid().nu(), // strides
             sizeof(T) * self.grid().nu() * self.grid().nv()});
      })
      .def("__repr__",
           [](const NXMClass &self) { return display_repr<T>(self); })
      .def("__str__",
           [](const NXMClass &self) { return display_repr<T>(self); })
      .def("__getitem__",
           [](const NXMClass &self, const MRI &ix) { return self[ix]; })
      .def("__setitem__",
           [](NXMClass &self, const MRI &ix, const T &val) { self[ix] = val; })
      .def("__getitem__",
           [](const NXMClass &self, const MRC &ix) { return self[ix]; })
      .def("__setitem__",
           [](NXMClass &self, const MRC &ix, const T &val) { self[ix] = val; })
      .def("get_data", [](const NXMClass &self,
                          const Coord_grid &pos) { return self.get_data(pos); })
      .def("set_data", [](NXMClass &self, const Coord_grid &pos,
                          const T &val) { self.set_data(pos, val); })
      .def("fill_map_with", [](NXMClass &self, const T &val) { self = val; })
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(
          "import_from_gemmi",
          [](NXMClass &self, const gemmi::Ccp4<float> &mapobj) {
            GEMMI::import_nxmap(self, mapobj);
          },
          py::arg("target"), "Import NXmap from gemmi.Ccp4Map.")
      .def(
          "export_to_gemmi",
          [](const NXMClass &self, gemmi::Ccp4<float> &mapobj,
             const Cell &unitcell) {
            GEMMI::export_nxmap(self, mapobj, unitcell);
          },
          py::arg("target"), py::arg("cell"),
          "Export NXmap from gemmi.Ccp4Map.");
}

void init_nxmap(py::module &m) {
  declare_map_base(m);
  // declare_map_reference_base(m);
  // declare_map_reference_index(m);
  // declare_map_reference_coord(m);
  declare_nxmap<int>(m, "int");
  declare_nxmap<ftype32>(m, "float");
  declare_nxmap<ftype64>(m, "double");
}