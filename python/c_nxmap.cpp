#include "helper_functions.h"
#include "type_conversions.h"
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
  nxmapbase.def("is_null", &NXmap_base::is_null)
      .def_property_readonly("grid", &NXmap_base::grid)
      .def_property_readonly("operator_orth_grid",
                             &NXmap_base::operator_orth_grid)
      .def_property_readonly("operator_grid_orth",
                             &NXmap_base::operator_grid_orth)
      .def("coord_orth", &NXmap_base::coord_orth, py::arg("coord_map"))
      .def("coord_map", &NXmap_base::coord_map, py::arg("coord_orth"))
      .def("in_map",
           (bool(NXmap_base::*)(const clipper::Coord_grid &) const) &
               NXmap_base::in_map,
           py::arg("pos"))
      .def("multiplicity", &NXmap_base::multiplicity,
           py::arg("coord_grid") = nullptr)
      .def("first", &NXmap_base::first)
      .def("first_coord", &NXmap_base::first_coord);

  using MRB = NXmap_base::Map_reference_base;
  py::class_<MRB>(nxmapbase, "NXmap_reference_base")
      .def_property_readonly("base_nxmap", &MRB::base_nxmap)
      .def_property_readonly("index", &MRB::index)
      .def("last", &MRB::last);

  using MRI = NXmap_base::Map_reference_index;
  py::class_<MRI, MRB>(nxmapbase, "NXmap_reference_index")
      .def_property(
          "coord", &MRI::coord,
          [](MRI &self, const Coord_grid &pos) -> void { self.set_coord(pos); })
      .def("coord_orth", &MRI::coord_orth)
      .def("next", [](MRI &self) -> void { self.next(); })
      .def("index_offset", &MRI::index_offset, py::arg("du"), py::arg("dv"),
           py::arg("dw"));

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
      .def(py::self -= py::self);
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