// bind maps stuff
#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/xmap.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// template <class T> std::string display_repr(const Xmap<T> xmapclass) {
template <class T> std::string display_repr(const Xmap<T> xmap) {
  if (strcmp(typeid(T).name(), "f") == 0)
    return "<clipper.Xmap_float: ASU grid " + xmap.grid_asu().format() + ">";
  else if (strcmp(typeid(T).name(), "d") == 0)
    return "<clipper.Xmap_double: ASU grid " + xmap.grid_asu().format() + ">";
  else if (strcmp(typeid(T).name(), "i") == 0)
    return "<clipper.Xmap_int: ASU grid " + xmap.grid_asu().format() + ">";
  else
    return "<clipper.Xmap>";
}

void declare_xmap_base(py::module &m) {
  py::class_<Xmap_base> xmap_base(m, "Xmap_base");
  xmap_base.def("is_null", &Xmap_base::is_null)
      .def_property_readonly("cell", &Xmap_base::cell)
      .def_property_readonly("spacegroup", &Xmap_base::spacegroup)
      .def_property_readonly("grid_sampling", &Xmap_base::grid_sampling)
      .def_property_readonly("grid_asu", &Xmap_base::grid_asu)
      .def("coord_of", &Xmap_base::coord_of)
      .def("index_of", &Xmap_base::index_of)
      .def("to_map_unit", &Xmap_base::to_map_unit)
      .def_property_readonly("operator_orth_grid",
                             &Xmap_base::operator_orth_grid)
      .def_property_readonly("operator_grid_orth",
                             &Xmap_base::operator_grid_orth)
      .def("coord_orth", &Xmap_base::coord_orth)
      .def("coord_map", &Xmap_base::coord_map)
      .def("in_map",
           (bool(Xmap_base::*)(const clipper::Coord_grid &) const) &
               Xmap_base::in_map,
           py::arg("Coord_grid") = nullptr)
      .def("multiplicity", &Xmap_base::multiplicity)
      .def("first", &Xmap_base::first)
      .def("first_coord", &Xmap_base::first_coord);

  using MRB = Xmap_base::Map_reference_base;
  py::class_<MRB>(xmap_base, "Xmap_reference_base")
      .def_property_readonly("base_xmap", &MRB::base_xmap)
      .def_property_readonly("index", &MRB::index)
      .def("last", &MRB::last);

  using MRI = Xmap_base::Map_reference_index;
  py::class_<MRI, Xmap_base::Map_reference_base>(xmap_base,
                                                 "Xmap_reference_index")
      .def_property(
          "coord", &MRI::coord,
          [](MRI &self, const Coord_grid &pos) -> void { self.set_coord(pos); })
      .def("coord_orth", &MRI::coord_orth)
      .def("next", [](MRI &self) -> void { self.next(); })
      .def("index_offset", &MRI::index_offset);

  using MRC = Xmap_base::Map_reference_coord;
  py::class_<MRC, Xmap_base::Map_reference_base>(xmap_base,
                                                 "Xmap_reference_coord")
      .def_property(
          "coord", &MRC::coord,
          [](MRC &self, const Coord_grid &pos) -> void { self.set_coord(pos); })
      .def("coord_orth", &MRC::coord_orth)
      .def_property_readonly("sym", &MRC::sym)
      .def("next", [](MRC &self) -> void { self.next(); })
      .def("next_u", [](MRC &self) -> void { self.next_u(); })
      .def("next_v", [](MRC &self) -> void { self.next_v(); })
      .def("next_w", [](MRC &self) -> void { self.next_w(); })
      .def("prev_u", [](MRC &self) -> void { self.prev_u(); })
      .def("prev_v", [](MRC &self) -> void { self.prev_v(); })
      .def("prev_w", [](MRC &self) -> void { self.prev_w(); });
}
// void declare_xmap_reference_base(py::module &m) {}
// void declare_xmap_reference_index(py::module &m) {}
//
// void declare_xmap_reference_coord(py::module &m) {}

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
      .def("get_data", [](const XMClass &self,
                          const Coord_grid &pos) { return self.get_data(pos); })
      .def("set_data", [](XMClass &self, const Coord_grid &pos,
                          const T &val) { self.set_data(pos, val); })
      .def("fill_map_with", [](XMClass &self, const T &val) { self = val; })
      .def(py::self += py::self)
      .def(py::self -= py::self);
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