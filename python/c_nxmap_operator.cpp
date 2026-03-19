// Wrapper for clipper nxmap operator
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>

using namespace clipper;

template<class Base, class T, class M>
void add_get_map_data_methods(nb::class_<Base> &pyclass, std::string funcname) {
  if constexpr (std::is_same<M, NXmap<T>>::value ) {
    pyclass.def( funcname.c_str(), [](Base &self, const M &map, const Coord_grid &c, int &order) {
      if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
      }
      switch (order) {
        case 0: return self.template nxmap_data<Interp_nearest, T>(map, c);
        case 1: return self.template nxmap_data<Interp_linear, T>(map, c);
        case 3: 
        default: return self.template nxmap_data<Interp_cubic, T>(map, c);
      }      
    }, nb::arg("nxmap"), nb::arg("coord_grid"), nb::arg("order")=3,
    "Get value of nxmap at xmap grid coord using fastest appropriate method with the given interpolation order. Default: 3 - cubic.");
  }
  if constexpr (std::is_same<M, Xmap<T>>::value ) {
    pyclass.def( funcname.c_str(), [](Base &self, const M &map, const Coord_grid &c, int &order) {
      if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
      }
      switch (order) {
        case 0: return self.template xmap_data<Interp_nearest, T>(map, c);
        case 1: return self.template xmap_data<Interp_linear, T>(map, c);
        case 3: 
        default: return self.template xmap_data<Interp_linear, T>(map, c);
      }
    }, nb::arg("xmap"), nb::arg("coord_grid"), nb::arg("order")=3,
    "Get value of xmap at nxmap grid coord using fastest appropriate method with the given interpolation order. Default: 3 - cubic.");
  }
}

template<class I, class T>
void xmap_interp_all_points(NX_operator &self, const Xmap<T> &map, NXmap<T> &outmap, const Coord_grid &cg0, const Coord_grid &cg1) {
  NXmap_base::Map_reference_coord i0, iu, iv, iw;
  i0 = NXmap_base::Map_reference_coord(outmap, cg0);
  for(iu = i0; iu.coord().u() <= cg1.u(); iu.next_u())
    for(iv = iu; iv.coord().v() <= cg1.v(); iv.next_v())
      for(iw = iv; iw.coord().w() <= cg1.w(); iw.next_w()) {
        outmap[iw] = self.template xmap_data<I, T>(map, iw.coord());
      }
}

template<class I, class T>
void nxmap_interp_all_points(NX_operator &self, const NXmap<T> &map, Xmap<T> &outmap, const Coord_grid &cg0, const Coord_grid &cg1) {
  Xmap_base::Map_reference_coord i0, iu, iv, iw;
  i0 = Xmap_base::Map_reference_coord(outmap, cg0);
  for(iu = i0; iu.coord().u() <= cg1.u(); iu.next_u())
    for(iv = iu; iv.coord().v() <= cg1.v(); iv.next_v())
      for(iw = iv; iw.coord().w() <= cg1.w(); iw.next_w()) {
        outmap[iw] = self.template nxmap_data<I, T>(map, iw.coord());
      }
}

template<class Base, class T>
void add_nxmap_interp_all(nb::class_<Base> &pyclass) {
  pyclass.def("nxmap_interp_all_points", [](Base &self, const NXmap<T> &nxmap, Xmap<T> &xmap, const Coord_grid &cg0, const Coord_grid &cg1, int &order) {
        if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
        }
        switch (order) {
          case 0: nxmap_interp_all_points<Interp_nearest, T>(self, nxmap, xmap, cg0, cg1);
                  break;
          case 1: nxmap_interp_all_points<Interp_linear, T>(self, nxmap, xmap, cg0, cg1);
                  break;
          case 3: nxmap_interp_all_points<Interp_cubic, T>(self, nxmap, xmap, cg0, cg1);
                  break;
        }
      }, nb::arg("nxmap"), nb::arg("result_xmap"), nb::arg("coord_grid_start"), nb::arg("coord_grid_end"), nb::arg("order")=3,
    "Perform n-th order interpolation for all points within given grid coordinates and store values in Xmap supplied.");
}

template<class Base, class T>
void add_xmap_interp_all(nb::class_<Base> &pyclass) {
  pyclass.def("xmap_interp_all_points", [](Base &self, const Xmap<T> &xmap, NXmap<T> &nxmap, const Coord_grid &cg0, const Coord_grid &cg1, int &order) {
        if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
        }
        switch (order) {
          case 0: xmap_interp_all_points<Interp_nearest, T>(self, xmap, nxmap, cg0, cg1);
                  break;
          case 1: xmap_interp_all_points<Interp_linear, T>(self, xmap, nxmap, cg0, cg1);
                  break;
          case 3: xmap_interp_all_points<Interp_cubic, T>(self, xmap, nxmap, cg0, cg1);
                  break;
        }
      }, nb::arg("xmap"), nb::arg("result_nxmap"), nb::arg("coord_grid_start"), nb::arg("coord_grid_end"), nb::arg("order")=3,
    "Perform n-th order interpolation for all points within given grid coordinates and store values in NXmap supplied.");
}

void declare_nxoperator(nb::module_ &m)
{
  nb::class_<NX_operator> nx_op(m, "NX_operator");
  nx_op
      .def(nb::init<>(), "Null constructor")
      .def(nb::init<const Xmap_base &, const NXmap_base &, const RTop_orth &>(),
           nb::arg("xmap"), nb::arg("nxmap"), nb::arg("rtop"),
           "Constructor: from Xmap, NXmap, and operator")
      .def(nb::init<const Cell &, const Grid_sampling &, const NXmap_base &, const RTop_orth &>(),
           nb::arg("cell"), nb::arg("grid"), nb::arg("nxmap"), nb::arg("rtop"),
           "Constructor: from cell, grid sampling, NXmap, and operator")
      .def("init", (void(NX_operator::*)(const Xmap_base &, const NXmap_base &, const RTop_orth &)) & NX_operator::init,
           nb::arg("xmap"), nb::arg("nxmap"), nb::arg("rtop"), "Initialiser: from Xmap, NXmap, and operator")
      .def("init", (void(NX_operator::*)(const Cell &, const Grid_sampling &, const NXmap_base &, const RTop_orth &)) & NX_operator::init,
           nb::arg("cell"), nb::arg("grid"), nb::arg("nxmap"), nb::arg("rtop"),
           "Initialiser: from cell, grid sampling, NXmap, and operator")
      .def("coord_map", &NX_operator::coord_map, nb::arg("coord_frac"), "Convert xtal frac coord to nxmap map coord")
      .def("coor_frac", &NX_operator::coord_frac, nb::arg("coord_map"), "Convert nxmap map coord to xtal frac coord")
      .def("is_null", &NX_operator::is_null, "Test if object has been initialised.")
      .doc() = "NX_operator: non-crystal map operator\n"
               "This class holds a reference to a non-crystal map frame from "
               "somewhere within a crystallographic map frame. In the general "
               "case, an orthogonal rotation-translation operator is provided "
               "which maps the orthogonal frame of the crystal space onto the "
               "orthogonal frame of the NXmap space.\n"
               "The object calculates and stores optimised transformations between "
               "the crystallographic frame (described either in fractional or grid "
               "coordinates), and the NXmap grid. Fast paths are generated "
               "automatically if the grids are related.";
  add_get_map_data_methods<NX_operator, float, NXmap<float>>(nx_op, "nxmap_data");
  add_get_map_data_methods<NX_operator, float, Xmap<float>>(nx_op, "xmap_data");
  add_get_map_data_methods<NX_operator, double, NXmap<double>>(nx_op, "nxmap_data");
  add_get_map_data_methods<NX_operator, double, Xmap<double>>(nx_op, "xmap_data");
  add_xmap_interp_all<NX_operator, float>(nx_op);
  add_nxmap_interp_all<NX_operator, float>(nx_op);
  add_xmap_interp_all<NX_operator, double>(nx_op);
  add_nxmap_interp_all<NX_operator, double>(nx_op);
}

template<class T>
void declare_nxmap_operator(nb::module_ &m, const std::string &name){
  using Class = NXmap_operator<T>;
  std::string PyClass = std::string("NXmap_operator_") + name;
  nb::class_<Class, NX_operator>(m, PyClass.c_str())
    .def(nb::init<>(), "Null constructor")
    .def(nb::init<const Xmap_base &, const NXmap<T> &, const RTop_orth &>(),
         nb::arg("xmap"), nb::arg("nxmap"), nb::arg("rtop"), "Constructor: from Xmap, NXmap, and operator.")
    .def(nb::init<const Cell &, const Grid_sampling &, const NXmap<T> &, const RTop_orth &>(),
         nb::arg("cell"), nb::arg("grid"), nb::arg("nxmap"), nb::arg("rtop"),
         "Constructor: from cell, grid sampling, NXmap, and operator.")
    .def("init", (void(Class::*)(const Xmap_base &, const NXmap<T> &, const RTop_orth &)) &  Class::init,
         nb::arg("xmap"), nb::arg("nxmap"), nb::arg("rtop"), "Initialiser: from Xmap, NXmap, and operator")
    .def("init", (void(Class::*)(const Cell &, const Grid_sampling &, const NXmap<T> &, const RTop_orth &)) & Class::init,
         nb::arg("cell"), nb::arg("grid"), nb::arg("nxma"), nb::arg("rtop"), 
         "Initialiser: from Xmap, NXmap, and operator")
    .def("nxmap", &Class::nxmap, "Get the target NXmap of this operator.")
    .def("nxmap_data", [](Class &self, const Coord_grid &c, int &order) {
      if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
      }
      switch (order) {
        case 0: return self.template nxmap_data<Interp_nearest>(c);
        case 1: return self.template nxmap_data<Interp_linear>(c);
        case 3: 
        default: return self.template nxmap_data<Interp_cubic>(c);
      }
      }, nb::arg("coord_grid"), nb::arg("order")=3,
      "Access NXmap directly from Xmap grid coord using fastest appropriate method with the given interpolation order. Default: 3 - cubic.")
    .def("nxmap_interp_all_points", [](Class &self, Xmap<T> &xmap, const Coord_grid &c0, const Coord_grid &c1, int &order) {
      if (order != 0 && order != 1 && order != 3) {
          throw std::invalid_argument("Invalid order value! (0 - nearest, 1 - linear, 3 - cubic)");
        }
      switch (order) {
        case 0: nxmap_interp_all_points<Interp_nearest, T>(self, self.nxmap(), xmap, c0, c1);
                break;
        case 1: nxmap_interp_all_points<Interp_linear, T>(self, self.nxmap(), xmap, c0, c1);
                break;
        case 3: nxmap_interp_all_points<Interp_cubic, T>(self, self.nxmap(), xmap, c0, c1);
                break;
      }
      }, nb::arg("result_xmap"), nb::arg("coord_grid_start"), nb::arg("coord_grid_end"), nb::arg("order")=3,
    "Perform n-th order interpolation for all points within given grid coordinates and store values on Xmap supplied.")
    .doc() = "NXmap_operator: non-crystal map operator referencing a particular NXmap.\n"
             "This class holds a reference to a non-crystal map object from "
             "somewhere within a crystallographic map frame. In the general "
             "case, an orthogonal rotation-translation operator is provided "
             "which maps the orthogonal frame of the crystal space onto the "
             "orthogonal frame of the NXmap space.\n"
             "The object calculates and stores optimised transformations between "
             "the crystallgoraphic frame (described either in fractional or grid "
             "coordinates), and the NXmap grid. Fast paths are generated "
             "automatically if the grids are related.\n"
             "Note: This object differes from NX_operator in that it keeps a "
             "reference to an individual NXmap, which may be used to access that "
             "object directly.";

}

void add_nxmap_operator(nb::module_ &m) {
  declare_nxoperator(m);
  declare_nxmap_operator<ftype32>(m, "float");
  declare_nxmap_operator<ftype64>(m, "double");
  declare_nxmap_operator<int>(m, "int");
}