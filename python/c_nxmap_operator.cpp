// Wrapper for clipper nxmap operator
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>

using namespace clipper;

template<class Base, class I, class T, class M>
void add_get_map_data_methods(nb::class_<Base> &pyclass, const char *map_type, const char *interp_data_type) {
  auto fn_name = map_type + std::string( "_data_" ) + interp_data_type;
  if (std::is_same<M, NXmap<T>>::value ) {
    pyclass.def( fn_name.c_str(), [](Base &self, const M& map, const Coord_grid& c) {
      return self.template nxmap_data<I, T>(map, c);
    }, nb::arg("nxmap"), nb::arg("coord_grid"),
    "Get value of nxmap at xmap grid coord using fastest appropriate method");
  }
  if (std::is_same<M, Xmap<T>>::value ) {
    pyclass.def( fn_name.c_str(), [](Base &self, const M& map, const Coord_grid& c) {
      return self.template xmap_data<I, T>(map, c);
    }, nb::arg("xmap"), nb::arg("coord_grid"),
    "Get value of xmap at nxmap grid coord using fastest appropriate method");
  }
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
      //.def("nxmap_data_linear_float", &NX_operator)
      .doc() = "NX_operator: non-crystal map operator\n"
               "This class holds a reference to a non-crystal map frame from "
               "somewhere within a crystallographic map frame. In the general "
               "case, an orthogonal rotation-translation operator is provided "
               "which maps the orthogonal frame of the crystal space onto the "
               "orthogonal frame of the NXmap space.\n"
               "The object calculates and stores optimised transformations between "
               "the crystallgoraphic frame (described either in fractional or grid "
               "coordinates), and the NXmap grid. Fast paths are generated "
               "automatically if the grids are related.";
  add_get_map_data_methods<NX_operator, Interp_linear, float, NXmap<float>>(nx_op, "nxmap", "interp_linear_float");
  add_get_map_data_methods<NX_operator, Interp_cubic, float, NXmap<float>>(nx_op, "nxmap", "interp_cubic_float");
  add_get_map_data_methods<NX_operator, Interp_nearest, float, NXmap<float>>(nx_op, "nxmap", "interp_nearest_float");
  add_get_map_data_methods<NX_operator, Interp_linear, float, Xmap<float>>(nx_op, "xmap", "interp_linear_float");
  add_get_map_data_methods<NX_operator, Interp_cubic, float, Xmap<float>>(nx_op, "xmap", "interp_cubic_float");
  add_get_map_data_methods<NX_operator, Interp_nearest, float, Xmap<float>>(nx_op, "xmap", "interp_nearest_float");
  
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
    .def("nxmap_data_interp_linear", [](Class &self, const Coord_grid &c) {
      return self.template nxmap_data<Interp_linear>(c);
    })
    .def("nxmap_data_interp_cubic", [](Class &self, const Coord_grid &c) {
      return self.template nxmap_data<Interp_cubic>(c);
    })
    .def("nxmap_data_interp_nearest", [](Class &self, const Coord_grid &c) {
      return self.template nxmap_data<Interp_nearest>(c);
    })
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
}