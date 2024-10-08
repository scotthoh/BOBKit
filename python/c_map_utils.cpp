// Wrapper for clipper map utils
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_map_stats(py::module &m) {
  py::class_<Map_stats>(m, "Map_stats")
      .def(py::init<>())
      .def(py::init<const Xmap<ftype32> &>(), "Constructor from Xmap")
      .def(py::init<const Xmap<ftype64> &>(), "Constructor from Xmap")
      .def(py::init<const NXmap<ftype32> &>(), "Constructor from NXmap")
      .def(py::init<const NXmap<ftype64> &>(), "Constructor from NXmap")
      .def_property_readonly("mean", &Map_stats::mean, "Return mean of map.")
      .def_property_readonly("std_dev", &Map_stats::std_dev,
                             "Return standard deviation of map.")
      .def_property_readonly("min", &Map_stats::min, "Return minimum of map.")
      .def_property_readonly("max", &Map_stats::max, "Return maximum of map.")
      .def_property_readonly("range", &Map_stats::range, "Return range.")
      .def("__repr__",
           [](const Map_stats &self) { return "<clipper.Map_stats class.>"; })
      .def(
          "format",
          [](const Map_stats &self) {
            return "Min: " + String(self.min(), 10, 4) +
                   ", Max: " + String(self.max(), 10, 4) +
                   ", Mean: " + String(self.mean(), 10, 4) +
                   ", SD: " + String(self.std_dev(), 10, 4);
          },
          "Return formatted string representation.")
      .def("__str__",
           [](const Map_stats &self) {
             return String(self.min(), 10, 4) + ", " +
                    String(self.max(), 10, 4) + ", " +
                    String(self.mean(), 10, 4) + ", " +
                    String(self.std_dev(), 10, 4);
           })
      .doc() = "Generic map statistics class.\nThis class is used "
               "to calculate and store the mean and standard deviation "
               "of a generic map object of scalar types (e.g. Xmap, "
               "NXmap). If the map contains NaN values, those points "
               "are excluded for the calculation. In the case of an Xmap, "
               "the proper multiplicty corrections are applied to give "
               "statistics for a whole unit cell.";
}

// void declare_map_index_sort(py::module &m) {
//   py::class_<Map_index_sort>(m, "Map_index_sort")
//       .def_static(
//           "sort_increasing",
//           static_cast<void (*)(const Xmap<ftype32> &, std::vector<int> &)>(
//               &Map_index_sort::sort_increasing))
//       .def_static();
// }

void init_map_utils(py::module &m) {
  declare_map_stats(m);
  // declare_map_index_sort(m);
}