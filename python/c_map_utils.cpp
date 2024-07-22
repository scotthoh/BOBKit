#include <clipper/clipper.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_map_stats(py::module &m) {
  py::class_<Map_stats>(m, "Map_stats")
      .def(py::init<>())
      .def(py::init<const Xmap<ftype32> &>())
      .def(py::init<const Xmap<ftype64> &>())
      .def(py::init<const NXmap<ftype32> &>())
      .def(py::init<const NXmap<ftype64> &>())
      .def_property_readonly("mean", &Map_stats::mean)
      .def_property_readonly("std_dev", &Map_stats::std_dev)
      .def_property_readonly("min", &Map_stats::min)
      .def_property_readonly("max", &Map_stats::max)
      .def_property_readonly("range", &Map_stats::range)
      .def("__repr__",
           [](const Map_stats &self) { return "<clipper.Map_stats class.>"; })
      .def("format",
           [](const Map_stats &self) {
             return "Min: " + String(self.min(), 10, 4) +
                    ", Max: " + String(self.max(), 10, 4) +
                    ", Mean: " + String(self.mean(), 10, 4) +
                    ", SD: " + String(self.std_dev(), 10, 4);
           })
      .def("__str__", [](const Map_stats &self) {
        return String(self.min(), 10, 4) + ", " + String(self.max(), 10, 4) +
               ", " + String(self.mean(), 10, 4) + ", " +
               String(self.std_dev(), 10, 4);
      });
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