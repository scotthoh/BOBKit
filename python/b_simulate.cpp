#include "buccaneer/simulate-lib.h"

#include "helper_functions.h"
#include "type_conversions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void declare_map_simulate(py::module &m) {
  py::class_<MapSimulate>(m, "MapSimulate")
      .def(py::init<int, int>(), py::arg("nresbins") = 100,
           py::arg("binmin") = 20)
      .def("__call__", &MapSimulate::operator(), py::arg("sim_f"),
           py::arg("sim_hl"), py::arg("ref_f"), py::arg("ref_hl"),
           py::arg("wrk_f"), py::arg("wrk_hl"))
      .def("__repr__", [](const MapSimulate &self) {
        return "<buccaneer.MapSimulate class.>";
      });
}