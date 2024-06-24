#include "buccaneer/buccaneer-grow.h"

#include "type_conversions.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_ca_grow(py::module &m) {
  py::class_<Ca_grow>(m, "Ca_grow")
      .def(py::init<int>(), py::arg("n_grow") = 25)
      .def("__call__", &Ca_grow::operator(), py::arg("mol"), py::arg("xmap"),
           py::arg("llktarget"))
      .def_static("grow", &Ca_grow::grow, py::arg("chain"), py::arg("xmap"),
                  py::arg("llktarget"), py::arg("rama1"), py::arg("rama2"),
                  py::arg("cutoff"), py::arg("ngrow"))
      .def_static("next_ca_group", &Ca_grow::next_ca_group, py::arg("chain"),
                  py::arg("xmap"), py::arg("llktarget"), py::arg("rama1"),
                  py::arg("rama2"))
      .def_static("prev_ca_group", &Ca_grow::prev_ca_group, py::arg("chain"),
                  py::arg("xmap"), py::arg("llktarget"), py::arg("rama1"),
                  py::arg("rama2"))
      .def_static("set_cpus", &Ca_grow::set_cpus, py::arg("ncpus"))
      .def("__repr__", [](const Ca_grow &self) {
        std::stringstream stream;
        stream << "<buccaneer.Ca_grow class>";
        return stream.str();
      });
}

void declare_grow_threaded(py::module &m) {
  py::class_<Grow_threaded>(m, "Grow_threaded")
      .def(py::init<>())
      .def(py::init<const std::vector<Ca_chain> &, const Xmap<float> &,
                    const LLK_map_target &, const double &, const int &>(),
           py::arg("chains"), py::arg("xmap"), py::arg("llktarget"),
           py::arg("cutoff"), py::arg("n_grow"))
      .def("grow", &Grow_threaded::grow, py::arg("chn"))
      .def_property_readonly("result", &Grow_threaded::result)
      .def("__call__", &Grow_threaded::operator(), py::arg("nthread") = 0)
      .def("merge", &Grow_threaded::merge, py::arg("other"))
      // inherited function/property
      .def_property_readonly("id", &Grow_threaded::id)
      .def("__repr__", [](const Grow_threaded &self) {
        std::stringstream stream;
        stream << "<buccaneer.Grow_threaded class, id=";
        stream << self.id() << ">";
        return stream.str();
      });
}

// Target_fn_refine_n_terminal_build defined in b_simplex.cpp
// to be within same scope as Target_fn_zero_order trampoline definition

void init_ca_grow(py::module &m) {
  declare_ca_grow(m);
  declare_grow_threaded(m);
}