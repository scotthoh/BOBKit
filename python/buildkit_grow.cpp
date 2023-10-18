#include "buccaneer/buccaneer-grow.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_ca_grow(py::module &m)
{
  py::class_<Ca_grow>(m, "Ca_grow")
      .def(py::init<int>(), py::arg("n_grow") = 25)
      .def("__call__", &Ca_grow::operator())
      .def_static("grow", &Ca_grow::grow)
      .def_static("next_ca_group", &Ca_grow::next_ca_group)
      .def_static("prev_ca_group", &Ca_grow::prev_ca_group)
      .def_static("set_cpus", &Ca_grow::set_cpus);
}

void declare_grow_threaded(py::module &m)
{
  py::class_<Grow_threaded>(m, "Grow_threaded")
      .def(py::init<>())
      .def(py::init<const std::vector<Ca_chain> &, const Xmap<float> &, const LLK_map_target &, const double &, const int &>())
      .def("grow", &Grow_threaded::grow)
      .def("result", &Grow_threaded::result)
      .def("__call__", &Grow_threaded::operator())
      .def("merge", &Grow_threaded::merge);
}

void declare_target_fn_refine_n_terminal_build(py::module &m)
{
  py::class_<Target_fn_refine_n_terminal_build>(m, "Target_fn_refine_n_terminal_build")
      .def(py::init<>())
      .def(py::init<const Xmap<float> &, const LLK_map_target &, const Ramachandran &, const Ramachandran &, const double &>())
      .def_property_readonly("num_params", &Target_fn_refine_n_terminal_build::num_params)
      .def("__call__", &Target_fn_refine_n_terminal_build::operator())
      .def("refine", &Target_fn_refine_n_terminal_build::refine);
}
void init_ca_grow(py::module &m)
{
  declare_ca_grow(m);
  declare_grow_threaded(m);
  declare_target_fn_refine_n_terminal_build(m);
}