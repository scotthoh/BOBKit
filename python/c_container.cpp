#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/container_hkl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void init_containers(py::module &m) {
  // clipper's generic container
  py::class_<Container>(m, "Container")
      .def(py::init<const String>(), py::arg("name") = "")
      .def(py::init<Container &, const String &>(), py::arg("parent"),
           py::arg("relpath"))
      .def("update", &Container::update)
      .def("path", &Container::path)
      .def("name", &Container::name)
      .def("set_name", &Container::set_name, py::arg("name"))
      .def("is_destroyed_with_parent", &Container::is_destroyed_with_parent)
      .def("set_destroyed_with_parent", &Container::set_destroyed_with_parent,
           py::arg("destroy") = true)
      .def("move", &Container::move, py::arg("path"))
      .def("has_parent", &Container::has_parent)
      .def("parent",
           (const Container &(Container::*)() const) & Container::parent)
      .def("parent", (Container & (Container::*)()) & Container::parent)
      .def("num_children", &Container::num_children)
      .def("child", (const Container &(Container::*)(const int &) const) &
                        Container::child)
      .def("child",
           (Container & (Container::*)(const int &)) & Container::child)
      .def("ultimate_parent", (const Container &(Container::*)() const) &
                                  Container::ultimate_parent)
      .def("ultimate_parent",
           (Container & (Container::*)()) & Container::ultimate_parent)
      .def("parent_prt", &Container::parent_ptr)
      .def("debug", &Container::debug);
  // parent_ptr, parent_of_type_ptr, find_path_ptr not bind

  // clipper's HKL list and indexing object container
  py::class_<CHKL_info, Container, HKL_info>(m, "CHKL_info")
      .def(py::init<const String, const Spacegroup &, const Cell &,
                    const Resolution &, const bool &>(),
           py::arg("name") = "", py::arg("spacegroup").none(true),
           py::arg("cell").none(true), py::arg("resolution").none(true),
           py::arg("generate") = false)
      .def(py::init<Container &, const String, const bool &>(),
           py::arg("parent"), py::arg("name") = "", py::arg("generate") = false)
      .def("init", &CHKL_info::init, py::arg("spacegroup"), py::arg("cell"),
           py::arg("resolution"), py::arg("generate") = false)
      .def("generate_hkl_list", &CHKL_info::generate_hkl_list)
      .def("update", &CHKL_info::update);
}