// Wrapper for clipper containers
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void init_containers(py::module &m) {
  // clipper's generic container
  py::class_<Container>(m, "Container")
      .def(py::init<const String>(), py::arg("name") = "",
           "Constructor: make null object or top object in a tree.")
      .def(py::init<Container &, const String &>(), py::arg("parent"),
           py::arg("relpath"),
           "Constructor: from any other member and a relative path.")
      .def("update", &Container::update, "Update hierarchical content.")
      .def("path", &Container::path, "Get the path of this tree object.")
      .def("name", &Container::name, "Get the name of this tree object.")
      .def("set_name", &Container::set_name, py::arg("name"),
           "Set the name of this tree object.")
      .def("is_destroyed_with_parent", &Container::is_destroyed_with_parent,
           "Is this object to be destroyed when parent is destroyed?")
      .def("set_destroyed_with_parent", &Container::set_destroyed_with_parent,
           py::arg("destroy") = true,
           "Set this object to be destroyed when parent is destroyed.")
      .def("move", &Container::move, py::arg("path"),
           "Move this object to somewhere else in the hierachy.")
      .def("has_parent", &Container::has_parent,
           "Test if this object has parent.")
      .def("parent",
           (const Container &(Container::*)() const) & Container::parent,
           "Get the parent of this object.")
      .def("parent", (Container & (Container::*)()) & Container::parent,
           "Get the parent of this object.")
      .def("num_children", &Container::num_children,
           "Return number of children.")
      .def("child",
           (const Container &(Container::*)(const int &) const) &
               Container::child,
           "Get the i'th child of this object.")
      .def("child",
           (Container & (Container::*)(const int &)) & Container::child,
           "Get the i'th child of this object.")
      .def("ultimate_parent",
           (const Container &(Container::*)() const) &
               Container::ultimate_parent,
           "Get the ultimate of this object - top of the tree.")
      .def("ultimate_parent",
           (Container & (Container::*)()) & Container::ultimate_parent,
           "Get the ultimate of this object - top of the tree.")
      .def("parent_prt", &Container::parent_ptr,
           "Get the parent of this object (NULL or failt)")
      .def("debug", &Container::debug, "Debug method.")
      .doc() =
      "Definition for a generic container Object\n"
      "Container is a definition for a generic container object with a "
      "name, parents, and children. Any object that wants to be part of the "
      "tree simply subclasses this class. The class also implements search "
      "and move objects. The tree is navigate using unix-like pathnames. A "
      "recursive update method can be overridden to update content after "
      "altering the hierarchy. "
      "The top container in a tree is created by passing Container() as its "
      "parent.";
  // parent_ptr, parent_of_type_ptr, find_path_ptr not bind

  // clipper's HKL list and indexing object container
  py::class_<CHKL_info, Container, HKL_info>(m, "CHKL_info")
      .def(py::init<const String, const Spacegroup &, const Cell &,
                    const Resolution &, const bool &>(),
           py::arg("name") = "", py::arg("spacegroup").none(true),
           py::arg("cell").none(true), py::arg("resolution").none(true),
           py::arg("generate") = false,
           "Constructor make null object or top object in tree.")
      .def(py::init<Container &, const String, const bool &>(),
           py::arg("parent"), py::arg("name") = "", py::arg("generate") = false,
           "Constructor: inherit spacegroup, cell and resolution.")
      .def("init", &CHKL_info::init, py::arg("spacegroup"), py::arg("cell"),
           py::arg("resolution"), py::arg("generate") = false,
           "Initialiser: supply or inherit spacegroup, cell and resolution.")
      .def("generate_hkl_list", &CHKL_info::generate_hkl_list,
           "Synthesize hkl list and update children.")
      .def("update", &CHKL_info::update, "Hierarchical update.")
      .doc() = "HKL list and indexing object container\n"
               "CHKL_info: This is the reflection list object for reflection "
               "data objects to reference in order to identify their data "
               "entries.";
}