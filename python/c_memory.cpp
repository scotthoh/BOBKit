#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/clipper_memory.h>
#include <clipper/core/clipper_sysdep.h>

namespace py = pybind11;
using namespace clipper;

template <class T>
void declare_property(py::module &m, const std::string &name) {
  using Class = Property<T>;
  std::string PyClass = std::string("Property_") + name;
  py::class_<Class, Property_base>(m, PyClass.c_str())
      .def(py::init<const T &>())
      .def("clone", &Class::clone)
      .def_property("value", &Class::value);
}

void declare_property_manager(py::module &m) {
  py::class_<PropertyManager>(m, "PropertyManager")
      .def(py::init<>())
      .def(py::init<const PropertyManager &>())
      .def("assign", &PropertyManager::operator=)
      .def("copy", &PropertyManager::copy)
      .def("set_property", &PropertyManager::set_property, py::arg("label"),
           py::arg("property"))
      .def("get_property", &PropertyManager::get_property, py::arg("label"))
      .def("exists_property", &PropertyManager::exists_property,
           py::arg("label"))
      .def("delete_property", &PropertyManager::delete_property);
}

void init_clipper_memory(py::module &m) {
  declare_property<std::string>(m, "string");
  declare_property<int>(m, "int");
  declare_property<float>(m, "float");
  declare_property<double>(m, "double");
  declare_property_manager(m);
}