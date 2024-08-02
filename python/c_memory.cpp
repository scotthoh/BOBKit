#include "buccaneer/buccaneer-sequence.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/clipper_memory.h>
#include <clipper/core/clipper_sysdep.h>
#include <pybind11/detail/common.h>

namespace py = pybind11;
using namespace clipper;

template <class T>
void declare_property(py::module &m, const std::string &name) {

  using Class = Property<T>;
  std::string PyClass = std::string("Property_") + name;
  py::class_<Class> property(m, PyClass.c_str());
  property.def(py::init<const T &>())
      .def("clone", &Class::clone)
      .def_property_readonly("value", &Class::value)
      .def("__repr__", [=](const Class &self) {
        return "<clipper.Property_" + name + ">";
      });
}

void declare_property_manager(py::module &m) {
  py::class_<PropertyManager>(m, "PropertyManager")
      .def(py::init<>())
      .def(py::init<const PropertyManager &>())
      .def("assign", &PropertyManager::operator=)
      .def("copy", &PropertyManager::copy)
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<std::string> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"))
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<int> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"))
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<float> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"))
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<double> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"))
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<Ca_sequence::Sequence_data> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"))
      //.def("set_property", & PropertyManager::set_property, py::arg("label"),
      // py::arg("property"))
      .def("get_property", &PropertyManager::get_property, py::arg("label"),
           py::return_value_policy::reference_internal)
      .def("exists_property", &PropertyManager::exists_property,
           py::arg("label"))
      .def("delete_property", &PropertyManager::delete_property)
      .def("__repr__", [](const PropertyManager &self) {
        return "<clipper.PropertyManager class.>";
      });
}

void init_clipper_memory(py::module &m) {
  declare_property<std::string>(m, "string");
  declare_property<int>(m, "int");
  declare_property<float>(m, "float");
  declare_property<double>(m, "double");
  declare_property<Ca_sequence::Sequence_data>(m, "sequence_data");
  declare_property_manager(m);
}