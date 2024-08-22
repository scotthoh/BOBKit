// Wrapper for clipper memory for Property
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-sequence.h"
#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/clipper_memory.h>
#include <clipper/core/clipper_sysdep.h>
#include <pybind11/detail/common.h>

namespace py = pybind11;
using namespace clipper;

// trampoline class for property base
class PyProperty_base : public Property_base {
public:
  /* Inherit the constructors. */
  using Property_base::Property_base;
  /* Trampoline (need one for each virtual function) */
  Property_base *clone() const override {
    PYBIND11_OVERRIDE_PURE(Property_base *, Property_base, clone);
  }
};

void declare_property_base(py::module &m) {
  py::class_<Property_base, PyProperty_base>(
      m, "Property_base", "Base class for properties of arbitrary types.")
      .def(py::init<>())
      .def("clone", &Property_base::clone, "Factory copy method.");
}

template <class T>
void declare_property(py::module &m, const std::string &name) {

  using Class = Property<T>;
  std::string PyClass = std::string("Property_") + name;
  py::class_<Class, Property_base> property(
      m, PyClass.c_str(), "Template for a property holding an arbitrary type.");
  property
      .def(py::init<const T &>(), py::arg("value"),
           "Constructor takes contents.")
      .def("clone", &Class::clone, "Copy method.")
      .def_property_readonly("value", &Class::value,
                             "Return value of contents.")
      .def("__repr__", [=](const Class &self) {
        return "<clipper.Property_" + name + ">";
      });
}

void declare_property_manager(py::module &m) {
  py::class_<PropertyManager>(m, "PropertyManager")
      .def(py::init<>())
      .def(py::init<const PropertyManager &>(), "Copy constructor")
      .def("assign", &PropertyManager::operator=, "Assign operator")
      .def("copy", &PropertyManager::copy, "Copy manager")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property_base &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled property to the list.")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<std::string> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled string property to the list.")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<int> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled int property to the list.")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<float> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled float property to the list.")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<double> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled double property to the list.")
      .def(
          "set_property",
          [](PropertyManager &self, const std::string label,
             const Property<Ca_sequence::Sequence_data> &property) {
            self.set_property(label, property);
          },
          py::arg("label"), py::arg("property"),
          "Add a labelled Sequence_data property to the list.")
      //.def("set_property", & PropertyManager::set_property, py::arg("label"),
      // py::arg("property"))
      .def("get_property", &PropertyManager::get_property, py::arg("label"),
           py::return_value_policy::reference_internal,
           "Get a labelled property from the list.")
      .def("exists_property", &PropertyManager::exists_property,
           py::arg("label"), "Test for property.")
      .def("delete_property", &PropertyManager::delete_property,
           py::arg("label"), "Delete property.")
      .def("__repr__",
           [](const PropertyManager &self) {
             return "<clipper.PropertyManager class.>";
           })
      .doc() =
      "Class for holding a list of labelled properties of arbitrary types.\n"
      "To add a property list to an object, derive it from this class, "
      "or include a member and mirror the methods. To add a property, "
      "simply call insert_property(label,property). Properties must be "
      "objects derived from clipper::Property_base. Usually, you can just "
      "use the template form, clipper::Property<T>.";
}

void init_clipper_memory(py::module &m) {
  // these are the ones used in Buccaneer and Clipper
  declare_property_base(m);
  declare_property<std::string>(m, "string");
  declare_property<int>(m, "int");
  declare_property<float>(m, "float");
  declare_property<double>(m, "double");
  declare_property<Ca_sequence::Sequence_data>(m, "sequence_data");
  declare_property_manager(m);
}