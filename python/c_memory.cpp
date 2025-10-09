// Nanobind bindings for clipper_memory; Property in MiniMol
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <buccaneer/buccaneer-sequence.h>
#include <clipper/core/clipper_memory.h>
#include <clipper/core/clipper_sysdep.h>
#include <nanobind/trampoline.h>

using namespace clipper;

// trampoline class for property base
class PyProperty_base : public Property_base {
public:
  /* Inherit the constructors. */
  // using Property_base::Property_base;
  NB_TRAMPOLINE( Property_base, 1 );

  /* Trampoline (need one for each virtual function) */
  Property_base *clone() const override { NB_OVERRIDE_PURE( clone ); }
};

void declare_property_base( nb::module_ &m ) {
  nb::class_<Property_base, PyProperty_base>( m, "Property_base", "Base class for properties of arbitrary types." )
      .def( nb::init<>() )
      .def( "clone", &Property_base::clone, "Factory copy method." );
}

template <class T> void declare_property( nb::module_ &m, const std::string &name ) {

  using Prop = Property<T>;
  nb::class_<Prop, Property_base> property( m, ( "Property_" + name ).c_str(),
                                            "Template for a property holding an arbitrary type." );
  property.def( nb::init<const T &>(), nb::arg( "value" ), "Constructor takes contents." )
      .def( "clone", &Prop::clone, "Copy method." )
      .def_prop_ro( "value", &Prop::value, "Return value of contents." )
      .def( "__repr__", [=]( const Prop &self ) { return "<clipper.Property_" + name + ">"; } );
}

// To include other types, the following example can be added to your bindings in another file.
// Replace double with the type and type name
// Also refer c_example_property_type.cpp
// void add_property_double(nb::module_& m) {
//    using Prop = Property<double>;
//    nb::class_<Prop, Property_base> property(
//        m, "Property_double", "Template for a property holding an arbitrary type.");
//    property
//        .def(nb::init<const double &>(), nb::arg("value"),
//             "Constructor takes contents.")
//        .def("clone", &Prop::clone, "Copy method.")
//        .def_prop_ro("value", &Prop::value,
//                               "Return value of contents.")
//        .def("__repr__", [=](const Prop &self) {
//          return "<clipper.Property_double>";
//        });
//

void declare_property_manager( nb::module_ &m ) {
  nb::class_<PropertyManager>( m, "PropertyManager" )
      .def( nb::init<>() )
      .def( nb::init<const PropertyManager &>(), "Copy constructor" )
      .def( "assign", &PropertyManager::operator=, "Assign operator" )
      .def( "copy", &PropertyManager::copy, "Copy manager" )
      .def(
          "set_property",
          []( PropertyManager &self, const std::string label, const Property_base &property ) {
            self.set_property( label, property );
          },
          nb::arg( "label" ), nb::arg( "property" ), "Add a labelled property to the list." )
      .def( "get_property", &PropertyManager::get_property, nb::arg( "label" ), nb::rv_policy::reference_internal,
            "Get a labelled property from the list." )
      .def( "exists_property", &PropertyManager::exists_property, nb::arg( "label" ), "Test for property." )
      .def( "delete_property", &PropertyManager::delete_property, nb::arg( "label" ), "Delete property." )
      .def( "__repr__", []( const PropertyManager &self ) { return "<clipper.PropertyManager class.>"; } )
      .doc() = "Class for holding a list of labelled properties of arbitrary types.\n"
               "To add a property list to an object, derive it from this class, "
               "or include a member and mirror the methods. To add a property, "
               "simply call insert_property(label,property). Properties must be "
               "objects derived from clipper::Property_base. Usually, you can just "
               "use the template form, clipper::Property<T>.";
}

void add_clipper_memory( nb::module_ &m ) {
  // These are the ones used in Clipper
  // Other types can be added by using the example binding above if needed.
  declare_property_base( m );
  declare_property<std::string>( m, "string" );
  declare_property<int>( m, "int" );
  declare_property<float>( m, "float" );
  declare_property<double>( m, "double" );
  declare_property<Ca_sequence::Sequence_data>(m, "sequence_data");
  declare_property_manager( m );
}