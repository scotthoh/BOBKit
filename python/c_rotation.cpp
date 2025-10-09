// Nanobind bindings for clipper rotation
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/core/rotation.h>
#include <nanobind/operators.h>

using namespace clipper;

// template loop to bind Rotation constructors from Euler code 0 to 24,
// and method to return euler angles with specified code
template <int N> void bind_euler_construct( nb::class_<Rotation> &pyclass ) {
  auto eulercode = std::string( "euler" ) + std::to_string( N );
  pyclass.def( nb::init<const Euler<N> &>(), nb::arg( "rot" ), "Constructor: copy/convert." )
      .def(
          eulercode.c_str(), []( const Rotation &self ) { return self.euler<N>(); },
          "Return Euler angles with specified integer convention code." );
  bind_euler_construct<N - 1>( pyclass );
}

// for N=-1, do nothing, stop recursion
template <> void bind_euler_construct<-1>( nb::class_<Rotation> &pyclass ) {}

template <int T> void declare_euler( nb::module_ &m, const std::string &name ) {
  using Class = Euler<T>;
  // auto pyclassname = std::string("Euler_") + std::to_string(T);
  // nb::class_<Class>(m, pyclassname.c_str())
  nb::class_<Class>( m, ( "Euler_" + name ).c_str() )
      .def( nb::init<>(), "Null constructor." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "alpha" ), nb::arg( "beta" ),
            nb::arg( "gamma" ), "Constructor: from specified angles." )
      .def( nb::init<const Rotation &>(), nb::arg( "rot" ), "Constructor: from rotation" )
      .def( "rotation", &Class::rotation, "Return rotation." )
      .def_prop_ro( "alpha", &Class::alpha, "Return alpha." )
      .def_prop_ro( "beta", &Class::beta, "Return beta." )
      .def_prop_ro( "gamma", &Class::gamma, "Return gamma." )
      .def( "format", &Class::format, "Return formatted string representation." )
      .doc() = "Euler angles class\nRotations are generally handled through the clipper::Rotation class."
               "This class only exists for conversion purposes. "
               "This particular class represents generic Euler angles. The  "
               "convention is selected from the 24 possible conventions according "
               "to the template parameter. The integer convention code is "
               "enumerated in the Rotation::EULERtype enumation in the form "
               "Rotation::EulerZYZr, Rotation::EulerXYZs etc., where the X/Y/Z  "
               "indicates the axes of rotation in order, and the r/s indicates  "
               "static or rotating axes. The type of an Euler class is also given "
               "as a prefix to the result of format(). ";
  declare_euler<T - 1>( m, std::to_string( T - 1 ) );
}

template <> void declare_euler<-1>( nb::module_ &m, const std::string &name ) {}

void declare_euler_ccp4( nb::module_ &m ) {
  nb::class_<Euler_ccp4>( m, "Euler_ccp4" )
      .def( nb::init<>(), "Null constructor." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "alpha" ), nb::arg( "beta" ),
            nb::arg( "gamma" ), "Constructor: from specified angles." )
      .def_prop_ro( "alpha", &Euler_ccp4::alpha, "Return alpha." )
      .def_prop_ro( "beta", &Euler_ccp4::beta, "Return beta." )
      .def_prop_ro( "gamma", &Euler_ccp4::gamma, "Return gamma." )
      .def( "format", &Euler_ccp4::format, "Return formatted string representation." )
      .doc() = "Euler_ccp4 angles class\nRotations are generally handled through the clipper::Rotation class. "
               "This class only exists for conversion purposes. "
               "This particular class represents Euler_ccp4 angles according to the "
               "CCP4 standard, i.e. "
               "- Rotation 1 (alpha) about K, "
               "- Rotation 2 (beta) about the new J, "
               "- Rotation 3 (gamma) about the new K. ";
}

void declare_polar_ccp4( nb::module_ &m ) {
  nb::class_<Polar_ccp4>( m, "Polar_ccp4" )
      .def( nb::init<>(), "Null constructor." )
      .def( nb::init<const double &, const double &, const double &>(), nb::arg( "omega" ), nb::arg( "phi" ),
            nb::arg( "kappa" ), "Constructor: from specified angles." )
      .def_prop_ro( "psi", &Polar_ccp4::psi, "Return psi." )
      .def_prop_ro( "omega", &Polar_ccp4::omega, "Return omega." )
      .def_prop_ro( "phi", &Polar_ccp4::phi, "Return phi." )
      .def_prop_ro( "kappa", &Polar_ccp4::kappa, "Return kappa." )
      .def( "format", &Polar_ccp4::format, "Return formatted string representation." )
      .doc() = "Polar_ccp4 angle class\n"
               "Rotations are generally handled through the clipper::Rotation class. "
               "This class only exists for conversion purposes. "
               "This particular class represents Polar_ccp4 angles according to the "
               "CCP4 standard, i.e. "
               "- omega gives inclination of rotation axis to K axis, "
               "- phi gives anticlockwise rotation from I to projection of "
               "rotation axis onto I-J plane, "
               "- kappa is the rotation about the rotation axis. ";
}

void declare_rotation( nb::module_ &m ) {
  nb::class_<Rotation> rotation( m, "Rotation" );

  nb::enum_<Rotation::EULERtype>( rotation, "EULERtype", "Enumeration of Euler conventions." )
      .value( "EulerXYZr", Rotation::EULERtype::EulerXYZr )
      .value( "EulerXYZs", Rotation::EULERtype::EulerXYZs )
      .value( "EulerXYXr", Rotation::EULERtype::EulerXYXr )
      .value( "EulerXYXs", Rotation::EULERtype::EulerXYXs )
      .value( "EulerXZXr", Rotation::EULERtype::EulerXZXr )
      .value( "EulerXZXs", Rotation::EULERtype::EulerXZXs )
      .value( "EulerXZYr", Rotation::EULERtype::EulerXZYr )
      .value( "EulerXZYs", Rotation::EULERtype::EulerXZYs )
      .value( "EulerYZXr", Rotation::EULERtype::EulerYZXr )
      .value( "EulerYZXs", Rotation::EULERtype::EulerYZXs )
      .value( "EulerYZYr", Rotation::EULERtype::EulerYZYr )
      .value( "EulerYZYs", Rotation::EULERtype::EulerYZYs )
      .value( "EulerYXYr", Rotation::EULERtype::EulerYXYr )
      .value( "EulerYXYs", Rotation::EULERtype::EulerYXYs )
      .value( "EulerYXZr", Rotation::EULERtype::EulerYXZr )
      .value( "EulerYXZs", Rotation::EULERtype::EulerYXZs )
      .value( "EulerZXYr", Rotation::EULERtype::EulerZXYr )
      .value( "EulerZXYs", Rotation::EULERtype::EulerZXYs )
      .value( "EulerZXZr", Rotation::EULERtype::EulerZXZr )
      .value( "EulerZXZs", Rotation::EULERtype::EulerZXZs )
      .value( "EulerZYZr", Rotation::EULERtype::EulerZYZr )
      .value( "EulerZYZs", Rotation::EULERtype::EulerZYZs )
      .value( "EulerZYXr", Rotation::EULERtype::EulerZYXr )
      .value( "EulerZYXs", Rotation::EULERtype::EulerZYXs )
      .export_values();

  bind_euler_construct<23>( rotation ); // bind all 24 templated constructors and methods in Rotation class
  rotation.def( nb::init<>(), "Null constructor" )
      .def( nb::init<const Euler_ccp4 &>(), nb::arg( "euler" ), "Constructor: from Euler_ccp4." )
      .def( nb::init<const Polar_ccp4 &>(), nb::arg( "polar" ), "Constructor: from Polar_ccp4." )
      .def( nb::init<const Mat33<> &>(), nb::arg( "matrix" ), "Constructor: from matrix." )
      .def( nb::init<const double &, const double &, const double &, const double &>(), nb::arg( "w" ), nb::arg( "x" ),
            nb::arg( "y" ), nb::arg( "z" ), "Constructor: from components." )
      .def_prop_ro( "w", &Rotation::w, "Return w component." )
      .def_prop_ro( "x", &Rotation::x, "Return x component." )
      .def_prop_ro( "y", &Rotation::y, "Return y component." )
      .def_prop_ro( "z", &Rotation::z, "Return z component." )
      //.def("euler", &Rotation::euler<0>, "Return Euler angles.")
      .def( "euler_ccp4", &Rotation::euler_ccp4, "Return Euler_ccp4 angles." )
      .def( "polar_ccp4", &Rotation::polar_ccp4, "Return Polar_ccp4 angles." )
      .def( "norm", &Rotation::norm, "Normalise this quaternion." )
      .def( "abs_angle", &Rotation::abs_angle, "Return positive magnitude of the angle of rotation." )
      .def( "matrix", &Rotation::matrix, "Return 3x3 matrix." )
      .def( "inverse", &Rotation::inverse, "Return inverse rotation." )
      .def_static( "zero", &Rotation::zero, "Return zero rotation." )
      .def_static( "null", &Rotation::null, "Return null rotation." )
      .def( "is_null", &Rotation::is_null, "Test if object has been initialised." )
      .def( nb::self * nb::self, "Combine two rotations")
      .def( "format", &Rotation::format, "Return formatted string representation." )
      .doc() = "Rotation class\n"
               "This class represents a rotation. The internal representation is as a unit quaternion, "
               "which is easily combined, inverted, or converted to or from other commonly used forms.";
}

void add_rotation( nb::module_ &m ) {
  // declare euler templates 0 to 23
  declare_euler<23>( m, "23" );
  // declare_euler<0>(m, "0"); declare_euler<1>(m, "1"); declare_euler<2>(m, "2"); declare_euler<3>(m, "3");
  // declare_euler<4>(m, "4"); declare_euler<5>(m, "5"); declare_euler<6>(m, "6"); declare_euler<7>(m, "7");
  // declare_euler<8>(m, "8"); declare_euler<9>(m, "9"); declare_euler<10>(m, "10"); declare_euler<11>(m, "11");
  // declare_euler<12>(m, "12"); declare_euler<13>(m, "13"); declare_euler<14>(m, "14"); declare_euler<15>(m, "15");
  // declare_euler<16>(m, "16"); declare_euler<17>(m, "17"); declare_euler<18>(m, "18"); declare_euler<19>(m, "19");
  // declare_euler<20>(m, "20"); declare_euler<21>(m, "21"); declare_euler<22>(m, "22"); declare_euler<23>(m, "23");
  declare_euler_ccp4( m );
  declare_polar_ccp4( m );
  declare_rotation( m );
}