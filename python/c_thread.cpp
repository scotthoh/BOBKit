// Nanobind bindings for clipper xmap
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/core/clipper_thread.h>
#include <nanobind/trampoline.h>

using namespace clipper;

// trampoline class for thread base
class Trampoline_Thread_base : public Thread_base {
 public:
  /* Inherit the constructors. */
  NB_TRAMPOLINE( Thread_base, 1 );
  /* Trampoline (need one for each virtual function) */
  void Run() override { NB_OVERRIDE_PURE( Run ); }
};

class PyThread_base : public Thread_base {
 public:
  using Thread_base::Run;
};

void add_thread_base( nb::module_ & m ) {
  nb::class_<Thread_base, Trampoline_Thread_base>( m, "Thread_base" )
      .def( nb::init<>() )
      .def( "run", &Thread_base::run )
      .def( "join", &Thread_base::join )
      .def_static( "lock", &Thread_base::lock )
      .def_static( "unlock", &Thread_base::unlock )
      .def( "id", &Thread_base::id )
      .def( "Run", &PyThread_base::Run )
      .doc() =
      "Thread base class: Override this to create new threads.\n"
      "Store data as members with accessors to set input and read output. "
      "Override the Run() method to do the actual work. ";
}
