// Wrapper for clipper_thread
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/core/clipper_thread.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

// trampoline class for thread base
class Trampoline_Thread_base : public Thread_base {
 public:
  /* Inherit the constructors. */
  using Thread_base::Thread_base;
  /* Trampoline (need one for each virtual function) */
  void Run() override { PYBIND11_OVERRIDE_PURE( void, Thread_base, Run ); }
};

class PyThread_base : public Thread_base {
 public:
  using Thread_base::Run;
};

void init_thread_base( py::module& m ) {
  py::class_<Thread_base, Trampoline_Thread_base>( m, "Thread_base" )
      .def( py::init<>() )
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
