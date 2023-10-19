// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "version.hpp" // for buildkit version
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
void init_clipper_types(py::module &m);
void init_minimol(py::module &m);
void init_ca_join(py::module &m);
void init_cell(py::module &m);
void init_coord_orth(py::module &m);
void init_spacegroup(py::module &m);
void init_maps(py::module &m);
void init_nxmap(py::module &m);
void init_map_io(py::module &m);
void init_ca_filter(py::module &m);
void init_buccaneer_lib(py::module &m);
void init_ca_prep(py::module &m);
void init_buccaneer_prot(py::module &m);
void init_ca_merge(py::module &m);
void init_ca_find(py::module &m);
void init_ca_grow(py::module &m);
void init_ca_link(py::module &m);
void init_ca_sequence(py::module &m);
void init_ca_correct(py::module &m);
void init_ca_ncsbuild(py::module &m);
void init_ca_prune(py::module &m);
void init_model_tidy(py::module &m);

PYBIND11_MODULE(buildkit, mbk)
{
  mbk.doc() = "Python bindings to Buccaneer and Clipper library for atomic model buildkit.";
  mbk.attr("__version__") = BUILDKIT_VERSION;

  py::register_exception_translator([](std::exception_ptr p)
                                    {
      try {
        if (p) std::rethrow_exception(p);
      } catch (const clipper::Message_fatal &e) {
        PyErr_SetString(PyExc_RuntimeError, e.text().c_str());
      } });

  auto mc = mbk.def_submodule("clipper");
  auto mb = mbk.def_submodule("buccaneer");
  init_clipper_types(mc);
  init_cell(mc);
  init_spacegroup(mc);
  init_coord_orth(mc);
  init_minimol(mc);
  init_maps(mc);
  init_nxmap(mc);

  init_map_io(mc);
  init_ca_join(mb);
  init_ca_filter(mb);
  init_buccaneer_lib(mb);
  init_ca_prep(mb);
  init_buccaneer_prot(mb);
  init_ca_merge(mb);
  init_ca_find(mb);
  init_ca_grow(mb);
  init_ca_link(mb);
  init_ca_sequence(mb);
  init_ca_correct(mb);
  init_ca_ncsbuild(mb);
  init_ca_prune(mb);
  init_model_tidy(mb);
  // init_clipper_types(mod);
  // init_cell(mod);
  // init_spacegroup(mod);
  // init_coord_orth(mod);
  // init_minimol(mod);
  // init_maps(mod);
  // init_nxmap(mod);
  // init_ca_join(mod);
  // init_map_io(mod);
  // init_ca_filter(mod);
  // init_buccaneer_lib(mod);
  // init_ca_prep(mod);
  // init_buccaneer_prot(mod);
  // init_ca_merge(mod);
  // init_ca_find(mod);
  // init_ca_grow(mod);
  // init_ca_link(mod);
  // init_ca_sequence(mod);
  // init_ca_correct(mod);
  // init_ca_ncsbuild(mod);
  // init_ca_prune(mod);
  // init_model_tidy(mod);
}