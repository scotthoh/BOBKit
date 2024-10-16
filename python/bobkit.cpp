// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "version.hpp" // for bobkit version
#include <clipper/clipper.h>
#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_clipper_types(py::module &m);
void init_ramachandran(py::module &m);
void init_clipper_memory(py::module &m);
void init_matomindex(py::module &m);
void init_minimol(py::module &m);
void init_ca_join(py::module &m);
void init_cell(py::module &m);
void init_coord_orth(py::module &m);
void init_atomlist(py::module &m);
void init_spacegroup(py::module &m);
void init_hklinfo(py::module &m);
void init_hkl_data(py::module &m);
void init_hkl_datatypes(py::module &m);
void init_sfscale(py::module &m);
void init_containers(py::module &m);
void init_maps(py::module &m);
void init_nxmap(py::module &m);
void init_map_io(py::module &m);
void init_symop(py::module &m);
void init_gemmi_structure(py::module &m);
void init_minimol_seq(py::module &m);
void init_clipper_stats(py::module &m);
void init_map_utils(py::module &m);
void init_clipper_util(py::module &m);
// void init_ccp4_mtz_io(py::module &m);
void init_clipper_tests(py::module &m);

void init_ca_build(py::module &m);
void init_ca_filter(py::module &m);
void init_buccaneer_lib(py::module &m);
void init_ca_prep(py::module &m);
void init_buccaneer_prot(py::module &m);
void init_protein_loop(py::module &m);
void init_ca_merge(py::module &m);
void init_ca_find(py::module &m);
void init_ca_grow(py::module &m);
void init_ca_link(py::module &m);
void init_ca_sequence(py::module &m);
void init_ca_correct(py::module &m);
void init_ca_ncsbuild(py::module &m);
void init_ca_prune(py::module &m);
void init_knownstructure(py::module &m);
void init_model_tidy(py::module &m);
void init_buccaneer_util(py::module &m);
void init_simplex_lib(py::module &m);
void init_map_simulate(py::module &m);
void init_proteindb(py::module &m);

PYBIND11_MODULE(bobkit, mbk) {
  mbk.doc() =
      "Python bindings to Buccaneer and Clipper library for biomacromolecular "
      "model building kit (BOBKIT).";
  mbk.attr("__version__") = BOBKIT_VERSION;
  // auto package = pybind11::module::import("gemmi");
  // auto module = package.attr("Structure");
  // mbk.add_object("Structure", module);

  py::register_exception_translator([](std::exception_ptr p)
                                    {
      try {
        if (p) std::rethrow_exception(p);
      } catch (const clipper::Message_fatal &e) {
        PyErr_SetString(PyExc_RuntimeError, e.text().c_str());
      } });

  mbk.def("long_running_func", []() {
    for (;;) {
      if (PyErr_CheckSignals() != 0)
        throw py::error_already_set();
      // Long running iteration
    }
  });

  auto mc = mbk.def_submodule("clipper");
  auto mb = mbk.def_submodule("buccaneer");
  auto mpdb = mbk.def_submodule("protein_db");

  init_clipper_types(mc);
  init_ramachandran(mc);
  init_clipper_memory(mc);
  init_cell(mc);
  init_symop(mc);
  init_spacegroup(mc);
  init_coord_orth(mc);
  init_matomindex(mc);
  init_atomlist(mc);
  init_hklinfo(mc);
  init_containers(mc);
  init_hkl_datatypes(mc);
  init_hkl_data(mc);
  init_sfscale(mc);

  init_minimol(mc);
  init_maps(mc);
  init_nxmap(mc);
  init_gemmi_structure(mc);
  init_minimol_seq(mc);
  init_map_io(mc);
  init_clipper_stats(mc);
  init_map_utils(mc);
  init_clipper_util(mc);
  // init_ccp4_mtz_io(mc);
  init_clipper_tests(mc);

  init_simplex_lib(mb);
  init_map_simulate(mb);
  init_ca_build(mb);
  init_ca_join(mb);
  init_ca_filter(mb);
  init_buccaneer_lib(mb);
  init_buccaneer_util(mb);
  init_ca_prep(mb);
  init_buccaneer_prot(mb);
  init_protein_loop(mb);
  init_ca_merge(mb);
  init_ca_find(mb);
  init_ca_grow(mb);
  init_ca_link(mb);
  init_ca_sequence(mb);
  init_ca_correct(mb);
  init_ca_ncsbuild(mb);
  init_ca_prune(mb);
  init_knownstructure(mb);
  init_model_tidy(mb);

  init_proteindb(mpdb);
}