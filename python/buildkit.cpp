// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

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

PYBIND11_MODULE(buildkit, mod)
{
  mod.doc() = "Python modules for Atomic Model BuildKit (C++) exposed via pybind11.";
  auto m_a = mod.def_submodule("clipper");
  auto m_b = mod.def_submodule("buccaneer");
  init_clipper_types(m_a);
  init_cell(m_a);
  init_spacegroup(m_a);
  init_coord_orth(m_a);
  init_minimol(m_a);
  init_maps(m_a);
  init_nxmap(m_a);

  init_map_io(m_a);
  init_ca_join(m_b);
  init_ca_filter(m_b);
  init_buccaneer_lib(m_b);
  init_ca_prep(m_b);
  init_buccaneer_prot(m_b);
  init_ca_merge(m_b);
  init_ca_find(m_b);
  init_ca_grow(m_b);
  init_ca_link(m_b);
  init_ca_sequence(m_b);
  init_ca_correct(m_b);
  init_ca_ncsbuild(m_b);
  init_ca_prune(m_b);
  init_model_tidy(m_b);
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