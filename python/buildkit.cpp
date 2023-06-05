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

PYBIND11_MODULE(buildkit, mod)
{
  mod.doc() = "Python wrapper for Buccaneer BuildKit (C++) exposed via pybind11.";
  init_clipper_types(mod);
  init_cell(mod);
  init_spacegroup(mod);
  init_coord_orth(mod);
  init_minimol(mod);
  init_maps(mod);
  init_nxmap(mod);
  init_ca_join(mod);
  init_map_io(mod);
  init_ca_filter(mod);
  init_buccaneer_lib(mod);
}