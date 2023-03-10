// Wrapper for buccaneer-join.cpp
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
void init_minimol(py::module &m);
void init_ca_join(py::module &m);
void init_cell(py::module &m);
void init_coord_orth(py::module &m);
void init_spacegroup(py::module &m);

PYBIND11_MODULE(buildkit, mod)
{
  mod.doc() = "Python wrapper for Buccaneer BuildKit (C++) exposed via pybind11.";
  init_cell(mod);
  init_spacegroup(mod);
  init_coord_orth(mod);
  init_minimol(mod);
  init_ca_join(mod);
}