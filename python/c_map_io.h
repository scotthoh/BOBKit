#ifndef BUILDKIT_MAP_IO_H_INCLUDED
#define BUILDKIT_MAP_IO_H_INCLUDED

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "type_conversions.h"
#include <clipper/clipper.h>
#include "helper_functions.h"

namespace buildkit
{

  template <class T>
  void numpy_to_nxmap(pybind11::array_t<T, 3> data, const clipper::Cell &cell, clipper::NXmap<T> &nxmap);
  template <class T>
  void numpy_to_xmap(pybind11::array_t<T, 3> data, const clipper::Cell &cell, clipper::Xmap<T> &xmap);

}

#endif