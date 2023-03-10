#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <clipper/clipper.h>

#include <vector>

namespace py = pybind11;
using namespace clipper;

template <typename T>
int normalise_index(int index, const T &container)
{
  if (index < 0)
    index += (int)container.size();
  if ((size_t)index >= container.size())
    throw pybind11::index_error();
  return index;
}

template <class C, typename T>
py::array_t<T> to_numpy_1d(const C &v, int n)
{
  py::array_t<T> ret(n);
  auto buf = ret.request();
  T *ptr = (T *)buf.ptr;
  for (int i = 0; i < n; ++i)
    ptr[i] = v[i];
  return ret;
}

template <class C, typename T>
void to_numpy_1d(const C &v, int n, py::array_t<T> target)
{
  auto buf = target.request();
  T *ptr = (T *)buf.ptr;
  for (int i = 0; i < n; ++i)
  {
    ptr[i] = v[i];
  }
}
