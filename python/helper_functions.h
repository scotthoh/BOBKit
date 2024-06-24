#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <clipper/clipper.h>

#include <vector>

namespace py = pybind11;
using namespace clipper;

template <typename T>
void delitem_at_index(T &container, pybind11::ssize_t index)
{
  container.erase(container.begin() + index);
}

template <typename Item>
void delitem_at_index(std::vector<Item> &items, pybind11::ssize_t index)
{
  items.erase(items.begin() + index);
}

template <typename T>
void delitem_range(T &container, pybind11::ssize_t start, pybind11::ssize_t end)
{
  container.erase(container.begin() + start, container.begin() + end);
}

template <typename Item>
void delitem_range(std::vector<Item> &items, pybind11::ssize_t start, pybind11::ssize_t end)
{
  items.erase(items.begin() + start, items.begin() + end);
}

template <typename Items>
void delitem_slice(Items &items, const pybind11::slice &slice)
{
  py::ssize_t start, stop, step, slice_len;
  if (!slice.compute((py::ssize_t)items.size(), &start, &stop, &step, &slice_len))
    throw py::error_already_set();
  if (step == 1)
  {
    delitem_range(items, start, start + slice_len);
  }
  else
  {
    for (int i = 0; i < slice_len; ++i)
      delitem_at_index(items, start + (step > 0 ? slice_len - 1 - i : i) * step);
  }
}

template <typename T, typename C>
C &add_item(T &container, C child, int pos)
{
  if ((size_t)pos > container.size()) // true for negative pos
    pos = (int)container.size();
  return *container.insert(container.begin() + pos, std::move(child));
}

template <typename P, typename C>
C &add_child(P &parent, C child, int pos)
{
  return add_item(parent.children(), std::move(child), pos);
}

template <typename P, typename C>
C &get_child(P &parent, int index)
{
  auto &children = parent.children();
  return children[normalise_index(index, children.size())];
}

template <typename P, typename C>
void set_child(P &parent, int index, C &child)
{
  auto &children = parent.children();
  children[normalise_index(index, children.size())] = child;
}

template <typename P>
void remove_child(P &parent, int index)
{
  auto &children = parent.children();
  children.erase(children.begin() + normalise_index(index, children.size()));
}

template <typename P>
void remove_children(P &parent, py::slice slice)
{
  delitem_slice(parent.children(), slice);
}

template <typename T>
int normalise_index(int index, const T &container_size)
{
  if (index < 0)
    index += (int)container_size;
  if ((size_t)index >= container_size)
    throw pybind11::index_error("Index out of bounds");
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

// template <class C, typename T>
// py::array_t<T> CoordList_to_numpy(C &v, size_t n, size_t m)
//{
//   // size_t l = v.size();
//   py::array_t<T, py::array::c_style> ret({n, m});
//   // T buf = ret.request();
//   // T *ptr = (T *)ret.request().ptr;
//   auto ra = ret.mutable_unchecked();
//   for (size_t i = 0; i < n; i++)
//   {
//     // for (size_t j = 0; j < n; j++)
//     //{
//     ra(i, 0) = v[i].x();
//     ra(i, 1) = v[i].y();
//     ra(i, 2) = v[i].z();
//     //}
//   }
//   return ret;
// }

static std::string triple(double x, double y, double z)
{
  using namespace std; // VS2015/17 doesn't like std::snprintf
  char buf[128];
  auto r = [](double d)
  { return std::fabs(d) >= 1e-15 ? d : 0; };
  snprintf(buf, 128, "%g, %g, %g", r(x), r(y), r(z));
  return std::string(buf);
}

// check that the given array matches expected dimensions, and throw an
// error if not. direction is true for incoming, false for outgoing.
template <typename T>
void check_array_shape(py::array_t<T> target, std::vector<int> dim, bool direction)
{
  py::buffer_info buf = target.request();
  bool fail = false;
  if ((size_t)buf.ndim != dim.size())
    fail = true;
  else
    for (ssize_t i = 0; i < buf.ndim; ++i)
      if (buf.shape[i] != dim[i])
        fail = true;
  if (fail)
  {
    auto shapet_txt = std::string("( ");
    auto shapeb_txt = std::string("( ");
    for (ssize_t i = 0; i < buf.ndim; ++i)
      shapet_txt += std::to_string(buf.shape[i]) + " ";
    for (size_t i = 0; i < dim.size(); ++i)
      shapeb_txt += std::to_string(dim[i]) + " ";
    shapet_txt += ")";
    shapeb_txt += ")";

    std::string message;

    auto msg = std::string("Array shape mismatch! ");
    if (direction)
      msg += "Input array shape is ";
    else
      msg += "Target array shape is ";
    msg += shapet_txt;
    msg += ", while expected shape is ";
    msg += shapeb_txt;
    msg += ".";
    throw std::runtime_error(msg.c_str());
  }
}

template <class C, typename T>
py::array_t<T> make_array_1d(const C &v, int n)
{
  py::array_t<T> arr(n);
  py::buffer_info buf = arr.request();
  T *ptr = (T *)buf.ptr;
  for (int i = 0; i < n; ++i)
    ptr[i] = v[i];
  return arr;
}

template <class C, typename T>
void array_as_numpy_1d(const C &v, int n, py::array_t<T> target)
{
  check_array_shape(target, {n}, true);
  auto buf = target.request();
  T *ptr = (T *)buf.ptr;
  for (int i = 0; i < n; ++i)
    ptr[i] = v[i];
}

template <class C, typename T>
void fill_array_1d(C &v, int n, py::array_t<T> arr)
{
  check_array_shape(arr, {n}, true);
  auto buf = arr.request();
  T *ptr = (T *)buf.ptr;
  for (int i = 0; i < n; ++i)
    v[i] = ptr[i];
}

template <class HKLdtype, class dtype>
py::array_t<dtype> hkl_data_export_numpy(const HKLdtype &self, const int &size)
{
  auto ret = py::array_t<dtype>(size);
  dtype *ptr = (dtype *)ret.request().ptr;
  self.data_export(ptr);
  return ret;
} // hkl_data_export_numpy

template <class HKLdtype, class dtype>
void hkl_data_import_numpy(HKLdtype &self, const int &size, py::array_t<dtype> vals)
{
  check_array_shape(vals, {size}, true);
  dtype *ptr = (dtype *)vals.request().ptr;
  self.data_import(ptr);
} // hkl_data_import_numpy

template <typename T>
std::unique_ptr<Vec3<T>> new_vec3_from_numpy(py::array_t<T> vals)
{
  check_array_shape(vals, {3}, true);
  return std::unique_ptr<Vec3<T>>(new Vec3<T>(vals.at(0), vals.at(1), vals.at(2)));
}

template <typename T>
std::unique_ptr<Mat33<T>> new_mat33_from_numpy(py::array_t<T> vals)
{
  check_array_shape(vals, {3, 3}, true);
  auto r = vals.template unchecked<2>();
  return std::unique_ptr<Mat33<T>>(new Mat33<T>(
      r(0, 0), r(0, 1), r(0, 2), r(1, 0), r(1, 1), r(1, 2), r(2, 0), r(2, 1), r(2, 2)));
}