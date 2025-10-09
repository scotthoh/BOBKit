#pragma once
#include "commons.h"
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/vector.h>

// for 1d arrays
template <typename T> using cpu_c_array = nb::ndarray<T, nb::shape<-1>, nb::device::cpu, nb::c_contig>;
// for 3d arrays
template <typename T> using np_cpu_c_3darray = nb::ndarray<T, nb::ndim<3>, nb::device::cpu, nb::c_contig>;
template <typename T> using np_cpu_f_3darray = nb::ndarray<T, nb::ndim<3>, nb::device::cpu, nb::f_contig>;
// for matrices
template <typename T, size_t M, size_t N>
using cpu_c_2darray = nb::ndarray<T, nb::shape<M, N>, nb::device::cpu, nb::c_contig>;

// make numpy array with given size/shape
template <typename T>
auto make_numpy_array(std::initializer_list<size_t> shape, std::initializer_list<int64_t> strides={}) {
  size_t array_size = 1;
  for ( size_t i : shape )
    array_size *= i;  
  // allocate memory region
  T* arr = new T[array_size];
  // Delete 'data' when the 'owner' capsule expires
  nb::capsule owner(arr, [](void *p) noexcept {
    delete [] static_cast<T*>(p);
  });
  return nb::ndarray<T, nb::numpy, nb::ndim<3>, nb::device::cpu, nb::c_contig>(arr, shape, owner, strides);
}

template <typename T>
auto make_numpy_array_norv(std::initializer_list<size_t> shape, std::initializer_list<int64_t> strides={}) {
  size_t array_size = 1;
  for ( size_t i : shape )
    array_size *= i;  
  // allocate memory region
  T* arr = new T(array_size);
  // Delete 'data' when the 'owner' capsule expires
  //nb::capsule owner(arr, [](void *p) noexcept {
  //  delete [] static_cast<T*>(p);
  //});
  return nb::ndarray<T, nb::numpy, nb::ndim<3>>(arr, shape, nullptr, strides);
  //ßreturn view.cast();
}

// 23 sept read in ndarray and turn into xmap/nxmap