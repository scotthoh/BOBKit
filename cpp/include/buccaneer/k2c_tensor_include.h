/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */

#ifndef K2C_TENSOR_INCLUDE_H_
#define K2C_TENSOR_INCLUDE_H_
#include <stdlib.h>
#define K2C_MAX_NDIM 5
struct k2c_tensor
{
    /** Pointer to array of tensor values flattened in row major order. */
    float * array;

    /** Rank of the tensor (number of dimensions). */
    size_t ndim;

    /** Number of elements in the tensor. */
    size_t numel;

    /** Array, size of the tensor in each dimension. */
    size_t shape[K2C_MAX_NDIM];
};
#endif /* K2C_TENSOR_INCLUDE_H_ */
