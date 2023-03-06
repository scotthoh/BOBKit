/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */

#ifndef K2C_HELPER_FUNCTIONS_H_
#define K2C_HELPER_FUNCTIONS_H_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "k2c_tensor_include.h"
class k2c_helper_functions {
public:
	k2c_helper_functions(){};
	virtual ~k2c_helper_functions();
	static void k2c_matmul(float * C, const float * A, const float * B, const size_t outrows,
	                const size_t outcols, const size_t innerdim);
	static void k2c_affine_matmul(float * C, const float * A, const float * B, const float * d,
	                       const size_t outrows,const size_t outcols, const size_t innerdim);


	static void k2c_dot(k2c_tensor* C, const k2c_tensor* A, const k2c_tensor* B, const size_t * axesA,
	             const size_t * axesB, const size_t naxes, const int normalize, float * fwork);

	static void k2c_bias_add(k2c_tensor* A, const k2c_tensor* b);
	static size_t k2c_sub2idx(const size_t * sub, const size_t * shape, const size_t ndim);
	static void k2c_idx2sub(const size_t idx, size_t * sub, const size_t * shape, const size_t ndim);
	static float* k2c_read_array(const char* layername,const size_t array_size  ,const char* filename);
	static float* k2c_zero_array(const size_t array_size);
};

#endif /* K2C_HELPER_FUNCTIONS_H_ */
