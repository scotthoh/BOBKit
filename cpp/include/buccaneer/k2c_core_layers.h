/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */

#ifndef K2C_CORE_LAYERS_H_
#define K2C_CORE_LAYERS_H_
#include "k2c_tensor_include.h"
#include "k2c_include.h"
class k2c_core_layers {
public:
	k2c_core_layers();
	virtual ~k2c_core_layers();
	static void k2c_dense(k2c_tensor* output, const k2c_tensor* input, const k2c_tensor* kernel,
	               const k2c_tensor* bias, k2c_activationType *activation, float * fwork);
};

#endif /* K2C_CORE_LAYERS_H_ */
