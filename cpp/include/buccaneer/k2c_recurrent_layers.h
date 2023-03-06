/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */

#ifndef K2C_RECURRENT_LAYERS_H_
#define K2C_RECURRENT_LAYERS_H_
#include "k2c_tensor_include.h"
#include "k2c_include.h"
class k2c_recurrent_layers {
public:
	k2c_recurrent_layers();
	virtual ~k2c_recurrent_layers();
	static void k2c_lstmcell(float * state, const float * input, const k2c_tensor* kernel,
	                  const k2c_tensor* recurrent_kernel, const k2c_tensor* bias, float * fwork,
	                  k2c_activationType *recurrent_activation,
	                  k2c_activationType *output_activation);

	static void k2c_lstm(k2c_tensor* output, const k2c_tensor* input, float * state,
	              const k2c_tensor* kernel, const k2c_tensor* recurrent_kernel,
	              const k2c_tensor* bias, float * fwork, const int go_backwards,
	              const int return_sequences, k2c_activationType *recurrent_activation,
	              k2c_activationType *output_activation);
};

#endif /* K2C_RECURRENT_LAYERS_H_ */
