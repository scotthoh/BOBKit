/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */

#ifndef K2C_ACTIVATIONS_H_
#define K2C_ACTIVATIONS_H_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
class k2c_activations {
public:
	k2c_activations();
	virtual ~k2c_activations();
	static void k2c_tanh_func(float * x, const size_t size);
	static void k2c_sigmoid_func(float * x, const size_t size);


};

#endif /* K2C_CORE_LAYERS_H_ */
