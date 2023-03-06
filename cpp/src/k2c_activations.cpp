/**
k2c_activations.c
This file is part of keras2c
Copyright 2020 Rory Conlin
Licensed under MIT License
https://github.com/f0uriest/keras2c

Converted to C++ by Emad Alharbi
University of York
 */


#include <math.h>
#include <stdio.h>
#include "k2c_include.h"
#include "k2c_activations.h"




/**
 * Tanh activation function.
 *   y = tanh(x)
 *
 * :param x: array of input values. Gets overwritten by output.
 * :param size: length of input array.
 */
void k2c_activations::k2c_tanh_func(float * x, const size_t size) {

    for (size_t i=0; i<size; ++i) {
        x[i] = tanh(x[i]);
    }
}



/**
 * Sigmoid activation function.
 *   y = 1/(1+exp(-x))
 *
 * :param x: array of input values. Gets overwritten by output.
 * :param size: length of input array.
 */
void k2c_activations::k2c_sigmoid_func(float * x, const size_t size) {

    for (size_t i=0; i < size; ++i) {
        x[i] = 1/(1+exp(-x[i]));
    }
}




