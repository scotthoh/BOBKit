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
#include <string.h>
#include "k2c_include.h"
#include "k2c_core_layers.h"
#include "k2c_helper_functions.h"
/**
 * Dense (fully connected) Layer.
 *
 * :param output: output tensor.
 * :param input: input tensor.
 * :param kernel: kernel tensor.
 * :param bias: bias tensor.
 * :param activation: activation function to apply to output.
 * :param fwork: array of working space, size(fwork) = size(input) + size(kernel)
 */
void k2c_core_layers::k2c_dense(k2c_tensor* output, const k2c_tensor* input, const k2c_tensor* kernel,
               const k2c_tensor* bias, k2c_activationType *activation, float * fwork) {

    if (input->ndim <=2) {
        size_t outrows;

        if (input->ndim>1) {
            outrows = input->shape[0];
        }
        else {
            outrows = 1;
        }
        const size_t outcols = kernel->shape[1];
        const size_t innerdim = kernel->shape[0];
        const size_t outsize = outrows*outcols;
        k2c_helper_functions::k2c_affine_matmul(output->array,input->array,kernel->array,bias->array,
                          outrows,outcols,innerdim);
        activation(output->array,outsize);
    }
    else {
        const size_t axesA[1] = {input->ndim-1};
        const size_t axesB[1] = {0};
        const size_t naxes = 1;
        const int normalize = 0;

        k2c_helper_functions::k2c_dot(output, input, kernel, axesA, axesB, naxes, normalize, fwork);
        k2c_helper_functions::k2c_bias_add(output, bias);
        activation(output->array, output->numel);
    }
}



