#ifndef BUILDKIT_MINIMOL_H_
#define BUILDKIT_MINIMOL_H_

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>

#include <thread>
#include <string>
// #include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>

#include "buccaneer-util.h"
// #include "toString.hpp"

// namespace buildkit
//{
//   class PyCMiniMol
//   {
//   public:
//     PyCMiniMol(){};
//     PyCMiniMol(std::string &filepath_to_structure, bool enable_messages);
//     void read_file(std::string &filepath_to_structure);
//     clipper::MiniMol *get_mmol(std::string &filepath_to_structure);
//     // void import_minimol(std::string &filepath_to_structure);
//     // void export_minimol();
//
//   private:
//     std::string input_model_file_;
//     clipper::MiniMol model_in_;
//     // clipper::MiniMol model_out_;
//     bool enable_messages_;
//   };
// }
#endif