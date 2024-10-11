// Wrapper for buccaneer-join
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-join.h"
#include <clipper/core/coords.h>
#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;
// PYBIND11_MAKE_OPAQUE(std::vector<int>, "IntList");
//  29May problem with array and vector not able to update values individually
//  only able to assign through the function with entire array/vector.
void init_ca_join(py::module &m) {
  py::class_<Ca_join> pyCaJoin(m, "Ca_join");

  pyCaJoin
      .def(py::init<double &, double &>(), py::arg("rad_merge") = 2.0,
           py::arg("rad_join") = 2.0)
      .def("__call__", &Ca_join::operator(), py::arg("mol"),
           "Build chains by merging and joining tri-residue fragments.")
      .def_static("join", &Ca_join::join, py::arg("mol"), py::arg("rmerge"),
                  py::arg("rjoin"), py::arg("com"))
      .def_static("join",
                  [](MiniMol &mol, const double &rmerg, const double &rjoin,
                     const std::array<ftype, 3> &com) -> bool {
                    return Ca_join::join(mol, rmerg, rjoin,
                                         Coord_orth(com[0], com[1], com[2]));
                  })
      .def("__repr__",
           [](const Ca_join &self) { return "<buccaneer.Ca_join class.>"; })
      .doc() =
      "Class for merging overlapped Ca chains and grouping by symmetry.";

  using NodeClass = Ca_join::Node;
  py::class_<NodeClass>(pyCaJoin, "Node")
      .def(py::init<>())
      //.def(py::init([](const float &score, const py::array_t<int> &pointers) {
      //       NodeClass *nptr = new Ca_join::Node();
      //       auto buf = pointers.request();
      //       if (buf.ndim != 1)
      //         throw std::runtime_error("Pointers must be a 1D array!");
      //       int *iptr = (int *)buf.ptr;
      //       int n = buf.shape[0];
      //       for (int i = 0; i < n; ++i) {
      //         nptr->ptrs.push_back(*iptr);
      //         iptr += 1;
      //       }
      //
      //       nptr->score = score;
      //       return std::unique_ptr<NodeClass>(nptr);
      //     }),
      //     py::arg("score"), py::arg("pointers"),
      //     "Constructor from score and array of pointers.")
      .def_readwrite("score", &NodeClass::score)
      //.def_readwrite("ptrs", &NodeClass::ptrs);
      .def(
          "get_ptrs",
          [](NodeClass &self) -> py::array_t<int> {
            return to_numpy_1d<std::vector<int>, int>(self.ptrs,
                                                      self.ptrs.size());
          },
          "Get the Node.ptrs as numpy arrays.")
      .def(
          "set_ptrs",
          [](NodeClass &self, const py::array_t<int> &a) {
            auto buf = a.request();
            self.ptrs.clear();
            if (buf.ndim != 1)
              throw std::runtime_error("Pointers must be a 1D array!");
            int *iptr = (int *)buf.ptr;
            int n = buf.shape[0];

            for (int i = 0; i < n; ++i) {
              self.ptrs.push_back(*iptr);
              iptr += 1;
            }
          },
          py::arg("array"),
          "Set the Node.ptrs from an array, Overrides "
          "existing data.")
      .def("__repr__",
           [](const NodeClass &self) {
             std::stringstream stream;
             stream << "<buccaneer.Ca_join.Node: score = ";
             stream << clipper::String(self.score, 10, 4) << " >";
             return stream.str();
           })
      .doc() = "Node class to store pointers and scores.";

  // py::bind_vector<std::vector<int>>(m, "VectorInt");

  using Res3 = Ca_join::Tri_residue;
  py::class_<Res3>(pyCaJoin, "Tri_residue")
      .def(py::init<>())
      .def_readwrite("flag", &Res3::flag)
      .def_readwrite("score", &Res3::score)
      .def("__repr__",
           [](const Res3 &self) {
             std::stringstream stream;
             stream << "<buccaneer.Ca_join.Tri_residue: ";
             stream << self.type[0] << "," << self.type[1] << ","
                    << self.type[2];
             stream << " with score = " << clipper::String(self.score, 6, 2)
                    << " >";
             return stream.str();
           })
      .def_property(
          "type",
          [](Res3 &self)
              -> std::vector<std::string> { // needs to set the buffer numpy
                                            // arrays and convert the strings
                                            // first
            std::vector<std::string> ret;
            // auto buf = ret.request();
            int n = 3;
            // auto *iptr = (std::string *)buf.ptr;
            for (int i = 0; i < n; ++i)
              ret.push_back(self.type[i]);
            return ret;
          },
          [](Res3 &self, const std::vector<std::string> &types) {
            for (int i = 0; i < 3; ++i)
              self.type[i] = types[i];
          })
      .def(
          "set_type",
          [](Res3 &self, const int pos, const std::string &restype) {
            self.type[pos] = clipper::String(restype);
          },
          "Set residue type with 3-letter code at given position.")
      .def(
          "set_res",
          [](Res3 &self, const int &res, const int &atm,
             const std::array<ftype, 3> &coord) {
            self.res[res][atm] =
                clipper::Coord_orth(coord[0], coord[1], coord[2]);
          },
          py::arg("res_index"), py::arg("atom_index"), py::arg("coords"),
          "Set coordinates for atom of given indices.")
      .def_property(
          "res",
          [](const Res3 &self) -> py::array_t<ftype> {
            // auto ibuf = indices.request();
            // if (ibuf.ndim != 1)
            //   throw std::runtime_error("Indices must be a 1D array!");
            int n = 3; // ibuf.shape[0];
            // int *iptr = (int *)ibuf.ptr;
            auto ret = py::array_t<ftype>({n, n, 3});
            ftype *rptr = (ftype *)ret.request().ptr;
            for (int i = 0; i < n; ++i) {
              for (int j = 0; j < n; ++j) {
                const auto &coord = self.res[i][j];
                for (int k = 0; k < 3; ++k)
                  *rptr++ = coord[k];
              }
            }
            return ret;
          },
          [](Res3 &self, const py::array_t<ftype>
                             coords) { // const int index, const int atom_index,
            // something wrong here not setting properly 26/june
            auto buf = coords.request();
            if (buf.ndim != 3 || buf.shape[1] != 3 || buf.shape[2] != 3)
              throw std::runtime_error("Array must be a 1D of length 3, "
                                       "and coords must be a 3 x 3 array.");

            ftype *ptr = (ftype *)buf.ptr;
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                // clipper::Coord_orth a(*ptr, *(ptr + 1), *(ptr + 2));
                self.res[i][j] =
                    clipper::Coord_orth(*ptr, *(ptr + 1), *(ptr + 2));
                // std::cout << a.format() << std::endl;
                ptr += 3;
              }
            }
          });
}