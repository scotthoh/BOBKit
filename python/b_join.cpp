// Nanobind bindings for buccaneer-join
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-join.h"
#include "commons.h"
#include "arrays.h"
#include <clipper/clipper.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/array.h>


using namespace clipper;

// PYBIND11_MAKE_OPAQUE(std::vector<int>, "IntList");
//  29May problem with array and vector not able to update values individually
//  only able to assign through the function with entire array/vector.
void add_ca_join(nb::module_ &m) {
  nb::class_<Ca_join> pyCaJoin(m, "Ca_join");

  pyCaJoin
      .def(nb::init<double &, double &>(), nb::arg("rad_merge") = 2.0,
           nb::arg("rad_join") = 2.0)
      .def("__call__", &Ca_join::operator(), nb::arg("mol"),
           "Build chains by merging and joining tri-residue fragments.")
      .def_static("join", &Ca_join::join, nb::arg("mol"), nb::arg("rmerge"),
                  nb::arg("rjoin"), nb::arg("com"))
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
  nb::class_<NodeClass>(pyCaJoin, "Node")
      .def(nb::init<>())
      //.def(nb::init([](const float &score, const nb::array_t<int> &pointers) {
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
      //     nb::arg("score"), nb::arg("pointers"),
      //     "Constructor from score and array of pointers.")
      .def_rw("score", &NodeClass::score)
      //.def_rw("ptrs", &NodeClass::ptrs);
      // 24th Sept NEeds changing
      .def(
          "get_ptrs",
          [](NodeClass &self) {
            auto v = new (std::vector<int>)(std::move(self.ptrs));
            nb::capsule owner(v, [](void *p) noexcept { delete static_cast<std::vector<int>*>(p); });
            return nb::ndarray<nb::numpy, int, nb::shape<-1>>(v->data(), {v->size()}, owner);
          },
          "Get the Node.ptrs as numpy arrays.")
      .def(
          "set_ptrs",
          [](NodeClass &self, const cpu_c_array<int> &a) {
            auto v = a.view();
            self.ptrs.clear();
            self.ptrs.resize(v.shape(0));
            if (v.ndim() != 1)
              throw std::runtime_error("Pointers must be a 1D array!");
            for (int i = 0; i < v.shape(0); ++i) {
              self.ptrs.push_back(v(i));
            }
          },
          nb::arg("array"),
          "Set the Node.ptrs from an array, Overrides existing data.")
      .def("__repr__",
           [](const NodeClass &self) {
             return "<buccaneer.Ca_join.Node: score = " + String(self.score, 10, 4) + " >"; 
           })
      .doc() = "Node class to store pointers and scores.";

  // nb::bind_vector<std::vector<int>>(m, "VectorInt");

  using Res3 = Ca_join::Tri_residue;
  nb::class_<Res3>(pyCaJoin, "Tri_residue")
      .def(nb::init<>())
      .def_rw("flag", &Res3::flag)
      .def_rw("score", &Res3::score)
      .def("__repr__",
           [](const Res3 &self)
           {
             return "<buccaneer.Ca_join.Tri_residue: " + self.type[0] +  "," + self.type[1] + "," +  self.type[2] + " with score = " + String(self.score, 6, 2) + " >";
           })
      //.def_rw("type", &Res3::type)
      .def_prop_rw(
          "type",
          [](Res3 &self)
          { // needs to set the buffer numpy
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
          [](Res3 &self, const std::vector<std::string> &types)
          {
            for (int i = 0; i < 3; ++i)
              self.type[i] = types[i];
          })
      .def(
          "set_type_at",
          [](Res3 &self, const int pos, const std::string &restype)
          {
            self.type[pos] = String(restype);
          },
          "Set residue type with 3-letter code at given position.")
      .def(
          "set_res_at",
          [](Res3 &self, const int &res, const int &atm,
             const std::array<ftype, 3> &coord)
          {
            self.res[res][atm] = Coord_orth(coord[0], coord[1], coord[2]);
          },
          nb::arg("res_index"), nb::arg("atom_index"), nb::arg("coords"),
          "Set coordinates for atom of given indices.")
      .def_prop_rw(
          "res",
          [](const Res3 &self)
          {
            // Tri-res, 3 atoms, coordinates-xyz
            auto np_array = make_numpy_array<ftype>({3,3,3});
            ftype* arr = np_array.data();
            for (size_t i = 0; i < 3; ++i)
              for (size_t j = 0; j < 3; ++j) { 
                const auto &coord = self.res[i][j];
                for (size_t k = 0; k < 3; ++k)
                  *arr++ = (ftype) coord[k];
              }
            return np_array;
            // auto ibuf = indices.request();
            // if (ibuf.ndim != 1)
            //   throw std::runtime_error("Indices must be a 1D array!");
            //int n = 3; // ibuf.shape[0];
            //// int *iptr = (int *)ibuf.ptr;
            //auto ret = nb::array_t<ftype>({n, n, 3});
            //ftype *rptr = (ftype *)ret.request().ptr;
            //for (int i = 0; i < n; ++i)
            //{
            //  for (int j = 0; j < n; ++j)
            //  {
            //    const auto &coord = self.res[i][j];
            //    for (int k = 0; k < 3; ++k)
            //      *rptr++ = coord[k];
            //  }
            //}
            //return ret;
          },
          [](Res3 &self, const np_cpu_c_3darray<ftype> &coords) { 
            // const int index, const int atom_index,
            // something wrong here not setting properly 26/june
            auto c = coords.view();
            if (c.ndim() != 3 || c.shape(0) !=3 || c.shape(1) != 3 || c.shape(2) != 3)
              throw std::runtime_error("Array shape must be of {3, 3, 3}!");
            for (size_t i = 0; i < 3; ++i)
              for (size_t j = 0; j < 3; ++j) {
                  Coord_orth co((ftype)c(i,j,0), (ftype)c(i,j,1), (ftype)c(i,j,2));
                  self.res[i][j] = co;
              }
          });
}