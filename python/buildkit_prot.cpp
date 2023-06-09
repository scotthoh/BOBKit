// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer-prot.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include "helper_functions.h"

PYBIND11_MAKE_OPAQUE(Ca_chain);

void declare_ca_group(py::module &m)
{
  py::class_<Ca_group>(m, "Ca_group")
      .def(py::init<>())
      .def(py::init<const clipper::Coord_orth &, const clipper::Coord_orth &, const clipper::Coord_orth &>())
      .def(py::init<const clipper::MResidue &>())
      .def_property_readonly("is_null", &Ca_group::is_null)
      .def_property_readonly("coord_n", &Ca_group::coord_n)
      .def_property_readonly("coord_ca", &Ca_group::coord_ca)
      .def_property_readonly("coord_c", &Ca_group::coord_c)
      .def_property_readonly("coord_cb", &Ca_group::coord_cb)
      .def("rtop_from_std_ori", &Ca_group::rtop_from_std_ori)
      .def("rtop_beta_carbon", &Ca_group::rtop_beta_carbon)
      .def("next_ca_group", &Ca_group::next_ca_group, py::arg("psi"), py::arg("phi"))
      .def("prev_ca_group", &Ca_group::prev_ca_group, py::arg("phi"), py::arg("psi"))
      .def_static("std_coord_ca", &Ca_group::std_coord_ca)
      .def_static("std_coord_c", &Ca_group::std_coord_c)
      .def_static("std_coord_n", &Ca_group::std_coord_n)
      .def_static("std_coord_cb", &Ca_group::std_coord_cb)
      .def_static("null", &Ca_group::null);
}

void declare_ca_chain(py::module &m)
{
  // py::bind_vector < std::vector < Ca_group>>(m, "DequeCagroup");
  py::class_<Ca_chain>(m, "Ca_chain")
      .def(py::init<>())
      .def("ramachandran_phi", &Ca_chain::ramachandran_phi)
      .def("ramachandran_psi", &Ca_chain::ramachandran_psi)
      .def("__len__", [](const Ca_chain &self)
           { return self.size(); });
}

void declare_pr_group(py::module &m)
{
  py::class_<Pr_group> pr_group(m, "Pr_group");

  py::enum_<Pr_group::TYPE>(pr_group, "TYPE")
      .value("CaCN", Pr_group::TYPE::CaCN)
      .value("CaCO", Pr_group::TYPE::CaCO);

  pr_group
      .def(py::init<>())
      .def(py::init<const clipper::Coord_orth &, const clipper::Coord_orth &,
                    const clipper::Coord_orth &, const Pr_group::TYPE &>())
      .def_property_readonly("coord_ca", &Pr_group::coord_ca)
      .def_property_readonly("coord_c", &Pr_group::coord_c)
      .def_property_readonly("coord_n_next", &Pr_group::coord_n_next)
      .def("coord_o", &Pr_group::coord_o)
      .def("coord_ca_next", &Pr_group::coord_ca_next)
      .def("rtop_from_std_ori", &Pr_group::rtop_from_std_ori)
      .def("next_pr_group", &Pr_group::next_pr_group)
      .def("prev_pr_group", &Pr_group::prev_pr_group);
}

void declare_proteinloop(py::module &m)
{
  py::class_<ProteinLoop>(m, "ProteinLoop")
      .def(py::init<int>(), py::arg("torsion_sampling") = 24)
      .def("Coord_O", &ProteinLoop::Coord_O)
      .def("Coord_Cb", &ProteinLoop::Coord_Cb)
      .def("rebuild5atoms", &ProteinLoop::rebuild5atoms)
      .def("rebuild5atoms", [](const ProteinLoop &self, std::array<float, 3> &c0,
                               std::array<float, 3> &n1, std::array<float, 3> &ca1,
                               std::array<float, 3> &ca3, std::array<float, 3> &c3,
                               std::array<float, 3> &n4)
           {
        clipper::Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]), CA1(ca1[0], ca1[1], ca1[2]);
        clipper::Coord_orth CA3(ca3[0], ca3[1], ca3[2]), C3(c3[0], c3[1], c3[2]), N4(n4[0], n4[1], n4[2]);
        return self.rebuild5atoms(C0, N1, CA1, CA3, C3, N4); })
      .def("rebuild8atoms", &ProteinLoop::rebuild8atoms)
      .def("rebuild8atoms", [](const ProteinLoop self, std::array<float, 3> &c0,
                               std::array<float, 3> &n1, std::array<float, 3> &ca1,
                               std::array<float, 3> &ca4, std::array<float, 3> &c4,
                               std::array<float, 3> &n5)
           {
        clipper::Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]), CA1(ca1[0], ca1[1], ca1[2]);
        clipper::Coord_orth CA4(ca4[0], ca4[1], ca4[2]), C4(c4[0], c4[1], c4[2]), N5(n5[0], n5[1], n5[2]);
        returb self.rebuild8atoms(C0, N1, CA1, CA4, C4, N5); });
  //.def("CoordList_to_numpy");
}

template <int N>
void declare_coordlist(py::module &m, const std::string &name)
{
  using Class = ProteinLoop::CoordList<N>;
  std::string PyClass = std::string("CoordList_") + name;
  py::class_<Class> coordlist(m, PyClass.c_str());
  coordlist
      .def(py::init<>())
      .def("__getitem__", [](Class &self, const int i)
           { return self[i]; })
      .def("__iter__", [](Class &self)
           { return py::make_iterator(&self[0], &self[N]); })
      .def("__len__", [](Class &self)
           { return N; });
}

void init_buccaneer_prot(py::module &m)
{
  declare_ca_group(m);
  declare_ca_chain(m);
  declare_pr_group(m);
  declare_proteinloop(m);
  declare_coordlist<5>(m, "5");
  declare_coordlist<8>(m, "8");
}