// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prot.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include "helper_functions.h"
#include "type_conversions.h"

void declare_ca_group(py::module &m)
{
  py::class_<Ca_group>(m, "Ca_group")
      .def(py::init<>())
      .def(py::init<const clipper::Coord_orth &, const clipper::Coord_orth &, const clipper::Coord_orth &>())
      // constructor from list/array of coordinates N, CA, C
      .def(py::init([](const std::array<ftype, 3> &n, const std::array<ftype, 3> &ca, const std::array<ftype, 3> &c)
                    {
        Coord_orth N(n[0], n[1], n[2]);
        Coord_orth CA(ca[0], ca[1], ca[2]);
        Coord_orth C(c[0], c[1], c[2]);
        return std::unique_ptr<Ca_group>(new Ca_group(N, CA, C)); }))
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
  // Have to bind some of the deque member functions.
  py::class_<Ca_chain>(m, "Ca_chain")
      .def(py::init<>())
      .def("ramachandran_phi", &Ca_chain::ramachandran_phi)
      .def("ramachandran_psi", &Ca_chain::ramachandran_psi)
      .def("__len__", [](const Ca_chain &self)
           { return self.size(); })
      .def("size", [](const Ca_chain &self)
           { return self.size(); })
      .def(
          "append", [](Ca_chain &self, Ca_group c)
          { self.push_back(c); },
          "Append item to the end.")
      .def(
          "append_left", [](Ca_chain &self, Ca_group c)
          { self.push_front(c); },
          "Append item to the beginning.")
      .def("__getitem__", [](Ca_chain &self, const int index)
           { return self.at(index); })
      .def("__setitem__", [](Ca_chain &self, Ca_group c, const int index)
           { self.at(index) = c; })
      .def(
          "__iter__", [](Ca_chain &self)
          { return py::make_iterator(self.begin(), self.end()); },
          py::keep_alive<0, 1>())
      .def(
          "__reversed__", [](Ca_chain &self)
          { return py::make_iterator(self.rbegin(), self.rend()); },
          py::keep_alive<0, 1>())
      .def(
          "delete_item", [](Ca_chain &self, const int &i)
          { delitem_at_index(self, i); },
          "Delete item at the given index.")
      .def(
          "insert_item", [](Ca_chain &self, Ca_group c, const int &i)
          { add_item(self, c, i); },
          "Insert item at the given index.")
      .def(
          "delete_range", [](Ca_chain &self, const int &start, const int &end)
          { delitem_range(self, start, end); },
          "Delete items at range of indices given [start, end).");
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
      // constructor from list/array of coordinates
      .def(py::init([](const std::array<ftype, 3> &ca, const std::array<ftype, 3> &c, const std::array<ftype, 3> &thirdatom, const Pr_group::TYPE &type)
                    {
                      Coord_orth CA(ca[0], ca[1], ca[2]);
                      Coord_orth C(c[0], c[1], c[2]);
                      Coord_orth other(thirdatom[0], thirdatom[1], thirdatom[2]);
                      return std::unique_ptr<Pr_group>(new Pr_group(CA, C, other, type)); }))
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
      .def("rebuild5atoms", [](const ProteinLoop &self, const std::array<float, 3> &c0,
                               const std::array<float, 3> &n1, const std::array<float, 3> &ca1,
                               const std::array<float, 3> &ca3, const std::array<float, 3> &c3,
                               const std::array<float, 3> &n4)
           {
        clipper::Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]), CA1(ca1[0], ca1[1], ca1[2]);
        clipper::Coord_orth CA3(ca3[0], ca3[1], ca3[2]), C3(c3[0], c3[1], c3[2]), N4(n4[0], n4[1], n4[2]);
        return self.rebuild5atoms(C0, N1, CA1, CA3, C3, N4); })
      .def("rebuild8atoms", &ProteinLoop::rebuild8atoms)
      .def("rebuild8atoms", [](const ProteinLoop self, const std::array<float, 3> &c0,
                               const std::array<float, 3> &n1, const std::array<float, 3> &ca1,
                               const std::array<float, 3> &ca4, const std::array<float, 3> &c4,
                               const std::array<float, 3> &n5)
           {
        clipper::Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]), CA1(ca1[0], ca1[1], ca1[2]);
        clipper::Coord_orth CA4(ca4[0], ca4[1], ca4[2]), C4(c4[0], c4[1], c4[2]), N5(n5[0], n5[1], n5[2]);
        return self.rebuild8atoms(C0, N1, CA1, CA4, C4, N5); });
  //.def("CoordList_to_numpy");
}

template <int N>
void declare_coordlist(py::module &m, const std::string &name)
{
  using Class = ProteinLoop::CoordList<N>;
  std::string PyClass = std::string("CoordList_") + name;
  py::class_<Class> coordlist(m, PyClass.c_str(), py::buffer_protocol());
  coordlist
      .def(py::init<>())
      .def_buffer([](Class &self) -> py::buffer_info
                  { return py::buffer_info(
                        &self[0][0],
                        {N, 3},
                        {sizeof(ftype) * N, sizeof(ftype)}); })
      .def("__getitem__", [](Class &self, const int i)
           { return self[i]; })
      .def("__iter__", [](Class &self)
           { return py::make_iterator(&self[0], &self[N]); })
      .def("__len__", [](const Class &self)
           { return N; })
      .def("size", [](const Class &self)
           { return N; });
}

void declare_proteintools(py::module &m)
{
  py::class_<ProteinTools>(m, "ProteinTools")
      .def(py::init<>())
      .def_static("residue_index", static_cast<int (*)(char)>(&ProteinTools::residue_index))
      .def_static("residue_index", static_cast<int (*)(String, bool)>(&ProteinTools::residue_index),
                  py::arg("code"), py::arg("translate") = true)
      //.def_static(
      //    "residue_index", [](std::string code, bool translate)
      //    { return ProteinTools::residue_index(code, translate); },
      //    py::arg("code"), py::arg("translate") = true)
      .def_static("residue_index_translate", &ProteinTools::residue_index_translate)
      .def_static("residue_index_1", &ProteinTools::residue_index_1, py::arg("code"), py::arg("translate") = true)
      .def_static("residue_index_3", &ProteinTools::residue_index_3, py::arg("code"), py::arg("translate") = true)
      .def_static("residue_code_1", &ProteinTools::residue_code_1)
      .def_static("residue_code_3", &ProteinTools::residue_code_3)
      .def_static("residue_code", &ProteinTools::residue_code, py::arg("code"), py::arg("translate") = true)
      .def_static("residue_codes", &ProteinTools ::residue_codes)
      .def_static("chain_sequence", &ProteinTools::chain_sequence)
      .def_static("chain_sequence_match", &ProteinTools::chain_sequence_match)
      .def_static("superpose", &ProteinTools::superpose)
      .def_static("chain_number", &ProteinTools::chain_number)
      .def_static("chain_label", &ProteinTools::chain_label) // might need manual since there is clipperMMDBManager::TYPE cifflag
      .def_static("get_usedlabels", &ProteinTools::get_usedlabels)
      .def_static("copy_residue_types", &ProteinTools::copy_residue_types)
      .def_static("globularise", static_cast<bool (*)(clipper::MiniMol &, const clipper::Coord_frac)>(&ProteinTools::globularise))
      .def_static("globularise", static_cast<bool (*)(clipper::MiniMol &)>(&ProteinTools::globularise))
      .def_static("symm_match", &ProteinTools::symm_match)
      .def_static("main_chain_densities", &ProteinTools::main_chain_densities, py::arg("mp"), py::arg("xmap"), py::arg("nsmooth") = 0)
      .def_static("main_chain_u_values", &ProteinTools::main_chain_u_values, py::arg("mp"), py::arg("nsmooth") = 0)
      .def_static("main_chain_u_mean", &ProteinTools::main_chain_u_mean)
      .def_static("split_chains_at_gap", &ProteinTools::split_chains_at_gap)
      .def_static("split_chains_at_unk", &ProteinTools::split_chains_at_unk)
      .def_static("tidy_peptide_bond", &ProteinTools::tidy_peptide_bond)
      .def_static("ca_chains", &ProteinTools::ca_chains)
      .def_static("insert_ca_chains", &ProteinTools::insert_ca_chains)
      .def_static("trim_to_protein", &ProteinTools::trim_to_protein)
      .def_static("is_protein", &ProteinTools::is_protein);
}

void init_buccaneer_prot(py::module &m)
{
  declare_ca_group(m);
  declare_ca_chain(m);
  declare_pr_group(m);
  declare_proteinloop(m);
  declare_coordlist<5>(m, "5");
  declare_coordlist<8>(m, "8");
  declare_proteintools(m);
}