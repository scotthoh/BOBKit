// Wrapper for buccaneer-prot protein_loop
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-prot.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

void declare_proteinloop(py::module &m) {
  py::class_<ProteinLoop>(m, "ProteinLoop")
      .def(py::init<int>(), py::arg("torsion_sampling") = 24, "Constructor.")
      .def("Coord_O", &ProteinLoop::Coord_O, py::arg("ca0"), py::arg("c0"),
           py::arg("n1"), "Return O from Ca, C, N")
      .def("Coord_Cb", &ProteinLoop::Coord_Cb, py::arg("n0"), py::arg("ca0"),
           py::arg("c0"), "Return C-beta from N, Ca, C")
      .def("rebuild5atoms", &ProteinLoop::rebuild5atoms, py::arg("c0"),
           py::arg("n1"), py::arg("ca1"), py::arg("ca3"), py::arg("c3"),
           py::arg("n4"), "Re-build 5 torsions worth of atoms.")
      .def(
          "rebuild5atoms",
          [](const ProteinLoop &self, const std::array<float, 3> &c0,
             const std::array<float, 3> &n1, const std::array<float, 3> &ca1,
             const std::array<float, 3> &ca3, const std::array<float, 3> &c3,
             const std::array<float, 3> &n4) {
            Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]),
                CA1(ca1[0], ca1[1], ca1[2]);
            Coord_orth CA3(ca3[0], ca3[1], ca3[2]), C3(c3[0], c3[1], c3[2]),
                N4(n4[0], n4[1], n4[2]);
            return self.rebuild5atoms(C0, N1, CA1, CA3, C3, N4);
          },
          py::arg("c0"), py::arg("n1"), py::arg("ca1"), py::arg("ca3"),
          py::arg("c3"), py::arg("n4"), "Re-build 5 torsions worth of atoms.")
      .def("rebuild8atoms", &ProteinLoop::rebuild8atoms, py::arg("c0"),
           py::arg("n1"), py::arg("ca1"), py::arg("ca4"), py::arg("c4"),
           py::arg("n5"), "Re-build 8 torsions worth of atoms.")
      .def(
          "rebuild8atoms",
          [](const ProteinLoop self, const std::array<float, 3> &c0,
             const std::array<float, 3> &n1, const std::array<float, 3> &ca1,
             const std::array<float, 3> &ca4, const std::array<float, 3> &c4,
             const std::array<float, 3> &n5) {
            Coord_orth C0(c0[0], c0[1], c0[2]), N1(n1[0], n1[1], n1[2]),
                CA1(ca1[0], ca1[1], ca1[2]);
            Coord_orth CA4(ca4[0], ca4[1], ca4[2]), C4(c4[0], c4[1], c4[2]),
                N5(n5[0], n5[1], n5[2]);
            return self.rebuild8atoms(C0, N1, CA1, CA4, C4, N5);
          },
          py::arg("c0"), py::arg("n1"), py::arg("ca1"), py::arg("ca4"),
          py::arg("c4"), py::arg("n5"), "Re-build 8 torsions worth of atoms.")
      .def("__repr__",
           [](const ProteinLoop &self) {
             return "<buccaneer.ProteinLoop builder class.>";
           })
      .doc() = "Protein loop builder class.\n"
               "Contains methods for rebuilding loops of various lengths, and "
               "for rebuilding a whole protein.";
  //.def("CoordList_to_numpy");
}

template <int N>
void declare_coordlist(py::module &m, const std::string &name) {
  using Class = ProteinLoop::CoordList<N>;
  std::string PyClass = std::string("CoordList_") + name;

  py::class_<Class> coordlist(m, PyClass.c_str(), py::buffer_protocol());
  coordlist.def(py::init<>())
      .def_buffer([](Class &self) -> py::buffer_info {
        return py::buffer_info(&self[0][0], {N, 3},
                               {sizeof(ftype) * N, sizeof(ftype)});
      })
      .def(
          "__getitem__", [](Class &self, const int i) { return self[i]; },
          "Get atom coordinates at given index.")
      .def(
          "__setitem__",
          [](Class &self, const int i, const Coord_orth &atom) {
            self[i] = atom;
          },
          "Set atom coordinates at given index.")
      .def("__iter__",
           [](Class &self) { return py::make_iterator(&self[0], &self[N]); })
      .def("__len__", [](const Class &self) { return N; })
      .def(
          "size", [](const Class &self) { return N; }, "Return size of list.")
      .def("__repr__",
           [](const Class &self) {
             return "<buccaneer.CoordList_" + String(N) + " class.>";
           })
      .doc() = "Class holding coordinates for a list of " + name + " atoms.";
}

void declare_proteintools(py::module &m) {
  py::class_<ProteinTools>(m, "ProteinTools")
      .def(py::init<>())
      .def_static("residue_index",
                  static_cast<int (*)(char)>(&ProteinTools::residue_index),
                  py::arg("c"), "Get index from 1-letter residue code.")
      .def_static(
          "residue_index",
          static_cast<int (*)(String, bool)>(&ProteinTools::residue_index),
          py::arg("code"), py::arg("translate") = true,
          "Get index from 3-letter residue code.")
      //.def_static(
      //    "residue_index", [](std::string code, bool translate)
      //    { return ProteinTools::residue_index(code, translate); },
      //    py::arg("code"), py::arg("translate") = true)
      .def_static("residue_index_translate",
                  &ProteinTools::residue_index_translate, py::arg("c"),
                  "Get index from residue code.")
      .def_static("residue_index_1", &ProteinTools::residue_index_1,
                  py::arg("code"), py::arg("translate") = true,
                  "Get index from 1-letter residue code.")
      .def_static("residue_index_3", &ProteinTools::residue_index_3,
                  py::arg("code"), py::arg("translate") = true,
                  "Get index from 3-letter residue code.")
      .def_static("residue_code_1", &ProteinTools::residue_code_1,
                  py::arg("index"), "Get 1-letter residue code from index.")
      .def_static("residue_code_3", &ProteinTools::residue_code_3,
                  py::arg("index"), "Get 3-letter residue code from index.")
      .def_static("residue_code", &ProteinTools::residue_code, py::arg("code"),
                  py::arg("translate") = true,
                  "Get residue code from given code.")
      .def_static("residue_codes", &ProteinTools ::residue_codes,
                  "Return residue codes.")
      .def_static("chain_sequence", &ProteinTools::chain_sequence,
                  py::arg("chn"), "Return sequence of given chain.")
      .def_static("chain_sequence_match", &ProteinTools::chain_sequence_match,
                  py::arg("chnseq"), py::arg("seq"), "Match chain sequence.")
      .def_static("superpose", &ProteinTools::superpose, py::arg("c1"),
                  py::arg("c2"), py::arg("rmsd"), py::arg("nmatch"),
                  py::arg("nmismatch"), "Superpose chains.")
      .def_static("chain_number", &ProteinTools::chain_number, py::arg("mol"),
                  "Number sequence in chain.")
      .def_static("chain_label", &ProteinTools::chain_label, py::arg("mol"),
                  py::arg("chainid_2char") = false, "Label chain ids.")
      .def_static(
          "get_usedlabels", &ProteinTools::get_usedlabels, py::arg("chainid"),
          py::arg("labels"),
          "Get a pair of indices of a chain id from a vector of labels.")
      .def_static("copy_residue_types", &ProteinTools::copy_residue_types,
                  py::arg("target"), py::arg("source"),
                  "Copy residue types from source.")
      .def_static("globularise",
                  static_cast<bool (*)(MiniMol &, const Coord_frac)>(
                      &ProteinTools::globularise),
                  py::arg("mol"), py::arg("com"),
                  "Globularise model with a given fractional coordinates "
                  "centre of mass.")
      .def_static("globularise",
                  static_cast<bool (*)(MiniMol &)>(&ProteinTools::globularise),
                  py::arg("mol"), "Find centre of mass and globularise model ")
      .def_static("symm_match", &ProteinTools::symm_match, py::arg("molwrk"),
                  py::arg("molref"),
                  "Perform symmetry match of a model to a "
                  "reference model.")
      .def_static("main_chain_densities", &ProteinTools::main_chain_densities,
                  py::arg("mp"), py::arg("xmap"), py::arg("nsmooth") = 0,
                  "Calculate density of all N, CA, and C atoms.")
      .def_static("main_chain_u_values", &ProteinTools::main_chain_u_values,
                  py::arg("mp"), py::arg("nsmooth") = 0,
                  "Calculate mean isotropic U values of each residue from N, "
                  "CA, and C atoms.")
      .def_static("main_chain_u_mean", &ProteinTools::main_chain_u_mean,
                  py::arg("mol"),
                  "Calculate mean isotropic U value of all N, CA, and C atoms. "
                  "If none exist, return default value of 0.5.")
      .def_static("split_chains_at_gap", &ProteinTools::split_chains_at_gap,
                  py::arg("mol"), "Separate chains at gaps.")
      .def_static("split_chains_at_unk", &ProteinTools::split_chains_at_unk,
                  py::arg("mol"), py::arg("xmap"), "Separate chains at UNK.")
      .def_static("tidy_peptide_bond", &ProteinTools::tidy_peptide_bond,
                  py::arg("mm1"), py::arg("mm2"),
                  "Tidy and rebuild peptide unit.")
      .def_static("ca_chains", &ProteinTools::ca_chains, py::arg("mol"),
                  "Return a list of Ca_chain.")
      .def_static("insert_ca_chains", &ProteinTools::insert_ca_chains,
                  py::arg("mol"), py::arg("chains"),
                  "Insert Ca_chains to model.")
      .def_static("trim_to_protein", &ProteinTools::trim_to_protein,
                  py::arg("mol"), "Trim to protein backbone.")
      .def_static("is_protein", &ProteinTools::is_protein, py::arg("res"),
                  "Return true if residue in amino acid.")
      .def("__repr__",
           [](const ProteinTools &self) {
             return "<buccaneer.ProteinTools class.>";
           })
      .doc() = "Useful tools for manipulating protein.";
}

void init_protein_loop(py::module &m) {
  declare_proteinloop(m);
  declare_coordlist<5>(m, "5");
  declare_coordlist<8>(m, "8");
  declare_proteintools(m);
}