// Nanobind bindings for buccaneer-prot proteinloop proteintools
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-prot.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>
#include <nanobind/make_iterator.h>

using namespace clipper;

void remove_sidechains( MiniMol &mol, const bool &keep_oxygen=true ) {
  MiniMol mold = mol;
  mol = MiniMol( mold.spacegroup(), mold.cell() );
  for ( int chn = 0; chn < mold.size(); chn++ ) {
    MPolymer mp;
    for ( int res = 0; res < mold[chn].size(); res++ ) {
      //MMonomer mm;
      clipper::String searchstr = "";
      int in = mold[chn][res].lookup( " N  ", clipper::MM::ANY );
      int ia = mold[chn][res].lookup( " CA ", clipper::MM::ANY );
      int ic = mold[chn][res].lookup( " C  ", clipper::MM::ANY );
      if ( in >= 0 && ia >= 0 && ic >= 0 ) {
        searchstr = "N,CA,C";
      //  mm.insert( mold[chn][res][in] );
      //  mm.insert( mold[chn][res][ia] );
      //  mm.insert( mold[chn][res][ic] );
      }
      if ( keep_oxygen ) {
        int io = mold[chn][res].lookup( " O  ", clipper::MM::ANY );
        if ( io >= 0 ) searchstr.append(",O"); //mm.insert( mold[chn][res][io] );
      }
      auto mm = mold[chn][res].select(searchstr);
      
      if ( mm.size() > 0 ) mp.insert( mm );
    }  
    if ( mp.size() > 0 ) mol.insert( mp );
  }
}

void declare_proteinloop(nb::module_ &m) {
  nb::class_<ProteinLoop>(m, "ProteinLoop")
      .def(nb::init<int>(), nb::arg("torsion_sampling") = 24, "Constructor.")
      .def("Coord_O", &ProteinLoop::Coord_O, nb::arg("ca0"), nb::arg("c0"),
           nb::arg("n1"), "Return O from Ca, C, N")
      .def("Coord_Cb", &ProteinLoop::Coord_Cb, nb::arg("n0"), nb::arg("ca0"),
           nb::arg("c0"), "Return C-beta from N, Ca, C")
      .def("rebuild5atoms", &ProteinLoop::rebuild5atoms, nb::arg("c0"),
           nb::arg("n1"), nb::arg("ca1"), nb::arg("ca3"), nb::arg("c3"),
           nb::arg("n4"), "Re-build 5 torsions worth of atoms.")
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
          nb::arg("c0"), nb::arg("n1"), nb::arg("ca1"), nb::arg("ca3"),
          nb::arg("c3"), nb::arg("n4"), "Re-build 5 torsions worth of atoms.")
      .def("rebuild8atoms", &ProteinLoop::rebuild8atoms, nb::arg("c0"),
           nb::arg("n1"), nb::arg("ca1"), nb::arg("ca4"), nb::arg("c4"),
           nb::arg("n5"), "Re-build 8 torsions worth of atoms.")
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
          nb::arg("c0"), nb::arg("n1"), nb::arg("ca1"), nb::arg("ca4"),
          nb::arg("c4"), nb::arg("n5"), "Re-build 8 torsions worth of atoms.")
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
void declare_coordlist(nb::module_ &m, const std::string &name) {
  using Class = ProteinLoop::CoordList<N>;
  std::string PyClass = std::string("CoordList_") + name;

  nb::class_<Class> coordlist(m, PyClass.c_str());
  coordlist.def(nb::init<>())
      .def_prop_rw("array", [](Class &self) {
        auto np_array = make_numpy_array<ftype>({N,3});
        ftype* arr = np_array.data();
        for (size_t i = 0; i < N; ++i) {
          const auto &coord = self[i];
          for (size_t j = 0; j < 3; ++j)
            *arr++ = (ftype) coord[j];
        }
        return np_array;
      },
        [](Class &self, const cpu_c_2darray<ftype, N, 3> &coords) { 
          // const int index, const int atom_index,
          // something wrong here not setting properly 26/june
          auto c = coords.view();
          if (c.ndim() != 2 || c.shape(0) != N || c.shape(1) !=3 )
            throw std::runtime_error("Array shape must be of {N, 3}!");
          for (size_t i = 0; i < N; ++i) {
                Coord_orth co((ftype)c(i,0), (ftype)c(i,1), (ftype)c(i,2));
                self[i] = co;
            }
      })
      //.def_buffer([](Class &self) -> nb::buffer_info {
      //  return nb::buffer_info(&self[0][0], {N, 3},
      //                         {sizeof(ftype) * N, sizeof(ftype)});
      //})
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
           [](Class &self) { return nb::make_iterator<nb::rv_policy::reference_internal>( nb::type<Class>(), "iterator", &self[0], &self[N]); }
          , nb::keep_alive<0, 1>())
      .def("__len__", [](const Class &self) { return N; })
      .def(
          "size", [](const Class &self) { return N; }, "Return size of list.")
      .def("__repr__",
           [](const Class &self) {
             return "<buccaneer.CoordList_" + String(N) + " class.>";
           })
      .doc() = "Class holding coordinates for a list of " + name + " atoms.";
}

void declare_proteintools(nb::module_ &m) {
  nb::class_<ProteinTools>(m, "ProteinTools")
      .def(nb::init<>())
      .def_static("residue_index",
                  static_cast<int (*)(char)>(&ProteinTools::residue_index),
                  nb::arg("c"), "Get index from 1-letter residue code.")
      .def_static(
          "residue_index",
          static_cast<int (*)(String, bool)>(&ProteinTools::residue_index),
          nb::arg("code"), nb::arg("translate") = true,
          "Get index from 3-letter residue code.")
      //.def_static(
      //    "residue_index", [](std::string code, bool translate)
      //    { return ProteinTools::residue_index(code, translate); },
      //    nb::arg("code"), nb::arg("translate") = true)
      .def_static("residue_index_translate",
                  &ProteinTools::residue_index_translate, nb::arg("c"),
                  "Get index from residue code.")
      .def_static("residue_index_1", &ProteinTools::residue_index_1,
                  nb::arg("code"), nb::arg("translate") = true,
                  "Get index from 1-letter residue code.")
      .def_static("residue_index_3", &ProteinTools::residue_index_3,
                  nb::arg("code"), nb::arg("translate") = true,
                  "Get index from 3-letter residue code.")
      .def_static("residue_code_1", &ProteinTools::residue_code_1,
                  nb::arg("index"), "Get 1-letter residue code from index.")
      .def_static("residue_code_3", &ProteinTools::residue_code_3,
                  nb::arg("index"), "Get 3-letter residue code from index.")
      .def_static("residue_code", &ProteinTools::residue_code, nb::arg("code"),
                  nb::arg("translate") = true,
                  "Get residue code from given code.")
      .def_static("residue_codes", &ProteinTools ::residue_codes,
                  "Return residue codes.")
      .def_static("chain_sequence", &ProteinTools::chain_sequence,
                  nb::arg("chn"), "Return sequence of given chain.")
      .def_static("chain_sequence_match", &ProteinTools::chain_sequence_match,
                  nb::arg("chnseq"), nb::arg("seq"), "Match chain sequence.")
      .def_static("superpose", &ProteinTools::superpose, nb::arg("c1"),
                  nb::arg("c2"), nb::arg("rmsd"), nb::arg("nmatch"),
                  nb::arg("nmismatch"), "Superpose chains.")
      .def_static("chain_number", &ProteinTools::chain_number, nb::arg("mol"),
                  "Number sequence in chain.")
      .def_static("chain_label", &ProteinTools::chain_label, nb::arg("mol"),
                  nb::arg("chainid_2char") = false, "Label chain ids.")
      .def_static(
          "get_usedlabels", &ProteinTools::get_usedlabels, nb::arg("chainid"),
          nb::arg("labels"),
          "Get a pair of indices of a chain id from a vector of labels.")
      .def_static("copy_residue_types", &ProteinTools::copy_residue_types,
                  nb::arg("target"), nb::arg("source"),
                  "Copy residue types from source.")
      .def_static("globularise",
                  static_cast<bool (*)(MiniMol &, const Coord_frac)>(
                      &ProteinTools::globularise),
                  nb::arg("mol"), nb::arg("com"),
                  "Globularise model with a given fractional coordinates "
                  "centre of mass.")
      .def_static("globularise",
                  static_cast<bool (*)(MiniMol &)>(&ProteinTools::globularise),
                  nb::arg("mol"), "Find centre of mass and globularise model ")
      .def_static("symm_match", &ProteinTools::symm_match, nb::arg("molwrk"),
                  nb::arg("molref"),
                  "Perform symmetry match of a model to a "
                  "reference model.")
      .def_static("main_chain_densities", &ProteinTools::main_chain_densities,
                  nb::arg("mp"), nb::arg("xmap"), nb::arg("nsmooth") = 0,
                  "Calculate density of all N, CA, and C atoms.")
      .def_static("main_chain_u_values", &ProteinTools::main_chain_u_values,
                  nb::arg("mp"), nb::arg("nsmooth") = 0,
                  "Calculate mean isotropic U values of each residue from N, "
                  "CA, and C atoms.")
      .def_static("main_chain_u_mean", &ProteinTools::main_chain_u_mean,
                  nb::arg("mol"),
                  "Calculate mean isotropic U value of all N, CA, and C atoms. "
                  "If none exist, return default value of 0.5.")
      .def_static("split_chains_at_gap", &ProteinTools::split_chains_at_gap,
                  nb::arg("mol"), "Separate chains at gaps.")
      .def_static("split_chains_at_unk", &ProteinTools::split_chains_at_unk,
                  nb::arg("mol"), nb::arg("xmap"), "Separate chains at UNK.")
      .def_static("tidy_peptide_bond", &ProteinTools::tidy_peptide_bond,
                  nb::arg("mm1"), nb::arg("mm2"),
                  "Tidy and rebuild peptide unit.")
      .def_static("ca_chains", &ProteinTools::ca_chains, nb::arg("mol"),
                  "Return a list of Ca_chain.")
      .def_static("insert_ca_chains", &ProteinTools::insert_ca_chains,
                  nb::arg("mol"), nb::arg("chains"),
                  "Insert Ca_chains to model.")
      .def_static("trim_to_protein", &ProteinTools::trim_to_protein,
                  nb::arg("mol"), "Trim to protein only.")
      .def_static("remove_sidechain", [](MiniMol &mol, const bool &keepoxy) { remove_sidechains(mol, keepoxy); },
                  nb::arg("mol"), nb::arg("keep_oxygen") = true,
                  "Remove sidechains, leave only atoms with labels "
                  "'N', 'Ca', 'C', or 'O' if keep_oxygen is True. Default: keep_oxygen = True")
      .def_static("is_protein", &ProteinTools::is_protein, nb::arg("res"),
                  "Return true if residue in amino acid.")
      .def("__repr__",
           [](const ProteinTools &self) {
             return "<buccaneer.ProteinTools class.>";
           })
      .doc() = "Useful tools for manipulating protein.";
}

void add_protein_loop(nb::module_ &m) {
  declare_proteinloop(m);
  declare_coordlist<5>(m, "5");
  declare_coordlist<8>(m, "8");
  declare_proteintools(m);
}