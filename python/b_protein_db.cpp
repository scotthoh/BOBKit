// Nanobind bindings for protein_db
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York


#include "buccaneer/protein_db.h"
#include "commons.h"
#include "arrays.h"
#include <clipper/minimol/minimol.h>

using namespace clipper;

void declare_residue(nb::module_ &m) {
  using PDBRes = ProteinDB::Residue;
  nb::class_<PDBRes> protdb_res(m, "Residue");

  nb::enum_<PDBRes::FLAG>(protdb_res, "FLAG")
      .value("NONE", PDBRes::FLAG::NONE)
      .value("NORMAL", PDBRes::FLAG::NORMAL)
      .value("CALPHA", PDBRes::FLAG::CALPHA)
      .export_values();

  protdb_res.def(nb::init<>())
      .def(nb::init<Coord_orth &, const String &>(), nb::arg("Calpha"),
           nb::arg("type"),
           "Constructor from C-alpha coordinates and residue type.")
      .def(nb::init<Coord_orth &, Coord_orth &, Coord_orth &, const String &>(),
           nb::arg("Cn"), nb::arg("Calpha"), nb::arg("Cc"), nb::arg("type"),
           "Constructor from main chain coordinates and residue type.")
      .def(nb::init<const MMonomer &>(), nb::arg("monomer"),
           "Constructor from MMonomer.")
      .def_prop_ro("coord_n", &PDBRes::coord_n,
                             "Get N atom coordinate.")
      .def_prop_ro("coord_ca", &PDBRes::coord_ca,
                             "Get C-alpha atom coordinate.")
      .def_prop_ro("coord_c", &PDBRes::coord_c,
                             "Get C atom coordinate.")
      .def_prop_ro("mmonomer", &PDBRes::mmonomer, "Build monomer.")
      .def("transform", &PDBRes::transform, nb::arg("rtop"),
           "Transform by rotation-translation operator.")
      .def("merge", &PDBRes::merge, nb::arg("res"), nb::arg("wn"),
           nb::arg("wa"), nb::arg("wc"),
           "Merge other residue coords with this one using given weights.")
      .def_prop_rw("type", &PDBRes::type, &PDBRes::set_type,
                    "Get and set 1-letter residue types.")
      .def_prop_rw("flag", &PDBRes::flag, &PDBRes::set_flag,
                    "Get and set flag.")
      //.def("set_type", &PDBRes::set_type)
      //.def("set_flag", &PDBRes::set_flag)
      .def("data_import", &PDBRes::data_import, nb::arg("d"),
           "Import from char array.")
      .def("data_export", [](const PDBRes &self) {
          char* cstr;
          self.data_export(cstr);
          return std::string(cstr);
      }, "Export to char array")
          
      //    &PDBRes::data_export, nb::arg("d"),
      //     "Export to char array")
      .def("is_null", &PDBRes::is_null, "Test for null.")
      .def_static("residue_type", &PDBRes::residue_type,
                  "Return 1-letter residue type from 1 or 3 letter char "
                  "string. Return space on error.")
      .def("__repr__",
           [](const PDBRes &self) {
             return "<buccaneer.ProteinDB.Residue: " +
                    clipper::String(self.type()) + " >";
           })
      .doc() =
      "Class for storing compact amino acid main chain information for use "
      "with Top500 DB. The same class is used for storing residues in the "
      "DB, and also for storing residues for searching against the DB. "
      "This class stores all the information concerning a single residue.";

  using ResTM = PDBRes::TypeMask;
  nb::class_<ResTM>(protdb_res, "TypeMask")
      .def(nb::init<>())
      .def(nb::init<const char>(), nb::arg("type_1_char"),
           "Initialise from a residue code, '?' for all.")
      .def_prop_ro(
          "mask", &ResTM::mask,
          "Return mask. Use '&' to test if a mask matches a type.")
      .def(
          "__or__",
          [](const ResTM &self, const ResTM &other) {
            return self.mask() | other.mask();
          },
          nb::is_operator())
      .def(
          "__and__",
          [](const ResTM &self, const ResTM &other) {
            return self.mask() & other.mask();
          },
          nb::is_operator())
      .def(
          "__not__", [](const ResTM &self) { return !self.mask(); },
          nb::is_operator())
      .doc() = "Class for describing a residue type mask, used to describe a "
               "list of allowed residue types.";
}

void declare_chain(nb::module_ &m) {
  using PDBChn = ProteinDB::Chain;
  nb::class_<PDBChn>(m, "Chain")
      .def(nb::init<>())
      .def("add_pdb", &PDBChn::add_pdb)
      .def("add_residue", &PDBChn::add_residue)
      .def("save_db", &PDBChn::save_db)
      .def("load_db", &PDBChn::load_db)
      .def("merge", &PDBChn::merge)
      .def("extract", &PDBChn::extract)
      .def("is_continuous", &PDBChn::is_continuous)
      .def("transform", &PDBChn::transform)
      .def("lsq_superpose",
           (void(PDBChn::*)(const PDBChn &)) & PDBChn::lsq_superpose,
           nb::arg("frag"))
      .def("lsq_superpose",
           (void(PDBChn::*)(const PDBChn &, const std::vector<double> &)) &
               PDBChn::lsq_superpose,
           nb::arg("frag"), nb::arg("weights"))
      .def("rmsd", (double(PDBChn::*)(const PDBChn &) const) & PDBChn::rmsd,
           nb::arg("other"))
      .def("rmsd",
           (double(PDBChn::*)(const PDBChn &, const std::vector<double> &)
                const) &
               PDBChn::rmsd,
           nb::arg("other"), nb::arg("weights"))
      .def(
          "__getitem__",
          [](PDBChn &self, const int &i) -> const ProteinDB::Residue & {
            return self[i];
          },
          nb::arg("i"), nb::rv_policy::reference_internal)
      .def(
          "__setitem__",
          [](PDBChn &self, const int &i, ProteinDB::Residue res) {
            self[i] = res;
          },
          nb::arg("i"), nb::arg("res"),
          nb::rv_policy::reference_internal)
      .def("size", &PDBChn::size)
      .def("__len__", &PDBChn::size)
      .def("debug", &PDBChn::debug)
      .def("__repr__",
           [](const PDBChn &self) {
             return "<buccaneer.ProteinDB.Chain with " +
                    clipper::String(int(self.size())) + " residues.>";
           })
      .doc() =
      "Class for storing compact amino acid main chain information for use "
      "with Top500 DB. The same class is used for storing residues in the "
      "DB, and also for storing residues for searching against the DB. "
      "This class stores a list of residues representing either a complete "
      "chain or chains (in the DB) or a complete search fragment to be "
      "searched against the DB.";

  using CDB = ProteinDB::ChainDB;
  nb::class_<CDB, PDBChn>(m, "ChainDB")
      .def(nb::init<>())
      .def(nb::init<const PDBChn &>(), "Constructor from Chain.")
      .def(nb::init<const String>(), nb::arg("file"),
           "Constructor from binary DB file name.")
      .def("init", &CDB::init, "Initialiser from binary DB file name.")
      .def("calc_distances", &CDB::calc_distances,
           "Calculate the distance matrix elements.")
      .def("score_distances",
           (double(CDB::*)(const CDB &, int) const) & CDB::score_distance,
           nb::arg("frag"), nb::arg("offset"),
           "Score a fragment against the nth fragment in the DB.")
      .def("score_distances",
           (double(CDB::*)(const CDB &, int, double) const) &
               CDB::score_distance,
           nb::arg("frag"), nb::arg("offset"), nb::arg("score_cutoff"),
           "Score a fragment against the nth fragment in the DB with cutoff.")
      .def("score_distances",
           (double(CDB::*)(const CDB &,
                           const std::vector<ProteinDB::Residue::TypeMask> &,
                           int, double) const) &
               CDB::score_distance,
           nb::arg("frag"), nb::arg("types"), nb::arg("offset"),
           nb::arg("score_cutoff"),
           "Score a fragment against the nth fragment in the DB with residue "
           "masks.")
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(const CDB &, int) const) &
               CDB::match_fragment_preliminary,
           nb::arg("fragdb"), nb::arg("nhit"),
           "Return list of tentative fragment offsets matching a given "
           "fragment.")
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int) const) &
               CDB::match_fragment_preliminary,
           nb::arg("fragdb"), nb::arg("type_mask"), nb::arg("nhit"),
           "Return list of tentative fragment offsets with residue masks.")
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(const CDB &, int, int) const) &
               CDB::match_fragment,
           nb::arg("fragdb"), nb::arg("nlsq"), nb::arg("nhit") = 0,
           "Return the best DB fragments matching a given fragment.")
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int, int) const) &
               CDB::match_fragment,
           nb::arg("fragdb"), nb::arg("types"), nb::arg("nlsq"),
           nb::arg("nhit") = 0,
           "Return the best DB fragments matching a given fragment with "
           "residue masks.")
      .def("__repr__",
           [](const PDBChn &self) {
             return "<buccaneer.ProteinDB.ChainDB class.>";
           })
      .doc() =
      "Class for storing compact amino acid main chain information for "
      "use with Top500 DB. The same class is used for storing residues in "
      "the DB, and also for storing residues for searching against the DB. "
      "This class is an extension of Chain, which also adds the fast "
      "distance matrix entries. These are filled out by calling the "
      "calc_distances() method after the chain has been completely "
      "assembled.";
}

void add_proteindb(nb::module_ &m) {
  declare_residue(m);
  declare_chain(m);
}