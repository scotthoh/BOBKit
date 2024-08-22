// Wrapper for protein_db
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/protein_db.h>
#include <clipper/minimol/minimol.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helper_functions.h"
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_residue(py::module &m) {
  using PDBRes = ProteinDB::Residue;
  py::class_<PDBRes> protdb_res(m, "Residue");

  py::enum_<PDBRes::FLAG>(protdb_res, "FLAG")
      .value("NONE", PDBRes::FLAG::NONE)
      .value("NORMAL", PDBRes::FLAG::NORMAL)
      .value("CALPHA", PDBRes::FLAG::CALPHA)
      .export_values();

  protdb_res.def(py::init<>())
      .def(py::init<Coord_orth &, const String &>(), py::arg("Calpha"),
           py::arg("type"),
           "Constructor from C-alpha coordinates and residue type.")
      .def(py::init<Coord_orth &, Coord_orth &, Coord_orth &, const String &>(),
           py::arg("Cn"), py::arg("Calpha"), py::arg("Cc"), py::arg("type"),
           "Constructor from main chain coordinates and residue type.")
      .def(py::init<const MMonomer &>(), py::arg("monomer"),
           "Constructor from MMonomer.")
      .def_property_readonly("coord_n", &PDBRes::coord_n,
                             "Get N atom coordinate.")
      .def_property_readonly("coord_ca", &PDBRes::coord_ca,
                             "Get C-alpha atom coordinate.")
      .def_property_readonly("coord_c", &PDBRes::coord_c,
                             "Get C atom coordinate.")
      .def_property_readonly("mmonomer", &PDBRes::mmonomer, "Build monomer.")
      .def("transform", &PDBRes::transform, py::arg("rtop"),
           "Transform by rotation-translation operator.")
      .def("merge", &PDBRes::merge, py::arg("res"), py::arg("wn"),
           py::arg("wa"), py::arg("wc"),
           "Merge other residue coords with this one using given weights.")
      .def_property("type", &PDBRes::type, &PDBRes::set_type,
                    "Get and set 1-letter residue types.")
      .def_property("flag", &PDBRes::flag, &PDBRes::set_flag,
                    "Get and set flag.")
      //.def("set_type", &PDBRes::set_type)
      //.def("set_flag", &PDBRes::set_flag)
      .def("data_import", &PDBRes::data_import, py::arg("d"),
           "Import from char array.")
      .def("data_export", &PDBRes::data_export, py::arg("d"),
           "Export to char array")
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
  py::class_<ResTM>(protdb_res, "TypeMask")
      .def(py::init<>())
      .def(py::init<const char>(), py::arg("type_1_char"),
           "Initialise from a residue code, '?' for all.")
      .def_property_readonly(
          "mask", &ResTM::mask,
          "Return mask. Use '&' to test if a mask matches a type.")
      .def(
          "__or__",
          [](const ResTM &self, const ResTM &other) {
            return self.mask() | other.mask();
          },
          py::is_operator())
      .def(
          "__and__",
          [](const ResTM &self, const ResTM &other) {
            return self.mask() & other.mask();
          },
          py::is_operator())
      .def(
          "__not__", [](const ResTM &self) { return !self.mask(); },
          py::is_operator())
      .doc() = "Class for describing a residue type mask, used to describe a "
               "list of allowed residue types.";
}

void declare_chain(py::module &m) {
  using PDBChn = ProteinDB::Chain;
  py::class_<PDBChn>(m, "Chain")
      .def(py::init<>())
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
           py::arg("frag"))
      .def("lsq_superpose",
           (void(PDBChn::*)(const PDBChn &, const std::vector<double> &)) &
               PDBChn::lsq_superpose,
           py::arg("frag"), py::arg("weights"))
      .def("rmsd", (double(PDBChn::*)(const PDBChn &) const) & PDBChn::rmsd,
           py::arg("other"))
      .def("rmsd",
           (double(PDBChn::*)(const PDBChn &, const std::vector<double> &)
                const) &
               PDBChn::rmsd,
           py::arg("other"), py::arg("weights"))
      .def(
          "__getitem__",
          [](PDBChn &self, const int &i) -> const ProteinDB::Residue & {
            return self[i];
          },
          py::arg("i"), py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](PDBChn &self, const int &i, ProteinDB::Residue res) {
            self[i] = res;
          },
          py::arg("i"), py::arg("res"),
          py::return_value_policy::reference_internal)
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
  py::class_<CDB, PDBChn>(m, "ChainDB")
      .def(py::init<>())
      .def(py::init<const PDBChn &>(), "Constructor from Chain.")
      .def(py::init<const String>(), py::arg("file"),
           "Constructor from binary DB file name.")
      .def("init", &CDB::init, "Initialiser from binary DB file name.")
      .def("calc_distances", &CDB::calc_distances,
           "Calculate the distance matrix elements.")
      .def("score_distances",
           (double(CDB::*)(const CDB &, int) const) & CDB::score_distance,
           py::arg("frag"), py::arg("offset"),
           "Score a fragment against the nth fragment in the DB.")
      .def("score_distances",
           (double(CDB::*)(const CDB &, int, double) const) &
               CDB::score_distance,
           py::arg("frag"), py::arg("offset"), py::arg("score_cutoff"),
           "Score a fragment against the nth fragment in the DB with cutoff.")
      .def("score_distances",
           (double(CDB::*)(const CDB &,
                           const std::vector<ProteinDB::Residue::TypeMask> &,
                           int, double) const) &
               CDB::score_distance,
           py::arg("frag"), py::arg("types"), py::arg("offset"),
           py::arg("score_cutoff"),
           "Score a fragment against the nth fragment in the DB with residue "
           "masks.")
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(const CDB &, int) const) &
               CDB::match_fragment_preliminary,
           py::arg("fragdb"), py::arg("nhit"),
           "Return list of tentative fragment offsets matching a given "
           "fragment.")
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int) const) &
               CDB::match_fragment_preliminary,
           py::arg("fragdb"), py::arg("type_mask"), py::arg("nhit"),
           "Return list of tentative fragment offsets with residue masks.")
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(const CDB &, int, int) const) &
               CDB::match_fragment,
           py::arg("fragdb"), py::arg("nlsq"), py::arg("nhit") = 0,
           "Return the best DB fragments matching a given fragment.")
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int, int) const) &
               CDB::match_fragment,
           py::arg("fragdb"), py::arg("types"), py::arg("nlsq"),
           py::arg("nhit") = 0,
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

void init_proteindb(py::module &m) {
  declare_residue(m);
  declare_chain(m);
}