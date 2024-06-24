#include "buccaneer/protein_db.h"

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/minimol/minimol.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace clipper;

void declare_residue(py::module &m) {
  using PDBRes = ProteinDB::Residue;
  py::class_<PDBRes> protdb_res(m, "Residue");

  py::enum_<PDBRes::FLAG>(protdb_res, "FLAG")
      .value("NONE", PDBRes::FLAG::NONE)
      .value("NORMAL", PDBRes::FLAG::NORMAL)
      .value("CALPHA", PDBRes::FLAG::CALPHA);

  protdb_res.def(py::init<>())
      .def(py::init<Coord_orth &, const String &>(), py::arg("Calpha"),
           py::arg("type"))
      .def(py::init<Coord_orth &, Coord_orth &, Coord_orth &, const String &>(),
           py::arg("Cn"), py::arg("Calpha"), py::arg("Cc"), py::arg("type"))
      .def(py::init<const MMonomer &>(), py::arg("monomer"))
      .def_property_readonly("coord_n", &PDBRes::coord_n)
      .def_property_readonly("coord_ca", &PDBRes::coord_ca)
      .def_property_readonly("coord_c", &PDBRes::coord_c)
      .def_property_readonly("mmonomer", &PDBRes::mmonomer)
      .def("transform", &PDBRes::transform)
      .def("merge", &PDBRes::merge)
      .def_property("type", &PDBRes::type, &PDBRes::set_type)
      .def_property("flag", &PDBRes::flag, &PDBRes::set_flag)
      //.def("set_type", &PDBRes::set_type)
      //.def("set_flag", &PDBRes::set_flag)
      .def("data_import", &PDBRes::data_import)
      .def("data_export", &PDBRes::data_export)
      .def("is_null", &PDBRes::is_null)
      .def_static("residue_type", &PDBRes::residue_type)
      .def("__repr__", [](const PDBRes &self) {
        return "<buccaneer.ProteinDB.Residue: " + clipper::String(self.type()) +
               " >";
      });

  using ResTM = PDBRes::TypeMask;
  py::class_<ResTM>(protdb_res, "TypeMask")
      .def(py::init<>())
      .def(py::init<const char>(), py::arg("type_1_char"))
      .def_property_readonly("mask", &ResTM::mask)
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
          py::is_operator());
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
      .def("__repr__", [](const PDBChn &self) {
        return "<buccaneer.ProteinDB.Chain with " +
               clipper::String(int(self.size())) + " residues.>";
      });

  using CDB = ProteinDB::ChainDB;
  py::class_<CDB, PDBChn>(m, "ChainDB")
      .def(py::init<>())
      .def(py::init<const PDBChn &>())
      .def(py::init<const String>(), py::arg("file"))
      .def("init", &CDB::init)
      .def("calc_distances", &CDB::calc_distances)
      .def("score_distances",
           (double(CDB::*)(const CDB &, int) const) & CDB::score_distance,
           py::arg("frag"), py::arg("offset"))
      .def("score_distances",
           (double(CDB::*)(const CDB &, int, double) const) &
               CDB::score_distance,
           py::arg("frag"), py::arg("offset"), py::arg("score_cutoff"))
      .def("score_distances",
           (double(CDB::*)(const CDB &,
                           const std::vector<ProteinDB::Residue::TypeMask> &,
                           int, double) const) &
               CDB::score_distance,
           py::arg("frag"), py::arg("type_mask"), py::arg("offset"),
           py::arg("score_cutoff"))
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(const CDB &, int) const) &
               CDB::match_fragment_preliminary,
           py::arg("fragdb"), py::arg("nhit"))
      .def("match_fragment_preliminary",
           (std::vector<int>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int) const) &
               CDB::match_fragment_preliminary,
           py::arg("fragdb"), py::arg("type_mask"), py::arg("nhit"))
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(const CDB &, int, int) const) &
               CDB::match_fragment,
           py::arg("fragdb"), py::arg("nlsq"), py::arg("nhit") = 0)
      .def("match_fragment",
           (std::vector<PDBChn>(CDB::*)(
               const CDB &, const std::vector<ProteinDB::Residue::TypeMask> &,
               int, int) const) &
               CDB::match_fragment,
           py::arg("fragdb"), py::arg("type_mask"), py::arg("nlsq"),
           py::arg("nhit") = 0)
      .def("__repr__", [](const PDBChn &self) {
        return "<buccaneer.ProteinDB.ChainDB class.>";
      });
}

void init_proteindb(py::module &m) {
  declare_residue(m);
  declare_chain(m);
}