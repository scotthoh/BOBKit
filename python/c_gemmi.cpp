// #include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper-minimol.h>
#include <clipper/core/clipper_types.h>
// #include <gemmi/symmetry.hpp>
#include <clipper/minimol/minimol_io_gemmi.h>
#include <gemmi/model.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

py::object structure =
    (py::object)py::module_::import("gemmi").attr("Structure");

void declare_gemmifile(py::module &m) {
  py::class_<GEMMIfile> gemmifile(m, "GEMMIfile");
  py::enum_<GEMMIfile::TYPE>(gemmifile, "TYPE")
      .value("DEFAULT", GEMMIfile::TYPE::Default)
      .value("PDB", GEMMIfile::TYPE::PDB)
      .value("CIF", GEMMIfile::TYPE::CIF)
      .value("Mmjson", GEMMIfile::TYPE::Mmjson)
      .export_values();

  gemmifile.def(py::init<>())
      .def(
          "read_file",
          [](GEMMIfile &self, const std::string &file, int max_line_length,
             bool split_chain_on_ter, bool force_label) {
            gemmi::PdbReadOptions opts;
            opts.max_line_length = max_line_length;
            opts.split_chain_on_ter = split_chain_on_ter;
            self.read_file(file, opts, force_label);
          },
          py::arg("file"), py::arg("max_line_length") = 0,
          py::arg("split_chain_on_ter") = false, py::arg("force_label") = false)
      .def(
          "write_pdb",
          [](GEMMIfile &self, const std::string &file, const bool &minimal,
             const bool &shortTer, const bool &linkr, const bool &shortchn) {
            QuickWriteOptions opts;
            opts.Minimal = minimal;
            opts.PdbShortTer = shortTer;
            opts.PdbShortChainNames = shortchn;
            opts.PdbLinkR = linkr;
            self.write_file(file, GEMMIfile::TYPE::PDB, opts);
          },
          py::arg("file"), py::arg("minimal") = false,
          py::arg("shortTer") = false, py::arg("linkr") = false,
          py::arg("short_chain_names") = false)
      .def(
          "write_cif",
          [](GEMMIfile &self, const std::string &file, const bool &allauth,
             const bool &update) {
            QuickWriteOptions opts;
            opts.CifAllAuth = allauth;
            opts.UpdateCifDoc = update;
            self.write_file(file, GEMMIfile::TYPE::CIF, opts);
          },
          py::arg("file"), py::arg("allauth") = false,
          py::arg("updatecif") = false)
      .def("import_minimol", &GEMMIfile::import_minimol, py::arg("minimol"),
           py::arg("model_num") = 1)
      .def("export_minimol", &GEMMIfile::export_minimol, py::arg("minimol"))
      .def_property("gemmi_structure", &GEMMIfile::get_gemmi_structure,
                    &GEMMIfile::set_gemmi_structure)
      .def("set_mmcif_output_groups", &GEMMIfile::set_mmcif_output_groups)
      .def("set_pdb_write_options", &GEMMIfile::set_pdb_write_options)
      .def_property("spacegroup", &GEMMIfile::spacegroup,
                    &GEMMIfile::set_spacegroup)
      .def_property("cell", &GEMMIfile::cell, &GEMMIfile::set_cell)
      .def("set_cif_style", &GEMMIfile::set_cif_style);
}
void declare_gemmi_structure(py::module &m) {
  // importing bindings for Structure from gemmi
  // caution: need to make sure gemmi is built with same pybind11 and stl
  // versions
  py::class_<GemmiStructure>(m, "gemmiStructure", structure)
      .def(py::init<>())
      .def(py::init<const gemmi::Structure &>())
      .def("spacegroup", &GemmiStructure::spacegroup)
      .def_static("to_minimol",
                  [](const gemmi::Structure &s, MiniMol &mol) {
                    GEMMIfile gfile;
                    gfile.set_gemmi_structure(s);
                    gfile.import_minimol(mol);
                    if (!mol.is_null())
                      return true;
                    else
                      return false;
                  })
      .def("export_minimol",
           [](const GemmiStructure &self, MiniMol &mol) {
             GEMMIfile gfile;
             gfile.set_gemmi_structure(self);
             gfile.import_minimol(mol);
             if (!mol.is_null())
               return true;
             else
               return false;
           })
      .def("import_minimol", [](GemmiStructure &self, MiniMol &mol) {
        GEMMIfile gfile;
        gfile.export_minimol(mol);
        self = gfile.get_gemmi_structure();
        if (!self.models.empty())
          return true;
        else
          return false;
      });
}

void declare_clipper_gemmi(py::module &m) {
  py::class_<GEMMI>(m, "clipper_gemmi")
      //.def_static("spacegroup",
      //            (const gemmi::SpaceGroup *(GEMMI::*)(const Spacegroup &))
      //            &
      //                GEMMI::spacegroup)
      .def_static(
          "to_gemmi_spacegroup",
          [](const Spacegroup &spgr) { return GEMMI::spacegroup(spgr); })
      .def_static("to_clipper_spacegroup", [](const gemmi::SpaceGroup &spgr) {
        return GEMMI::spacegroup(spgr);
      });
}

void init_gemmi_structure(py::module &m) {
  declare_gemmi_structure(m);
  declare_gemmifile(m);
  declare_clipper_gemmi(m);
}