// Wrapper for clipper gemmi
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper-minimol-gemmi.h>
#include <clipper/core/clipper_types.h>
#include <clipper/gemmi/clipper_gemmi_model.h>
#include <clipper/minimol/minimol_io_gemmi.h>
#include <gemmi/model.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

// Importing bindings for Structure from gemmi
// caution: need to make sure gemmi is built with same pybind11 and stl
// versions
nb::object structure =
    (nb::object)nb::module_::import("gemmi").attr("Structure");

void declare_gemmifile(nb::module &m) {
  nb::class_<GEMMIfile> gemmifile(m, "GEMMIfile");
  nb::enum_<GEMMIfile::TYPE>(gemmifile, "TYPE", "File formats.")
      .value("DEFAULT", GEMMIfile::TYPE::Default)
      .value("PDB", GEMMIfile::TYPE::PDB)
      .value("CIF", GEMMIfile::TYPE::CIF)
      .value("Mmjson", GEMMIfile::TYPE::Mmjson)
      .export_values();

  gemmifile.def(nb::init<>())
      .def(
          "read_file",
          [](GEMMIfile &self, const std::string &file, int max_line_length,
             bool split_chain_on_ter, bool force_label) {
            gemmi::PdbReadOptions opts;
            opts.max_line_length = max_line_length;
            opts.split_chain_on_ter = split_chain_on_ter;
            self.read_file(file, opts, force_label);
          },
          nb::arg("file"), nb::arg("max_line_length") = 0,
          nb::arg("split_chain_on_ter") = false, nb::arg("force_label") = false,
          "Read structure from file.")
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
          nb::arg("file"), nb::arg("minimal") = false,
          nb::arg("shortTer") = false, nb::arg("linkr") = false,
          nb::arg("short_chain_names") = false,
          "Write structure as PDB format.")
      .def(
          "write_cif",
          [](GEMMIfile &self, const std::string &file, const bool &allauth,
             const bool &update) {
            QuickWriteOptions opts;
            opts.CifAllAuth = allauth;
            opts.UpdateCifDoc = update;
            self.write_file(file, GEMMIfile::TYPE::CIF, opts);
          },
          nb::arg("file"), nb::arg("allauth") = false,
          nb::arg("updatecif") = false, "Write structure as CIF format.")
      .def("import_minimol", &GEMMIfile::import_minimol, nb::arg("minimol"),
           nb::arg("model_num") = 1, "Import MiniMol from Gemmi Structure.")
      .def("export_minimol", &GEMMIfile::export_minimol, nb::arg("minimol"),
           "Export MiniMol to Gemmi Structure.")
      .def_property("gemmi_structure", &GEMMIfile::get_gemmi_structure,
                    &GEMMIfile::set_gemmi_structure,
                    "Get/set structure if not read from file.")
      .def("set_mmcif_output_groups", &GEMMIfile::set_mmcif_output_groups,
           "Set mmcif output groups.")
      .def("set_pdb_write_options", &GEMMIfile::set_pdb_write_options,
           "Set pdb write options.")
      .def_property("spacegroup", &GEMMIfile::spacegroup,
                    &GEMMIfile::set_spacegroup, "Get/set spacegroup.")
      .def_property("cell", &GEMMIfile::cell, &GEMMIfile::set_cell,
                    "Get/set cell from/to structure.")
      .def("set_cif_style", &GEMMIfile::set_cif_style, "Set cif output style.")
      .doc() = "GEMMI file object for MiniMol i/o.\nThis object is an i/o "
               "object for MiniMol, representing an interface between MiniMol"
               "and PDB or CIF files. The model is read using GEMMI implemented"
               "methods and then imported to MiniMol object and vice versa."
               "Works fine for gemmi==v0.6.4.";
}
void declare_gemmi_structure(nb::module &m) {
  nb::class_<GemmiStructure>(m, "gemmiStructure", structure)
      .def(nb::init<>())
      .def(nb::init<const gemmi::Structure &>(),
           "Constructor copy/convert from gemmi.Structure.")
      .def_property("clipper_spacegroup", &GemmiStructure::spacegroup,
                    &GemmiStructure::set_spacegroup,
                    "Get/set spacegroup (clipper format).")
      .def_property("clipper_cell", &GemmiStructure::get_cell,
                    &GemmiStructure::set_cell, "Get/set cell (clipper format).")
      .def("get_experimental_method", &GemmiStructure::get_exptlmethod,
           "Get experimental method from records.")
      .def_property("origx_as_rtop", &GemmiStructure::get_transform,
                    &GemmiStructure::set_transform,
                    "Get/set origx transformation matrix as/from RTop_orth.")
      .def_static("get_id_str", &GemmiStructure::GetID_str,
                  nb::arg("model_name"), nb::arg("cra"), nb::arg("entity"),
                  "Get IDs for model, chain, residue, or atom.")
      .def_static(
          "to_minimol",
          [](const gemmi::Structure &s, MiniMol &mol) {
            GEMMIfile gfile;
            gfile.set_gemmi_structure(s);
            gfile.import_minimol(mol);
            if (!mol.is_null())
              return true;
            else
              return false;
          },
          "Static function to import structure as MiniMol.")
      .def(
          "export_minimol",
          [](const GemmiStructure &self, MiniMol &mol) {
            GEMMIfile gfile;
            gfile.set_gemmi_structure(self);
            gfile.import_minimol(mol);
            if (!mol.is_null())
              return true;
            else
              return false;
          },
          "Export structure to MiniMol.")
      .def(
          "import_minimol",
          [](GemmiStructure &self, MiniMol &mol) {
            GEMMIfile gfile;
            gfile.export_minimol(mol);
            self = gfile.get_gemmi_structure();
            if (!self.models.empty())
              return true;
            else
              return false;
          },
          "Import structure from MiniMol.")
      .doc() = "Gemmi Structure object wrapper.\nThis class is a "
               "trivial derivation of the corresponding Gemmi, "
               "providing access in terms of Clipper types. Thus, "
               "when you need such access, simply cast the Gemmi object "
               "reference to this type to access additional functions.";
}

void init_gemmi_structure(nb::module &m) {
  declare_gemmi_structure(m);
  declare_gemmifile(m);
}