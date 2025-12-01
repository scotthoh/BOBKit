// Wrapper for clipper gemmi
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper-minimol-gemmi.h>
#include <clipper/core/clipper_types.h>

#include <gemmi/model.hpp>

using namespace clipper;

// Importing bindings for Structure from gemmi
// caution: need to make sure gemmi is built with same pybind11 and stl
// versions
//nb::object structure = nb::module_::import_("gemmi").attr("Structure");
//
void declare_gemmi_structure(nb::module_ &m) {
  nb::class_<GemmiStructure>(m, "gemmiStructure") //, nb::module_::import_("gemmi").attr("Structure"))
      .def(nb::init<>())
      .def(nb::init<const gemmi::Structure &>(),
           "Constructor copy/convert from gemmi.Structure.")
      .def_prop_rw("clipper_spacegroup", &GemmiStructure::spacegroup,
                    &GemmiStructure::set_spacegroup,
                    "Get/set spacegroup (clipper format).")
      .def_prop_rw("clipper_cell", &GemmiStructure::get_cell,
                    &GemmiStructure::set_cell, "Get/set cell (clipper format).")
      .def("get_experimental_method", &GemmiStructure::get_exptlmethod,
           "Get experimental method from records.")
      .def_prop_rw("origx_as_rtop", &GemmiStructure::get_transform,
                    &GemmiStructure::set_transform,
                    "Get/set origx transformation matrix as/from RTop_orth.")
      .def_static("get_id_str", 
                  static_cast<String ( * )( const int &, const gemmi::CRA &, const String &)>(& GemmiStructure::GetID_str),
                  nb::arg("model_num"), nb::arg("cra"), nb::arg("entity"),
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

void init_gemmi_structure(nb::module_ &m) {
  declare_gemmi_structure(m);
  //declare_gemmifile(m);
}