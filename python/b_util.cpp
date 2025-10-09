// Nanobind bindings for buccaneer-util
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-util.h"
#include "commons.h"
#include "arrays.h"
#include <clipper/minimol/minimol_seq.h>

using namespace clipper;

void declare_buccaneer_log(nb::module_ &m) {
  nb::class_<BuccaneerLog>(m, "Log", "Class for logging and simple profiling.")
      .def(nb::init<String &>(), nb::arg("title"))
      .def(
          "log", nb::overload_cast<const String &>(&BuccaneerLog::log),
          //[](BuccaneerLog &self, const String &id) { self.log(id); },
          nb::arg("id"), "Start time log for specified id.")
      .def(
          "log", nb::overload_cast<const String &, const MiniMol &, bool>(&BuccaneerLog::log),
          //[](BuccaneerLog &self, const String &id,
          //   const MiniMol &mol,
          //   bool view) { self.log(id, mol, view); },
          nb::arg("id"), nb::arg("mol"), nb::arg("view") = false,
          "Start time log for specified id and model.")
      .def(
          "summary",
          [](BuccaneerLog &self, const MiniMol &mol,
             const MiniMol &mr,
             const MMoleculeSequence &seq) { return self.log(mol, mr, seq); },
          nb::arg("mol"), nb::arg("mol_mr"), nb::arg("seq"),
          "Return summary of model stats as string.")
      .def("xml", &BuccaneerLog::xml, nb::arg("xml"), "Write XML output.")
      .def("profile", &BuccaneerLog::profile,
           "Outputs time profile for each id logged.")
      .def("evaluate", &BuccaneerLog::evaluate, nb::arg("mol"),
           nb::arg("mol_mr"), nb::arg("seq"), "Evaluate model.");
}

void add_buccaneer_util(nb::module_ &m) {
  m.def( "set_reference", &BuccaneerUtil::set_reference, nb::arg( "mtz" ), nb::arg( "pdb" ),
         "Set reference MTZ and PDB, only work when CCP4's CLIBD "
         "environment path is set." );
  // declare_buccaneer_util( m );
  declare_buccaneer_log(m);
}
