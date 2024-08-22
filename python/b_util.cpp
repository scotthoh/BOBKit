// Wrapper for buccaneer-util
// Author: S.W.Hoh
// (C) 2023 -
// York Structural Biology Laboratory
// The University of York

#include <buccaneer/buccaneer-util.h>
#include <clipper/minimol/minimol_seq.h>
#include <pybind11/pybind11.h>

#include "helper_functions.h"
#include "type_conversions.h"

void declare_buccaneer_util(py::module &m) {
  py::class_<BuccaneerUtil>(m, "Util", "Utility functions.")
      .def_static("set_reference", &BuccaneerUtil::set_reference,
                  py::arg("mtz"), py::arg("pdb"),
                  "Set reference MTZ and PDB, only work when CCP4's CLIBD "
                  "environment path is set.")
      .def_static("read_model", &BuccaneerUtil::read_model, py::arg("mol"),
                  py::arg("file"), py::arg("verbose") = true,
                  "Read model and export to minimol.");
}

void declare_buccaneer_log(py::module &m) {
  py::class_<BuccaneerLog>(m, "Log", "Class for logging and simple profiling.")
      .def(py::init<clipper::String &>(), py::arg("title"))
      .def(
          "log",
          [](BuccaneerLog &self, const clipper::String &id) { self.log(id); },
          py::arg("id"), "Start time log for specified id.")
      .def(
          "log",
          [](BuccaneerLog &self, const clipper::String &id,
             const clipper::MiniMol &mol,
             bool view) { self.log(id, mol, view); },
          py::arg("id"), py::arg("mol"), py::arg("view") = false,
          "Start time log for specified id and model.")
      .def(
          "summary",
          [](BuccaneerLog &self, const clipper::MiniMol &mol,
             const clipper::MiniMol &mr,
             const MMoleculeSequence &seq) { return self.log(mol, mr, seq); },
          py::arg("mol"), py::arg("mol_mr"), py::arg("seq"),
          "Return summary of model stats as string.")
      .def("xml", &BuccaneerLog::xml, py::arg("xml"), "Write XML output.")
      .def("profile", &BuccaneerLog::profile,
           "Outputs time profile for each id logged.")
      .def("evaluate", &BuccaneerLog::evaluate, py::arg("mol"),
           py::arg("mol_mr"), py::arg("seq"), "Evaluate model.");
}

void init_buccaneer_util(py::module &m) {
  declare_buccaneer_util(m);
  declare_buccaneer_log(m);
}
