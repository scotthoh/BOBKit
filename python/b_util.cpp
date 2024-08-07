// Author: S.W.Hoh
// (C) 2023 -
// York Structural Biology Laboratory
// The University of York

#include "helper_functions.h"
#include "type_conversions.h"
#include <buccaneer/buccaneer-util.h>
#include <clipper/minimol/minimol_seq.h>
#include <pybind11/pybind11.h>

void declare_buccaneer_util(py::module &m) {
  py::class_<BuccaneerUtil>(m, "Util")
      .def_static("set_reference", &BuccaneerUtil::set_reference,
                  py::arg("mtz"), py::arg("pdb"))
      .def_static("read_model", &BuccaneerUtil::read_model, py::arg("mol"),
                  py::arg("file"), py::arg("verbose"));
}

void declare_buccaneer_log(py::module &m) {
  py::class_<BuccaneerLog>(m, "Log")
      .def(py::init<clipper::String &>(), py::arg("title"))
      //.def("log",
      //     static_cast<void (BuccaneerLog::*)(const clipper::String &)>(
      //         &BuccaneerLog::log),
      //     py::arg("id"))
      //.def("log",
      //     static_cast<void (BuccaneerLog::*)(const clipper::String &,
      //                                        const clipper::MiniMol &,
      //                                        bool)>(
      //         &BuccaneerLog::log),
      //     py::arg("id"), py::arg("mol"), py::arg("view"))
      .def(
          "log",
          [](BuccaneerLog &self, const clipper::String &id) { self.log(id); },
          py::arg("id"))
      .def(
          "log",
          [](BuccaneerLog &self, const clipper::String &id,
             const clipper::MiniMol &mol,
             bool view) { self.log(id, mol, view); },
          py::arg("id"), py::arg("mol"), py::arg("view") = false)
      .def(
          "summary",
          [](BuccaneerLog &self, const clipper::MiniMol &mol,
             const clipper::MiniMol &mr,
             const MMoleculeSequence &seq) { return self.log(mol, mr, seq); },
          py::arg("mol"), py::arg("mol_mr"), py::arg("seq"))
      .def("xml", &BuccaneerLog::xml, py::arg("xml"))
      .def("profile", &BuccaneerLog::profile)
      .def("evaluate", &BuccaneerLog::evaluate, py::arg("mol"),
           py::arg("mol_mr"), py::arg("seq"));
}

void init_buccaneer_util(py::module &m) {
  declare_buccaneer_util(m);
  declare_buccaneer_log(m);
}
