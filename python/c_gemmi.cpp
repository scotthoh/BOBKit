// #include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper-minimol.h>
#include <clipper/core/clipper_types.h>
// #include <gemmi/symmetry.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

py::object structure =
    (py::object)py::module_::import("gemmi").attr("Structure");

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
      //            (const gemmi::SpaceGroup *(GEMMI::*)(const Spacegroup &)) &
      //                GEMMI::spacegroup)
      .def_static(
          "to_gemmi_spacegroup",
          [](const Spacegroup &spgr) { return GEMMI::spacegroup(spgr); })
      .def_static("to_clipper_spacegroup", [](const gemmi::SpaceGroup &spgr) {
        return GEMMI::spacegroup(spgr);
      });
}

void init_gemmi_structure(py::module &m) {
  // declare_gemmi_structure(m);
  declare_clipper_gemmi(m);
}