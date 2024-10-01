// Wrapper for clipper Ramachandran
// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace clipper;

void init_ramachandran(py::module &m) {
  py::class_<Ramachandran> ramachandran(m, "Ramachandran");

  py::enum_<Ramachandran::TYPE>(ramachandran, "TYPE",
                                "Enumeration of built-in Ramachandran tables.")
      .value("Gly", Ramachandran::TYPE::Gly)
      .value("Pro", Ramachandran::TYPE::Pro)
      .value("NonGlyPro", Ramachandran::TYPE::NonGlyPro)
      .value("NonGly", Ramachandran::TYPE::NonGly)
      .value("All", Ramachandran::TYPE::All)
      .value("Gly5", Ramachandran::TYPE::Gly5)
      .value("Pro5", Ramachandran::TYPE::Pro5)
      .value("NonGlyPro5", Ramachandran::TYPE::NonGlyPro5)
      .value("NonGly5", Ramachandran::TYPE::NonGly5)
      .value("All5", Ramachandran::TYPE::All5)
      .value("All2", Ramachandran::TYPE::All2)
      .value("Gly2", Ramachandran::TYPE::Gly2)
      .value("Pro2", Ramachandran::TYPE::Pro2)
      .value("PrePro2", Ramachandran::TYPE::PrePro2)
      .value("IleVal2", Ramachandran::TYPE::IleVal2)
      .value("NoGPIVpreP2", Ramachandran::TYPE::NoGPIVpreP2)
      .export_values();

  ramachandran.def(py::init<>(), "Null constructor.")
      .def(py::init<Ramachandran::TYPE>(), py::arg("type"),
           "Constructor from standard plot.")
      .def("init", &Ramachandran::init, py::arg("type"),
           "Initialise from standard plot.")
      .def("set_thresholds", &Ramachandran::set_thresholds,
           py::arg("prob_favoured") = 0.01, py::arg("prob_allowed") = 0.0005,
           "Change thresholds to difference values.")
      .def("probability", &Ramachandran::probability, py::arg("phi"),
           py::arg("psi"), "Get probability for a particular pair of angles.")
      .def("favoured", &Ramachandran::favored, py::arg("phi"), py::arg("psi"),
           "Test if a pair of angles are in the favoured region.")
      .def("favored", &Ramachandran::favored, py::arg("phi"), py::arg("psi"),
           "Test if a pair of angles are in the favoured region.")
      .def("allowed", &Ramachandran::allowed, py::arg("phi"), py::arg("psi"),
           "Test if a pair of angles are in an allowed region.")
      .def("__repr__",
           [](const Ramachandran &self) {
             return "<clipper.Ramachandran plot class.>";
           })
      .doc() = "Ramachandran plot class.\nThis class provides a reference "
               "Ramachandran plot for Gly, Pro, other, and combinations "
               "of those types of residues. The source data comes from the "
               "best residues from the 'top500' best-determined structures "
               "list of D. C. and J. S. Richardson, "
               "http://kinemage.biochem.duke.edu/index.html \n"
               "The Ramachandran plot is normalised in inverse radians squared,"
               "so the mean value of a probability is 1/(2 pi)<sup>2</sup>.";
}
