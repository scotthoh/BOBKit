// Nanobind bindings for simplex-lib
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-find.h"
#include "buccaneer/buccaneer-grow.h"
#include "buccaneer/buccaneer-build.h"
#include "buccaneer/simplex-lib.h"
#include "commons.h"
#include "arrays.h"
#include <nanobind/operators.h>
#include <nanobind/trampoline.h>

using namespace clipper;

void declare_target_functions(nb::module_ &m) {
  // Define trampoline class for abstract class
  class PyTarget_fn_order_zero : public Target_fn_order_zero {
  public:
    NB_TRAMPOLINE( Target_fn_order_zero, 2 );

    // Trampoline
    int num_params() const override { NB_OVERRIDE_PURE( num_params ); }
    double operator()(const std::vector<double> &args) const override {
      NB_OVERRIDE_PURE_NAME( "__call__", operator(), args );
      // "function name", function, arguments
    }
  };
  // Define as such to be used by private inherited class, need it so
  // in python other functions that require Target_fn_order_zero
  // as argument knows the inherited class is from the abstract class
  // https://github.com/pybind/pybind11/issues/1383
  auto clsBase =
      nb::class_<Target_fn_order_zero, PyTarget_fn_order_zero>(
          m, "Target_fn_order_zero",
          "Abstract base class for zero-th order function.")
          .def(nb::init<>())
          .def("num_params", &Target_fn_order_zero::num_params)
          .def("__call__", &Target_fn_order_zero::operator(), nb::arg("args"));

  // Both the inherited class defined here to be at the same
  // scope as the base asbstract class
  using refine_nterm = Target_fn_refine_n_terminal_build;
  nb::class_<refine_nterm>(m, "Target_fn_refine_n_terminal_build", clsBase)
      .def(nb::init<>())
      .def(nb::init<const Xmap<float> &, const LLK_map_target &,
                    const Ramachandran &, const Ramachandran &,
                    const double &>(),
           nb::arg("xmap"), nb::arg("llktarget"), nb::arg("rama1"),
           nb::arg("rama2"), nb::arg("rot_step"))
      .def("num_params", &refine_nterm::num_params,
           "Return number of parameters.")
      .def("__call__", &refine_nterm::operator(), nb::arg("args"),
           "Evaluate target function for EulerXYZr offset from rotation.")
      .def("refine", &refine_nterm::refine, nb::arg("chain"), nb::arg("args"),
           "Refine rotation.")
      .def("__repr__",
           [](const refine_nterm &self) {
             return "<buccaneer.Target_fn_refine_n_terminal_build class.>";
           })
      .doc() = "Class for refining grown Ca groups.";

  using refine_llk = Target_fn_refine_llk_map_target;
  nb::class_<refine_llk>(m, "Target_fn_refine_llk_map_target", clsBase)
      .def(nb::init<>())
      .def(nb::init<const Xmap<float> &, const LLK_map_target &, const double &,
                    const double &>(),
           nb::arg("xmap"), nb::arg("llktarget"), nb::arg("rot_step"),
           nb::arg("trn_step"),
           "Constructor with density map, LLK target, rotation and "
           "translation steps.")
      .def("num_params", &refine_llk::num_params,
           "Return number of parameters.")
      .def("__call__", nb::overload_cast<const RTop_orth &>(& refine_llk::operator(), nb::const_),
           nb::arg("rtop"), "Evaluate target function for given rotation.")
      .def("__call__",
           nb::overload_cast<const std::vector<double> &>(& refine_llk::operator(), nb::const_),
           nb::arg("args"), "Evaluate target function for EulerXYZr offset from rotation.")
      .def("rtop_orth", &refine_llk::rtop_orth, nb::arg("args"),
           "Convert parameters to rotation.")
      .def("refine", &refine_llk::refine, nb::arg("rtop"), "Refine rotation.")
      .def("__repr__",
           [](const refine_llk &self) {
             return "<buccaneer.Target_fn_refine_llk_map_target class.>";
           })
      .doc() = "Class for refining Ca groups.";

//  using refine_fragment = Target_fn_refine_amino_acid_fragment;
//  nb::class_<refine_fragment>(m, "Target_fn_refine_amino_acid_fragment", clsBase)
//      .def(nb::init<>())
//      .def(nb::init<const Xmap<float> &, const double &, const double &>(),
//           nb::arg("xmap"), nb::arg("rot_step"), nb::arg("trn_step"))
//      .def("num_params", &refine_fragment::num_params,
//           "Return number of parameters.")
//      .def("__call__", &refine_fragment::operator(), nb::arg("args"))
//      .def("rtop_orth", &refine_fragment::rtop_orth, nb::arg("args"))
//      .def("refine", &refine_fragment::refine, nb::arg("residue"))
//      .def("__repr__", [](const refine_fragment &self) {
//        return "<buccaneer.Target_fn_refine_amino_acid_fragment class.>";
//      })
//      .doc() = "Class for refining built amino acid fragment";
}

void declare_optimiser_simplex(nb::module_ &m) {
  nb::class_<Optimiser_simplex> opt_simp(m, "Optimiser_simplex",
                                         "Simplex optimiser.");

  nb::enum_<Optimiser_simplex::TYPE>(opt_simp, "TYPE", "Optimiser type.")
      .value("NORMAL", Optimiser_simplex::TYPE::NORMAL)
      .value("GRADIENT", Optimiser_simplex::TYPE::GRADIENT)
      .export_values();

  opt_simp
      .def(nb::init<double, int, Optimiser_simplex::TYPE>(),
           nb::arg("tolerance") = 0.001, nb::arg("max_cycles") = 50,
           nb::arg("type") = Optimiser_simplex::TYPE::NORMAL,
           "Constructor with tolerance, max cycle and optimiser type.")
      .def("__call__", &Optimiser_simplex::operator(), nb::arg("target_fn"),
           nb::arg("args"), "Run optimiser.")
      .def("debug", &Optimiser_simplex::debug, "Turn on debug mode.")
      .def("__repr__", [](const Optimiser_simplex &self) {
        return "<buccaneer.Optimiser_simplex class.>";
      });
}

void add_simplex_lib(nb::module_ &m) {
  declare_target_functions(m);
  declare_optimiser_simplex(m);
}