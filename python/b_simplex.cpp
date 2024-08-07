#include <buccaneer/buccaneer-find.h>
#include <buccaneer/buccaneer-grow.h>
#include <buccaneer/simplex-lib.h>

#include "helper_functions.h"
#include "type_conversions.h"
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void declare_target_functions(py::module &m) {
  // Define trampoline class for abstract class
  class PyTarget_fn_order_zero : public Target_fn_order_zero {
  public:
    // inherit constructors
    using Target_fn_order_zero::Target_fn_order_zero;
    // Trampoline
    int num_params() const override {
      PYBIND11_OVERRIDE_PURE(int,                  // return type
                             Target_fn_order_zero, // Parent class
                             num_params,           // function name
      );
    }
    double operator()(const std::vector<double> &args) const override {
      PYBIND11_OVERRIDE_PURE_NAME(double, Target_fn_order_zero,
                                  "__call__", operator(),
                                  args // arguments
      );
    }
  };
  // Define as such to be used by private inherited class, need it so
  // in python other functions that require Target_fn_order_zero
  // as argument knows the inherited class is from the abstract class
  // https://github.com/pybind/pybind11/issues/1383
  auto clsBase =
      py::class_<Target_fn_order_zero, PyTarget_fn_order_zero>(
          m, "Target_fn_order_zero")
          .def(py::init<>())
          .def("num_params", &Target_fn_order_zero::num_params)
          .def("__call__", &Target_fn_order_zero::operator(), py::arg("args"));

  // Both the inherited class defined here to be at the same
  // scope as the base asbstract class
  using refine_nterm = Target_fn_refine_n_terminal_build;
  py::class_<refine_nterm>(m, "Target_fn_refine_n_terminal_build", clsBase)
      .def(py::init<>())
      .def(py::init<const Xmap<float> &, const LLK_map_target &,
                    const Ramachandran &, const Ramachandran &,
                    const double &>(),
           py::arg("xmap"), py::arg("llktarget"), py::arg("rama1"),
           py::arg("rama2"), py::arg("rot_step"))
      .def("num_params", &refine_nterm::num_params)
      .def("__call__", &refine_nterm::operator(), py::arg("args"))
      .def("refine", &refine_nterm::refine, py::arg("chain"), py::arg("args"))
      .def("__repr__", [](const refine_nterm &self) {
        return "<buccaneer.Target_fn_refine_n_terminal_build class.>";
      });

  using refine_llk = Target_fn_refine_llk_map_target;
  py::class_<refine_llk>(m, "Target_fn_refine_llk_map_target", clsBase)
      .def(py::init<>())
      .def(py::init<const Xmap<float> &, const LLK_map_target &, const double &,
                    const double &>(),
           py::arg("xmap"), py::arg("llktarget"), py::arg("rot_step"),
           py::arg("trn_step"))
      .def("num_params", &refine_llk::num_params)
      .def("__call__",
           (double(refine_llk::*)(const RTop_orth &) const) &
               refine_llk::operator(),
           py::arg("rtop"))
      .def("__call__",
           (double(refine_llk::*)(const std::vector<double> &) const) &
               refine_llk::operator(),
           py::arg("args"))
      .def("rtop_orth", &refine_llk::rtop_orth, py::arg("args"))
      .def("refine", &refine_llk::refine, py::arg("rtop"))
      .def("__repr__", [](const refine_llk &self) {
        return "<buccaneer.Target_fn_refine_llk_map_target class.>";
      });
}

void declare_optimiser_simplex(py::module &m) {
  py::class_<Optimiser_simplex> opt_simp(m, "Optimiser_simplex");

  py::enum_<Optimiser_simplex::TYPE>(opt_simp, "TYPE")
      .value("NORMAL", Optimiser_simplex::TYPE::NORMAL)
      .value("GRADIENT", Optimiser_simplex::TYPE::GRADIENT);

  opt_simp
      .def(py::init<double, int, Optimiser_simplex::TYPE>(),
           py::arg("tolerance") = 0.001, py::arg("max_cycles") = 50,
           py::arg("type") = Optimiser_simplex::TYPE::NORMAL)
      .def("__call__", &Optimiser_simplex::operator(), py::arg("target_fn"),
           py::arg("args"))
      .def("debug", &Optimiser_simplex::debug)
      .def("__repr__", [](const Optimiser_simplex &self) {
        return "<buccaneer.Optimiser_simplex class.>";
      });
}

void init_simplex_lib(py::module &m) {
  declare_target_functions(m);
  declare_optimiser_simplex(m);
}