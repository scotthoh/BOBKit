// Author: S.W.Hoh
// (C)2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-lib.h"

#include "helper_functions.h"
#include "type_conversions.h"
#include <clipper/core/coords.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// PYBIND11_MAKE_OPAQUE(std::vector<LLK_map_target>)

template <class T>
void declare_score_list(py::module &m, const std::string name) {
  using SLClass = Score_list<T>;
  auto PyClass = std::string("Score_list_" + name);
  py::class_<SLClass> score_list(m, PyClass.c_str());
  score_list.def(py::init<>())
      .def(py::init<const int &>(), py::arg("nmax"))
      .def("init", &SLClass::init, py::arg("nmax"))
      .def("addable", &SLClass::addable, py::arg("score"))
      .def("add", &SLClass::add, py::arg("score"), py::arg("data"))
      .def("del", &SLClass::del, py::arg("pos"))
      .def("score", &SLClass::score, py::arg("pos"))
      .def("__getitem__",
           [](const SLClass &self, const int &i) {
             return self[normalise_index(i, self.size())];
           })
      .def("size", &SLClass::size)
      .def("__len__", &SLClass::size)
      .def("__repr__", [](const SLClass &self) {
        return "<Buccaneer.Score_list with " + clipper::String(self.size()) +
               " entries.>";
      });
}

void declare_LLK_map_target(py::module &m) {
  py::class_<LLK_map_target> llkmaptgt(m, "LLK_map_target");

  py::enum_<LLK_map_target::TYPE>(llkmaptgt, "TYPE")
      .value("NORMAL", LLK_map_target::TYPE::NORMAL)
      .value("CORREL", LLK_map_target::TYPE::CORREL);

  llkmaptgt.def(py::init<>())
      .def(py::init<const ftype &, const ftype &, LLK_map_target::TYPE>(),
           py::arg("rad"), py::arg("sampling"),
           py::arg("type") = LLK_map_target::TYPE::NORMAL)
      .def("init", &LLK_map_target::init, py::arg("rad"), py::arg("sampling"),
           py::arg("type") = LLK_map_target::TYPE::NORMAL)
      .def("prep_llk", &LLK_map_target::prep_llk)
      .def("accumulate", &LLK_map_target::accumulate, py::arg("xmap"),
           py::arg("rtop"))
      .def("set_scale", &LLK_map_target::set_scale, py::arg("scale"),
           py::arg("offset"))
      .def("rtop_list", &LLK_map_target::rtop_list, py::arg("spacegroup"),
           py::arg("step"))
      .def("search", &LLK_map_target::search, py::arg("resultscr"),
           py::arg("resultrot"), py::arg("resulttrn"), py::arg("xmap"),
           py::arg("rtops"))
      .def("llk_approx", &LLK_map_target::llk_approx, py::arg("xmap"),
           py::arg("rtop"))
      .def("llk", &LLK_map_target::llk, py::arg("xmap"), py::arg("rtop"))
      .def("format", &LLK_map_target::format)
      .def("__repr__",
           [](LLK_map_target &self) {
             std::string msg = "<Buccaneer.LLK_map_target with ";
             if (!self.llk_target().is_null())
               return msg + clipper::String(self.num_samples()) +
                      " map samples.>";
             else
               return msg + "0 map samples.>";
           })
      .def("prep_llk_distribution", &LLK_map_target::prep_llk_distribution,
           py::arg("xmap"))
      .def("llk_distribution", &LLK_map_target::prep_llk_distribution,
           py::arg("ordinal"))
      .def_property_readonly("sampled", &LLK_map_target::sampled)

      // maybe this is correct binding
      // auto x = static_cast< const Dummy& (Alpha::*)() const>(&Alpha::force);
      // auto y = static_cast< Dummy& (Alpha::*)()>(&Alpha::force);
      .def_property(
          "llk_target",
          static_cast<const clipper::NXmap<float> &(LLK_map_target::*)() const>(
              &LLK_map_target::llk_target),
          static_cast<clipper::NXmap<float> &(LLK_map_target::*)()>(
              &LLK_map_target::llk_target))
      //.def_property(
      //    "llk_target",
      //    [](const LLK_map_target &self) -> const NXmap<float> & {
      //      return self.llk_target();
      //    },
      //    [](LLK_map_target &self, const NXmap<float> &target_map) {
      //      self.llk_target() = target_map;
      //    })
      .def_property(
          "llk_weight",
          static_cast<const clipper::NXmap<float> &(LLK_map_target::*)() const>(
              &LLK_map_target::llk_weight),
          static_cast<clipper::NXmap<float> &(LLK_map_target::*)()>(
              &LLK_map_target::llk_weight))
      .def_property_readonly("num_samples", &LLK_map_target::num_samples);
  // py::bind_vector<std::vector<LLK_map_target>>(m, "LLK_map_target_List");
}

void declare_sampled(py::module &m) {
  using SampledClass = LLK_map_target::Sampled;
  py::class_<SampledClass>(m, "Sampled")
      .def(py::init<>())
      .def("insert", &SampledClass::insert, py::arg("coord"), py::arg("tgt"),
           py::arg("wgt"))
      .def("llk", &SampledClass::llk, py::arg("xmap"), py::arg("rtop"))
      .def("correl", &SampledClass::correl, py::arg("xmap"), py::arg("rtop"))
      .def("target",
           (ftype(SampledClass::*)(const clipper::Xmap<float> &,
                                   const clipper::RTop_orth &) const) &
               SampledClass::target,
           py::arg("xmap"), py::arg("rtop"))
      .def("size", &SampledClass::size)
      .def("__len__", &SampledClass::size)
      .def("__repr__",
           [](const SampledClass &self) {
             return "<buccaneer.Sampled with " + clipper::String(self.size()) +
                    " sampled values.>";
           })
      .def("set_type", &SampledClass::set_type, py::arg("type"))
      .def("coord_orth", &SampledClass::coord_orth, py::arg("pos"))
      .def("target", (ftype(SampledClass::*)(int) const) & SampledClass::target,
           py::arg("pos"))
      .def("weight", &SampledClass::weight, py::arg("pos"));
}

void init_buccaneer_lib(py::module &m) {
  declare_score_list<clipper::String>(m, "String");
  declare_score_list<clipper::RTop_orth>(m, "RTop_orth");
  declare_score_list<py::array_t<ftype>>(m, "float_array");
  declare_LLK_map_target(m);
  declare_sampled(m);
}