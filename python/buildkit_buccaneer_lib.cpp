// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer-lib.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "helper_functions.h"

template <class T>
void declare_score_list(py::module &m, const std::string name)
{
  using SLClass = Score_list<T>;
  auto PyClass = std::string("Score_list_" + name);
  py::class_<SLClass> score_list(m, PyClass.c_str());
  score_list
      .def(py::init<>())
      .def(py::init<const int &>())
      .def("init", &SLClass::init)
      .def("addable", &SLClass::addable)
      .def("add", &SLClass::add)
      .def("del", &SLClass::del)
      .def("score", &SLClass::score)
      .def("__getitem__", [](const SLClass &self, const int &i)
           { return self[normalise_index(i, self.size())]; })
      .def("size", &SLClass::size)
      .def("__repr__", [](const SLClass &self)
           {  std::stringstream stream;
              stream << "<Buccaneer.Score_list with ";
              stream << self.size() << " entries.>";
              return stream.str(); })
      .def("__str__", [](const SLClass &self)
           {  std::stringstream stream;
              stream << "<Buccaneer.Score_list with ";
              stream << self.size() << " entries.>";
              return stream.str(); });
}

void declare_LLK_map_target(py::module &m)
{
  py::class_<LLK_map_target> llkmaptgt(m, "LLK_map_target");

  py::enum_<LLK_map_target::TYPE>(llkmaptgt, "TYPE")
      .value("NORMAL", LLK_map_target::TYPE::NORMAL)
      .value("CORREL", LLK_map_target::TYPE::CORREL);

  llkmaptgt
      .def(py::init<>())
      .def(py::init<const ftype &, const ftype &, LLK_map_target::TYPE>())
      .def("init", &LLK_map_target::init)
      .def("prep_llk", &LLK_map_target::prep_llk)
      .def("accumulate", &LLK_map_target::accumulate)
      .def("set_scale", &LLK_map_target::set_scale)
      .def("rtop_list", &LLK_map_target::rtop_list)
      .def("search", &LLK_map_target::search)
      .def("llk_approx", &LLK_map_target::llk_approx)
      .def("llk", &LLK_map_target::llk)
      .def("format", &LLK_map_target::format)

      .def("prep_llk_distribution", &LLK_map_target::prep_llk_distribution)
      .def("llk_distribution", &LLK_map_target::prep_llk_distribution)
      .def_property_readonly("sampled", &LLK_map_target::sampled)
      .def_property(
          "llk_target",
          [](const LLK_map_target &self) -> const NXmap<float> &
          { return self.llk_target(); },
          [](LLK_map_target &self, const NXmap<float> &target_map)
          { self.llk_target() = target_map; })
      .def_property(
          "llk_weight",
          [](const LLK_map_target &self) -> const NXmap<float> &
          { return self.llk_weight(); },
          [](LLK_map_target &self, const NXmap<float> &target_map)
          { self.llk_weight() = target_map; })
      .def_property_readonly("num_sampled", &LLK_map_target::num_samples);
}

void declare_sampled(py::module &m)
{
  using SampledClass = LLK_map_target::Sampled;
  py::class_<SampledClass>(m, "Sampled")
      .def(py::init<>())
      .def("insert", &SampledClass::insert)
      .def("llk", &SampledClass::llk)
      .def("correl", &SampledClass::correl)
      .def("target", (ftype(SampledClass::*)(const clipper::Xmap<float> &, const clipper::RTop_orth &) const) & SampledClass::target)
      .def("size", &SampledClass::size)
      .def("set_type", &SampledClass::set_type)
      .def("coord_orth", &SampledClass::coord_orth)
      .def("target", (ftype(SampledClass::*)(int) const) & SampledClass::target)
      .def("weight", &SampledClass::weight);
}

void init_buccaneer_lib(py::module &m)
{
  declare_score_list<ftype64>(m, "double");
  declare_LLK_map_target(m);
  declare_sampled(m);
}