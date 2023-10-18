// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-find.h"

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "type_conversions.h"

namespace py = pybind11;
using namespace clipper;

void declare_search_result(py::module &m)
{
  py::class_<SearchResult>(m, "SearchResult")
      .def(py::init<>())
      .def(py::init([](const ftype32 &score, const int &rot, const int &trn)
                    {
                    //SearchResult *result{score, rot, trn} = ;
                    return std::unique_ptr<SearchResult>(new SearchResult({score, rot, trn})); }),
           py::arg("score"), py::arg("rot"), py::arg("trn"))
      .def_readwrite("score", &SearchResult::score)
      .def_readwrite("rot", &SearchResult::rot)
      .def_readwrite("trn", &SearchResult::trn)
      .def(py::self < py::self)
      .def("__str__", [](const SearchResult self)
           { return String(self.score, 6, 4); })
      .def("__repr__", [](const SearchResult self)
           { return "{Score = " + String(self.score, 6, 4) + ", Rot = " + String(self.rot) + ", Trn = " + String(self.trn) + "}"; });
}

void declare_ca_find(py::module &m)
{
  py::class_<Ca_find> ca_find(m, "Ca_find");

  py::enum_<Ca_find::TYPE>(ca_find, "TYPE")
      .value("LIKELIHOOD", Ca_find::TYPE::LIKELIHOOD)
      .value("SECSTRUC", Ca_find::TYPE::SECSTRUC);

  ca_find
      .def(py::init<int, double>(), py::arg("n_find") = 500, py::arg("resol") = 1.0)
      .def("__call__", &Ca_find::operator(),
           py::arg("mol"), py::arg("knownstruc"), py::arg("xmap"),
           py::arg("llktarget"), py::arg("type") = Ca_find::TYPE::LIKELIHOOD, py::arg("modelindex") = 0)
      .def_static("set_cpus", &Ca_find::set_cpus);
}

void declare_search_threaded(py::module &m)
{
  py::class_<Search_threaded>(m, "Search_threaded")
      .def(py::init<>())
      .def(py::init<const Xmap<int> &, const FFFear_fft<float> &, const LLK_map_target &, const std::vector<RTop_orth> &, const int>())
      .def("set_range", &Search_threaded::set_range)
      .def("search", &Search_threaded::search)
      .def_property_readonly("results", &Search_threaded::results)
      .def("__call__", &Search_threaded::operator())
      .def("merge", &Search_threaded::merge);
}

void declare_ssfind(py::module &m)
{
  py::class_<SSfind> ssfind(m, "SSfind");

  py::enum_<SSfind::SSTYPE>(ssfind, "TYPE")
      .value("ALPHA2", SSfind::SSTYPE::ALPHA2)
      .value("ALPHA3", SSfind::SSTYPE::ALPHA3)
      .value("ALPHA4", SSfind::SSTYPE::ALPHA4)
      .value("BETA2", SSfind::SSTYPE::BETA2)
      .value("BETA3", SSfind::SSTYPE::BETA3)
      .value("BETA4", SSfind::SSTYPE::BETA4);

  ssfind
      .def(py::init<>())
      .def("prep_xmap", &SSfind::prep_xmap)
      .def("prep_search", (void(SSfind::*)(const Xmap<float> &)) & SSfind::prep_search)
      .def("prep_search", (void(SSfind::*)(const Xmap<float> &, const double, const double, const Coord_orth)) & SSfind::prep_search)
      .def("search", &SSfind::search);

  using Class = SSfind::Target;
  py::class_<Class> target(m, "Target");
  target
      .def(py::init<SSfind::SSTYPE, int>())
      .def_property_readonly("target_coords", &Class::target_coords)
      .def_property_readonly("calpha_coords", &Class::calpha_coords);
}

void declare_target_fn_refine_llk_map_target(py::module &m)
{
  py::class_<Target_fn_refine_llk_map_target>(m, "Target_fn_refine_llk_map_target")
      .def(py::init<>())
      .def(py::init<const Xmap<float> &, const LLK_map_target &, const double &, const double &>())
      .def_property_readonly("num_params", &Target_fn_refine_llk_map_target::num_params)
      .def("__call__", (double(Target_fn_refine_llk_map_target::*)(const RTop_orth &) const) & Target_fn_refine_llk_map_target::operator())
      .def("__call__", (double(Target_fn_refine_llk_map_target::*)(const std::vector<double> &) const) & Target_fn_refine_llk_map_target::operator())
      .def("rtop_orth", &Target_fn_refine_llk_map_target::rtop_orth)
      .def("refine", &Target_fn_refine_llk_map_target::refine);
}

void init_ca_find(py::module &m)
{
  declare_search_result(m);
  declare_ca_find(m);
  declare_search_threaded(m);
  declare_ssfind(m); // weird Target(ALPHA2, 4) results.
  declare_target_fn_refine_llk_map_target(m);
}