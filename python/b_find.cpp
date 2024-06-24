// Author: S.W.Hoh
// 2023 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-find.h"

#include "type_conversions.h"
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void declare_search_result(py::module &m) {
  py::class_<SearchResult>(m, "SearchResult")
      .def(py::init<>())
      .def(py::init([](const ftype32 &score, const int &rot, const int &trn) {
             // SearchResult *result{score, rot, trn} = ;
             return std::unique_ptr<SearchResult>(
                 new SearchResult({score, rot, trn}));
           }),
           py::arg("score"), py::arg("rot_ind"), py::arg("trn_ind"))
      .def_readwrite("score", &SearchResult::score)
      .def_readwrite("rot", &SearchResult::rot)
      .def_readwrite("trn", &SearchResult::trn)
      .def("__lt__", &SearchResult::operator<, py::is_operator())
      .def("__str__",
           [](const SearchResult self) {
             return (String(self.score, 6, 4) + "," + String(self.rot) + "," +
                     String(self.trn));
           })
      .def("__repr__", [](const SearchResult self) {
        return "<buccaneer.SearchResult {Score = " + String(self.score, 6, 4) +
               ", Rot_index = " + String(self.rot) +
               ", Trn_index = " + String(self.trn) + "}>";
      });
}

void declare_ca_find(py::module &m) {
  py::class_<Ca_find> ca_find(m, "Ca_find");

  py::enum_<Ca_find::TYPE>(ca_find, "TYPE")
      .value("LIKELIHOOD", Ca_find::TYPE::LIKELIHOOD)
      .value("SECSTRUC", Ca_find::TYPE::SECSTRUC);

  ca_find
      .def(py::init<int, double>(), py::arg("n_find") = 500,
           py::arg("resol") = 1.0)
      .def("__call__", &Ca_find::operator(), py::arg("mol"),
           py::arg("knownstruc"), py::arg("xmap"), py::arg("llktarget"),
           py::arg("type") = Ca_find::TYPE::LIKELIHOOD,
           py::arg("modelindex") = 0)
      .def_static("set_cpus", &Ca_find::set_cpus, py::arg("ncpus"))
      .def("__repr__", [](const Ca_find &self) {
        std::stringstream stream;
        stream << "<buccaneer.Ca_find class>";
        return stream.str();
      });
}

void declare_search_threaded(py::module &m) {
  py::class_<Search_threaded>(m, "Search_threaded")
      .def(py::init<>())
      .def(py::init<const Xmap<int> &, const FFFear_fft<float> &,
                    const LLK_map_target &, const std::vector<RTop_orth> &,
                    const int>(),
           py::arg("xlookp1"), py::arg("srch"), py::arg("llktarget"),
           py::arg("RToperators"), py::arg("lresult"))
      .def("set_range", &Search_threaded::set_range, py::arg("n1"),
           py::arg("n2"))
      .def("search", &Search_threaded::search, py::arg("op"))
      .def_property_readonly("results", &Search_threaded::results)
      .def("__call__", &Search_threaded::operator(), py::arg("nthread") = 0)
      .def("merge", &Search_threaded::merge, py::arg("other"))
      .def("__repr__",
           [](const Search_threaded &self) {
             std::stringstream stream;
             stream << "<buccaneer.Search_threaded class>";
             return stream.str();
           })
      // inherited function/property
      .def_property_readonly("id", &Search_threaded::id);
}

void declare_ssfind(py::module &m) {
  py::class_<SSfind> ssfind(m, "SSfind");

  py::enum_<SSfind::SSTYPE>(ssfind, "TYPE")
      .value("ALPHA2", SSfind::SSTYPE::ALPHA2)
      .value("ALPHA3", SSfind::SSTYPE::ALPHA3)
      .value("ALPHA4", SSfind::SSTYPE::ALPHA4)
      .value("BETA2", SSfind::SSTYPE::BETA2)
      .value("BETA3", SSfind::SSTYPE::BETA3)
      .value("BETA4", SSfind::SSTYPE::BETA4);

  ssfind.def(py::init<>())
      .def("prep_xmap", &SSfind::prep_xmap, py::arg("xmap"), py::arg("radius"))
      .def("prep_search",
           (void(SSfind::*)(const Xmap<float> &)) & SSfind::prep_search,
           py::arg("xmap"))
      .def("prep_search",
           (void(SSfind::*)(const Xmap<float> &, const double, const double,
                            const Coord_orth)) &
               SSfind::prep_search,
           py::arg("xmap"), py::arg("rhocut"), py::arg("radcut"),
           py::arg("centre"))
      .def("search", &SSfind::search, py::arg("target_coords"), py::arg("op"),
           py::arg("rhocut"), py::arg("frccut") = 0.0)
      .def("__repr__", [](const SSfind &self) {
        std::stringstream stream;
        stream << "<buccaneer.SSfind class>";
        return stream.str();
      });

  using Class = SSfind::Target;
  py::class_<Class> target(ssfind, "Target");
  target
      .def(py::init<SSfind::SSTYPE, int>(), py::arg("type"),
           py::arg("num_residues"))
      .def_property_readonly("target_coords", &Class::target_coords)
      .def_property_readonly("calpha_coords", &Class::calpha_coords)
      .def("__repr__", [](Class &self) {
        std::stringstream stream;
        stream << "<buccaneer.SSfind.Target with backbone coordinates for ";
        stream << self.calpha_coords().size() << " residues.>";
        return stream.str();
      });
}

// Target_fn_refine_llk_map_target defined in b_simplex.cpp
// to be within same scope as Target_fn_zero_order trampoline definition

void init_ca_find(py::module &m) {
  declare_search_result(m);
  declare_ca_find(m);
  declare_search_threaded(m);
  declare_ssfind(m); // weird Target(ALPHA2, 4) results.
                     // declare_target_fn_refine_llk_map_target(m);
}