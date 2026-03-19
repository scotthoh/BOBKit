// Nanobind bindings for buccaneer-lib
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-lib.h"
#include "buccaneer/buccaneer-grow.h"
#include "commons.h"
#include "arrays.h"
#include <clipper/core/coords.h>

using namespace clipper;

template <class T>
void declare_score_list(nb::module_ &m, const std::string name) {
  using SLClass = Score_list<T>;
  auto PyClass = std::string("Score_list_" + name);
  nb::class_<SLClass> score_list(m, PyClass.c_str());
  score_list.def(nb::init<>())
      .def(nb::init<const int &>(), nb::arg("nmax"),
           "Constructor taking maximum size of list.")
      .def("init", &SLClass::init, nb::arg("nmax"),
           "Initialiser taking maximum size of list")
      .def("addable", &SLClass::addable, nb::arg("score"),
           "Is given score good enough to add to list?")
      .def("add", &SLClass::add, nb::arg("score"), nb::arg("data"),
           "Add a score to list, if it is good enough.")
      .def("delete", &SLClass::del, nb::arg("pos"), "Delete score from list.")
      .def("score", &SLClass::score, nb::arg("pos"), "Access score.")
      .def(
          "__getitem__",
          [](const SLClass &self, const int &i) {
            return self[normalise_index(i, self.size())];
          },
          "Access list.")
      .def("size", &SLClass::size, "List size.")
      .def("__len__", &SLClass::size)
      .def("__repr__",
           [](const SLClass &self) {
             return "<Buccaneer.Score_list with " +
                    clipper::String(self.size()) + " entries.>";
           })
      .doc() =
      "Store list class. \nThe class keeps a list of the best fits to some target function, \
      scoring both the score and an object of a user defined type representing the \
      parameters whuich gave that score.On construction, a maximum size is specified. \
      Only the best N matches will be stored.";
}

void declare_LLK_map_target(nb::module_ &m) {
  nb::class_<LLK_map_target> llkmaptgt(m, "LLK_map_target");

  nb::enum_<LLK_map_target::TYPE>(llkmaptgt, "TYPE")
      .value("NORMAL", LLK_map_target::TYPE::NORMAL)
      .value("CORREL", LLK_map_target::TYPE::CORREL);

  llkmaptgt.def(nb::init<>())
      .def(nb::init<const ftype &, const ftype &, LLK_map_target::TYPE>(),
           nb::arg("rad"), nb::arg("sampling"),
           nb::arg("type") = LLK_map_target::TYPE::NORMAL,
           "Constructor: provide radius and sampling in A for LLK target.")
      .def("init", &LLK_map_target::init, nb::arg("rad"), nb::arg("sampling"),
           nb::arg("type") = LLK_map_target::TYPE::NORMAL,
           "Initialiser: provide radius and sampling in A for LLK target.")
      .def("prep_llk", &LLK_map_target::prep_llk,
           "Prepare LLK target after accumulating density or loading map.")
      .def("accumulate", &LLK_map_target::accumulate, nb::arg("xmap"),
           nb::arg("rtop"),
           "Accumulate density statistics from a sample orientation in a "
           "sample map.")
      .def("set_scale", &LLK_map_target::set_scale, nb::arg("scale"),
           nb::arg("offset"), "Set the scaling for llk and llk_approx.")
      .def("rtop_list", &LLK_map_target::rtop_list, nb::arg("spacegroup"),
           nb::arg("step"), "Return a list of RTops.")
      .def("search", &LLK_map_target::search, nb::arg("resultscr"),
           nb::arg("resultrot"), nb::arg("resulttrn"), nb::arg("xmap"),
           nb::arg("rtops"),
           "Perform 6-d (RT) search for LLK target in given map with prior.")
      .def("llk_approx", &LLK_map_target::llk_approx, nb::arg("xmap"),
           nb::arg("rtop"),
           "Calculate fast approx to LLK for given orientation.")
      .def("llk", &LLK_map_target::llk, nb::arg("xmap"), nb::arg("rtop"),
           "Calculate full LLK for given orientation.")
      .def("format", &LLK_map_target::format,
           "Output formatted representation.")
      .def("__repr__",
           [](LLK_map_target &self) {
             return "<Buccaneer.LLK_map_target>";
           })
      .def("prep_llk_distribution", &LLK_map_target::prep_llk_distribution,
           nb::arg("xmap"),
           "Generate a distribution of LLK values for a given map.")
      .def("llk_distribution", &LLK_map_target::prep_llk_distribution,
           nb::arg("ordinal"),
           "Return llk value by position in cumulative distribution.")
      .def("sampled", &LLK_map_target::sampled,
                             "Access to fine sampled target function.")

      // maybe this is correct binding
      // auto x = static_cast< const Dummy& (Alpha::*)() const>(&Alpha::force);
      // auto y = static_cast< Dummy& (Alpha::*)()>(&Alpha::force);
      .def_prop_rw(
          "llk_target",
          nb::overload_cast<>(&LLK_map_target::llk_target, nb::const_),
          nb::overload_cast<>(&LLK_map_target::llk_target),
          //static_cast<const clipper::NXmap<float> &(LLK_map_target::*)() const>(
          //    &LLK_map_target::llk_target),
          //static_cast<clipper::NXmap<float> &(LLK_map_target::*)()>(
          //    &LLK_map_target::llk_target),
          "Access to LLK target for load/save")
      .def_prop_rw(
          "llk_weight",
          nb::overload_cast<>(&LLK_map_target::llk_weight, nb::const_),
          nb::overload_cast<>(&LLK_map_target::llk_weight),
          //static_cast<const clipper::NXmap<float> &(LLK_map_target::*)() const>(
          //    &LLK_map_target::llk_weight),
          //static_cast<clipper::NXmap<float> &(LLK_map_target::*)()>(
          //    &LLK_map_target::llk_weight),
          "Access to LLK weight for load/save.")
      .def_prop_ro("num_samples", &LLK_map_target::num_samples,
                             "Access to number of samples.")
      .doc() =
      "Log-likelihood map matching target.\n"
      "This class is used in determining the log likelihood fit of some "
      "desired density features from some region of a noisy electron density "
      "map. It contains methods to accumulate the log likelihood target from a "
      "number of sample density regions from a map with similar noise levels "
      "to the target map; methods for a FFFear 6-dimensional "
      "(rotation/orientation) search of the target map; and methods for "
      "testing individual sample positions and orientations.\n"
      "Note the results from the 6-d search and the fast and full LLK "
      "calculations are on different scales and so cannot be compared "
      "directly.";

  using SampledClass = LLK_map_target::Sampled;
  nb::class_<SampledClass>(llkmaptgt, "Sampled")
      .def(nb::init<>())
      .def("insert", &SampledClass::insert, nb::arg("coord"), nb::arg("tgt"),
           nb::arg("wgt"), "Insert values.")
      .def("llk", &SampledClass::llk, nb::arg("xmap"), nb::arg("rtop"),
           "Evaluate llk.")
      .def("correl", &SampledClass::correl, nb::arg("xmap"), nb::arg("rtop"),
           "Evaluate correlation.")
      .def("target", nb::overload_cast<const clipper::Xmap<float> &, const clipper::RTop_orth &>(&SampledClass::target, nb::const_),
           nb::arg("xmap"), nb::arg("rtop"), "Evaluate target function.")
      .def("size", &SampledClass::size, "Get number of samples.")
      .def("__len__", &SampledClass::size)
      .def("__repr__",
           [](const SampledClass &self) {
             return "<buccaneer.Sampled with " + clipper::String(self.size()) +
                    " sampled values.>";
           })
      .def("set_type", &SampledClass::set_type, nb::arg("type"), "Set type.")
      .def("coord_orth", &SampledClass::coord_orth, nb::arg("pos"),
           "Get target coordinates at given index.")
      .def("target", nb::overload_cast<int>(& SampledClass::target, nb::const_),
           nb::arg("pos"), "Get target at given index.")
      .def("weight", &SampledClass::weight, nb::arg("pos"),
           "Get weights at given index.")
      .doc() = "Class to hold sampled values from LLK map target.";
}

void add_buccaneer_lib(nb::module_ &m) {
  declare_score_list<clipper::String>(m, "String");
  declare_score_list<clipper::RTop_orth>(m, "RTop_orth");
  // another way for similar entries like Rama_ang in buccaneer-grow
  declare_score_list<std::vector<ftype>>(m, "FloatList");
  declare_LLK_map_target(m);
}