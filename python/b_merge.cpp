// Nanobind bindings for buccaneer-merge
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-merge.h"
#include "commons.h"
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>

using namespace clipper;

void add_ca_merge(nb::module_ &m) {

  nb::class_<Ca_merge> camerge(m, "Ca_merge");
  camerge.def(nb::init<double>(), nb::arg("reliability") = 0.5)
      .def("__call__", &Ca_merge::operator(), nb::arg("mol"), nb::arg("xmap"),
           nb::arg("llktgt"), nb::arg("seq"), "Merge model.")
      .def_static("merge_mr", &Ca_merge::merge_mr, nb::arg("mol"),
                  nb::arg("mol_mr"), nb::arg("sigcut"), nb::arg("nseed"),
                  nb::arg("mr_filter"), nb::arg("mr_seed"),
                  "Merge with molecular replacement model.")
      .def_static( "merge", [](clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq, double reliability, const nb::object &pystream ) -> bool {
        Ca_merge camerge(reliability);
        std::string m = " C-alphas before model merge : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
        if (pystream.is_valid()) to_pystream(m, pystream);   
        else std::cout << m;
        bool success = camerge(mol, xmap, llktarget, seq);
        m = " C-alphas after model merge  : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
        if (pystream.is_valid()) to_pystream(m, pystream);   
        else std::cout << m;
        return success;
        }, nb::arg("mol"), nb::arg("xmap"), nb::arg("llktgt"), nb::arg("seq"), nb::arg("reliability")=0.5, nb::arg("stdout")=nullptr,
        "Static function to merge model."  )
      .def("__repr__",
           [](const Ca_merge &self) { return "<buccaneer.Ca_merge class.>"; })
      .doc() = "Class for augmenting model with MERGE model.";
}