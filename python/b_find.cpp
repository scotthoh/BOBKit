// Nanobind bindings for buccaneer-find
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "buccaneer/buccaneer-find.h"
#include "commons.h"
#include <nanobind/stl/vector.h>


using namespace clipper;

void declare_search_result(nb::module_ &m) {
  nb::class_<SearchResult>(m, "SearchResult")
      .def(nb::init<>())
      .def("__init__", [](SearchResult *sr, const ftype32 &score, const int &rot, const int &trn) {
             // SearchResult *result{score, rot, trn} = ;
             new ( sr ) SearchResult({score, rot, trn});
             //return std::unique_ptr<SearchResult>(
             //    new SearchResult({score, rot, trn}));
           },
           nb::arg("score"), nb::arg("rot_ind"), nb::arg("trn_ind"))
      .def_rw("score", &SearchResult::score)
      .def_rw("rot", &SearchResult::rot)
      .def_rw("trn", &SearchResult::trn)
      .def("__lt__", &SearchResult::operator<, nb::is_operator())
      .def("__str__",
           [](const SearchResult self) {
             return (String(self.score, 6, 4) + "," + String(self.rot) + "," +
                     String(self.trn));
           })
      .def("__repr__",
           [](const SearchResult self) {
             return "<buccaneer.SearchResult {Score = " +
                    String(self.score, 6, 4) +
                    ", Rot_index = " + String(self.rot) +
                    ", Trn_index = " + String(self.trn) + "}>";
           })
      .doc() = "Results class.";
}

void declare_ca_find(nb::module_ &m) {
  nb::class_<Ca_find> ca_find(m, "Ca_find",
                              "Class for finding Ca's from density.");

  nb::enum_<Ca_find::TYPE>( ca_find, "TYPE", "Find methods." )
      .value( "LIKELIHOOD", Ca_find::TYPE::LIKELIHOOD )
      .value( "SECSTRUC", Ca_find::TYPE::SECSTRUC )
      //.value( "CENTROIDS", Ca_find::TYPE::CENTROIDS )
      .export_values();

  ca_find.def( nb::init<int, double>(), nb::arg( "n_find" ) = 500, nb::arg( "resol" ) = 1.0 )
      .def( "__call__", &Ca_find::operator(), nb::arg( "mol" ), nb::arg( "knownstruc" ), nb::arg( "xmap" ),
            nb::arg( "llktarget" ), nb::arg( "type" ) = Ca_find::TYPE::LIKELIHOOD, nb::arg( "modelindex" ) = 0,
            "Find Ca using density." )
      .def_static(
          "find",
          []( MiniMol &mol, const KnownStructure &knownstruc, const clipper::Xmap<float> &xmap,
              const LLK_map_target &llktarget, const Ca_find::TYPE &type, int &modelind, int &nfind, double &resol, int &cpus,
              const nb::object &pystream ) -> bool {
            Ca_find cafind( nfind );
            cafind.set_cpus( cpus );
            bool success = cafind( mol, knownstruc, xmap, llktarget, type, modelind );
            std::string m = " C-alphas after finding    : " + clipper::String( int( mol.select( "*/*/CA" ).atom_list().size() ), 7 ) + "\n";
            if ( pystream.is_valid() ) to_pystream(m, pystream);
            else  std::cout << m;
            return success;
          },
          nb::arg( "mol" ), nb::arg( "knownstruc" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
          nb::arg( "type" ) = Ca_find::TYPE::LIKELIHOOD, nb::arg( "modelindex" ) = 0, nb::arg( "nfind" ) = 500,
          nb::arg( "resol" ) = 1.0, nb::arg( "ncpus" ) = 1, nb::arg( "stdout" ) = nullptr,
          "Static function to find Ca using density, with an option to print summary." )
      //.def(
      //    "__call__",
      //    []( Ca_find& self, clipper::MiniMol& mol, const KnownStructure& knownstruc,
      //        const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget,
      //        const Ca_find::TYPE type, const int modelindex ) {
      //      return self( mol, knownstruc, xmap, llktarget, type, modelindex );
      //    },
      //    nb::arg( "mol" ), nb::arg( "knownstruc" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
      //    nb::arg( "type" ) = Ca_find::TYPE::LIKELIHOOD, nb::arg( "modelindex" ) = 0,
      //    "Find Ca using density." )
      //.def(
      //    "__call__",
      //    []( Ca_find& self, clipper::MiniMol& mol, const KnownStructure& knownstruc,
      //        const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget,
      //        const std::vector<Coord_orth>& aa_instance, const Ca_find::TYPE type,
      //        const int modelindex ) {
      //      return self( mol, knownstruc, xmap, llktarget, aa_instance, type, modelindex );
      //    },
      //    nb::arg( "mol" ), nb::arg( "knownstruc" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
      //    nb::arg( "centroids" ), nb::arg( "type" ) = Ca_find::TYPE::LIKELIHOOD,
      //    nb::arg( "modelindex" ) = 0, "Find Ca using centroids and density." )
      .def_static( "set_cpus", &Ca_find::set_cpus, nb::arg( "ncpus" ), "Set number of cpu threads to use." )
      .def( "set_starting_instance_coords", &Ca_find::set_starting_instance_coords,
            //( void( Ca_find::* )( const std::vector<clipper::Coord_orth>&, const Xmap<float>&,
            //      const LLK_map_target& llktgt, const Ca_find::TYPE type ) ) &
            //    Ca_find::set_starting_instance_coords,
            nb::arg( "aa_instance" ), nb::arg( "xmap" ), nb::arg( "llktarget" ),
            nb::arg( "type" ) = Ca_find::TYPE::LIKELIHOOD, nb::arg( "refine_coords" ) = false,
            "Set starting instance coordinates from a list of orthogonal coordinates of amino acid "
            "instances." )
      .def( "get_initial_results", [](Ca_find &self, const clipper::Xmap<float> &xmap) {
        std::vector<clipper::Coord_orth> co;
        self.get_initial_results(co, xmap);
        return co;
      }
      )
      .def( "__repr__", []( const Ca_find &self ) { return "<buccaneer.Ca_find class>"; } );
}

void declare_search_threaded(nb::module_ &m) {
  nb::class_<Search_threaded>(m, "Search_threaded",
                              "Class with thread methods to search Ca groups.")
      .def(nb::init<>())
      .def(nb::init<const Xmap<int> &, const FFFear_fft<float> &,
                    const LLK_map_target &, const std::vector<RTop_orth> &,
                    const int>(),
           nb::arg("xlookp1"), nb::arg("srch"), nb::arg("llktarget"),
           nb::arg("RToperators"), nb::arg("lresult"))
      .def("set_range", &Search_threaded::set_range, nb::arg("n1"),
           nb::arg("n2"), "Set search range.")
      .def("search", &Search_threaded::search, nb::arg("op"),
           "Search Ca groups.")
      .def_prop_ro("results", &Search_threaded::results,
                             "Return search results.")
      .def("__call__", &Search_threaded::operator(), nb::arg("nthread") = 0,
           "Run single or multi-threaded.")
      .def("merge", &Search_threaded::merge, nb::arg("other"),
           "Merge results from multiple threads.")
      .def("__repr__",
           [](const Search_threaded &self) {
             return "<buccaneer.Search_threaded class>";
           })
      // inherited function/property
      .def_prop_ro("id", &Search_threaded::id, "Return thread id.");
}

void declare_ssfind(nb::module_ &m) {
  nb::class_<SSfind> ssfind(
      m, "SSfind",
      "Class for fast secondary structure finding (alternative to fffear).");

  nb::enum_<SSfind::SSTYPE>(ssfind, "TYPE", "Secondary structure type.")
      .value("ALPHA2", SSfind::SSTYPE::ALPHA2)
      .value("ALPHA3", SSfind::SSTYPE::ALPHA3)
      .value("ALPHA4", SSfind::SSTYPE::ALPHA4)
      .value("BETA2", SSfind::SSTYPE::BETA2)
      .value("BETA3", SSfind::SSTYPE::BETA3)
      .value("BETA4", SSfind::SSTYPE::BETA4)
      .export_values();

  ssfind.def( nb::init<>() )
      .def( "prep_xmap", &SSfind::prep_xmap, nb::arg( "xmap" ), nb::arg( "radius" ),
            "Prepare target map." )
      .def( "prep_search", nb::overload_cast<const Xmap<float>&> ( & SSfind::prep_search ),
          //( void( SSfind::* )( const Xmap<float>& ) ) & SSfind::prep_search,
            nb::arg( "xmap" ), "Prepare search with given map." )
      .def( "prep_search",
            nb::overload_cast<const Xmap<float>&, const double, const double, const Coord_orth>( &
                SSfind::prep_search ),
            nb::arg( "xmap" ), nb::arg( "rhocut" ), nb::arg( "radcut" ), nb::arg( "centre" ),
            "Prepare search with given map, density and radius cutoff, centre "
            "coordinates." )
      .def( "prep_search",
            nb::overload_cast<const Xmap<float>&, const std::vector<clipper::Coord_grid>&> ( &
                SSfind::prep_search ),
            nb::arg( "xmap" ), nb::arg( "centroids" ),
            "Prepare search with given map, density and radius cutoff, centre "
            "coordinates." )
      .def( "search", &SSfind::search, nb::arg( "target_coords" ), nb::arg( "op" ),
            nb::arg( "rhocut" ), nb::arg( "frccut" ) = 0.0, "Search secondary structure elements." )
      .def( "__repr__", []( const SSfind& self ) {
        return "<buccaneer.SSfind class>";
      } );

  using Class = SSfind::Target;
  nb::class_<Class> target(ssfind, "Target",
                           "Class to hold target coordinates.");
  target
      .def(nb::init<SSfind::SSTYPE, int>(), nb::arg("type"),
           nb::arg("num_residues"),
           "Constructor with secondary structure type and number of residues.")
      .def_prop_ro("target_coords", &Class::target_coords,
                             "Return list of target coordinates pairs")
      .def_prop_ro("calpha_coords", &Class::calpha_coords,
                             "Return list of C-alpha coordinates")
      .def("__repr__", [](Class &self) {
        return "<buccaneer.SSfind.Target with backbone coordinates for " + String(int(self.calpha_coords().size())) + " residues.>";
      });
}

void declare_search_op_aa_instance( nb::module_ &m ){
  nb::class_<Search_op_aa_instance_threaded>(m, "Search_op_aa_instance_threaded",
                              "Class for searching RTop for Ca positions from centroids.")
      .def(nb::init<>())
      .def(nb::init<const Xmap<float> &, const std::vector<Coord_grid>&,
                    const FFFear_fft<float> &, const LLK_map_target &,
                    const std::vector<RTop_orth> &, const int>(),
           nb::arg("xmap"), nb::arg("aa_instances"), nb::arg("srch"), nb::arg("llktarget"),
           nb::arg("RToperators"), nb::arg("lresult"))
      .def("set_range", &Search_op_aa_instance_threaded::set_range, nb::arg("n1"),
           nb::arg("n2"), "Set search range.")
      .def("search_op", &Search_op_aa_instance_threaded::search_op, nb::arg("op"),
           "Search RTop for amino acid instances.")
      .def_prop_ro("results", &Search_op_aa_instance_threaded::results,
                             "Return search results.")
      .def("__call__", &Search_op_aa_instance_threaded::operator(), nb::arg("nthread") = 0,
           "Run single or multi-threaded.")
      .def("merge", &Search_op_aa_instance_threaded::merge, nb::arg("other"),
           "Merge results from multiple threads.")
      .def("__repr__",
           [](const Search_op_aa_instance_threaded &self) {
             return "<buccaneer.Search_op_aa_instance_threaded class>";
           })
      // inherited function/property
      .def_prop_ro("id", &Search_op_aa_instance_threaded::id, "Return thread id.");
}

// Target_fn_refine_llk_map_target defined in b_simplex.cpp
// to be within same scope as Target_fn_zero_order trampoline definition

void add_ca_find(nb::module_ &m) {
  declare_search_result(m);
  declare_ca_find(m);
  declare_search_threaded(m);
  declare_ssfind(m); // weird Target(ALPHA2, 4) results.
  declare_search_op_aa_instance(m);
}