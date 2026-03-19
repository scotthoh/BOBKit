// Nanobind bindings for clipper minimol io
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include "c_minimol_io.h"
// #include <clipper/clipper-minimol.h>
#include <clipper/clipper-minimol-gemmi.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>

using namespace clipper;

void declare_minimol_io( nb::module_ &m ) {
  nb::class_<QuickWriteOptions>( m, "QuickWriteOptions" )
      .def( nb::init<>() )
      .def_rw( "minimal", &QuickWriteOptions::Minimal, "Disable many records (HEADER, TITLE, ...)" )
      //.def_rw( "only_headers", &QuickWriteOptions::only_headers, "Write only headers" )
      .def_rw( "pdb_short_ter", &QuickWriteOptions::PdbShortTer, "Write PDB TER records without numbers" )
      .def_rw( "pdb_linkr", &QuickWriteOptions::PdbLinkR, "Use non-standard Refmac LINKR records instead of LINK" )
      .def_rw( "pdb_short_chain_names", &QuickWriteOptions::PdbShortChainNames, "Shorten chain names for chain name size" )
      .def_rw( "cif_all_auth", &QuickWriteOptions::CifAllAuth,
               "Include _atom_site.auth_atom_id and auth_comp_id in mmcif" )
      .def_rw( "update_cif_doc", &QuickWriteOptions::UpdateCifDoc, "Update existing mmcif blocks from read in document" );

  nb::class_<GEMMIfile> gfile( m, "GEMMIfile" );
  nb::enum_<GEMMIfile::TYPE>( gfile, "TYPE" )
      .value( "default", GEMMIfile::Default )
      .value( "pdb", GEMMIfile::PDB )
      .value( "cif", GEMMIfile::CIF )
      //.value( "Mmjson", GEMMIfile::Mmjson )
      .export_values();

  gfile.def( nb::init<>() )
      .def(
          "read_file",
          []( GEMMIfile &self, const std::string fname, bool force_label, int max_line_length, bool split_chain_on_ter,
              bool skip_remarks ) {
            gemmi::PdbReadOptions opt{ max_line_length, split_chain_on_ter, skip_remarks };
            self.read_file( fname, opt, force_label );
          },
          nb::arg( "filename" ), nb::arg( "force_label" ) = false, nb::arg( "max_line_length" ) = 0,
          nb::arg( "split_chain_on_ter" ) = false, nb::arg( "skip_remarks" ) = false, "Read structure file." )
      .def(
          "write_pdb",
          []( GEMMIfile &self, const std::string fname, QuickWriteOptions &opts ) {
            self.write_file( fname, GEMMIfile::TYPE::PDB, opts );
          },
          nb::arg( "filename" ), nb::arg( "options" ).sig( "clipper.QuickWriteOptions()" ) = QuickWriteOptions(),
          "Write PDB file. To set more write options use GEMMIfile.set_pdb_write_options()" )
      .def(
          "write_cif",
          []( GEMMIfile &self, const std::string fname, QuickWriteOptions &opts ) {
            self.write_file( fname, GEMMIfile::TYPE::CIF, opts );
          },
          nb::arg( "filename" ), nb::arg( "options" ).sig( "clipper.QuickWriteOptions()" ) = QuickWriteOptions(),
          "Write CIF file. To set mmcif output groups, use set_mmcif_output_groups()" )
      .def( "import_minimol", &GEMMIfile::import_minimol, nb::arg( "minimol" ), nb::arg( "model_num" ) = 1,
            "Import MiniMol from GEMMI Structure." )
      .def( "export_minimol", &GEMMIfile::export_minimol, nb::arg( "minimol" ), "Export MiniMol to GEMMI Structure." )
      .def_static( "to_minimol", [](const gemmi::Structure &s, MiniMol &mol) {
            GEMMIfile gfile;
            gfile.set_gemmi_structure(s);
            gfile.import_minimol(mol);
            if (!mol.is_null())
              return true;
            else
              return false;
          },
          "Static function to import structure as MiniMol.")
      // clang-format off
      GETSET( gemmi_structure, GEMMIfile, GemmiStructure, get_gemmi_structure, set_gemmi_structure,
              "Get structure held in GEMMIfile", "Set structure if not read from file" )
      GETSET( spacegroup, GEMMIfile, Spacegroup, spacegroup, set_spacegroup, "Get spacegroup in CLIPPER format",
              "Set spacegroup" )
      GETSET( cell, GEMMIfile, Cell, cell, set_cell, "Get cell in CLIPPER format", "Set cell" )
      // clang-format on
      //.def_prop_rw(
      //    "gemmi_structure", &GEMMIfile::get_gemmi_structure,
      //    []( GEMMIfile &self, const GemmiStructure &st ) { self.set_gemmi_structure( st ); },
      //    nb::for_getter( "Get structure held in GEMMIfile." ), nb::for_setter( nb::arg( "structure " ) ),
      //    nb::for_setter( "Set structure if not read from file." ) )
      // might need to put all args in the []()
      .def(
          "set_mmcif_output_groups",
          []( GEMMIfile &self, bool all, const nb::kwargs &kwargs ) {
            gemmi::MmcifOutputGroups opts( all );
            // Expand call per field
            for (auto item : kwargs) {
              bool matched = false;
#define X(field) SETOPTS(opts, item, field, matched);
              GROUPS_FIELDS( X )
              #undef X
              if (!matched) {
                std::string _msg = "set_mmcif_output_groups: Unknown keyword argument '" + nb::cast<std::string>(item.first) + "'";
                clipper::Message::message(clipper::Message_warn(_msg.c_str()));
                //PyErr_WarnEx(PyExc_UserWarning, _msg.c_str(), 1);
              }
            }
            //#define X( field ) SETOPTS( opts, kwargs, field );
            //GROUPS_FIELDS( X )
            //#undef X
            self.set_mmcif_output_groups( opts );
          },
          nb::arg( "all" ), nb::arg( "kwargs" ),
          nb::sig( "def set_mmcif_output_groups(self, all: bool, atoms: bool, block_name: bool, entry: bool, "
                   "database_status: bool, author: bool, cell: bool, symmetry: bool, "
                   "entity: bool, entity_poly: bool, struct_ref: bool, chem_comp: bool, "
                   "exptl: bool, diffrn: bool, reflns: bool, refine: bool, "
                   "title_keywords: bool, ncs: bool, struct_asym: bool, origx: bool, "
                   "struct_conf: bool, struct_sheet: bool, struct_biol: bool, assembly: bool, "
                   "conn: bool, cis: bool, modres: bool, scale: bool, atom_type: bool, "
                   "entity_poly_seq: bool, tls: bool, software: bool, group_pdb: bool, "
                   "auth_all: bool /) -> None" ),
          // nb::kw_only(), nb::arg( "atoms" ), nb::arg( "block_name" ), nb::arg( "entry" ),
          // nb::arg( "database_status" ), nb::arg( "author" ), nb::arg( "cell" ), nb::arg( "symmetry" ),
          // nb::arg( "entity" ), nb::arg( "entity_poly" ), nb::arg( "struct_ref" ), nb::arg( "chem_comp" ),
          // nb::arg( "exptl" ), nb::arg( "diffrn" ), nb::arg( "reflns" ), nb::arg( "refine" ),
          // nb::arg( "title_keywords" ), nb::arg( "ncs" ), nb::arg( "struct_asym" ), nb::arg( "origx" ),
          // nb::arg( "struct_conf" ), nb::arg( "struct_sheet" ), nb::arg( "struct_biol" ), nb::arg( "assembly" ),
          // nb::arg( "conn" ), nb::arg( "cis" ), nb::arg( "modres" ), nb::arg( "scale" ), nb::arg( "atom_type" ),
          // nb::arg( "entity_poly_seq" ), nb::arg( "tls" ), nb::arg( "software" ), nb::arg( "group_pdb" ),
          // nb::arg( "auth_all" ),
          "Set MMcifOutputGroups to control what will be included in the control. Default is all. Refer "
          "MmcifOutputGroups in GEMMI" )
      .def(
          "set_pdb_write_options",
          []( GEMMIfile &self, const QuickWriteOptions &quick_opts, const nb::kwargs &kwargs ) {
            gemmi::PdbWriteOptions opts;  
            if ( quick_opts.Minimal )
              opts = gemmi::PdbWriteOptions::minimal();
            //else if ( quick_opts.only_headers )
            //  opts = gemmi::PdbWriteOptions::headers_only();
            // EXPAND call per field
            for (auto item : kwargs) {
              bool matched = false;
              #define X(field) SETOPTS(opts, item, field, matched);
              PDB_WRITE_OPTS( X )
              #undef X
              if (!matched) {
                std::string _msg = "set_pdb_write_options: Unknown keyword argument '" + nb::cast<std::string>(item.first) + "'";
                clipper::Message::message(clipper::Message_warn(_msg.c_str()));
                //PyErr_WarnEx(PyExc_UserWarning, _msg.c_str(), 1);
              }
            }
            //#define X( field ) SETOPTS( opts, kwargs, field );
            //PDB_WRITE_OPTS( X )
            //#undef X
            self.set_pdb_write_options( opts );
          },
          nb::arg( "quick_opts" ).sig( "clipper.QuickWriteOptions()" ) = QuickWriteOptions(), nb::arg( "kwargs" ),
          nb::sig( "def set_pdb_write_options(self, quick_opts: QuickWriteOptsions = QuickWriteOptions(), "
                   "minimal_file: bool, atom_records: bool, seqres_records: bool, ssbond_records: bool, "
                   "link_records: bool, cispep_records: bool, cryst1_record: bool, ter_records: bool, "
                   "conect_records: bool, end_record: bool, numbered_ter: bool, ter_ignores_type: bool, "
                   "use_linkr: bool, preserve_serial: bool /) -> None" ),
          // nb::kw_only(), nb::arg( "minimal_file" ) = false,
          // nb::arg( "atom_records" ) = true, nb::arg( "seqres_records" ) = true, nb::arg( "ssbond_records" ) = true,
          // nb::arg( "link_records" ) = true, nb::arg( "cispep_records" ) = true, nb::arg( "cryst1_record" ) = true,
          // nb::arg( "ter_records" ) = true, nb::arg( "conect_records" ) = false, nb::arg( "end_record" ) = true,
          // nb::arg( "numbered_ter" ) = true, nb::arg( "ter_ignores_type" ) = false, nb::arg( "use_linkr" ) = false,
          // nb::arg( "preserve_serial" ) = false,
          "Set PDB write options. Refer PdbWriteOptions in GEMMI" )
      .def(
          "set_cif_style",
          []( GEMMIfile &self, const bool &pairs, const bool &compact, const bool &mhash, const int &align_pairs,
              const int &align_loops ) {
            gemmi::cif::WriteOptions opts;
            opts.prefer_pairs = pairs;
            opts.compact = compact;
            opts.misuse_hash = mhash;
            opts.align_pairs = align_pairs;
            opts.align_loops = align_loops;
            self.set_cif_style( opts );
          },
          nb::arg( "prefer_pairs" ) = false, nb::arg( "compact" ) = false, nb::arg( "misuse_hash" ) = false,
          nb::arg( "align_pair" ) = 0, nb::arg( "align_loops" ) = 0,
          "Set cif output options. Refer CIF::WriteOptions in GEMMI" );
}

// void declare_gemmi_structure( nb::module_ &m ) {
//   nb::class_<GemmiStructure>( m, "GemmiStructure" )
//       .def( nb::init<>(), "Null constructor" )
//       .def( nb::init<const gemmi::Structure &>(), nb::arg( "structure" ), "Constructor from GEMMI structure" )
//       .def_static(
//           "get_id_str",
//           []( const String &model_num, const gemmi::CRA &cra, const std::string &entity ) {
//             return GemmiStructure::GetID_str( model_num, cra, entity );
//           },
//           nb::arg( "modelnum" ), nb::arg( "cra" ), nb::arg( "entity" ),
//           "Return a string ID for atom, residue, chain or model. e.g. /1/C/LYS.ins/NZ[N]:altloc" )
//       .def_static(
//           "get_id_str",
//           []( const int &num, const gemmi::CRA &cra, const std::string &entity ) {
//             return GemmiStructure::GetID_str( String( num ), cra, entity );
//           },
//           nb::arg( "modelnum" ), nb::arg( "cra" ), nb::arg( "entity" ),
//           "Return a string ID for atom, residue, chain or model. e.g. /1/C/LYS.ins/NZ[N]:altloc" )
//       // clang-format off
//       GETSET( spacegroup, GemmiStructure, Spacegroup, spacegroup, set_spacegroup,
//               "Get spacegroup in CLIPPER format", "Set spacegroup" )
//       GETSET( cell, GemmiStructure, Cell, get_cell, set_cell, "Get cell in CLIPPER format", "Set cell" )
//       GETSET( resolution, GemmiStructure, Resolution, get_resolution, set_resolution,
//               "Get resolution in CLIPPER format", "Set resolution" )
//       GETSET( transform, GemmiStructure, RTop_orth, get_transform, set_transform,
//               "Get transformation/orgix matrix in CLIPPER format", "Set transform" )
//       // clang-format on
//       .def( "get_exptlmethod", &GemmiStructure::get_exptlmethod, "Get experimental method" );
// 
//   nb::class_<GemmiAtom_list, Atom_list>( m, "GemmiAtom_list" )
//       .def( nb::init<gemmi::CraProxy>(), "Constructor from GEMMI CraProxy" );
// }

void add_minimol_io_gemmi( nb::module_ &m ) {
  declare_minimol_io( m );
  //declare_gemmi_structure( m );
}