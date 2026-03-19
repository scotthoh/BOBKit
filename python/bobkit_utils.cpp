// Nanobind bindings of some useful methods
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/minimol/minimol_seq.h>
#include "commons.h"
#include "arrays.h"
#include "buccaneer/buccaneer-util.h"

using namespace clipper;

void add_utils( nb::module_ & m ) {
  // nb::class_<BuccaneerUtil>( m, "Util", "Utility functions." )
  //.def_static("set_reference", &BuccaneerUtil::set_reference,
  //             nb::arg("mtz"), nb::arg("pdb"),
  //             "Set reference MTZ and PDB, only work when CCP4's CLIBD "
  //             "environment path is set.")
  //.def_static("read_model", &BuccaneerUtil::read_model, nb::arg("mol"),
  //             nb::arg("file"), nb::arg("verbose") = true,
  //             "Read model and export to minimol.")
  m.def(
      "read_structure",
      []( const std::string& fpath, bool enable_messages ) {
        if ( fpath == "undefined" ) {
          throw std::invalid_argument( "No path/filename provided for input model! Aborting..." );
        }

        clipper::MiniMol mmol;
        BuccaneerUtil::read_model( mmol, fpath, enable_messages );
        // MiniMol *pymmol = new MiniMol(mmol);
        return mmol;
      },
      // need to see how to read in spacegroup/cell
      // maybe should update clipper to exchange with gemmi
      // return std::unique_ptr<MiniMol>(new MiniMol(mmol)); }, // pymmol;
      // },
      nb::arg( "filepath" ) = "undefined", nb::arg( "enable_user_messages" ) = true,
      "Reads a coordinate file and return a structure in MiniMol." );
  m.def(
      "read_structure",
      []( MiniMol& mmol, const std::string& fpath, bool enable_messages ) {
        if ( fpath == "undefined" ) {
          throw std::invalid_argument( "No path/filename provided for input model! Aborting..." );
        }
        BuccaneerUtil::read_model( mmol, fpath, enable_messages );
        return ( mmol.model().size() > 0 );
        // MiniMol *pymmol = new MiniMol(mmol);
      },
      // return mmol; },
      //  need to see how to read in spacegroup/cell
      //  maybe should update clipper to exchange with gemmi
      //  return std::unique_ptr<MiniMol>(new MiniMol(mmol)); }, // pymmol;
      //  },
      nb::arg( "minimol" ), nb::arg( "filepath" ) = "undefined",
      nb::arg( "enable_user_messages" ) = true,
      "Reads a coordinate file and store structure into MiniMol" );
  m.def(
      "write_structure",
      []( MiniMol& mmol, const std::string& fpath, bool cif_format ) {
        clipper::GEMMIfile gfile;
        gfile.export_minimol( mmol );
        if ( fpath == "undefined" ) {
          throw std::invalid_argument( "No output path/filename provided! Aborting..." );
        }
        std::string filename = fpath.substr( 0, fpath.rfind( "." ) + 1 );
        if ( cif_format )
          gfile.write_file( filename + "cif", clipper::GEMMIfile::CIF );
        else
          gfile.write_file( filename + "pdb" );
      },
      nb::arg( "minimol" ), nb::arg( "filepath" ) = "undefined", nb::arg( "cif_format" ) = true );
}