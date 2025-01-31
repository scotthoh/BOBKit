// Wrapper for buccaneer-util
// Author: S.W.Hoh
// (C) 2023 -
// York Structural Biology Laboratory
// The University of York

#include <clipper/minimol/minimol_seq.h>
#include <pybind11/pybind11.h>

#include "buccaneer/buccaneer-util.h"
#include "helper_functions.h"
#include "type_conversions.h"

void test_array( const std::vector<double>& a ) { std::cout << &a << std::endl; }

void declare_test_array( py::module& m ) { m.def( "test_array", &test_array, py::arg( "array" ) ); }

void init_utils( py::module& m ) {
  // py::class_<BuccaneerUtil>( m, "Util", "Utility functions." )
  //.def_static("set_reference", &BuccaneerUtil::set_reference,
  //             py::arg("mtz"), py::arg("pdb"),
  //             "Set reference MTZ and PDB, only work when CCP4's CLIBD "
  //             "environment path is set.")
  //.def_static("read_model", &BuccaneerUtil::read_model, py::arg("mol"),
  //             py::arg("file"), py::arg("verbose") = true,
  //             "Read model and export to minimol.")
  declare_test_array( m );
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
      py::arg( "filepath" ) = "undefined", py::arg( "enable_user_messages" ) = true,
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
      py::arg( "minimol" ), py::arg( "filepath" ) = "undefined",
      py::arg( "enable_user_messages" ) = true,
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
      py::arg( "minimol" ), py::arg( "filepath" ) = "undefined", py::arg( "cif_format" ) = true );
}