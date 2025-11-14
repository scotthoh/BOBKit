// Nanobind bindings for clipper edcalc
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>

using namespace clipper;

void add_cifdata_io( nb::module_ &m ) {
  nb::class_<CIFfile>( m, "CIFfile" )
      .def( nb::init<>(), "Constructor: does nothing" )
      .def( "open_read", &CIFfile::open_read, nb::arg("filename"), "Open a file for read access" )
      .def( "close_read", &CIFfile::close_read, "Close a file after reading." )
      .def_prop_ro( "spacegroup", &CIFfile::spacegroup, "Get file spacegroup." )
      .def_prop_ro( "cell", &CIFfile::cell, "Get file cell." )
      .def( "resolution", ( const Resolution &( CIFfile::* )() const ) & CIFfile::resolution,
                    "Get file resolution." )
      .def( "resolution", ( Resolution ( CIFfile::* )( const Cell & ) const ) & CIFfile::resolution,
            nb::arg( "cell" ),
            "Get resolution using supplied Cell object."
            "Resolution determined by the most extreme reflection in the file." )
      .def_prop_ro( "hkl_sampling", &CIFfile::hkl_sampling )
      .def( "import_hkl_info", &CIFfile::import_hkl_info, nb::arg( "hklinfo" ),
            "Import list of reflections from CIF file into given HKL_info object. If resolution is not found, "
            "then resolution will be determined from the input hkl data." )
      .def( "import_hkl_data", ( void ( CIFfile::* )( HKL_data_base & ) ) & CIFfile::import_hkl_data,
            "Mark a hkl_data for import from CIF file.")
      .def( "import_hkl_data", ( void ( CIFfile::* )( HKL_data_base &, const std::vector<String> &, const String & ) ) & CIFfile::import_hkl_data,
            nb::arg("hkldata"), nb::arg("column_names"), nb::arg("blockname")="*",
            "Mark a hkl_data with the given column name and optionally block name for import from CIF file")
      .def( "contains_phases_p", &CIFfile::contains_phases_p, "Check if file contains phases")
      .def( "contains_phases", &CIFfile::contains_phases,
            nb::arg("column_name"), nb::arg("verbose")=false,
            "Check all blocks if phases with given column name exists.")
      .def( "contains_phases_column", ( bool ( CIFfile::* )( const String &, const int &, const bool & ) const ) & CIFfile::contains_phases_column,
            nb::arg("column_name"), nb::arg("block_index")=-1, nb::arg("verbose")=false,
            "Check all blocks if phases with given column name and block index exists.")
      .def( "contains_phases", ( bool ( CIFfile::* )( const String &, const String &, const bool & ) const ) & CIFfile::contains_phases_column,
            nb::arg("column_name"), nb::arg("blockname")="*", nb::arg("verbose")=false,
            "Check all blocks if phases with given column name and blockname exists.")
      .doc() = "CIF import/export parent class for clipper objects.\n"
               "This is the import class which can be linked to an cif data "
               "file and be used to transfer data into a Clipper data structure.\n"
               "It is currently a read-only class.";
}