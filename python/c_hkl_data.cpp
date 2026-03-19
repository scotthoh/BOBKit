// Nanobind bindings for clipper hkl data
// Author: S.W.Hoh
// 2025 -
// York Structural Biology Laboratory
// The University of York

#include "commons.h"
#include <clipper/clipper-gemmi.h>
#include <gemmi/mtz.hpp>
#include "arrays.h"
#include <clipper/clipper.h>
//#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/unique_ptr.h>
//#include <nanobind/stl/vector.h>
//#include <nanobind/stl/array.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/trampoline.h>

using namespace clipper;

template <class C> void catch_null( const C &c ) {
  if ( c.is_null() )
    throw std::length_error( "Array is not initialised!" );
}

template <class C1, class C2> void catch_mismatched_lengths( const C1 &c1, const C2 &c2 ) {
  if ( c1.is_null() || c1.data_size() != c2.data_size() )
    throw std::out_of_range( "Array sizes must match!" );
}

template <class C1, class C2, class C3> void catch_mismatched_lengths( const C1 &c1, const C2 &c2, const C3 &c3 ) {
  if ( c1.is_null() || !( ( c1.data_size() == c2.data_size() ) && ( c2.data_size() == c3.data_size() ) ) )
    throw std::out_of_range( "Array sizes must match!" );
}

template <class Derived, class dtype> auto hkldatatype_to_array( const Derived &self, const HKL &hkl ) {
  dtype *data = new dtype[self.data_size()];
  self.data_export( hkl, data );
  // Delete 'data' when the 'owner' capsule expires
  nb::capsule owner( data, []( void *p ) noexcept { delete[] ( dtype * )p; } );
  return nb::ndarray<nb::numpy, dtype, nb::ndim<1>>( data, { ( size_t )self.data_size() }, owner ).cast();
}

// trampoline class for HKL data base
class PyHKLdata_base : public HKL_data_base {
public:
  /* Inherit the constructors. */
  // using Property_base::Property_base;
  NB_TRAMPOLINE( HKL_data_base, 9 );

  /* Trampoline (need one for each virtual function) */
  void update() override { NB_OVERRIDE_PURE( update ); }

  String type() const override { NB_OVERRIDE_PURE( type ); }

  bool missing( const int &index ) const override { NB_OVERRIDE_PURE( missing, index ); }

  void set_null( const int &index ) override { NB_OVERRIDE_PURE( set_null, index ); }

  int data_size() const override { NB_OVERRIDE_PURE( data_size ); }

  String data_names() const override { NB_OVERRIDE_PURE( data_names ); }

  void data_export( const HKL &hkl, xtype array[] ) const override { NB_OVERRIDE_PURE( data_export, hkl, array ); }

  void data_import( const HKL &hkl, const xtype array[] ) override { NB_OVERRIDE_PURE( data_import, hkl, array ); }

  void mask( const HKL_data_base &msk ) override { NB_OVERRIDE_PURE( mask, msk ); }
};

template <class C, class B> void add_scale_methods( nb::class_<HKL_data<C>, B> &pyclass ) {
  using Class = HKL_data<C>;
  pyclass
      .def(
          "scale",
          []( Class &self, const int &index, const double &factor ) {
            catch_null( self );
            if ( !self[index].missing() ) {
              self[index].scale( factor );
            } else {
              throw std::out_of_range( "No equivalent HKL has been indexed for this dataset!" );
            }
          },
          "Scale data at index with given factor" )
      .def(
          "scale_all",
          []( Class &self, const double &factor ) {
            catch_null( self );
            for ( auto ih = self.first(); !ih.last(); ih.next() ) {
              if ( !self[ih].missing() )
                self[ih].scale( factor );
            }
          },
          "Scale all data with given factor" );
}

template <class C, class B>
void add_hkl_data_base_methods( nb::class_<HKL_data<C>, B> &pyclass, std::string &docstring ) {
  using Class = HKL_data<C>;
  docstring.append( "An actual hkl_data object, containing actual data of "
                    "type T. This implements the generic interface, and "
                    "in addition provides type-specific access functions." );
  pyclass.def( nb::init<>() )
      .def( nb::init<const HKL_info &>(), nb::arg( "hklinfo" ), "Constructor from parent HKL_info." )
      .def( nb::init<const HKL_info &, const Cell &>(), nb::arg( "hklinfo" ), nb::arg( "cell" ),
            "Constructor from parent HKL_info and cell." )
      .def( nb::init<const Spacegroup &, const Cell &, const HKL_sampling &>(), nb::arg( "spacegroup" ),
            nb::arg( "cell" ), nb::arg( "hkl_sampling" ), "Constructor from spacegroup, cell and HKL_sampling." )
      .def( nb::init<const HKL_data_base &>(), nb::arg( "hkldata" ), "Constructor from another HKL_data object." )
      .def( "init", ( void ( Class::* )( const HKL_info &, const Cell & ) )&Class::init, nb::arg( "hklinfo" ),
            nb::arg( "cell" ), "Initialiser from parent hkl_info and cell." )
      .def( "init", ( void ( Class::* )( const Spacegroup &, const Cell &, const HKL_sampling & ) )&Class::init,
            nb::arg( "spacegroup" ), nb::arg( "cell" ), nb::arg( "hkl_sampling" ),
            "Initialiser from spacegroup, cell, and hkl_sampling." )
      .def( "init", ( void ( Class::* )( const HKL_data_base & ) )&Class::init, nb::arg( "hkldata" ),
            "Initialiaser from another HKL_data object." )
      .def( "update", &Class::update, "Update: synchronise info with parent HKL_info." )
      .def_prop_ro( "type", &Class::type,
                    "Get data type (a list of names corresponding to the im/export "
                    "values)." )
      .def( "missing", &Class::missing, "Check if a data entry in the list is marked as \'missing\'." )
      .def( "set_null", &Class::set_null, "Set data entry in the list to its null value." )
      .def_prop_ro( "data_size", &Class::data_size, "Return number of data elements in this type." )
      .def_prop_ro( "data_names", &Class::data_names, "Return names of data elements in this type" )
      //.def("get_data", &Class::get_data,
      //     "Get data by HKL_reference_coord (returns false if no equivalent "
      //     "hkl).")
      //.def("set_data", &Class::set_data,
      //     "Set data by HKL_reference_coord (returns false if no equivalent "
      //     "hkl).")
      .def(
          "data_export",
          []( const Class &self, const HKL &hkl ) { //}, std::vector<xtype> &vals) {
            return hkldatatype_to_array<Class, xtype>( self, hkl );
          },
          nb::rv_policy::reference_internal, "Export data at given HKL as array (for I/O)" )
      .def( "data_import", []( Class &self, const HKL &hkl, const std::vector<xtype> &vals) {
        self.data_import(hkl, vals.data());
      }, "Import data for given HKL from array (for I/O)" )
      .def_prop_ro( "hkl_array", []( const Class &self ) {
        auto arr = make_numpy_array_2d<int>({(size_t)self.num_obs(), 3});
        size_t count = 0;
        for (auto ih = self.first_data(); !ih.last(); self.next_data(ih)) {
          arr(count,0) = ih.hkl().h();
          arr(count,1) = ih.hkl().k();
          arr(count,2) = ih.hkl().l();
          count++;
        }
        return arr;
      }, nb::rv_policy::automatic)
      .def_prop_ro( "array", []( const Class &self ) {
          auto arr = make_numpy_array_2d<xtype>({(size_t)self.num_obs(), (size_t)self.data_size()});
          //xtype* arr = new xtype[self.num_obs() * self.data_size()];
          size_t count = 0;
          std::vector<xtype> v(self.data_size());
          for (auto ih = self.first_data(); !ih.last(); self.next_data(ih)) {
            self.data_export(ih.hkl(), v.data());
            for (auto j =0; j<self.data_size(); j++) {
              arr(count, j) = v[j];
            }
            count++;
          }
          return arr;
      }, "Export all non-missing data to numpy array.", nb::rv_policy::automatic)
      .def( "data_import_hkls", [](Class &self, const std::vector<std::array<int,3>> &hkls,const nb::ndarray<nb::numpy, xtype, nb::ndim<2>> &vals) {
        auto dv = vals.view();
        if ( dv.ndim() != 2 ) throw std::runtime_error("Data array must be 2-Dimension");
        if ( hkls.size() != dv.shape(0) ) throw std::runtime_error("HKL array has different length from data array.");
        for ( size_t i = 0; i <  dv.shape(0); i++ ) {
          self.data_import(HKL(hkls[i][0], hkls[i][1], hkls[i][2]), &dv(i,0));
        }
      }, nb::arg("hkls"), nb::arg("data"), "Import data from array for a list of HKL indices.")
      .def( "mask", &Class::mask, "Mask the data by marking any data missing in \'mask\' as missing." )
      .def(
          "__getitem__",
          []( const Class &self, const HKL_info::HKL_reference_index &i ) {
            catch_null( self );
            return self[i];
          },
          "Data accessor by HKL_reference_index.", nb::rv_policy::reference_internal )
      .def(
          "__setitem__",
          []( Class &self, const HKL_info::HKL_reference_index &i, const C &data ) {
            catch_null( self );
            self[i] = data;
          },
          "Data writer byb HKL_reference_index" )
      .def(
          "__getitem__",
          []( const Class &self, const HKL_info::HKL_reference_coord &ih ) {
            catch_null( self );
            C data;
            if ( self.get_data( ih, data ) )
              return data;
            throw std::out_of_range( "No data equivalent to that HKL!" );
          },
          "Data accessor by HKL_reference_coord." )
      .def(
          "__setitem__",
          []( Class &self, const HKL_info::HKL_reference_coord &ih, const C &data ) {
            if ( !self.set_data( ih, data ) )
              throw std::out_of_range( "No equivalent HKL has been indexed for this dataset!" );
          },
          "Data writer by HKL_reference_coord." )

      .def(
          "__getitem__",
          []( const Class &self, const int &index ) {
            catch_null( self );
            return self[index];
          },
          "Data accessor by index." )
      .def(
          "__setitem__",
          []( Class &self, const int &index, const C &data ) {
            catch_null( self );
            self[index] = data;
          },
          "Data writer by index." )
      .def(
          "__getitem__",
          []( const Class &self, const HKL &hkl ) {
            catch_null( self );
            C data;
            if ( self.get_data( hkl, data ) )
              return data;
            throw std::out_of_range( "No data equivalent to that HKL!" );
          },
          "Data accessor by hkl." )
      .def(
          "__setitem__",
          []( Class &self, const HKL &hkl, const C &data ) {
            if ( !self.set_data( hkl, data ) )
              throw std::out_of_range( "No equivalent HKL has been indexed for this dataset!" );
          },
          "Data writer by hkl." )
      .def(
          "copy_from", []( Class &self, const Class &other ) { self = other; }, "Copy from another hkl_data object." )
      .def(
          "set_all_values_to", []( Class &self, const C &value ) { self = value; }, "Set all values to given value." )
      .def(
          "restrict_to",
          [](Class &self, const HKL_info &hklinfo) {
            std::unique_ptr<Class> hkldata(new Class(hklinfo));
            auto ih = hklinfo.first();
            while (!ih.last()) {
              hkldata->set_data(ih.hkl(), self[ih.hkl()]);
              ih.next();
            }
            return hkldata.release();
            //return std::unique_ptr<Class>(new_array);
          }, nb::arg("hklinfo"),
          "Returns a new HKL_data object containing only those reflections in "
          "the given HKL_info.")
      .def( nb::self & nb::self )
      .def( nb::self | nb::self )
      .def( nb::self ^ nb::self )
      .def( !nb::self )
      // from clipper-gemmi
      .def( "import_from_gemmi", []( Class &self, const gemmi::Mtz &mtzobj, const std::string mtzpath, const bool &legacy) {
        String colpath(mtzpath);
        if (legacy){
          if (colpath.find("/") == String::npos && colpath.find("[") == String::npos) {
            colpath = "/*/*/[" + colpath + "]";
          }
        }
        GEMMI::import_hkl_data(self, mtzobj, colpath);
        }, nb::arg("mtz"), nb::arg("column_names"), nb::arg("legacy") = false,
        "Import HKL_data of specified columns from gemmi.Mtz object."
      )
      .doc() = docstring.c_str();
}

void declare_hkl_data_base( nb::module_ &m ) {
  nb::class_<HKL_data_base, PyHKLdata_base> /*, std::unique_ptr<HKL_data_base, nb::deleter<HKL_data_base>>>*/
      hkldatabase( m, "_HKL_data_base" );
  hkldatabase.def( "is_null", &HKL_data_base::is_null, "Test if object has been initialised." )
      .def_prop_ro( "base_hkl_info", &HKL_data_base::base_hkl_info, "Get the parent HKL_info object.",
                    nb::rv_policy::reference_internal )
      .def_prop_ro( "base_cell", &HKL_data_base::base_cell, "Get the parent cell.", nb::rv_policy::reference_internal )
      .def_prop_ro( "spacegroup", &HKL_data_base::spacegroup, "Get spacegroup." )
      .def_prop_ro( "cell", &HKL_data_base::cell, "Get cell." )
      .def_prop_ro( "resolution", &HKL_data_base::resolution, "Get resolution." )
      .def_prop_ro( "hkl_sampling", &HKL_data_base::hkl_sampling, "Get HKL_sampling" )
      .def_prop_ro( "hkl_info", &HKL_data_base::hkl_info, "Get HKL_info object.", nb::rv_policy::reference_internal )
      .def( "invresolsq", &HKL_data_base::invresolsq, "Get resolution by reflection index (based on true cell)." )
      .def( "invresolsq_range", &HKL_data_base::invresolsq_range,
            "Get resolution limits of the list (based on true cell and missing data)." )
      .def( "num_obs", &HKL_data_base::num_obs, "Get number of observations in this list (based on missing data)." )
      .def( "first", &HKL_data_base::first, "Return HKL_reference_index pointing to first reflection." )
      .def( "first_data", &HKL_data_base::first_data, "Return HKL_reference_index pointing to first non-missing data" )
      .def(
          "next_data", []( const HKL_data_base &self, HKL_info::HKL_reference_index &ih ) { self.next_data( ih ); },
          "increment HKL_reference_index to next non-missing data." )
      .def( "debug", &HKL_data_base::debug, "Output debugging details." )
      .doc() = "HKL_data_base.\n"
               "This is the virtual base for the typed hkl_data objects. "
               "It exists to guarantee and interface by which data can "
               "be managed without knowledge of the specific data type.";
}

template <class T> void declare_hkl_data_i_sigi( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<I_sigI<T>>;
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_I_sigI_" + std::string( dtype ) ).c_str() );
  std::string docstring = "HKL_data<I_sigI<>>.\n";
  add_hkl_data_base_methods<I_sigI<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<I_sigI<T>, HKL_data_base>( pyclass );
  pyclass
      // compute methods from hkl_compute.h
      .def(
          "compute_scale_u_iso_isigi",
          []( Class &self, const T &scale, const T &u_value, const Class &isigi ) {
            catch_null( self );
            self.compute( isigi, Compute_scale_u_iso<I_sigI<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and isotropic U to any scalable datatype." )
      .def(
          "compute_scale_u_aniso_isigi",
          []( Class &self, const T &scale, const U_aniso_orth &u_value, const Class &isigi ) {
            catch_null( self );
            self.compute( isigi, Compute_scale_u_aniso<I_sigI<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and anisotropic U to any scalable datatype." );
} // declare_hkl_data_isigi

template <class T> void declare_hkl_data_i_sigi_ano( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<I_sigI_ano<T>>;
  std::string docstring = "HKL_data<I_sigI_ano<>>.\n";
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_I_sigI_ano_" + std::string( dtype ) ).c_str() );
  add_hkl_data_base_methods<I_sigI_ano<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<I_sigI_ano<T>, HKL_data_base>( pyclass );
}

template <class T> void declare_hkl_data_f_sigf( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<F_sigF<T>>;
  std::string docstring = "HKL_data<F_sigF<>>.\n";
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_F_sigF_" + std::string( dtype ) ).c_str() );
  add_hkl_data_base_methods<F_sigF<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<F_sigF<T>, HKL_data_base>( pyclass );
  pyclass
      // compute methods from hkl_compute.h
      .def(
          "compute_mean_from_fano",
          []( Class &self, const HKL_data<F_sigF_ano<T>> &fano ) {
            catch_null( self );
            self.compute( fano, Compute_mean_fsigf_from_fsigfano<T>() );
          },
          nb::arg( "fanom" ), "Compute from F_sigF_anom to F_sigF (mean structure factor)." )
      .def(
          "compute_diff_from_fano",
          []( Class &self, const HKL_data<F_sigF_ano<T>> &fano ) {
            catch_null( self );
            self.compute( fano, Compute_diff_fsigf_from_fsigfano<T>() );
          },
          nb::arg( "fanom" ), "Compute from F_sigF_anom to F_sigF (difference structure factor)." )
      .def(
          "compute_scale_u_iso_fsigf",
          []( Class &self, const T &scale, const T &u_value, const Class &fsigf ) {
            catch_null( self );
            self.compute( fsigf, Compute_scale_u_iso<F_sigF<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and isotropic U to any scalable datatype." )
      .def(
          "compute_scale_u_aniso_fsigf",
          []( Class &self, const T &scale, const U_aniso_orth &u_value, const Class &fsigf ) {
            catch_null( self );
            self.compute( fsigf, Compute_scale_u_aniso<F_sigF<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and anisotropic U to any scalable datatype." );
} // declare_hkl_data_fsigf

template <class T> void declare_hkl_data_f_sigf_ano( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<F_sigF_ano<T>>;
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_F_sigF_ano_" + std::string( dtype ) ).c_str() );
  std::string docstring = "HKL_data<F_sigF_ano<>>.\n";
  add_hkl_data_base_methods<F_sigF_ano<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<F_sigF_ano<T>, HKL_data_base>( pyclass );
  pyclass
      // compute methods from hkl_compute.h
      .def(
          "compute_scale_u_iso_fsigfano",
          []( Class &self, const T &scale, const T &u_value, const Class &fsigfano ) {
            catch_null( self );
            self.compute( fsigfano, Compute_scale_u_iso<F_sigF_ano<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and isotropic U to any scalable datatype." )
      .def(
          "compute_scale_u_aniso_fsigfano",
          []( Class &self, const T &scale, const U_aniso_orth &u_value, const Class &fsigfano ) {
            catch_null( self );
            self.compute( fsigfano, Compute_scale_u_aniso<F_sigF_ano<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and anisotropic U to any scalable datatype." );
} // declare_hkl_data_fsigf_ano

template <class T> void declare_hkl_data_e_sige( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<E_sigE<T>>;
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_E_sigE_" + std::string( dtype ) ).c_str() );
  std::string docstring = "HKL_data<E_sigE<>>.\n";
  add_hkl_data_base_methods<E_sigE<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<E_sigE<T>, HKL_data_base>( pyclass );
  // auto pyclass = declare_HKL_data<E_sigE<T>>(m, class_str, docstring);
  pyclass
      // compute methods from hkl_compute.h
      .def(
          "compute_from_fsigf",
          []( Class &self, const HKL_data<F_sigF<T>> &fsigf ) {
            catch_null( self );
            self.compute( fsigf, Compute_EsigE_from_FsigF<T>() );
          },
          nb::arg( "data" ), "Compute from F_sigF to E_sigE" )
      // extra useful methods carried over from SWIG wrappings
      .def(
          "scale_by_sqrt_resolution",
          []( Class &self, const ResolutionFn &escale ) {
            catch_null( self );
            for ( clipper::HKL_data_base::HKL_reference_index ih = self.first(); !ih.last(); ih.next() )
              if ( !self[ih].missing() )
                self[ih].scale( sqrt( escale.f( ih ) ) );
          },
          nb::arg( "resfn" ), "Apply scale (sqrt resolution)." )
      .def(
          "scale_by_resolution",
          []( Class &self, const ResolutionFn &escale ) {
            catch_null( self );
            for ( clipper::HKL_data_base::HKL_reference_index ih = self.first(); !ih.last(); ih.next() )
              if ( !self[ih].missing() )
                self[ih].scale( escale.f( ih ) );
          },
          nb::arg( "resfn" ), "Apply scale (resolution)." );
} // declare_hkl_data_e_sige

template <class T> void declare_hkl_data_f_phi( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<F_phi<T>>;
  std::string docstring = "HKL_data<F_phi<>>.\n";
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_F_phi_" + std::string( dtype ) ).c_str() );
  add_hkl_data_base_methods<F_phi<T>, HKL_data_base>( pyclass, docstring );
  add_scale_methods<F_phi<T>, HKL_data_base>( pyclass );
  pyclass
      .def(
          "compute_neg",
          []( Class &self, const Class &other ) {
            // catch_mismatched_lengths(self, other);
            catch_null( self );
            self.compute( other, Compute_neg_fphi<T>() );
          },
          nb::arg( "fphi" ), "Negate F_phi (i.e. advance phase by pi)." )
      .def(
          "compute_add_fphi",
          []( Class &self, const Class &fphi1, const Class &fphi2 ) {
            catch_null( self );
            self.compute( fphi1, fphi2, Compute_add_fphi<T>() );
          },
          nb::arg( "fphi1" ), nb::arg( "fphi2" ), "Add two F_phi datalists." )
      .def(
          "compute_sub_fphi",
          []( Class &self, const Class &fphi1, const Class &fphi2 ) {
            catch_null( self );
            self.compute( fphi1, fphi2, Compute_sub_fphi<T>() );
          },
          nb::arg( "fphi1" ), nb::arg( "fphi2" ), "Subtract two F_phi datalists." )
      .def(
          "compute_from_fsigf_phifom",
          []( Class &self, const HKL_data<F_sigF<T>> &fsigf, const HKL_data<Phi_fom<T>> &phifom ) {
            catch_null( self );
            self.compute( fsigf, phifom, Compute_fphi_from_fsigf_phifom<T>() );
          },
          nb::arg( "fsigf" ), nb::arg( "phifom" ), "Compute from F_sigF+Phi_fom to F_phi." )
      .def(
          "compute_scale_u_iso_fphi",
          []( Class &self, const T &scale, const T &u_value, const HKL_data<F_phi<T>> &fphi ) {
            catch_null( self );
            self.compute( fphi, Compute_scale_u_iso<F_phi<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and isotropic U to any scalable datatype." )
      .def(
          "compute_scale_u_aniso_fphi",
          []( Class &self, const T &scale, const U_aniso_orth &u_value, const HKL_data<F_phi<T>> &fphi ) {
            catch_null( self );
            self.compute( fphi, Compute_scale_u_aniso<F_phi<T>>( scale, u_value ) );
          },
          nb::arg( "scale" ), nb::arg( "u" ), nb::arg( "data" ),
          "Apply scale and anisotropic U to any scalable datatype." )
      // from hkl_operators.h
      .def( nb::self + nb::self )
      .def( nb::self - nb::self )
      .def( nb::self * float() )
      .def( float() * nb::self )
      .def( -nb::self );
} // declare_hkl_data_f_phi

template <class T> void declare_hkl_data_phi_fom( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<Phi_fom<T>>;
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_Phi_fom_" + std::string( dtype ) ).c_str() );
  std::string docstring = "HKL_data<Phi_fom<>>.\n";
  add_hkl_data_base_methods<Phi_fom<T>, HKL_data_base>( pyclass, docstring );
  pyclass
      // compute methods from hkl_compute.h
      .def(
          "compute_from_abcd",
          []( Class &self, const HKL_data<ABCD<T>> &abcd ) {
            catch_null( self );
            self.compute( abcd, Compute_phifom_from_abcd<T>() );
          },
          nb::arg( "abcd" ),
          "Compute from ABCD to Phi_fom by phase integration (loses "
          "bimodality)." );
} // declare_hkl_data_phi_fom

template <class T> void declare_hkl_data_abcd( nb::module_ &m, const char *dtype ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<ABCD<T>>;
  nb::class_<Class, HKL_data_base> pyclass( m, ( "HKL_data_ABCD_" + std::string( dtype ) ).c_str() );
  std::string docstring = "HKL_data<ABCD<>>.\n";
  add_hkl_data_base_methods<ABCD<T>, HKL_data_base>( pyclass, docstring );
  pyclass
      // compute method from hkl_compute.h
      .def(
          "compute_from_phi_fom",
          []( Class &self, const HKL_data<Phi_fom<T>> &phiw ) {
            catch_null( self );
            self.compute( phiw, Compute_abcd_from_phifom<T>() );
          },
          nb::arg( "phifom" ), "Compute from Phi_fom to ABCD ( C = D = 0 )." )
      .def(
          "compute_add_abcd",
          []( Class &self, const Class &abcd1, const Class &abcd2 ) {
            catch_null( self );
            self.compute( abcd1, abcd2, Compute_add_abcd<T>() );
          },
          nb::arg( "abcd1" ), nb::arg( "abcd2" ), "Add two ABCD datalists." )
      // from hkl_operators.h
      .def( nb::self + nb::self );
  ;
}

void declare_hkl_data_flag( nb::module_ &m ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<Flag>;
  nb::class_<Class, HKL_data_base> pyclass( m, "HKL_data_Flag" );
  std::string docstring = "HKL_data<Flag>.\n";
  add_hkl_data_base_methods<Flag, HKL_data_base>( pyclass, docstring );
  pyclass
      // from hkl_operators.h
      .def( nb::self == int() )
      .def( nb::self != int() )
      .def( nb::self >= int() )
      .def( nb::self <= int() )
      .def( nb::self > int() )
      .def( nb::self < int() );
}

void declare_hkl_data_flag_bool( nb::module_ &m ) {
  using namespace clipper::datatypes;
  using Class = HKL_data<Flag_bool>;
  nb::class_<Class, HKL_data_base> pyclass( m, "HKL_data_Flag_bool" );
  std::string docstring = "HKL_data<Flag_bool>.\n";
  add_hkl_data_base_methods<Flag_bool, HKL_data_base>( pyclass, docstring );
}

void init_hkl_data( nb::module_ &m ) { // nb::module_ &m32, nb::module_ &m64) {
  declare_hkl_data_base( m );
  {
    using namespace clipper::datatypes;
    // Non-floating-point datatypes go in the main module
    declare_hkl_data_flag( m );
    declare_hkl_data_flag_bool( m );
  }

  {
    using namespace clipper::data32;
    // 32-bit floating point datatypes
    const char *suffix = "float";
    typedef ftype32 dtype;
    auto module = m; // m32
    declare_hkl_data_i_sigi<dtype>( module, suffix );
    declare_hkl_data_i_sigi_ano<dtype>( module, suffix );
    declare_hkl_data_f_sigf<dtype>( module, suffix );
    declare_hkl_data_f_sigf_ano<dtype>( module, suffix );
    declare_hkl_data_e_sige<dtype>( module, suffix );
    declare_hkl_data_f_phi<dtype>( module, suffix );
    declare_hkl_data_phi_fom<dtype>( module, suffix );
    declare_hkl_data_abcd<dtype>( module, suffix );
  }

  {
    using namespace clipper::data64;
    // 64-bit floating point datatypes
    const char *suffix = "double";
    typedef ftype64 dtype;
    auto module = m; // m64;
    declare_hkl_data_i_sigi<dtype>( module, suffix );
    declare_hkl_data_i_sigi_ano<dtype>( module, suffix );
    declare_hkl_data_f_sigf<dtype>( module, suffix );
    declare_hkl_data_f_sigf_ano<dtype>( module, suffix );
    declare_hkl_data_e_sige<dtype>( module, suffix );
    declare_hkl_data_f_phi<dtype>( module, suffix );
    declare_hkl_data_phi_fom<dtype>( module, suffix );
    declare_hkl_data_abcd<dtype>( module, suffix );
  }
}