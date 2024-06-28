#include "type_conversions.h"
#include <clipper/clipper-gemmi.h>
#include <clipper/clipper.h>
#include <clipper/core/hkl_info.h>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace clipper;

void init_hklinfo(py::module &m) {
  py::class_<HKL_info> hklinfo(m, "HKL_info");
  hklinfo.def(py::init<>())
      .def(py::init<const Spacegroup &, const Cell &, const Resolution &,
                    const bool &>(),
           py::arg("spacegroup"), py::arg("cell"), py::arg("resolution"),
           py::arg("generate") = false)
      .def("init",
           (void(HKL_info::*)(const Spacegroup &, const Cell &,
                              const Resolution &, const bool &)) &
               HKL_info::init,
           py::arg("spacegroup"), py::arg("cell"), py::arg("resolution"),
           py::arg("generate") = false)
      .def("init",
           (void(HKL_info::*)(const Spacegroup &, const Cell &,
                              const HKL_sampling &, const bool &)) &
               HKL_info::init,
           py::arg("spacegroup"), py::arg("cell"), py::arg("hkl_sampling"),
           py::arg("generate") = true)
      // from gemmi types
      .def(
          "init",
          [](HKL_info &self, const gemmi::SpaceGroup &sg,
             const gemmi::UnitCell &uc, const double dmin, const double tol,
             const bool &generate) {
            self.init(GEMMI::spacegroup(sg), GEMMI::cell(uc),
                      Resolution(dmin - tol), generate);
          },
          py::arg("spacegroup"), py::arg("cell"), py::arg("dmin"),
          py::arg("tol") = 1.e-8, py::arg("generate") = false)
      .def("is_null", &HKL_info::is_null)
      .def_property_readonly("cell", &HKL_info::cell)
      .def_property_readonly("spacegroup", &HKL_info::spacegroup)
      .def_property_readonly("hkl_sampling", &HKL_info::hkl_sampling)
      .def_property_readonly("resolution", &HKL_info::resolution)

      // from gemmi::Mtz
      .def_static(
          "from_gemmi_mtz",
          [](const gemmi::Mtz &mtzobj, const double tol, const bool &generate) {
            return GEMMI::as_HKL_info(mtzobj, tol, generate);
          },
          py::arg("mtz"), py::arg("tol") = 1.e-8, py::arg("generate") = false)
      //.def_static("from_cell_spacegroup")
      .def("generate_hkl_list", &HKL_info::generate_hkl_list)
      //.def("add_hkl_list", &HKL_info::add_hkl_list, py::arg("hkls"))

      //.def("add_hkl_list",
      //     (void(HKL_info::*)(const std::vector<HKL> &add)) &
      //         HKL_info::add_hkl_list,
      //     py::arg("hkls"))
      .def("add_hkl_list",
           [](HKL_info &self, const py::array_t<int> &hkl) {
             std::vector<HKL> hkl_list;
             auto hbuf = hkl.request();
             if (hbuf.shape[1] != 3)
               throw std::logic_error(
                   "HKL list does not have the expected width!");
             auto l = hbuf.shape[0];
             int *hptr = (int *)hbuf.ptr;
             for (int i = 0; i < l; ++i) {
               hkl_list.emplace_back(HKL(*hptr, *(hptr + 1), *(hptr + 2)));
               hptr += 3;
             }
             // check hkl_list size and array size are the same
             if (hkl_list.size() != l)
              throw std::logic_error(
                  "Error in adding HKL list, length is not the same as
                  input.");
            else
              self.add_hkl_list(hkl_list);
           })
      .def("num_reflections", &HKL_info::num_reflections)
      .def("hkl_of", &HKL_info::hkl_of, py::arg("index"))
      .def("index_of", &HKL_info::index_of, py::arg("hkl"))
      .def("invresolsq", &HKL_info::invresolsq, py::arg("index"))
      .def_property_readonly("invresolsq_range", &HKL_info::invresolsq_range)
      .def("hkl_class", &HKL_info::hkl_class, py::arg("index"))
      .def("find_sym", &HKL_info::find_sym, py::arg("hkl"), py::arg("sym"),
           py::arg("friedel"))
      .def("first", &HKL_info::first)
      .def("debug", &HKL_info::debug)
      .def("__repr__", [](const HKL_info &self) {
        return "<clipper.HKL_info with spacegroup " +
               self.spacegroup().symbol_hm() + ", " +
               clipper::String(self.num_reflections()) + " reflections>";
      });

  using HKLB = HKL_info::HKL_reference_base;
  py::class_<HKLB>(hklinfo, "HKL_reference_base")
      .def("base_hkl_info", &HKLB::base_hkl_info)
      .def("index", &HKLB::index)
      .def("invresolsq",
           (ftype(HKLB::*)(const HKL_data_base &) const) & HKLB::invresolsq,
           py::arg("hkldata"))
      .def("invresolsq", (ftype(HKLB::*)() const) & HKLB::invresolsq)
      .def("last", &HKLB::last);

  using HKLI = HKL_info::HKL_reference_index;
  py::class_<HKLI, HKLB>(hklinfo, "HKL_reference_index")
      .def(py::init<>())
      .def(py::init<const HKL_info &, const int &>(), py::arg("hklinfo"),
           py::arg("index"))
      .def("hkl", &HKLI::hkl)
      .def("hkl_class", &HKLI::hkl_class)
      // note from Tristan: avoid creating new Python objects when incrementing
      .def("next", [](HKLI &self) { self.next(); });

  using HKLC = HKL_info::HKL_reference_coord;
  py::class_<HKLC, HKLB>(hklinfo, "HKL_reference_coord")
      .def(py::init<>())
      .def(py::init<const HKL_info &, const HKL &>(), py::arg("hklinfo"),
           py::arg("hkl"))
      .def_property("hlk", &HKLC::hkl,
                    [](HKLC &self, const HKL &hkl) { self.set_hkl(hkl); })
      .def_property_readonly("sym", &HKLC::sym)
      .def_property_readonly("friedel", &HKLC::friedel)
      // note from Tristan: avoid creating new Python objects when incrementing
      .def("next", [](HKLC &self) { self.next(); })
      .def("next_h", [](HKLC &self) { self.next_h(); })
      .def("next_k", [](HKLC &self) { self.next_k(); })
      .def("next_l", [](HKLC &self) { self.next_l(); })
      .def("prev_h", [](HKLC &self) { self.prev_h(); })
      .def("prev_k", [](HKLC &self) { self.prev_k(); })
      .def("prev_l", [](HKLC &self) { self.prev_l(); });
}