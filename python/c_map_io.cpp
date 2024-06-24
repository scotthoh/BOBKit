#include "c_map_io.h"

#define DBG std:: << "[" << __FUNCTION__ << "] - "

// functions numpy arrays to nxmap/xmap
template <class T>
void buildkit::numpy_to_nxmap(pybind11::array_t<T, 3> data,
                              const clipper::Cell &cell,
                              clipper::NXmap<T> &nxmap) {
  auto r = data.unchecked();
  clipper::Grid_sampling grid_sam(r.shape(0), r.shape(1), r.shape(2));
  clipper::Grid_range grid_map(
      clipper::Coord_grid(0, 0, 0),
      clipper::Coord_grid(r.shape(0) - 1, r.shape(1) - 1, r.shape(2) - 1));

  nxmap.init(cell, grid_sam, grid_map);
  // set map data to nxmap
  int g[3], gfms[3];
  gfms[0] = r.shape(0) - 1;
  gfms[1] = r.shape(1) - 1;
  gfms[2] = r.shape(2) - 1;
  NXmap_base::Map_reference_coord x(nxmap);
  for (g[2] = 0; g[2] <= gfms[2]; g[2]++)
    for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
      for (g[0] = 0; g[0] <= gfms[0]; g[0]++) {
        x.set_coord(clipper::Coord_grid(g[0], g[1], g[2]));
        nxmap[x] = r(g[0], g[1], g[2]);
      }
  std::cout << nxmap.grid().format() << std::endl;
}

template <class T>
void buildkit::numpy_to_xmap(pybind11::array_t<T, 3> data,
                             const clipper::Cell &cell,
                             clipper::Xmap<T> &xmap) {
  auto r = data.unchecked();
  clipper::Grid_sampling grid_sam(r.shape(0), r.shape(1), r.shape(2));
  clipper::Grid_range grid_map(
      clipper::Coord_grid(0, 0, 0),
      clipper::Coord_grid(r.shape(0) - 1, r.shape(1) - 1, r.shape(2) - 1));
  clipper::Spacegroup s = clipper::Spacegroup::p1();
  xmap.init(s, cell, grid_sam);
  // set map data to nxmap
  int g[3], gfms[3];
  gfms[0] = r.shape(0) - 1;
  gfms[1] = r.shape(1) - 1;
  gfms[2] = r.shape(2) - 1;
  Xmap_base::Map_reference_coord x(xmap);
  for (g[2] = 0; g[2] <= gfms[2]; g[2]++)
    for (g[1] = 0; g[1] <= gfms[1]; g[1]++)
      for (g[0] = 0; g[0] <= gfms[0]; g[0]++) {
        x.set_coord(clipper::Coord_grid(g[0], g[1], g[2]));
        xmap[x] = r(g[0], g[1], g[2]);
      }
  std::cout << xmap.grid_sampling().format() << std::endl;
}
// pybind11

namespace py = pybind11;
namespace bk = buildkit;
using namespace clipper;

template <class T> void declare_numpy_to_map(py::module &m) {
  m.def("numpy_to_nxmap", &bk::numpy_to_nxmap<T>);
  m.def("numpy_to_xmap", &bk::numpy_to_xmap<T>);
}

void init_map_io(py::module &m) {
  declare_numpy_to_map<ftype32>(m);
  declare_numpy_to_map<ftype64>(m);
}