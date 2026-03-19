"""
pytest for xmap bindings
"""

import pathlib
import pytest
import math
import bobkit.clipper as clipper
import bobkit.buccaneer as buccaneer
from bobkit.util import read_structure
import numpy as np
import gemmi

class TestNXmap:
    def test_nxmap_array(self):
        """Test if the array returned from Xmap is correctly filled
        """
        #gmap = gemmi.read_ccp4_map("/Users/swh514/Work/data/emd_3488/emd_3488_trimmed.mrc")
        gmap = gemmi.read_ccp4_map("/Users/swh514/Work/data/emd_0406/emd_0406.map")
        nxmap = clipper.NXmap_float()
        cell = clipper.Cell(gmap.grid.unit_cell.parameters)
        print("nxm2")
        nxm2 = clipper.NXmap_float(gmap.grid.array, cell)
        print("nxmap")
        nxmap.import_from_gemmi(gmap)
        grid_samp = clipper.Grid_sampling(nxmap.grid.nu, nxmap.grid.nv, nxmap.grid.nw)
        print("gmap array")
        garr = gmap.grid.array
        print("nxmap array")
        nxm_arr = nxmap.array
        print("assertions")
        for u in range(nxmap.grid.nu):
            for v in range(nxmap.grid.nv):
                for w in range(nxmap.grid.nw):
                    cg = clipper.Coord_grid(u, v, w)
                    cg = cg.unit(grid_samp)
                    #assert math.isclose(nxmap.get_data(cg), nxm_arr[u,v,w], abs_tol= 1.e-9)
                    assert nxmap.get_data(cg) == nxm_arr[u,v,w]
                    assert nxmap.get_data(cg) == gmap.grid.get_value(u,v,w)
                    assert nxmap.get_data(cg) == nxm2.get_data(cg)
                    assert nxm_arr[u,v,w] == garr[u,v,w]

if __name__ == "__main__":
    t = TestNXmap()
    t.test_nxmap_array()