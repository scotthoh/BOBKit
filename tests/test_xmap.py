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
import time


class TestXmap:
    def test_xmap_array(self):
        """Test if the array returned from Xmap is correctly filled
        """
        grid = clipper.Grid_sampling(10, 10, 10)
        a = np.zeros((grid.nu, grid.nv, grid.nw), dtype=np.float32)
        # fill a
        count = 0
        axis_pos = [0, 1, 2]
        cell = clipper.Cell([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
        sg = clipper.Spacegroup(4)
        xm1 = clipper.Xmap_float(sg, cell, grid)
        ix = xm1.first_coord()
        while not ix.last():
            cg = ix.coord
            a[cg.u, cg.v, cg.w] = count
            xm1[ix] = count
            count += 1
            ix.next()
        xmap = clipper.Xmap_float(a, sg, cell)
        ix = xmap.first_coord()
        while not ix.last():
            cg = ix.coord
            print(ix.coord, xmap[ix], xm1.get_data(ix.coord), a[cg.u,cg.v,cg.w])
            assert xmap[ix] == a[cg.u,cg.v,cg.w]
            ix.next()
        xm_array = xmap.array
        ix = xmap.first_coord()
        while not ix.last():
            cg = ix.coord
            print(cg, xmap[ix], xm_array[cg.u, cg.v, cg.w])
            ix.next()

        #for s in range(1000):
        #    if s < 10:
        #        ix.next_u()
        #        ix.next_w()
        #    if s >= 10 and s < 700:
        #        ix.next_v()
        #        ix.prev_w()
        #    if s >= 700:
        #        ix.next_u()
        #        ix.prev_v()
        #        ix.next_w()
        #    print(xmap[ix])
        #    cg = xmap.to_map_unit(ix.coord)
        #    ind, sym = xmap.find_sym(cg)
        #    cgs = xmap.coord_of(ind)
        #    print(cg, xmap[ix], a[cgs.u, cgs.v, cgs.w]) #, xm_array[cg.u, cg.v, cg.w] )
        #    assert xmap[ix] == a[cgs.u, cgs.v, cgs.w]
        #    #assert xmap[ix] == xm_array[cg.u, cg.v, cg.w]


    def r(self):
        """Test if the array returned from Xmap is correctly filled
        """
        grid = clipper.Grid_sampling(72, 72, 72)
        a = np.zeros((grid.nu, grid.nv, grid.nw))
        # fill a
        count = 0
        axis_pos = [0, 1, 2]
        cell = clipper.Cell([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
        sg = clipper.Spacegroup(2)
        xm1 = clipper.Xmap_float(sg, cell, grid)
        ix = xm1.first()
        while not ix.last():
            cg = ix.coord
            a[cg.u, cg.v, cg.w] = count
            count += 1
            ix.next()
        #for i in range(72):
        #    for j in range(72):
        #        for k in range(72):
        #            a[i, j, k] = count
        #            count += 1

        xmap = clipper.Xmap_float(a, sg, cell)
        #ix = xmap.map_reference_coord(clipper.Coord_grid(0, 0, 0))
        ix = xmap.first()
        count = 0
        #xm_array = xmap.array
        while not ix.last():
            print(xmap[ix], )
        #for s in range(1000):
        #    if s < 10:
        #        ix.next_u()
        #        ix.next_w()
        #    if s >= 10 and s < 700:
        #        ix.next_v()
        #        ix.prev_w()
        #    if s >= 700:
        #        ix.next_u()
        #        ix.prev_v()
        #        ix.next_w()
        #    print(xmap[ix])
            cg = ix.coord.unit(xmap.grid_sampling)
            #ind, sym = xmap.find_sym(cg)
            #cgs = xmap.coord_of(ind)
            #print(ind, cgs)
            assert xmap[ix] == a[cg.u, cg.v, cg.w]
            print(cg, xmap[ix], a[cg.u, cg.v, cg.w]) #, xm_array[cg.u, cg.v, cg.w] )
            ix.next()

    def from_file(self):
        gm = gemmi.read_ccp4_map("/Users/swh514/Work/data/7ins.ccp4")
        xm1 = clipper.Xmap_float()
        xm1.import_from_gemmi(gm)
        ix = xm1.first()
        gma = gm.grid.array
        print("test0")
        while not ix.last():
            cg = ix.coord.unit(xm1.grid_sampling)
            #print(cg, xm1[ix], gm.grid.get_value(cg.u,cg.v,cg.w), gma[cg.u,cg.v,cg.w])
            assert xm1[ix]-gma[cg.u,cg.v,cg.w] < 1.e-5
            ix.next()
        print("test1")
        xm2 = clipper.Xmap_float(gma, clipper.Spacegroup(gm.grid.spacegroup.number), clipper.Cell(gm.grid.unit_cell.parameters))#, [0,1,2])
        ix = xm2.first()
        while not ix.last():
            cg = ix.coord.unit(xm2.grid_sampling)
            #print(cg, xm2[ix], gma[cg.u, cg.v, cg.w])
            assert xm2[ix]-gma[cg.u,cg.v,cg.w]<1.e-5
            ix.next()
        print("test2")
        print(f"gma array shape : {gma.shape}")
        xm2a = xm2.array
        print(f"xm2s array shape  : {xm2a.shape}")
        print(xm2a.shape)
        ix = xm2.first_coord()
        #if xm2.is_null(): print("null")
        for i in range(xm2a.shape[0]):
            for j in range(xm2a.shape[1]):
                for k in range(xm2a.shape[2]):
                    ix = xm2.map_reference_coord(clipper.Coord_grid(i,j,k))
                    print(f"{i}, {j}, {k} : {xm2[ix]}, {xm2a[i,j,k]}, {gma[i,j,k]}")
                    assert xm2a[i,j,k]-gma[i,j,k] < 1.e-5
                    assert xm2a[i,j,k]-xm2[ix] < 1.e-5
        print("Done")
        # not working array all corrupted 

    def test_gemmi(self):
        #gm = gemmi.FloatGrid(10,10,10)
        grid = clipper.Grid_sampling(10,10,10)
        a = np.zeros((10,10,10), dtype=np.float32)
        # fill a
        count = 0
        axis_pos = [0, 1, 2]
        cell = clipper.Cell([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
        sg = clipper.Spacegroup(4)
        xm1 = clipper.Xmap_float(sg, cell, grid)
        ix = xm1.first_coord()
        while not ix.last():
            cg = ix.coord
            a[cg.u, cg.v, cg.w] = count
            xm1[ix] = count
            count += 1
            ix.next()
        fg = gemmi.FloatGrid(a, gemmi.UnitCell(20,20,20,90,90,90),gemmi.SpaceGroup(4))
        arr = fg.array
        ix = xm1.first_coord()
        while not ix.last():
            cg = ix.coord
            print(cg, arr[cg.u,cg.v,cg.w])
            ix.next()

    def test_xmap_array_time(self):
        gm = gemmi.read_ccp4_map("/Users/swh514/Work/data/7ins.ccp4")
        xm1 = clipper.Xmap_float()
        xm1.import_from_gemmi(gm)
        start = time.perf_counter_ns()
        xm1_arr = xm1.array
        end = time.perf_counter_ns()
        elapsed_s = (end - start) / 1000000000
        print(f"array export elapsed time : {elapsed_s:.6f} s")
        print(f"xm2s array shape  : {xm1_arr.shape}")
        ix = xm1.first_coord()
        while not ix.last():
            cg = ix.coord
            #print(xm1[ix], xm1_arr[cg.u,cg.v,cg.w])
            assert xm1[ix] == xm1_arr[cg.u, cg.v, cg.w]
            ix.next()
        
        print("Done")

if __name__ == "__main__":
    t = TestXmap()
    #t.from_file()
    t.test_xmap_array_time()
    #t.test_xmap_array()
    #t.test_gemmi()