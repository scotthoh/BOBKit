"""
unit testing running implementation of Buccaneer in python
"""
import pathlib

import pytest
import math
import bobkit.clipper as clipper
import bobkit.buccaneer as buccaneer


@pytest.fixture
def cell_descr_instance():
    return clipper.Cell_descr(63.0, 63.0, 63.0, 90.0, 90.0, 90.0)


@pytest.fixture
def cell_instance(cell_descr_instance):
    return clipper.Cell(cell_descr_instance)


@pytest.fixture
def get_find_cif_instance(request, cell_instance):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mmol = clipper.MiniMol(clipper.Spacegroup.p1(), cell_instance)
    flag = clipper.read_structure(str(testdir / "build1_find.cif"), mmol, False)
    return flag, mmol


@pytest.fixture
def get_grow_cif_instance(request, cell_instance):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mmol = clipper.MiniMol(clipper.Spacegroup.p1(), cell_instance)
    flag = clipper.read_structure(str(testdir / "build1_find_grow.mmcif"), mmol, False)
    return flag, mmol

# find, grow, join, merge, link, sequence, correct, filter, ncs, prune, rebuild
# class _TestCaGrow:
#    def _test_check_find_points(self, get_find_cif_instance):
#        flag, mmol = get_find_cif_instance
#        atmlist = mmol.model().atom_list()
#
#        assert flag
#        assert len(atmlist) == 297


class TestCaJoin:
    def test_run_ca_join(self, get_grow_cif_instance):
        flag, mmol = get_grow_cif_instance
        cajoin = buccaneer.Ca_join(2.0, 2.0)
        
        assert len(mmol.model()) == 99
        cajoin(mmol)
        assert len(mmol.model()) == 20
