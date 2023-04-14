"""
unit testing
"""
import os
import sys
import pathlib

import pytest
import math
import buildkit


@pytest.fixture
def cell_descr_instance():
    return buildkit.Cell_descr(105.0, 105.0, 105.0, 90.0, 90.0, 90.0)


@pytest.fixture
def cell_instance(cell_descr_instance):
    return buildkit.Cell(cell_descr_instance)


@pytest.fixture
def read_structure_instance(request, cell_instance):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mmol = buildkit.MiniMol(buildkit.Spacegroup.p1(), cell_instance)
    flag = buildkit.read_structure(str(testdir / "pdb5ni1_cryst1.pdb"), mmol, False)
    return flag, mmol


@pytest.fixture
def model_instance(read_structure_instance):
    return read_structure_instance[1].model()


@pytest.fixture
def chain_instance(model_instance):
    return model_instance[0]


@pytest.fixture
def residue_instance(chain_instance):
    return chain_instance[0]


@pytest.fixture
def atom_instance(residue_instance):
    return residue_instance[0]


class TestReadStructure:
    def test_read_structure(self, read_structure_instance):
        flag, mmol = read_structure_instance

        assert flag
        assert len(mmol) == 4


class TestMiniMol:
    def test___init__(self):
        mmol = buildkit.MiniMol()
        assert mmol.is_null

    def test_empty_minimol(
        self,
        cell_descr_instance,
        cell_instance,
    ):
        mmol = buildkit.MiniMol(buildkit.Spacegroup.p1(), cell_instance)
        assert mmol.cell.equals(cell_instance)
        assert math.isclose(mmol.cell.a, cell_descr_instance.a, rel_tol=1e-6)
        assert math.isclose(mmol.cell.b, cell_descr_instance.b, rel_tol=1e-6)
        assert math.isclose(mmol.cell.c, cell_descr_instance.c, rel_tol=1e-6)
        assert math.isclose(
            mmol.cell.alpha_deg, cell_descr_instance.alpha_deg, rel_tol=1e-6
        )
        assert math.isclose(
            mmol.cell.beta_deg, cell_descr_instance.beta_deg, rel_tol=1e-6
        )
        assert math.isclose(
            mmol.cell.gamma_deg, cell_descr_instance.gamma_deg, rel_tol=1e-6
        )
        assert mmol.is_empty()
        assert not mmol.is_null

        exp = " ".join(repr(mmol).split())
        assert exp == "<clipper.MiniMol containing model with 0 chain(s)>"


class TestModel:
    def test_model(self, model_instance):
        model1 = model_instance
        assert model1.size() == 4
        assert model1[0]
        assert model1["A"]

        chn_A = model1.find("A")
        assert chn_A.id == "A"

        exp = " ".join(repr(model1).split())
        assert exp == "<clipper.MModel containing 4 chain(s)>"

        clone_model = model1.clone()
        assert clone_model.size() == 4
        assert clone_model[0].id == "A"
        clone_model[0] = clone_model[3]
        assert clone_model[0].id == "D"
