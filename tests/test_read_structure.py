"""
unit testing
"""
import os
import sys
import pathlib

import pytest
import math
import buildkit.clipper as clipper

@pytest.fixture
def cell_descr_instance():
    return clipper.Cell_descr(105.0, 105.0, 105.0, 90.0, 90.0, 90.0)


@pytest.fixture
def cell_instance(cell_descr_instance):
    return clipper.Cell(cell_descr_instance)


@pytest.fixture
def read_structure_instance(request, cell_instance):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mmol = clipper.MiniMol(clipper.Spacegroup.p1(), cell_instance)
    flag = clipper.read_structure(str(testdir / "pdb5ni1_cryst1.pdb"), mmol, False)
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
        mmol = clipper.MiniMol()
        assert mmol.is_null

    def test_empty_minimol(
        self,
        cell_descr_instance,
        cell_instance,
    ):
        mmol = clipper.MiniMol(clipper.Spacegroup.p1(), cell_instance)
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

    def test_minimol(self, read_structure_instance):
        flag, mmol = read_structure_instance
        assert flag
        exp = " ".join(repr(mmol).split())
        assert exp == "<clipper.MiniMol containing model with 4 chain(s)>"
        assert len(mmol) == 4
        assert str(mmol.spacegroup) == str(clipper.Spacegroup.p1())
        mmol_clone = mmol.clone()
        assert len(mmol_clone) == len(mmol)
        assert mmol_clone.model()[0].id == mmol.model()[0].id
        mmol2 = clipper.MiniMol()
        mmol2.model(mmol.model())
        assert mmol2.model().size() == 4
        assert mmol2.model()[3].id == "D"
        assert not mmol.is_empty()


class TestModel:
    def test_empty_model(self):
        model1 = clipper.MModel()
        assert len(model1) == 0
        exp = " ".join(repr(model1).split())
        assert exp == "<clipper.MModel containing 0 chain(s)>"

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


class TestChain:
    def test_empty_chain(self):
        chn = clipper.MChain()
        exp = " ".join(repr(chn).split())
        assert chn.size() == 0
        assert exp == "<clipper.MChain containing 0 residue(s)>"

    def test_chain(self, chain_instance):
        assert chain_instance.id == "A"
        assert chain_instance.size() == 149
        assert len(chain_instance) == 149
        exp = " ".join(repr(chain_instance).split())
        assert exp == "<clipper.MChain A containing 149 residue(s)>"
        assert chain_instance[0].id == "1"
        chn_clone = chain_instance.clone()
        assert chn_clone.id == "A"
        assert chn_clone["1"].type == "VAL"
        chn_clone[0] = chn_clone[1]
        assert chn_clone[0].id == chn_clone[1].id
        assert chn_clone[0].type == chn_clone[1].type
        res = chn_clone.find("6")
        assert res.id == "6"
        chn_tmp = clipper.MChain()
        assert chn_tmp.size() == 0
        chn_tmp.copy_from(chain_instance, clipper.COPY.COPY_MPC)
        assert chn_tmp.size() == 149
        assert chn_tmp[1].type == "LEU"
        assert chn_tmp[1].seqnum == 2


class TestResidue:
    def test_empty_residue(self):
        res = clipper.MResidue()
        exp = " ".join(repr(res).split())
        assert res.size() == 0
        assert exp == "<clipper.MResidue () containing 0 atom(s)>"

    def test_residue(self, residue_instance, chain_instance):
        assert clipper.MResidue.id_match(
            residue_instance.id, "1", clipper.MODE.UNIQUE
        )
        assert residue_instance.type == "VAL"
        assert residue_instance.seqnum == 1
        assert residue_instance.size() == 7
        assert len(residue_instance) == 7
        exp = " ".join(repr(residue_instance).split())
        assert exp == "<clipper.MResidue 1(VAL) containing 7 atom(s)>"
        assert residue_instance[0].element == "N"
        res_clone = residue_instance.clone()
        res_clone[0] = residue_instance[-1]
        assert res_clone[0].id == "CG2"
        assert res_clone[0].element == "C"
        atm_CG = res_clone.find("CG2")
        assert atm_CG.name == "CG2"
        res_tmp = clipper.MResidue()
        assert res_tmp.size() == 0
        res_tmp.copy_from(residue_instance, clipper.COPY.COPY_MP)
        assert res_tmp.type == "VAL"
        res_tmp.insert(residue_instance[0], pos=-1)
        res_tmp.insert(residue_instance[1], pos=-1)
        res_tmp.insert(residue_instance[2], pos=-1)
        assert res_tmp.size() == 3
        res_inclusive = residue_instance & res_tmp
        assert res_inclusive.size() == 3
        res_either = residue_instance | res_tmp
        assert res_either.size() == 7
        # UTILITY
        res_tmp.build_carbonyl_oxygen()
        assert math.isclose(res_tmp["O"].pos.x, 44.2493, rel_tol=1e-6)
        assert math.isclose(res_tmp["O"].pos.y, 53.5123, rel_tol=1e-6)
        assert math.isclose(res_tmp["O"].pos.z, 68.0153, rel_tol=1e-6)
        res_tmp.build_carbonyl_oxygen(chain_instance[1])
        assert math.isclose(res_tmp["O"].pos.x, 45.4563, rel_tol=1e-6)
        assert math.isclose(res_tmp["O"].pos.y, 52.3961, rel_tol=1e-6)
        assert math.isclose(res_tmp["O"].pos.z, 66.849, rel_tol=1e-6)
        assert res_tmp.number_of_rotamers() == 3
        assert res_tmp.number_of_rotamers(t=clipper.TYPE.Dunbrack) == 3

        assert clipper.MResidue.protein_peptide_bond(
            residue_instance, chain_instance[1]
        )
        assert not clipper.MResidue.protein_peptide_bond(
            residue_instance, chain_instance[2]
        )
        phi = clipper.MResidue.protein_ramachandran_phi(
            residue_instance, chain_instance[1]
        )
        assert math.isclose(phi, -2.898123, rel_tol=1e-6)
        psi = clipper.MResidue.protein_ramachandran_psi(
            residue_instance, chain_instance[1]
        )
        assert math.isclose(psi, -1.266930, rel_tol=1e-6)

        res_tmp2 = res_tmp.clone()
        res_tmp.build_sidechain_numbered_rotamer(0, t=clipper.TYPE.Dunbrack)
        res_tmp2.build_sidechain_numbered_rotamer(0)

        assert math.isclose(res_tmp["CB"].pos.x, 47.0931, rel_tol=1e-6)
        assert math.isclose(res_tmp["CB"].pos.y, 54.6763, rel_tol=1e-6)
        assert math.isclose(res_tmp["CB"].pos.z, 68.9786, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG1"].pos.x, 47.7043, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG1"].pos.y, 53.3674, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG1"].pos.z, 69.4531, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG2"].pos.x, 48.1477, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG2"].pos.y, 55.7693, rel_tol=1e-6)
        assert math.isclose(res_tmp["CG2"].pos.z, 68.9065, rel_tol=1e-6)

        assert math.isclose(res_tmp2["CB"].pos.x, 47.0931, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CB"].pos.y, 54.6763, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CB"].pos.z, 68.9786, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG1"].pos.x, 47.6851, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG1"].pos.y, 53.3616, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG1"].pos.z, 69.4642, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG2"].pos.x, 48.1634, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG2"].pos.y, 55.7541, rel_tol=1e-6)
        assert math.isclose(res_tmp2["CG2"].pos.z, 68.9008, rel_tol=1e-6)


class TestAtom:
    def test_empty_atom(self):
        atm = clipper.MAtom()
        assert atm.id == ""
        assert atm.name == ""
        assert atm.element == ""

    def test_atom(self, atom_instance):
        assert atom_instance.id == "N"
        assert atom_instance.name == "N"
        assert atom_instance.occupancy == 1.00
        assert atom_instance.b_iso == 80.64
        atom_instance.b_iso = 90.56
        assert atom_instance.b_iso == 90.56
        coord_pos1 = clipper.Coord_orth(45.716, 55.727, 67.167)
        assert atom_instance.pos.x == coord_pos1.x
        assert atom_instance.pos.y == coord_pos1.y
        assert atom_instance.pos.z == coord_pos1.z
        exp = " ".join(repr(atom_instance).split())
        assert exp == "<clipper.MAtom N at (45.716, 55.727, 67.167)>"
        atm_tmp = clipper.MAtom()
        assert atm_tmp.id == ""
        atm_tmp.copy_from(atom_instance, mode=clipper.COPY.COPY_MP)
        assert clipper.MAtom.id_match(
            atm_tmp.id, atom_instance.id, mode=clipper.MODE.UNIQUE
        )
