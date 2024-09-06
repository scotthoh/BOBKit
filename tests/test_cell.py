import pytest
from bobkit import clipper
import gemmi
import math


@pytest.fixture
def cell_descr_instance():
    return clipper.Cell_descr(100.0, 105.0, 125.0, 90.0, 90.0, 90.0)


@pytest.fixture
def cell_instance(cell_descr_instance):
    return clipper.Cell(cell_descr_instance)


@pytest.fixture
def gemmi_cell_instance():
    return gemmi.UnitCell(26.42, 30.72, 33.01, 88.319, 108.095, 112.075)


class TestCellDescription:
    def test_cell_descr(self, cell_descr_instance):
        cdes = cell_descr_instance
        assert cdes.a == 100.0
        assert cdes.b == 105.0
        assert cdes.c == 125.0
        assert cdes.alpha_deg == 90.0
        assert cdes.alpha == clipper.Util.d2rad(90.0)


class TestCell:
    def test_cell(self, cell_instance, gemmi_cell_instance):
        cell = cell_instance
        cdes = clipper.Cell_descr(105.0, 105.0, 105.0, 90, 90, 90)
        assert not cell.equals(clipper.Cell(cdes))
        assert cell.a == 100.0
        assert cell.a_star == (
            cell.b * cell.c * math.sin(cell.alpha) / (cell.a * cell.b * cell.c)
        )
        assert cell.alpha_star == (
            math.acos(
                (
                    math.cos(cell.gamma) * math.cos(cell.beta)
                    - math.cos(cell.alpha)  # noqa 501
                )  # noqa 501
            )
            / (math.sin(cell.beta) * math.sin(cell.gamma))
        )
        cell.init(cdes)
        assert cell.a == 105.0
        assert cell.alpha_deg == 90.0
        assert cell.a_star == (
            cell.b * cell.c * math.sin(cell.alpha) / (cell.a * cell.b * cell.c)
        )
        assert cell.volume == (105.0 * 105.0 * 105.0)
        assert isinstance(cell.description, clipper.Cell_descr)
        cell_tmp = clipper.Cell(cdes)
        assert cell.equals(cell_tmp)
        orth_mat = cell.matrix_orth
        frac_mat = cell.matrix_frac
        assert orth_mat.get(0, 0) == 105.0
        assert orth_mat.get(1, 1) == 105.0
        assert orth_mat.get(2, 2) == 105.0
        g_cell = clipper.Cell.to_gemmi_cell(cell)
        assert g_cell.a == 105.0
        assert g_cell.alpha == 90.0

    def test_triclinic_cell(self, gemmi_cell_instance):
        c_temp = clipper.Cell.from_gemmi_cell(gemmi_cell_instance)
        assert c_temp.a == gemmi_cell_instance.a
        assert c_temp.b == gemmi_cell_instance.b
        assert c_temp.c == gemmi_cell_instance.c
        real_mat = "m00=698.016 m11=943.718 m22=1089.66 m01=-610.048 m02=-541.752 m12=59.4949"  # noqa 501
        assert c_temp.metric_real.format() == real_mat
