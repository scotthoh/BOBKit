import pytest
import bobkit.clipper as clipper
import numpy as np


@pytest.fixture(scope="function")
def reso_obj(request):
    if request.param == "":
        return clipper.Resolution()
    else:
        return clipper.Resolution(float(request.param))


@pytest.fixture(scope="function")
def hkl_class_obj(request):
    spg, hkl = request.param
    return clipper.HKL_class(spg, hkl)


@pytest.fixture(scope="function")
def rtop_orth_obj(request):
    if request.param == "":
        return clipper.RTop_orth()
    else:
        rot, vec = request.param
        return clipper.RTop_orth(rot, vec)


class TestResolution:
    @pytest.mark.parametrize(
        "reso_obj, expected",
        [
            pytest.param("", True, id="resolution_null"),
            pytest.param(2.0, False, id="resolution_2.0"),
        ],
        indirect=["reso_obj"],
    )
    def test_resolution_null(self, reso_obj, expected):
        assert reso_obj.is_null() == expected

    @pytest.mark.parametrize(
        "reso_obj, limit, invresolsq",
        [
            pytest.param(2.0, 2.0, 1.0 / (2.0 * 2.0), id="resol_limit_2.0"),
            pytest.param(3.4, 3.4, 1.0 / (3.4 * 3.4), id="resol_limit_3.4"),
        ],
        indirect=["reso_obj"],
    )
    def test_resolution_limit(self, reso_obj, limit, invresolsq):
        assert reso_obj.limit() == limit
        assert reso_obj.invresolsq_limit() == invresolsq


class TestHKLclass:
    @pytest.mark.parametrize(
        "hkl_class_obj, expected",
        [
            pytest.param(
                (clipper.Spacegroup.p1(), [0, 1, 0]),
                (1.0, 1.0, 66.758844, False, False),
                id="hklclass_p1_010",
            ),
            pytest.param(
                (clipper.Spacegroup(5), clipper.HKL(2, 3, 1)),
                (0.0, 0.0, 0.0, True, True),
                id="hklclass_c2_231",
            ),
        ],
        indirect=["hkl_class_obj"],
    )
    def test_hklclass(self, hkl_class_obj, expected):
        epsilon, epsilonc, allowed, is_centric, is_sys_abs = expected
        assert hkl_class_obj.epsilon == epsilon
        assert (hkl_class_obj.allowed() - allowed) < 1.0e-6
        assert hkl_class_obj.epsilonc() == epsilonc
        assert hkl_class_obj.is_centric() == is_centric
        assert hkl_class_obj.is_sys_abs() == is_sys_abs


class TestRToporth:
    @pytest.mark.parametrize("rtop_orth_obj", [pytest.param((""))], indirect=["rtop_orth_obj"])
    def test_null_constructor(self, rtop_orth_obj):
        # null constructor RTop, values not null
        assert not rtop_orth_obj.is_null()
        rtop_null = rtop_orth_obj.null()
        assert rtop_null.is_null()

    @pytest.mark.parametrize(
        "rtop_orth_obj, expected",
        [
            #pytest.param(
            #    (clipper.Mat33_double.null(), clipper.Vec3_double.null()),
            #    (clipper.Mat33_double.null(), clipper.Vec3_double.null()),
            #),
            pytest.param(
                (
                    clipper.Mat33_double(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0),
                    clipper.Vec3_double(11.0, 22.0, 33.0),
                ),
                (
                    clipper.Mat33_double(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0),
                    clipper.Vec3_double(11.0, 22.0, 33.0),
                ),
            ),
        ],
        indirect=["rtop_orth_obj"],
    )
    def test_constructors(self, rtop_orth_obj, expected):
        rtop = rtop_orth_obj
        rot, vec = expected
        assert rtop.rot.get(0, 0) == rot.get(0, 0)
        assert rtop.trn[0] == vec[0]

# class TestHKLClass:
#    @pytest.mark.parametrize
#    def test_construct()

# class TestHKLClass:


#class TestRTopOrth:
#    r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
#    t = [0.0, 0.0, 0.0]
#    rt = clipper.RTop_orth(r, t)
#    rot = rt.rot
#    trn = rt.trn
#    assert trn.equals(clipper.Vec3_double(t), tol=0.1e-6)
#    assert rot.equals(clipper.Mat33_double(r), tol=0.1e-6)
#    inv = rt.inverse()
#    transform = rt.to_gemmi_transform(inv)
#    assert inv.rot.get(0, 0) == transform.mat.row_copy(0)[0]
#    assert inv.rot.get(1, 1) == transform.mat.row_copy(1)[1]
#    assert inv.rot.get(2, 2) == transform.mat.row_copy(2)[2]
#    assert inv.trn[0] == transform.vec.x
#    assert inv.trn[1] == transform.vec.y
#    assert inv.trn[2] == transform.vec.z
#