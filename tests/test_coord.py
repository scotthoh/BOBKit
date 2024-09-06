import pytest
from bobkit import clipper
import math


class TestResolution:
    def test_resolution(self):
        res = clipper.Resolution()
        assert res.is_null()
        res = clipper.Resolution(2.0)
        assert res.limit() == 2.0
        assert res.invresolsq_limit() == 1.0 / (2.0 * 2.0)


# class TestHKLClass:


class TestRTopOrth:
    r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    t = [0.0, 0.0, 0.0]
    rt = clipper.RTop_orth(r, t)
    rot = rt.rot()
    trn = rt.trn()
    assert trn.equals(clipper.Vec3_double(t), tol=0.1e-6)
    assert rot.equals(clipper.Mat33_double(r), tol=0.1e-6)
    inv = rt.inverse()
    transform = rt.to_gemmi_transform(inv)
    assert inv.rot().get(0, 0) == transform.mat.row_copy(0)[0]
    assert inv.rot().get(1, 1) == transform.mat.row_copy(1)[1]
    assert inv.rot().get(2, 2) == transform.mat.row_copy(2)[2]
    assert inv.trn()[0] == transform.vec.x
    assert inv.trn()[1] == transform.vec.y
    assert inv.trn()[2] == transform.vec.z
