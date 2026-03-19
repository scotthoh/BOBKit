import faulthandler
import pytest
import bobkit.clipper as clipper  # noqa: E402
import bobkit.data as data
import numpy as np
faulthandler.enable()


@pytest.fixture(scope="module")
def core_object():
    return Clipper_Test_Core()


def test_test_core(core_object):
    core_object.run_test()


class Clipper_Test_Core():
    '''
        Reimplementation of clipper's test_core in python
        Just to test if things work the same
    '''
    def __init__(self):
        self.error_count = 0
        self.data = data.Test_data()
        self.fsig = self.data.hkl_data_f_sigf
        self.abcd = self.data.hkl_data_abcd
        self.spgr = self.fsig.hkl_info.spacegroup
        self.cell = self.fsig.hkl_info.cell
        self.reso = self.fsig.hkl_info.resolution

    def run_test(self):
        self.test_nanfuncs()
        self.test_hkldata()
        self.test_resolord()
        self.test_spacegroup()
        self.test_atomshapefn()
        return (self.error_count == 0)

    def compare_values(self, id: str, val1: float, val2: float, tol: float = 0.01):
        equal = False
        if np.fabs(val1 - val2) < tol:
            equal = True
        else:
            print(f"Test: {id}, Result: {(val2-tol):.6f} < {val1:.6f} < {(val2+tol):.6f}")
        if not equal:
            self.error_count += 1
        return equal

    def test_nanfuncs(self):
        # test Nan functions
        assert clipper.Util.is_nan(float('NaN'))

        for i in np.arange(-30.0, 30.0, 1.0, dtype=float):
            assert not clipper.Util.is_nan(np.exp(i))

    def test_hkldata(self):
        # HRI = clipper.HKL_info.HKL_reference_index
        # need to check compute
        # use hkl_data_ABCD_float compute ABCD from phifom
        ##abcd1 = clipper.HKL_data_ABCD_float(abcd.hkl_info())
        ##abcd1.compute_from_phifom()
        ##ih = abcd.first()
        ##if ih != abcd.last():
        ##    cap = clipper.

        # test map calculation
        pw = clipper.HKL_data_Phi_fom_float(self.fsig.hkl_info)
        fp1 = clipper.HKL_data_F_phi_float(self.fsig.hkl_info)
        fp2 = clipper.HKL_data_F_phi_float(self.fsig.hkl_info)
        pw.compute_from_abcd(self.abcd)
        fp1.compute_from_fsigf_phifom(self.fsig, pw)
        grid = clipper.Grid_sampling(self.spgr, self.cell, self.reso)
        xmap = clipper.Xmap_float(self.spgr, self.cell, grid)
        xmap.fft_from(fp1)
        xmap.fft_to(fp2)
        ih = fp1.first()
        while not ih.last():
            ab1 = complex(0.0, 0.0)
            ab2 = complex(0.0, 0.0)
            if not fp1[ih].missing():
                ab1 = fp1[ih].complex
            if not fp2[ih].missing():
                ab2 = fp2[ih].complex
            self.compare_values("FFT-A", ab1.real, ab2.real, 0.01)
            self.compare_values("FFT-B", ab1.imag, ab2.imag, 0.01)
            ih.next()

    def test_resolord(self):
        # test resolution ordinals
        resols = []
        resinv = []
        ih = self.fsig.first()
        while not ih.last():
            resols.append(ih.invresolsq())
            resinv.append(ih.invresolsq())
            ih.next()
        ordinal = clipper.Generic_ordinal()
        ordinal.init(resols)
        ordinv = clipper.Generic_ordinal()
        ordinv.init(resinv)
        ordinv.invert()
        resord = clipper.Resolution_ordinal(self.fsig.hkl_info, 1.0)
        # resord.init(fsig.hkl_info(), 1.0)
        sorted(resols)
        for i in range(0, len(resols)):
            self.compare_values(
                "RESORD",
                resols[i],
                ordinv.ordinal(resord.ordinal(resols[i])),
                0.001,
            )

        # test resolution functions
        param = [1.0] * 10
        basisfn = clipper.BasisFn_spline(self.fsig, 10)
        targetfn = clipper.TargetFn_meanFnth_F_sigF_float(self.fsig, 2.0)
        rfn = clipper.ResolutionFn(self.fsig.hkl_info, basisfn, targetfn, param)
        self.compare_values("RESFN0", rfn.params[0], 229690.9746, 0.1)
        self.compare_values("RESFN1", rfn.params[1], 216481.7609, 0.1)
        self.compare_values("RESFN2", rfn.params[2],  78484.9498, 0.1)
        self.compare_values("RESFN3", rfn.params[3], 148774.2654, 0.1)
        self.compare_values("RESFN4", rfn.params[4],  69255.6000, 0.1)
        self.compare_values("RESFN5", rfn.params[5], 143032.5088, 0.1)
        self.compare_values("RESFN6", rfn.params[6], 110371.3524, 0.1)
        self.compare_values("RESFN7", rfn.params[7], 108711.3487, 0.1)
        self.compare_values("RESFN8", rfn.params[8], 150487.5496, 0.1)
        self.compare_values("RESFN9", rfn.params[9], 141713.7420, 0.1)

    def test_atomshapefn(self):
        # test atom shape function derivatives
        d = 0.001
        clipper.ScatteringFactors.select_SFType(
            clipper.ScatteringFactorsType.SF_WAASMAIER_KIRFEL
        )
        co = clipper.Coord_orth(1.0, 2.0, 3.0)
        sf = clipper.AtomShapeFn(co, "N", 0.25, 1.0)
        sfx = clipper.AtomShapeFn(co + clipper.Coord_orth(d, 0.0, 0.0), "N", 0.25, 1.0)
        sfy = clipper.AtomShapeFn(co + clipper.Coord_orth(0.0, d, 0.0), "N", 0.25, 1.0)
        sfz = clipper.AtomShapeFn(co + clipper.Coord_orth(0.0, 0.0, d), "N", 0.25, 1.0)
        sfo = clipper.AtomShapeFn(co, "N", 0.25, 1.0 + d)
        sfu = clipper.AtomShapeFn(co, "N", 0.251, 1.0)
        sfx2 = clipper.AtomShapeFn(
            co + clipper.Coord_orth(2.0 * d, 0.0, 0.0), "N", 0.25, 1.0
        )
        sfy2 = clipper.AtomShapeFn(
            co + clipper.Coord_orth(0.0, 2.0 * d, 0.0), "N", 0.25, 1.0
        )
        sfz2 = clipper.AtomShapeFn(
            co + clipper.Coord_orth(0.0, 0.0, 2.0 * d), "N", 0.25, 1.0
        )
        sfo2 = clipper.AtomShapeFn(co, "N", 0.25, 1.0 + d + d)
        sfu2 = clipper.AtomShapeFn(co, "N", 0.252, 1.0)
        sfxy = clipper.AtomShapeFn(co + clipper.Coord_orth(d, d, 0.0), "N", 0.25, 1.0)
        sfyz = clipper.AtomShapeFn(co + clipper.Coord_orth(0.0, d, d), "N", 0.25, 1.0)
        sfzx = clipper.AtomShapeFn(co + clipper.Coord_orth(d, 0.0, d), "N", 0.25, 1.0)
        uai = clipper.U_aniso_orth(0.25, 0.25, 0.25, 0.0, 0.0, 0.0)
        sfuai = clipper.AtomShapeFn(co, "N", uai, 1.0)
        params = []
        params.append(clipper.AtomShapeFn.X)
        params.append(clipper.AtomShapeFn.Y)
        params.append(clipper.AtomShapeFn.Z)
        params.append(clipper.AtomShapeFn.Occ)
        params.append(clipper.AtomShapeFn.Uiso)
        params.append(clipper.AtomShapeFn.U11)
        params.append(clipper.AtomShapeFn.U22)
        params.append(clipper.AtomShapeFn.U33)
        params.append(clipper.AtomShapeFn.U12)
        params.append(clipper.AtomShapeFn.U13)
        params.append(clipper.AtomShapeFn.U23)
        sf.agarwal_params = params

        for i in range(0, 100):
            c2 = co + clipper.Coord_orth(0.11 * (i % 5-1.9),
                                         0.13 * (i % 7-2.8),
                                         0.15 * (i % 9-3.7))
            curv = clipper.Matrix_double(len(sf.agarwal_params), len(sf.agarwal_params))
            self.compare_values("ATOMSF-A", sfuai.rho(c2), sf.rho(c2), 1.0e-8)
            rho, grad = sf.rho_curv(c2, curv)
            self.compare_values("ATOMSF-G", (sfx.rho(c2)-sf.rho(c2))/d, grad[0], 0.01)
            self.compare_values("ATOMSF-G", (sfy.rho(c2)-sf.rho(c2))/d, grad[1], 0.01)
            self.compare_values("ATOMSF-G", (sfz.rho(c2)-sf.rho(c2))/d, grad[2], 0.01)
            self.compare_values("ATOMSF-G", (sfo.rho(c2)-sf.rho(c2))/d, grad[3], 0.01)
            self.compare_values("ATOMSF-G", (sfu.rho(c2)-sf.rho(c2))/d, grad[4], 0.05)
            self.compare_values("ATOMSF-C", (sfx2.rho(c2)-2*sfx.rho(c2)+sf.rho(c2))/(d * d), curv.get(0, 0), 0.1)
            self.compare_values("ATOMSF-C", (sfy2.rho(c2)-2*sfy.rho(c2)+sf.rho(c2))/(d * d), curv.get(1, 1), 0.1)
            self.compare_values("ATOMSF-C", (sfz2.rho(c2)-2*sfz.rho(c2)+sf.rho(c2))/(d * d), curv.get(2, 2), 0.1)
            self.compare_values("ATOMSF-C", (sfxy.rho(c2)-sfx.rho(c2)-sfy.rho(c2) + sf.rho(c2))/(d*d), curv.get(0, 1), 0.1)
            self.compare_values("ATOMSF-C", (sfyz.rho(c2)-sfy.rho(c2)-sfz.rho(c2) + sf.rho(c2))/(d*d), curv.get(1, 2), 0.1)
            self.compare_values("ATOMSF-C", (sfzx.rho(c2)-sfz.rho(c2)-sfx.rho(c2) + sf.rho(c2))/(d*d), curv.get(2, 0), 0.1)

        for j in range(0, 20):
            x = 0.19 * (j % 5 - 1.9)
            y = 0.15 * (j % 7 - 2.8)
            z = 0.13 * (j % 9 - 3.7)
            ua = clipper.U_aniso_orth(
                x * x + 0.2, y * y + 0.2, z * z + 0.2, y * z, z * x, x * y
            )
            ua00 = clipper.U_aniso_orth(
                ua.m00 + d, ua.m11, ua.m22, ua.m01, ua.m02, ua.m12
            )
            ua11 = clipper.U_aniso_orth(
                ua.m00, ua.m11 + d, ua.m22, ua.m01, ua.m02, ua.m12
            )
            ua22 = clipper.U_aniso_orth(
                ua.m00, ua.m11, ua.m22 + d, ua.m01, ua.m02, ua.m12
            )
            ua01 = clipper.U_aniso_orth(
                ua.m00, ua.m11, ua.m22, ua.m01 + d, ua.m02, ua.m12
            )
            ua02 = clipper.U_aniso_orth(
                ua.m00, ua.m11, ua.m22, ua.m01, ua.m02 + d, ua.m12
            )
            ua12 = clipper.U_aniso_orth(
                ua.m00, ua.m11, ua.m22, ua.m01, ua.m02, ua.m12 + d
            )
            sfua = clipper.AtomShapeFn(co, "N", ua, 1.0)
            sfuax = clipper.AtomShapeFn(
                co + clipper.Coord_orth(d, 0.0, 0.0), "N", ua, 1.0
            )
            sfuay = clipper.AtomShapeFn(
                co + clipper.Coord_orth(0.0, d, 0.0), "N", ua, 1.0
            )
            sfuaz = clipper.AtomShapeFn(
                co + clipper.Coord_orth(0.0, 0.0, d), "N", ua, 1.0
            )
            sfuao = clipper.AtomShapeFn(co, "N", ua, 1.0 + d)
            sfua00 = clipper.AtomShapeFn(co, "N", ua00, 1.0)
            sfua11 = clipper.AtomShapeFn(co, "N", ua11, 1.0)
            sfua22 = clipper.AtomShapeFn(co, "N", ua22, 1.0)
            sfua01 = clipper.AtomShapeFn(co, "N", ua01, 1.0)
            sfua02 = clipper.AtomShapeFn(co, "N", ua02, 1.0)
            sfua12 = clipper.AtomShapeFn(co, "N", ua12, 1.0)
            sfua.agarwal_params = params
            for i in range(0, 50):
                c2 = co + clipper.Coord_orth(
                    0.11 * (i % 5 - 1.9), 0.13 * (i % 7 - 2.8), 0.15 * (i % 9 - 3.7)
                )
                # c = clipper.Matrix_double(11, 11)
                rho, g = sfua.rho_grad(c2)
                self.compare_values(
                    "ATOMSF-AG", (sfuax.rho(c2) - sfua.rho(c2)) / d, g[0], 0.01
                )
                self.compare_values(
                    "ATOMSF-AG", (sfuay.rho(c2) - sfua.rho(c2)) / d, g[1], 0.01
                )
                self.compare_values(
                    "ATOMSF-AG", (sfuaz.rho(c2) - sfua.rho(c2)) / d, g[2], 0.01
                )
                self.compare_values(
                    "ATOMSF-AG", (sfuao.rho(c2) - sfua.rho(c2)) / d, g[3], 0.01
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua00.rho(c2) - sfua.rho(c2)) / d, g[5], 0.05
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua11.rho(c2) - sfua.rho(c2)) / d, g[6], 0.05
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua22.rho(c2) - sfua.rho(c2)) / d, g[7], 0.05
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua01.rho(c2) - sfua.rho(c2)) / d, g[8], 0.05
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua02.rho(c2) - sfua.rho(c2)) / d, g[9], 0.05
                )
                self.compare_values(
                    "ATOMSF-AG", (sfua12.rho(c2) - sfua.rho(c2)) / d, g[10], 0.05
                )

    def test_spacegroup(self):
        # test spacegroups
        pgs = [
            "-P 1",
            "-P 2",
            "-P 2y",
            "-P 2x",
            '-P 2"',
            '-P 2y"',
            '-P 2x"',
            "-P 2'",
            "-P2y'",
            "-P 2x'",
            "-P 2 2",
            '-P 2 2"',
            '-P 2 2"(y,z,x)',
            '-P 2 2"(z,x,y)',
            "-P3",
            "-P 3 (y,z,x)",
            "-P 3 (z,x,y)",
            "-P 3 (-x,y,z)",
            "-P 3 (y,z,-x)",
            "-P 3 (z,-x,y)",
            "-P 3*",
            "-P 3* (-x,y,z)",
            "-P 3* (x,-y,z)",
            "-P 3* (x,y,-z)",
            "-P 3 2",
            "-P 3 2 (y,z,x)",
            "-P 3 2 (z,x,y)",
            "-P 3* 2",
            "-P 3* 2 (-x,y,z)",
            "-P 3* 2 (x,-y,z)",
            "-P 3* 2 (-x,-y,z)",
            '-P 3 2"',
            '-P 3 2"(z,x,y)',
            '-P 3 2"(y,z,x)',
            '-P 3 2"(-x,y,z)',
            '-P 3 2"(z,-x,y)',
            '-P 3 2"(y,z,-x)',
            "-P 4",
            "-P 4 (y,z,x)",
            "-P 4 (z,x,y)",
            "-P 4 2",
            "-P 4 2 (y,z,x)",
            "-P 4 2 (z,x,y)",
            "-P 6",
            "-P 6 (y,z,x)",
            "-P 6 (z,x,y)",
            "-P 6 2",
            "-P 6 2 (y,z,x)",
            "-P 6 2 (z,x,y)",
            "-P 2 2 3",
            "-P 4 2 3",
        ]
        hallsymbols = []
        for sgd in data.spacegroup_datatable():
            hallsymbols.append(sgd.hall)
        for i in range(0, len(pgs)):
            hallsymbols.append(pgs[i])
        cellc = clipper.Cell(clipper.Cell_descr(37, 37, 37))
        cellha = clipper.Cell(clipper.Cell_descr(37, 37, 37, 120, 90, 90))
        cellhb = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 120, 90))
        cellhc = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 90, 120))
        cellha1 = clipper.Cell(clipper.Cell_descr(37, 37, 37, 60, 90, 90))
        cellhb1 = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 60, 90))
        cellhc1 = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 90, 60))
        for s in range(0, len(hallsymbols)):
            try:
                symbol = hallsymbols[s]
                sg = clipper.Spacegroup(
                    clipper.Spgr_descr(symbol, clipper.Spgr_descr.Hall)
                )
                # identify trigonal/hexagonal groups
                cell = clipper.Cell(cellc.description)
                for sym in range(1, sg.num_symops):
                    if (sg.symop(sym).rot.get(1, 1) * sg.symop(sym).rot.get(1, 2) == -1) or (
                        sg.symop(sym).rot.get(2, 1) * sg.symop(sym).rot.get(2, 2) == -1
                    ):
                        cell = cellha
                    if (sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 2) == -1) or (
                        sg.symop(sym).rot.get(2, 0) * sg.symop(sym).rot.get(2, 2) == -1
                    ):
                        cell = cellhb
                    if (sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 1) == -1) or (
                        sg.symop(sym).rot.get(1, 0) * sg.symop(sym).rot.get(1, 1) == -1
                    ):
                        cell = cellhc
                    if (sg.symop(sym).rot.get(1, 1) * sg.symop(sym).rot.get(1, 2) == 1) or (
                        sg.symop(sym).rot.get(2, 1) * sg.symop(sym).rot.get(2, 2) == 1
                    ):
                        cell = cellha1
                    if (sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 2) == 1) or (
                        sg.symop(sym).rot.get(2, 0) * sg.symop(sym).rot.get(2, 2) == 1
                    ):
                        cell = cellhb1
                    if (sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 1) == 1) or (
                        sg.symop(sym).rot.get(1, 0) * sg.symop(sym).rot.get(1, 1) == 1
                    ):
                        cell = cellhc1
                for i in range(0, 100):
                    rfl = clipper.HKL(i % 5-2, i % 7-3, i % 9-4)
                    s0 = rfl.invresolsq(cell)
                    for sym in range(1, sg.num_symops):
                        s1 = rfl.transform(sg.symop(sym)).invresolsq(cell)
                        self.compare_values("SG_"+symbol, s0, s1, 1.0e-12)
                grid = clipper.Grid_sampling(sg, cell, clipper.Resolution(4.0))
                xmap = clipper.Xmap_float(sg, cell, grid)
            except Exception:
                self.compare_values("SG_"+symbol, sg.spacegroup_number(), -1)

        # test rotation
        for x in np.arange(-1.0, 1.0, 0.02, dtype=float):
            for y in np.arange(-1.0, 1.0, 0.02, dtype=float):
                for z in np.arange(-1.0, 1.0, 0.02, dtype=float):
                    s = x*x + y*y + z*z
                    if s < 1.0:
                        w = np.sqrt(1.0 - s)
                        rot = clipper.Rotation(w, x, y, z)
                        rotinv = rot.inverse()
                        r1 = clipper.Rotation(rot.matrix())
                        r2 = clipper.Rotation(rot.polar_ccp4())
                        self.compare_values(
                            "ROT/MAT " + rot.format(),
                            (rotinv * r1).abs_angle(),
                            0.0,
                            1.0e-6,
                        )
                        self.compare_values(
                            "ROT/POL " + rot.format(),
                            (rotinv * r2).abs_angle(),
                            0.0,
                            1.0e-6,
                        )


if __name__ == "__main__":
    from datetime import datetime
    start = datetime.now()
    runtest = Clipper_Test_Core()
    print("Test core python implementations : ", end="")
    if runtest.run_test():
        print("Pass\n")
    else:
        print("Fail\n")
    end = datetime.now()
    print(f"Test time : {(end-start).microseconds/1e6:.4f} secs")