import faulthandler
import pytest
import pathlib
import os
import bobkit.clipper as clipper
import bobkit.data as data
import numpy as np
from datetime import datetime

try:
    from test_data.contrib_test_data import contrib_vals, contrib_tols
except ImportError:
    print("Did not managed to get contrib data")
faulthandler.enable()


@pytest.fixture(scope="module")
def data_fsig_obj():
    return data.Test_data().hkl_data_f_sigf


@pytest.fixture(scope="module")
def data_abcd_obj():
    return data.Test_data().hkl_data_abcd


@pytest.fixture(scope="module")
def data_atomlist():
    return data.Test_data().atom_list


@pytest.fixture(scope="module")
def data_contribs(request):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    #cwd = os.path.dirname(__file__)
    contrib_datafile = str(testdir / "contrib_test_data.py")  # os.path.join(cwd, "test_data/contrib_test_data.py")
    globals_dict = {}
    with open(contrib_datafile, "r", encoding="utf-8") as f:
        src = f.read()
    exec(compile(src, "err_data_contrib.txt", "exec"), globals_dict)
    contrib_vals = globals_dict["contrib_vals"]
    contrib_tols = globals_dict["contrib_tols"]
    return (contrib_vals, contrib_tols)


@pytest.fixture(scope="module")
def datacontrib_vals(data_contribs):
    return data_contribs[0]


@pytest.fixture(scope="module")
def datacontrib_tols(data_contribs):
    return data_contribs[1]


class Test_Contrib():
    '''
        Reimplementation of clipper's test_contrib in python
        Just to test if things work the same
    '''
    def test_load_data(
        self,
        data_fsig_obj,
        data_abcd_obj,
        data_atomlist,
        data_contribs
    ):
        contrib_vals, contrib_tols = data_contribs
        assert len(contrib_vals) > 0
        assert len(contrib_tols) > 0
        fsig = data_fsig_obj
        abcd = data_abcd_obj
        atoms = data_atomlist
        spgr = fsig.hkl_info.spacegroup
        cell = fsig.hkl_info.cell
        reso = fsig.hkl_info.resolution
        assert spgr.symbol_hall() == "P 2ac 2ab"
        assert cell.a == 64.8970
        assert cell.b == 78.3230
        assert cell.c == 38.7920
        reso = fsig.hkl_info.resolution
        assert reso.limit() == (5.0 - 1.0e-3)
        assert fsig.num_obs() > 0
        assert abcd.num_obs() > 0
        assert len(atoms) > 0

    def run_test(self):
        ''' For when using python3 test_clipper_contrib.py
        '''
        self.error_count = 0
        self.count = 0
        fsig = data.Test_data().hkl_data_f_sigf
        abcd = data.Test_data().hkl_data_abcd
        atomlist = data.Test_data().atom_list
        self.test_sfcalc_objects()
        self.test_sfweight(fsig, atomlist, contrib_vals, contrib_tols)
        self.test_map_filters(fsig, abcd)
        self.test_convolution_fffear(fsig, abcd)
        self.test_anisotropic_scaling(fsig)
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

    def compare_value(self, id: str, val1):
        equal = False

        if self.count < len(self.contrib_vals):
            v2 = self.contrib_vals[self.count]
            tol = self.contrib_tols[self.count]
            if np.fabs(val1 - v2) < tol:
                equal = True
            else:
                print(f"Test: {id} {self.count}, Result: {(v2-tol):.6f} < {val1:.6f} < {(v2+tol):.6f}")
        elif self.count == len(self.contrib_vals):
            print("Self test: out of data")
        self.count += 1
        if not equal:
            self.error_count += 1
        return equal

    def test_sfcalc_objects(self):
        # select spacegroup
        hallsymbols = []
        count = 0
        for sgd in data.spacegroup_datatable():
            if count > 0: count += 1
            if count % 10 == 0:
                hallsymbols.append(sgd.hall)
            count += 10

        # build model
        atoms = clipper.Atom_list()
        atom = clipper.Atom.null()
        atom.occupancy = 1.0
        atom.u_iso = 0.5
        atom.element = "C"
        atom.pos = clipper.Coord_orth(12, 8, 5)
        atoms.append(atom)
        atom.element = "N"
        atom.pos = clipper.Coord_orth(11, 6, 4)
        atoms.append(atom)
        atom.element = "O"
        atom.pos = clipper.Coord_orth(13, 5, 5)
        atoms.append(atom)
        # calc cell
        cellc = clipper.Cell(clipper.Cell_descr(37, 37, 37))
        cellha = clipper.Cell(clipper.Cell_descr(37, 37, 37, 120, 90, 90))
        cellhb = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 120, 90))
        cellhc = clipper.Cell(clipper.Cell_descr(37, 37, 37, 90, 90, 120))
        for s in range(0, len(hallsymbols)):
            try:
                symbol = hallsymbols[s]
                sg = clipper.Spacegroup(clipper.Spgr_descr(symbol, clipper.Spgr_descr.Hall))
                # identify trifonal/hexagonal groups
                cg = clipper.Cell(cellc.description)
                for sym in range(1, sg.num_symops):
                    if ((sg.symop(sym).rot.get(1, 1) * sg.symop(sym).rot.get(1, 2) == -1) or
                       (sg.symop(sym).rot.get(2, 1) * sg.symop(sym).rot.get(2, 2) == -1)):
                        cg = clipper.Cell(cellha.description)
                    if ((sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 2) == -1) or
                       (sg.symop(sym).rot.get(2, 0) * sg.symop(sym).rot.get(2, 2) == -1)):
                        cg = clipper.Cell(cellhb.description)
                    if ((sg.symop(sym).rot.get(0, 0) * sg.symop(sym).rot.get(0, 1) == -1) or
                       (sg.symop(sym).rot.get(1, 0) * sg.symop(sym).rot.get(1, 1) == -1)):
                        cg = clipper.Cell(cellhc.description)
                hklinfo = clipper.HKL_info(sg, cg, 5.0, True)
                fp1 = clipper.HKL_data_F_phi_float(hklinfo)
                fp2 = clipper.HKL_data_F_phi_float(hklinfo)
                clipper.SFcalc_iso_sum_float(fp1, atoms)
                clipper.SFcalc_iso_fft_float(fp2, atoms, 2.5, 2.5, 0.25)
                # need to bind fftmap.h if want to try extra fft tests
                tol = 0.005 * fp1[clipper.HKL(0, 0, 0)].f
                ih = fp1.first()
                if not ih.last():
                    ab1 = complex(0.0, 0.0)
                    ab2 = complex(0.0, 0.0)
                    if not fp1[ih].missing():
                        ab1 = fp1[ih].complex
                    if not fp2[ih].missing():
                        ab2 = fp2[ih].complex
                    self.compare_values("SF-A", ab1.real, ab2.real, tol)
                    self.compare_values("SF-B", ab1.imag, ab2.imag, tol)
            except Exception:
                self.compare_values("SFSG "+symbol, sg.spacegroup_number(), -1)

    def test_sfweight(self, data_fsig_obj, data_atomlist, datacontrib_vals, datacontrib_tols):
        # test sfcalc_obs and sfweight
        # sf calc obj
        self.contrib_vals = datacontrib_vals
        self.contrib_tols = datacontrib_tols
        fsig = data_fsig_obj
        atomlist = data_atomlist
        self.count = 0
        fcal = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        sfcb = clipper.SFcalc_obs_bulk_float()
        sfcb(fcal, fsig, atomlist)
        # sfcalc results
        ih = fcal.first()
        while not ih.last():
            if ih.index % 20 == 0:
                ab = fcal[ih].complex
                self.compare_value("SFO-A", ab.real)
                self.compare_value("SFO-B", ab.imag)
            ih.next()

        fb1 = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        fb2 = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        fd1 = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        fd2 = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        phiw1 = clipper.HKL_data_Phi_fom_float(fsig.hkl_info)
        phiw2 = clipper.HKL_data_Phi_fom_float(fsig.hkl_info)
        flag = clipper.HKL_data_Flag(fsig.hkl_info)
        abcd = clipper.HKL_data_ABCD_float(fsig.hkl_info)
        abcd2 = clipper.HKL_data_ABCD_float(fsig.hkl_info)
        abcd.set_all_values_to(clipper.ABCD_float(0., 0., 0., 0.,))
        ih = flag.first()
        while not ih.last():
            if not fsig[ih].missing():
                flag[ih].flag = clipper.SFweight_spline_float.BOTH.value
            else:
                flag[ih].flag = clipper.SFweight_spline_float.NONE.value
            ih.next()

        sfw = clipper.SFweight_spline_float(600, 12)
        fl = sfw(fb1, fd1, phiw1, fsig, fcal, flag)
        if not fl:
            self.compare_values("SFW-FAIL", 1, 0)

        # sfweight results
        ih = phiw1.first()
        while not ih.last():
            if not fsig[ih].missing():
                if ih.index % 20 == 0:
                    ab_b = fb1[ih].complex
                    ab_d = fd1[ih].complex
                    self.compare_value("SFWB_A", ab_b.real)
                    self.compare_value("SFWB_B", ab_b.imag)
                    self.compare_value("SFWD_A", ab_d.real)
                    self.compare_value("SFWD_B", ab_d.imag)
                    self.compare_value("SFW-W", phiw1[ih].fom)
            ih.next()

        # sfweight-hl results
        fsig0 = clipper.HKL_data_F_sigF_float(fsig.hkl_info)
        ih = fsig.first()
        while not ih.last():
            if not fsig[ih].missing():
                fsig0[ih] = clipper.F_sigF_float(fsig[ih].f, 0.0)
            ih.next()
        sfw1 = clipper.SFweight_spline_float(600, 12)
        sfw2 = clipper.SFweight_spline_float(600, 12)
        fl1 = sfw1(fb1, fd1, phiw1, fsig0, fcal, flag)
        if not fl1:
            self.compare_values("SFW-FAIL1", 1, 0)
        fl2 = sfw2(fb2, fd2, phiw2, abcd2, fsig0, abcd, fcal, flag)
        if not fl2:
            self.compare_values("SFW-FAIL2", 1, 0)
        params_e1 = sfw1.params_error
        params_s1 = sfw1.params_scale
        params_e2 = sfw2.params_error
        params_s2 = sfw2.params_scale
        for i in range(0, len(params_e1)):
            self.compare_values("SFWHL-E", params_e1[i], params_e2[i], 0.025)
            self.compare_values("SFWHL-S", params_s1[i], params_s2[i], 0.025)
        ih = fsig0.first()
        while not ih.last():
            if not fsig0[ih].missing():
                r00 = clipper.SFweight_spline_float.TargetResult()
                r01 = clipper.SFweight_spline_float.TargetResult()
                r10 = clipper.SFweight_spline_float.TargetResult()
                r11 = clipper.SFweight_spline_float.TargetResult()
                rhl = clipper.SFweight_spline_float.TargetResult()
                for s in np.arange(0.20, 1.01, 0.2, dtype=float):
                    for p in np.arange(0.20, 1.01, 0.2, dtype=float):
                        w = p * fsig0[ih].f
                        w = w * w
                        d = 0.000001
                        r00 = sfw1.targetfn(ih.hkl_class(), fsig0[ih], fcal[ih], s, w)
                        rhl = sfw1.targethl(ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w)
                        self.compare_values("SFW-TGT-CMP", r00.r, rhl.r, 1.0)
                        r00 = sfw1.targetfn(ih.hkl_class(), fsig0[ih], fcal[ih], s, w)
                        r01 = sfw1.targetfn(
                            ih.hkl_class(), fsig0[ih], fcal[ih], s, w + d
                        )
                        r10 = sfw1.targetfn(
                            ih.hkl_class(), fsig0[ih], fcal[ih], s + d, w
                        )
                        r11 = sfw1.targetfn(
                            ih.hkl_class(), fsig0[ih], fcal[ih], s + d, w + d
                        )
                        self.compare_values("SFW-FN-DW", (r01.r-r00.r)/d, r00.dw, 0.02)
                        self.compare_values("SFW-FN-DS", (r10.r-r00.r)/d, r00.ds, 0.02)
                        self.compare_values("SFW-FN-DWW", (r01.dw-r00.dw)/d, r00.dww, 0.02)
                        self.compare_values("SFW-FN-DSS", (r10.ds-r00.ds)/d, r00.dss, 0.02)
                        r00 = sfw1.targethl(ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w )
                        r01 = sfw1.targethl(ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w+d )
                        r10 = sfw1.targethl(ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s+d, w )
                        r11 = sfw1.targethl(ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s+d, w+d )
                        self.compare_values("SFW-HL-DW ", (r01.r-r00.r)/d, r00.dw, 0.02 )
                        self.compare_values("SFW-HL-DS ", (r10.r-r00.r)/d, r00.ds, 0.02 )
                        self.compare_values("SFW-HL-DWW", (r01.dw-r00.dw)/d, r00.dww, 0.02 )
                        self.compare_values("SFW-HL-DSS", (r10.ds-r00.ds)/d, r00.dss, 0.02 )
            ih.next()

    def test_map_filters(self, data_fsig_obj, data_abcd_obj, debug: bool = True):
        fsig = data_fsig_obj
        abcd = data_abcd_obj
        pw = clipper.HKL_data_Phi_fom_float(fsig.hkl_info)
        fp = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        pw.compute_from_abcd(abcd)
        fp.compute_from_fsigf_phifom(fsig, pw)
        grid = clipper.Grid_sampling(fsig.hkl_info.spacegroup, fsig.hkl_info.cell, fsig.hkl_info.resolution, 2.5)
        xmap = clipper.Xmap_float(fsig.hkl_info.spacegroup, fsig.hkl_info.cell, grid)
        step = clipper.MapFilterFn_step(2.5)
        fltr1 = clipper.MapFilter_slow_float(
            step, 1.0, clipper.MapFilter_slow_float.Relative
        )
        fltr2 = clipper.MapFilter_fft_float(
            step, 1.0, clipper.MapFilter_fft_float.Relative
        )
        f1 = clipper.Xmap_float()
        f2 = clipper.Xmap_float()
        start = datetime.now()
        fltr1(f1, xmap)
        end = datetime.now()
        start2 = datetime.now()
        fltr2(f2, xmap)
        end2 = datetime.now()
        if debug:
            print(f"Slow filter : {(end-start).microseconds} microsecs")
            print(f"FFT filter : {(end2-start2).microseconds} microsecs")

    def test_convolution_fffear(self, data_fsig_obj, data_abcd_obj):
        fsig = data_fsig_obj
        abcd = data_abcd_obj
        pw = clipper.HKL_data_Phi_fom_float(fsig.hkl_info)
        fp = clipper.HKL_data_F_phi_float(fsig.hkl_info)
        pw.compute_from_abcd(abcd)
        fp.compute_from_fsigf_phifom(fsig, pw)
        grid = clipper.Grid_sampling(fsig.hkl_info.spacegroup, fsig.hkl_info.cell, fsig.hkl_info.resolution, 2.5)
        xmap = clipper.Xmap_float(fsig.hkl_info.spacegroup, fsig.hkl_info.cell, grid)
        xmap.fft_from(fp)
        r1 = clipper.Xmap_float(clipper.Spacegroup.p1(), xmap.cell, xmap.grid_sampling)
        r2 = clipper.Xmap_float(clipper.Spacegroup.p1(), xmap.cell, xmap.grid_sampling)
        irad = 3
        tg = clipper.Grid_range(clipper.Coord_grid(-irad, -irad, -irad),
                                clipper.Coord_grid(irad, irad, irad))
        target = clipper.NXmap_float(xmap.cell, xmap.grid_sampling, tg)
        weight = clipper.NXmap_float(xmap.cell, xmap.grid_sampling, tg)
        target.fill_map_with(0.)
        weight.fill_map_with(0.)
        c = tg.min
        while not c.last(tg):
            if c * c <= 5:
                target.set_data(c - tg.min, xmap.get_data(c))
                weight.set_data(c - tg.min, 1.0)
            c.next(tg)

        # convolution search not exposed to nanobind for now
        #conv1 = clipper.Convolution_search_slow_float(xmap)
        # conv2 = clipper.Convolution_search_fft_float(xmap)
        # #conv1(r1, target)
        # conv2(r2, target)
        # #ix = r1.first()
        # #while not ix.last():
        # #    self.compare_values("CONVOL", r1[ix], r2[ix], 0.0001)
        # #    ix.next()
        # #srch1 = clipper.FFFear_slow_float(xmap)
        # srch2 = clipper.FFFear_fft_float(xmap)
        # #srch1(r1, target, weight)
        # srch2(r2, target, weight)
        # ix = r1.first()
        # with open("fffear_output_test.log", "w") as wopen:
        #     while not ix.last():
        #     #   self.compare_values("FFFEAR", r1[ix], r2[ix], 0.001)
        #         wopen.write(f"{r2[ix]:.3f}/n")
        #         ix.next()

    def test_anisotropic_scaling(self, data_fsig_obj):
        # something wrong here not passing test
        # expand to P1
        fsig = data_fsig_obj
        spgrp1 = clipper.Spacegroup.p1()
        hkl1 = clipper.HKL_info(spgrp1, fsig.hkl_info.cell, fsig.hkl_info.resolution, True)
        fs = clipper.HKL_data_F_sigF_float(hkl1)
        # print(fsig.hkl_info)
        # print(hkl1)
        ih = hkl1.first()
        while not ih.last():
            try:
                fs[ih] = fsig[ih.hkl()]
            # print(f"{fs[ih].f}, {fs[ih].sigf}, {fsig[ih.hkl()].f}, {fsig[ih.hkl()].sigf}")
            except IndexError:
                fs[ih].set_null()
            ih.next()
        # make data objects

        fp = clipper.HKL_data_F_phi_float(fs.hkl_info)
        # print(fp.hkl_info)
        fs1 = clipper.HKL_data_F_sigF_float()
        fs1.copy_from(fs)
        # print(fs1.hkl_info)
        fs2 = clipper.HKL_data_F_sigF_float()
        fs2.copy_from(fs)
        # print(fs2.hkl_info)
        u_ref = clipper.U_aniso_orth(0.10, 0.13, 0.17, -0.02, 0.03, -0.04)
        # simulate aniso data
        ih = fs.first()
        # found problem, the return value not refererence for __getitem__
        while not ih.last():
            if not fs[ih].missing():
                c = ih.hkl().coord_reci_orth(fsig.hkl_info.cell)
                s = np.exp(clipper.Util.twopi2() * u_ref.quad_form(c))
                f1tmp = fs[ih]
                f1tmp.scale(s)
                fs1[ih] = f1tmp
                f2tmp = fs2[ih]
                f2tmp.scale(1.0/s)
                fs2[ih] = f2tmp
                # fs1[ih].scale(s)
                # fs2[ih].scale(1.0/s)
                fp[ih] = clipper.F_phi_float(fs[ih].f, 0.0)
                # print(f"{fs1[ih].f:6f}, {fs1[ih].sigf:6f} | {fs2[ih].f:6f}, {fs2[ih].sigf:6f} | {fp[ih].f:6f} | {s:6f} | {c.format()}")
            ih.next()
        # now attempt scaling
        sfscl = clipper.SFscale_aniso_float()
        sfscl(fs1, fp)
        u_wrk1 = sfscl.u_aniso_orth(clipper.SFscale_aniso_float.F)
        sfscl(fp, fs2)
        u_wrk2 = sfscl.u_aniso_orth(clipper.SFscale_aniso_float.F)
        self.compare_values("ANISO-O-00", u_ref.m00, u_wrk1.m00, 1.0e-6)
        self.compare_values("ANISO-O-11", u_ref.m11, u_wrk1.m11, 1.0e-6)
        self.compare_values("ANISO-O-22", u_ref.m22, u_wrk1.m22, 1.0e-6)
        self.compare_values("ANISO-O-01", u_ref.m01, u_wrk1.m01, 1.0e-6)
        self.compare_values("ANISO-O-02", u_ref.m02, u_wrk1.m02, 1.0e-6)
        self.compare_values("ANISO-O-12", u_ref.m12, u_wrk1.m12, 1.0e-6)
        self.compare_values("ANISO-C-00", u_ref.m00, u_wrk2.m00, 1.0e-6)
        self.compare_values("ANISO-C-11", u_ref.m11, u_wrk2.m11, 1.0e-6)
        self.compare_values("ANISO-C-22", u_ref.m22, u_wrk2.m22, 1.0e-6)
        self.compare_values("ANISO-C-01", u_ref.m01, u_wrk2.m01, 1.0e-6)
        self.compare_values("ANISO-C-02", u_ref.m02, u_wrk2.m02, 1.0e-6)
        self.compare_values("ANISO-C-12", u_ref.m12, u_wrk2.m12, 1.0e-6)


if __name__ == "__main__":
    from datetime import datetime
    start0 = datetime.now()
    runtest = Test_Contrib()
    print("Test contrib python implementations : ", end="")
    if runtest.run_test():
        print("Pass\n")
    else:
        print("Fail\n")
    end0 = datetime.now()
    print(f"Test time : {(end0-start0).microseconds/1e6:.4f} secs")
