"""pytest for hkldata bindings
"""
import pathlib
import pytest
import bobkit.clipper as clipper
import gemmi
import cmath
import numpy as np

@pytest.fixture
def read_mtzfile(request):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mtz = gemmi.read_mtz_file(str(testdir / "7ins.mtz"))
    return mtz


class TestHKLData:
    def test___init__(self):
        fphif = clipper.HKL_data_F_phi_float()
        fphid = clipper.HKL_data_F_phi_double()
        abcdf = clipper.HKL_data_ABCD_float()
        abcdd = clipper.HKL_data_ABCD_double()
        fsigf = clipper.HKL_data_F_sigF_float()
        fsigfd = clipper.HKL_data_F_sigF_double()
        fsigfano = clipper.HKL_data_F_sigF_ano_float()
        fsigfanod = clipper.HKL_data_F_sigF_ano_double()
        isig = clipper.HKL_data_I_sigI_float()
        isigd = clipper.HKL_data_I_sigI_double()
        flag = clipper.HKL_data_Flag()
        flagb = clipper.HKL_data_Flag_bool()
        isigano = clipper.HKL_data_I_sigI_ano_float()
        isiganod = clipper.HKL_data_I_sigI_ano_double()
        phifom = clipper.HKL_data_Phi_fom_float()
        phifomd = clipper.HKL_data_Phi_fom_double()

    def test_hkldata(self, read_mtzfile):
        mtz = read_mtzfile
        miller = mtz.make_miller_array()
        print(f" miller size : {miller.shape}, {mtz.nreflections}")
        hklinfo = clipper.HKL_info.from_gemmi_mtz(mtz, 1e-7)
        #hkls1 = hklinfo.hkl_array
        print(f"num obs : {hklinfo.num_reflections()}" )
        #11788, [28  7 17], [28  8  0]
        print(f"cell {mtz.cell}, {hklinfo.cell}")
        hkl0 = clipper.HKL(28, 7, 17)
        resol = hkl0.invresolsq(hklinfo.cell)
        resog = clipper.Resolution(mtz.resolution_high()).invresolsq_limit()
        print(f"invresol {resol}, {resog}")
        #for i in range(0, len(miller)):
        #    print(f"{i}, {miller[i]}, {hkls1[i]}")
            #assert np.all(miller[i] == hkls1[i])
        fphidata = clipper.HKL_data_F_phi_float(hklinfo)
        data = mtz.get_f_phi("FWT", "PHWT", False)
        for d in data:
            hkl = clipper.HKL(d.hkl)
            fphidata[hkl] = clipper.F_phi_float(d.value)
        hkls = fphidata.hkl_array
        arr = fphidata.array
        print(len(arr), len(data))
        assert len(arr) == len(data)
        for i in range(0, len(data)):
            hkl = clipper.HKL(data[i].hkl)
            assert (fphidata[hkl].f - abs(data[i].value)) < 5.e-5
            assert fphidata[hkl].phi - cmath.phase(data[i].value) < 5.e-5
            if not clipper.Util.is_nan(arr[i][0]):
                assert arr[i][0] - abs(data[i].value) < 5.e-5
                assert arr[i][1] - cmath.phase(data[i].value) < 5.e-5
            assert np.all(hkls[i] == hkl)

        hklsamp = clipper.HKL_sampling(hklinfo.cell, clipper.Resolution(4.0))
        hklinfo2 = clipper.HKL_info()
        hklinfo2.init(hklinfo.spacegroup, hklinfo.cell, hklsamp)
        fpdata2 = fphidata.restrict_to(hklinfo2)
        ih = fpdata2.first_data()
        while not ih.last():
            v = clipper.F_phi_float()
            try:
                v = fphidata[ih.hkl()]
            except IndexError:
                print("Index Error")
            if not v.missing():
                assert fpdata2[ih].f == fphidata[ih.hkl()].f
                assert fpdata2[ih].phi == fphidata[ih.hkl()].phi
            fpdata2.next_data(ih)
        # for non missing data
        ih = fpdata2.first_data()
        while not ih.last():
            assert fpdata2[ih.hkl()].f == fphidata[ih.hkl()].f
            assert fpdata2[ih.hkl()].phi == fphidata[ih.hkl()].phi
            fpdata2.next_data(ih)

        # data export and import from list/vector/array
        ih = fpdata2.first_data()
        arr2 = fpdata2.data_export(ih.hkl())
        assert arr2[0] == fpdata2[ih].f
        assert arr2[1] == fpdata2[ih].phi
        fpdata3 = clipper.HKL_data_F_phi_float(hklinfo2)
        ih = fpdata3.first_data()
        b = []
        while not ih.last():
            a3 = np.array(fpdata2.data_export(ih.hkl())) #, dtype=np.float64)
            fpdata3.data_import(ih.hkl(), a3) #a3.tolist())
            b.append(a3)
            #assert arr2[0] == fpdata3[ih].f
            #assert arr2[1] == fpdata3[ih].phi
            assert fpdata2[ih].f == fpdata3[ih].f
        fsig = clipper.HKL_data_F_sigF_float(hklinfo)
        hkls = fphidata.hkl_array
        #print(f"len hkls : {len(hkls)}; {hklinfo.num_reflections()}")
#
        #for i in range(0, len(hkls)):
        #    hkli = data[i].hkl
        #    assert np.all(hkls[i] == hkli)
        #print(type(hkls), type(hkls[0]))
        #ba = np.array(b)
        #print(type(np.array(b)))
        #a = []
        #for i in hkls:
        #    a.append(clipper.HKL(i))
        ##hlist = hkls.tolist()
        #print(type(a), type(a[0]))
        #print(arr.size(), )
        fsig.data_import_hkls(hkls, arr)
        #
        #for i in range(0, hkls.shape[0]):
        #    assert np.all(hkls[i] == hklinfo2.hkl_of(i))
        #    assert np.all(hkls[i] == hklinfo2.hkl_of(i).array)
        #fsig.data_import_hkls(hkls.tolist(), a3)
        ih = fsig.first_data()
        while not ih.last():
            assert fphidata[ih].f == fsig[ih].f
            #assert arr[1] == fsig[ih].phi
            fsig.next_data(ih)
            print(f'{fsig[ih].f}, {fpdata2[ih].f}')



if __name__ == "__main__":
    mtzfile = "/Users/swh514/Work/data/7ins.mtz"
    m = gemmi.read_mtz_file(mtzfile)
    hklinfo = clipper.HKL_info.from_gemmi_mtz(m)
    fphidata = clipper.HKL_data_F_phi_float(hklinfo)
    data = m.get_f_phi("FWT", "PHWT")
    for d in data:
        hkl = clipper.HKL(d.hkl)
        fphidata[hkl] = clipper.F_phi_float(d.value)

    arr = fphidata.array
    ih = fphidata.first()
    ih.next()
    a = fphidata.data_export(ih.hkl())
    #print(ih.hkl(), a)
    #print(arr[])
    print(fphidata.num_obs())
    print(arr.shape)
    for i in range(0, fphidata.num_obs()):
        hkl = clipper.HKL(data[i].hkl)
        print(i, arr[i], fphidata[hkl].f, fphidata[hkl].phi, abs(data[i].value), cmath.phase(data[i].value))
    
    hklsamp = clipper.HKL_sampling(hklinfo.cell, clipper.Resolution(4.0))
    hklinfo2 = clipper.HKL_info()
    hklinfo2.init(hklinfo.spacegroup, hklinfo.cell, hklsamp)
    fpdata2 = fphidata.restrict_to(hklinfo2)
    ih = hklinfo2.first()
    while not ih.last():
        ih1 = clipper.HKL_info.HKL_reference_index(hklinfo, ih.index)
        v = clipper.F_phi_float()
        try:
            v = fphidata[ih.hkl()]
        except IndexError:
            print(ih.hkl().format(), fpdata2[ih].f, fpdata2[ih].phi, "nan, nan")
        if not v.missing():
            print(ih.hkl().format(), fpdata2[ih].f, fpdata2[ih].phi, fphidata[ih.hkl()].f, fphidata[ih.hkl()].phi)
        ih.next()
