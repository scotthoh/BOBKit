"""
Python script to run pybind11 version of Buccaneer
Author: S.W.Hoh, University of York, 2024
"""

# import faulthandler
import bobkit.clipper as clipper
import bobkit.buccaneer as buccaneer
import bobkit.util as util
import numpy as np
import sys
import gemmi
from .logger import set_logger
from .set_parameters import BuccaneerParams, BucArgParse
import datetime
from sklearn.cluster import DBSCAN

# faulthandler.enable()


def get_aa_pred(
    datapath: str, cell: clipper.Cell, grid_sampling: clipper.Grid_sampling
):
    aa_pred = np.load(f"{datapath}/inst_pred.npy")
    for i in range(0, len(aa_pred)):
        aa_pred[i] = np.swapaxes(aa_pred[i], 0, 2)
    nxs = -aa_pred[0].shape[0] // 2
    nys = -aa_pred[0].shape[1] // 2
    nzs = -aa_pred[0].shape[2] // 2
    ncorrect = 1
    spacing = np.array(
        [
            cell.a / aa_pred[0].shape[0],
            cell.b / aa_pred[0].shape[1],
            cell.c / aa_pred[0].shape[2],
        ]
    )
    corrections = [spacing[0], spacing[1], spacing[2], nxs, nys, nzs, ncorrect]
    # here 27 jan 2025


def get_coordinates_from_predicted_instance(
    datapath: str,
    cell: clipper.Cell,
    grid_sampling: clipper.Grid_sampling,
    fix_axis_positions=False,
    fix_origin=True,
    mapin_path="NONE",
    write_npy=False,
    return_map_index=True,
):
    offsets = np.load(f"{datapath}/inst_pred.npy")
    density = np.load(f"{datapath}/density.npy")
    nxs = 0
    nys = 0
    nzs = 0
    if mapin_path != "NONE":
        gmap = gemmi.read_ccp4_map(mapin_path)
        if fix_axis_positions and (gmap.grid.axis_order != gemmi.AxisOrder.XYZ):
            axis_pos = np.array(gmap.axis_positions())
            offsets = np.swapaxes(offsets, axis_pos[0], np.where(axis_pos == 0)[0][0])
            offsets = np.swapaxes(offsets, axis_pos[1], np.where(axis_pos == 1)[0][0])
            density = np.swapaxes(density, axis_pos[0], np.where(axis_pos == 0)[0][0])
            density = np.swapaxes(density, axis_pos[1], np.where(axis_pos == 1)[0][0])
            gmap.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
        if fix_origin:
            nxs = gmap.header_i32(5)
            nys = gmap.header_i32(6)
            nzs = gmap.header_i32(7)
            ncorrect = 0
        spacing = np.array(
            [
                cell.a / gmap.header_i32(8),
                cell.b / gmap.header_i32(9),
                cell.c / gmap.header_i32(10),
            ]
        )
    else:
        if fix_axis_positions:  # will swap the x and z axes from the ML npy output
            offsets = np.swapaxes(offsets, 1, 3)
            density = np.swapaxes(density, 0, 2)
        if fix_origin:
            nxs = -density.shape[0] // 2
            nys = -density.shape[1] // 2
            nzs = -density.shape[2] // 2
            ncorrect = 1
        spacing = np.array(
            [
                cell.a / density.shape[0],
                cell.b / density.shape[1],
                cell.c / density.shape[2],
            ]
        )
    if fix_origin:
        xyz = np.mgrid[
            nxs : density.shape[0] + (nxs + ncorrect),
            nys : density.shape[1] + (nys + ncorrect),
            nzs : density.shape[2] + (nzs + ncorrect),
        ].astype(np.float64)
    else:
        xyz = np.mgrid[
            nxs : density.shape[0], nys : density.shape[1], nzs : density.shape[2]
        ].astype(np.float64)

    print(f"xyz shape : {xyz.shape}")
    print(f"offsets shape : {offsets.shape}")
    print(f"density shape : {density.shape}")
    print(f"grid sampling : {grid_sampling.format()}")
    print(f"spacing: {spacing}")
    if not return_map_index:
        for i in (0, 1, 2):
            offsets[i] = offsets[i] * spacing[i]
            xyz[i] = xyz[i] * spacing[i]
    points = (xyz + offsets) * (density > 0)
    points = points[:, density > 0].T
    points_int = np.floor(points).astype(int)
    uniq, indices, count = np.unique(
        points_int, return_inverse=True, return_counts=True, axis=0
    )
    counts = count[indices]
    aa_instance_coordinates = []
    if return_map_index:
        if write_npy:
            np.save("aa_inst_mean.npy", uniq, allow_pickle=False)
        for point_ind in uniq:
            aa_instance_coordinates.append(clipper.Coord_map(point_ind))
    else:
        counts = count[indices]
        # use dbscan to group nearby points together
        db = DBSCAN(eps=np.max(np.sqrt(3 * (spacing * spacing))), min_samples=10).fit(
            points
        )
        labels = db.labels_
        groups = {}
        weights = {}
        for label in np.unique(labels):
            if label != -1:
                groups[label] = points[labels == label]
                weights[label] = counts[labels == label]

        tmp = np.zeros((len(groups), 3))
        count = 0
        for label, group in groups.items():
            sum_x = 0.0
            sum_y = 0.0
            sum_z = 0.0
            w = 0.0
            for p in zip(group, weights[label]):
                sum_x += p[0][0] * p[1]
                sum_y += p[0][1] * p[1]
                sum_z += p[0][2] * p[1]
                w += p[1]
            tmp[count] = np.array([sum_x / w, sum_y / w, sum_z / w])
            aa_instance_coordinates.append(clipper.Coord_orth(tmp[count]))
            count += 1
        if write_npy:
            np.save("aa_inst_mean.npy", tmp, allow_pickle=False)
    corrections = [spacing[0], spacing[1], spacing[2], nxs, nys, nzs, ncorrect]
    return aa_instance_coordinates, corrections  # spacing


def sequence_predictions_to_map(
    datapath: str,
    mapin_path: str = "NONE",
    fix_axis_positions: bool = False,
    fix_origin: bool = True,
):
    aa_pred = np.load(f"{datapath}/pred.npy")
    nxs = 0
    nys = 0
    nzs = 0
    if mapin_path != "NONE":
        gmap = gemmi.read_ccp4_map(mapin_path)
        if fix_axis_positions and (gmap.grid.axis_order != gemmi.AxisOrder.XYZ):
            axis_pos = np.array(gmap.axis_positions())
            aa_pred = np.swapaxes(aa_pred, axis_pos[0], np.where(axis_pos == 0)[0][0])
            aa_pred = np.swapaxes(aa_pred, axis_pos[1], np.where(axis_pos == 1)[0][0])
            gmap.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
        if fix_origin:
            nxs = gmap.header_i32(5)
            nys = gmap.header_i32(6)
            nzs = gmap.header_i32(7)
            ncorrect = 0
    else:
        if fix_axis_positions:  # will swap the x and z axes from the ML npy output
            aa_pred = np.swapaxes(aa_pred, 1, 3)
        if fix_origin:
            nxs = -aa_pred[0].shape[0] // 2
            nys = -aa_pred[0].shape[1] // 2
            nzs = -aa_pred[0].shape[2] // 2
            ncorrect = 1

    pred = np.zeros(aa_pred[0].shape, dtype=np.float32)
    for i in range(0, aa_pred.shape[0] - 1):
        pred += float(i + 1.0) * (aa_pred[i] > 0.8)
    #### here 16 Dec 2024
    xmap = clipper.Xmap_float()
    xmap.init(
        clipper.Spacegroup.p1(),
    )
    xmap.import_from_gemmi(gmap)

    if fix_origin:
        xyz = np.mgrid[
            nxs : aa_pred[0].shape[0] + (nxs + ncorrect),
            nys : aa_pred[0].shape[1] + (nys + ncorrect),
            nzs : aa_pred[0].shape[2] + (nzs + ncorrect),
        ].astype(np.float64)
    else:
        xyz = np.mgrid[
            nxs : aa_pred[0].shape[0],
            nys : aa_pred[0].shape[1],
            nzs : aa_pred[0].shape[2],
        ].astype(np.float64)


class Buccaneer:
    """Class for running buccaneer."""

    def __init__(self, buccaneer_args: BuccaneerParams):
        self.buccaneer_args = buccaneer_args
        self.log = buccaneer.Log(buccaneer_args.title)

    def _setup(self):
        buccaneer.Ca_prep.set_cpus(self.buccaneer_args.ncpu)
        buccaneer.Ca_find.set_cpus(self.buccaneer_args.ncpu)
        buccaneer.Ca_grow.set_cpus(self.buccaneer_args.ncpu)
        buccaneer.Ca_sequence.set_cpus(self.buccaneer_args.ncpu)
        buccaneer.Ca_sequence.set_semet(self.buccaneer_args.semet)
        self.buccaneer_args.check_steps()
        self.buccaneer_args.set_refcol_fo("FP.F_sigF.F,FP.F_sigF.sigF")
        self.buccaneer_args.set_refcol_hl(
            "FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D"
        )  # noqa 501

    def _print_steps_Casummaries(self, step, ca_num):
        """Print summaries from each step.

        Args:
            step (str): Step name
            ca_num (int): Number of C-alphas/residues.
        """
        msg = "C-alphas " + step
        print(" {0: <25} : {1: >7}".format(msg, ca_num))
        sys.stdout.flush()

    def run(self):
        """Run setup and Buccaneer"""
        self._setup()
        self.run_buccaneer()

    def run_buccaneer(self):
        """Run Buccaneer"""
        self.proteintools = buccaneer.ProteinTools()
        args = self.buccaneer_args
        if args.findtype == buccaneer.Ca_find.TYPE.SECSTRUC:
            print("Fast mode selected.\n")
        if args.correl:
            print("Correlation mode selected.\n")
        sys.stdout.flush()
        mtz = gemmi.read_mtz_file(args.mtzin_ref)
        mtzwrk = gemmi.read_mtz_file(args.mtzin)
        # resolution as set by clipper when using clipper::MTZfile
        res_ref = max(0.9999 / np.sqrt(mtz.max_1_d2), args.inresol)
        res_wrk = max(0.9999 / np.sqrt(mtzwrk.max_1_d2), args.inresol)
        resol = clipper.Resolution(max(res_ref, res_wrk))
        if args.verbose > 8:
            print("\nResolutions read : ")
            print(
                "Ref: {0:.5f}, Wrk: {1:.5f}, Input: {2:.5f}, "
                "Chosen: {3:.5f}\n".format(
                    res_ref,
                    res_wrk,
                    args.inresol,
                    resol.limit(),
                )
            )
            sys.stdout.flush()
        hkls_ref = clipper.HKL_info()

        hkls_ref.init(
            clipper.Spacegroup.from_gemmi_spacegroup(mtz.spacegroup),
            clipper.Cell.from_gemmi_cell(mtz.cell),
            resol,
            True,
        )
        sys.stdout.flush()

        ref_f = clipper.HKL_data_F_sigF_float(hkls_ref)
        ref_hl = clipper.HKL_data_ABCD_float(hkls_ref)
        ref_f.import_from_gemmi(mtz, args.ipcol_ref_fo, True)
        ref_hl.import_from_gemmi(mtz, args.ipcol_ref_hl, True)
        hkls_wrk = clipper.HKL_info()
        hkls_wrk.init(
            clipper.Spacegroup.from_gemmi_spacegroup(mtzwrk.spacegroup),
            clipper.Cell.from_gemmi_cell(mtzwrk.cell),
            resol,
            True,
        )
        if args.verbose > 8:
            print(
                f"\nmtzfile num_reflections : {mtzwrk.nreflections}",
            )
            print(
                f"hkl info num_reflections : {hkls_wrk.num_reflections()}\n",
            )
        sys.stdout.flush()
        wrk_f = clipper.HKL_data_F_sigF_float(hkls_wrk)
        wrk_hl = clipper.HKL_data_ABCD_float(hkls_wrk)
        wrk_pw = clipper.HKL_data_Phi_fom_float(hkls_wrk)
        wrk_fp = clipper.HKL_data_F_phi_float(hkls_wrk)
        flag = clipper.HKL_data_Flag(hkls_wrk)
        wrk_f.import_from_gemmi(mtzwrk, args.ipcol_fo, True)

        if args.ipcol_hl != "NONE":
            wrk_hl.import_from_gemmi(mtzwrk, args.ipcol_hl, True)
        if args.ipcol_phifom != "NONE":
            wrk_pw.import_from_gemmi(mtzwrk, args.ipcol_phifom, True)
        if args.ipcol_fc != "NONE":
            wrk_fp.import_from_gemmi(mtzwrk, args.ipcol_fc, True)
        if args.ipcol_free != "NONE":
            flag.import_from_gemmi(mtzwrk, args.ipcol_free, True)
        # aniso correction
        uaniso = clipper.U_aniso_orth(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if args.doanis:
            sfscl = clipper.SFscale_aniso_float(
                3.0, clipper.SFscale_aniso_float.SHARPEN
            )
            sfscl(wrk_f, -1.0, 12)  # this is the default for sfscal(wrk_f)
            uaniso = sfscl.u_aniso_orth(clipper.SFscale_aniso_float.F)
            # scale map coeffs
            if args.ipcol_fc != "NONE":
                wrk_fp.compute_scale_u_aniso_fphi(1.0, -uaniso, wrk_fp)
            print("Applying anisotropy correction:\n")
            print(uaniso.format())
            print("\n")
        # apply freer-flag
        wrk_fwrk = clipper.HKL_data_F_sigF_float()
        wrk_fwrk.copy_from(wrk_f)
        ih = hkls_wrk.first()
        sys.stdout.flush()
        while not ih.last():
            if flag[ih].flag == args.freerindex:
                wrk_fwrk[ih] = clipper.HKL_data_F_sigF_float()
            ih.next()
        # fill in hl
        if args.ipcol_hl == "NONE":
            wrk_hl.compute_from_phi_fom(wrk_pw)
        # get reference model
        cspg = hkls_wrk.spacegroup
        wrkcell = hkls_wrk.cell
        mol_ref = clipper.MiniMol()
        mol_wrk = clipper.MiniMol(cspg, wrkcell)
        mol_mr = clipper.MiniMol(cspg, wrkcell)
        mol_seq = clipper.MiniMol(cspg, wrkcell)

        gfile_ref = clipper.GEMMIfile()
        gfile_ref.read_file(args.pdbin_ref)
        gfile_ref.import_minimol(mol_ref)
        if args.pdbin != "NONE":
            util.read_structure(mol_wrk, args.pdbin)
        if mol_seq.size() > 0:
            buccaneer.Ca_sequence.set_prior_model(mol_seq)
        mol_wrk_in = mol_wrk.copy()  # store a copy of the input model
        knownstruc = buccaneer.KnownStructure(
            mol_wrk_in,
            args.known_ids,
            args.nprad,
        )
        if args.verbose >= 1:
            knownstruc.debug()
        # get work sequence
        seq_wrk = clipper.MMoleculeSequence()
        if args.ipseq_wrk != "NONE":
            seqf_wrk = clipper.SEQfile()
            seqf_wrk.read_file(args.ipseq_wrk)
            seqf_wrk.import_molecule_sequence(seq_wrk)

        # check input files match
        if not args.find and mol_wrk.is_null():
            ValueError("Missing work model!")
        if args.seqnc and seq_wrk.is_null():
            ValueError("Missing work sequence!")
        sys.stdout.flush()
        # initial map simulation
        sim_f = clipper.HKL_data_F_sigF_float(hkls_ref)
        sim_hl = clipper.HKL_data_ABCD_float(hkls_ref)
        mapsim = buccaneer.MapSimulate(100, 20)
        mapsim(sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl)

        # make llk target objects
        llktgt = buccaneer.LLK_map_target()
        llkcls = buccaneer.LLK_TargetList(20)

        # STAGE 1. Calculate target from reference data
        ref_fp = clipper.HKL_data_F_phi_float(hkls_ref)
        ref_pw = clipper.HKL_data_Phi_fom_float(hkls_ref)
        ref_pw.compute_from_abcd(sim_hl)
        ref_fp.compute_from_fsigf_phifom(sim_f, ref_pw)
        grid = clipper.Grid_sampling(
            hkls_ref.spacegroup, hkls_ref.cell, hkls_ref.resolution
        )
        xref = clipper.Xmap_float(hkls_ref.spacegroup, hkls_ref.cell, grid)
        xref.fft_from(ref_fp)
        # ccp41 = gemmi.Ccp4Map()
        # xref.export_to_gemmi(ccp41)
        # ccp41.write_ccp4_map("scaled_refmap.map")

        # prepare llk targets
        sys.stdout.flush()
        caprep = buccaneer.Ca_prep(
            args.main_tgt_rad,
            args.side_tgt_rad,
            args.rama_fltr,
            args.correl,
            args.seqnc,
            args.verbose > 3,
        )
        caprep(llktgt, llkcls, mol_ref, xref)
        sys.stdout.flush()
        self.log.log("PREP")

        # STAGE 2: Apply target to work data
        print("apply target to work data\n")
        wrk_pw.compute_from_abcd(wrk_hl)
        if args.ipcol_fc == "NONE":
            wrk_fp.compute_from_fsigf_phifom(wrk_fwrk, wrk_pw)
        gridxwrk = clipper.Grid_sampling(cspg, wrkcell, hkls_wrk.resolution)
        xwrk = clipper.Xmap_float(cspg, wrkcell, gridxwrk)
        xwrk.fft_from(wrk_fp)
        # write out work map (optional)
        if args.mapout != "NONE":
            ccp4 = gemmi.Ccp4Map()
            ccp4.grid = gemmi.FloatGrid(
                np.zeros(
                    (
                        xwrk.grid_sampling.nu,
                        xwrk.grid_sampling.nv,
                        xwrk.grid_sampling.nw,
                    ),
                    dtype=np.float32,
                )
            )
            ccp4.grid.unit_cell = clipper.Cell.to_gemmi_cell(wrkcell)
            ccp4.grid.spacegroup = clipper.Spacegroup.to_gemmi_spacegroup(cspg)
            ccp4.update_ccp4_header()
            xwrk.export_to_gemmi(ccp4)
            ccp4.write_ccp4_map(args.mapout)
        # initial number of fragments/residues to find
        vol = xwrk.cell.volume / float(xwrk.spacegroup.num_symops())

        nres = int(vol / 320.0)  # 320A^3/residue on average (inc solvent)
        args.nfrag = min(args.nfrag, int((args.nfragr * nres) / 100))
        sys.stdout.flush()
        # generate llk distribution of target values for cutoff
        llktgt.prep_llk_distribution(xwrk)
        buccaneer.ProteinTools.split_chains_at_gap(mol_wrk)
        sys.stdout.flush()
        # prepare search target
        cafind = buccaneer.Ca_find(args.nfrag, resol.limit())
        if args.aa_instance_directory != "NONE":
            aa_instance_coords, spacing = get_coordinates_from_predicted_instance(
                args.aa_instance_directory,
                xwrk.cell,
                xwrk.grid_sampling,
                fix_axis_positions=True,
                mapin_path=args.mapin,
                write_npy=True,
                return_map_index=False,
            )
            cafind.set_starting_centroid_coords(aa_instance_coords)
            aa_pred = np.load(f"{args.aa_instance_directory}/pred.npy")
            caseq_ml = buccaneer.Ca_sequence_ml(aa_pred, spacing, correl=args.correl)

        # merge multi model results
        if args.merge:
            camerge = buccaneer.Ca_merge(args.seq_rel)
            print(
                " C-alphas before model merge: {0}".format(
                    mol_wrk.select("*/*/CA").atom_list().size()
                )
            )
            camerge(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
            print(
                " C-alphas after model merge: {0}".format(
                    mol_wrk.select("*/*/CA").atom_list().size()
                )
            )
            sys.stdout.flush()
            self.log.log("MRG ", mol_wrk, args.verbose > 9)
            sys.stdout.flush()
        # filter input model
        if args.model_filter:
            self._print_steps_Casummaries(
                "before model_filter",
                mol_wrk.select("*/*/CA").atom_list().size(),
            )
            buccaneer.Ca_filter(mol_wrk, xwrk, args.model_filter_sig)
            self._print_steps_Casummaries(
                "after model filter",
                mol_wrk.select("*/*/CA").atom_list().size(),
            )
            sys.stdout.flush()
            self.log.log("FLT ", mol_wrk, args.verbose > 9)
            sys.stdout.flush()
        # augment input model with mr model
        if args.mr_model:
            result = buccaneer.Ca_merge.merge_mr(
                mol_wrk,
                mol_mr,
                args.mr_filter_sig,
                3,
                args.mr_model_filter,
                args.mr_model_seed,
            )
            print(
                " MR residues input: {0}, after filter: {1}, "
                "after seeding: {2}".format(result[0], result[1], result[2])
            )
            sys.stdout.flush()

        # trim input model to only protein residues
        buccaneer.ProteinTools.trim_to_protein(mol_wrk)

        # model building loop
        for cyc in range(0, args.ncyc):
            print(f"\nCycle: {cyc+1}\n")
            history = ""
            sys.stdout.flush()
            # find C-alphas by slow likelihood search
            if args.find:
                cafind(
                    mol_wrk,
                    knownstruc,
                    xwrk,
                    llktgt,
                    args.findtype,
                    args.modelindex,
                )
                self._print_steps_Casummaries(
                    "after finding",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("FIND", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.grow:
                cagrow = buccaneer.Ca_grow(25)
                cagrow(mol_wrk, xwrk, llktgt)
                self._print_steps_Casummaries(
                    "after growing",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("GROW", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.join:
                cajoin = buccaneer.Ca_join(2.0, 2.0)
                cajoin(mol_wrk)
                self._print_steps_Casummaries(
                    "after joining",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("JOIN", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.link:
                calnk = buccaneer.Ca_link(10.0, 24)
                calnk(mol_wrk, xwrk, llktgt)
                self._print_steps_Casummaries("linked", calnk.num_linked)
                sys.stdout.flush()
                self.log.log("LINK", mol_wrk, args.verbose > 9)
                sys.stdout.flush()
                # util.write_structure(mol_wrk, "linked.pdb", cif_format=False)

            if args.seqnc:
                # caseq_ml = buccaneer.Ca_sequence_ml(aa_pred, corrections)
                if args.aa_instance_directory != "NONE":
                    caseq_ml(mol_wrk)
                    if not caseq_ml.check_is_seqprob_set():
                        print("SEQPROB not set!")
                buccaneer.Ca_sequence.set_use_ml_sequence_probability(True)
                caseq = buccaneer.Ca_sequence(args.seq_rel)
                # caseq.set_use_set_use_ml_sequence_probability(True)
                caseq(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                self._print_steps_Casummaries(
                    "sequenced", caseq.num_sequenced()
                )  # noqa E501
                history = caseq.format()
                sys.stdout.flush()
                self.log.log("SEQU", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.corct:
                cacor = buccaneer.Ca_correct(12)
                cacor(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                self._print_steps_Casummaries("corrected", cacor.num_corrected)
                sys.stdout.flush()
                self.log.log("CORR", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.filtr:
                cafiltr = buccaneer.Ca_filter(1.0)
                cafiltr(mol_wrk, xwrk)
                self._print_steps_Casummaries(
                    "after filtering",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("FILT", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.ncsbd:
                cancsbuild = buccaneer.Ca_ncsbuild(args.seq_rel, 1.0, 12)
                cancsbuild(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                self._print_steps_Casummaries(
                    "after NCS build",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("NCSB", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.prune:
                caprune = buccaneer.Ca_prune(3.0)
                caprune(mol_wrk, xwrk)
                self._print_steps_Casummaries(
                    "after pruning",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("PRUN", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            if args.build:
                cabuild = buccaneer.Ca_build(args.newrestype)
                cabuild(mol_wrk, xwrk)
                self._print_steps_Casummaries(
                    "after rebuilding",
                    mol_wrk.select("*/*/CA").atom_list().size(),
                )
                sys.stdout.flush()
                self.log.log("REBU", mol_wrk, args.verbose > 9)
                sys.stdout.flush()

            knownstruc.prune(mol_wrk)
            buccaneer.ProteinTools.split_chains_at_gap(mol_wrk)
            buccaneer.ProteinTools.chain_number(mol_wrk)
            buccaneer.ProteinTools.chain_label(mol_wrk, args.chainid_2char)
            self.log.log("TDYI", mol_wrk, args.verbose > 9)
            if args.verbose > 7:
                print(history)
            sys.stdout.flush()
            msg = self.log.summary(mol_wrk, mol_mr, seq_wrk)
            print(f"\nInternal cycle {cyc+1:3d}\n {msg}")
            sys.stdout.flush()
            # file output
            if args.xmlout != "NONE":
                self.log.xml(args.xmlout)
            sys.stdout.flush()
            # output intermediates (optional)
            if args.optemp:
                gfile = clipper.GEMMIfile()
                gfile.export_minimol(mol_wrk)
                gfile.write_pdb(args.get_outfile_name(cyc + 1))
                if args.write_cif:
                    gfile.write_cif(args.get_outfile_name(cyc + 1, True))
            # end cycle
        # move model
        if args.fixpos:
            buccaneer.ProteinTools.symm_match(mol_wrk, mol_wrk_in)
        # model tidy
        if args.tidy:
            mtidy = buccaneer.ModelTidy(
                1.0,
                12,
                args.newrestype,
                args.verbose > 6,
            )
            chk = mtidy.tidy(mol_wrk, mol_mr, seq_wrk)
            if not chk:
                print("ModelTidy error")
            sys.stdout.flush()
            self.log.log("TIDY", mol_wrk, args.verbose > 9)
            sys.stdout.flush()
        # assign default B-factors to missing values
        default_u_iso = buccaneer.ProteinTools.main_chain_u_mean(mol_wrk_in)
        for c in mol_wrk:
            for r in c:
                for a in r:
                    if clipper.Util.is_nan(a.u_iso):
                        a.u_iso = default_u_iso

        if args.newresname != "NONE":
            for c in mol_wrk:
                for r in c:
                    if buccaneer.ProteinTools.residue_index_3(r.type) < 0:
                        r.type = args.newresname

        # add known structure from input model
        known_copy = knownstruc.copy_to(mol_wrk)
        if not known_copy:
            print("$TEXT:Warning: $$ $$\n")
            print(
                "WARNING: chain ID clash between known-structure and pdbin-mr."
                "Chains renamed.\n$$"
            )
        sys.stdout.flush()
        # label unlabelled chains
        buccaneer.ProteinTools.chain_label(mol_wrk, args.chainid_2char)
        sys.stdout.flush()
        self.log.log("LABE", mol_wrk, args.verbose > 9)
        sys.stdout.flush()

        # write answers
        gfile = clipper.GEMMIfile()
        gfile.export_minimol(mol_wrk)
        gfile.write_pdb(args.get_outfile_name())
        if args.write_cif:
            gfile.write_cif(args.get_outfile_name(cifout=True))

        if args.verbose > 8:
            self.log.profile()
        # end buccaneer


def main(args=None):
    """Example buccaneer run script"""
    print(f"\n### Script started : {datetime.datetime.now()} ###\n")
    raw_args = args or sys.argv[1:]
    parser = BucArgParse(sys.argv[0])
    parsed_args = parser.parse_args(raw_args)
    set_logger("pybuc.log", "w", stdout_level=1, stderr_level=1)
    parser.print_command_with_args()
    buc_params = BuccaneerParams(parsed_args)
    # if "CLIB" not in os.environ:
    # setting EM reference data
    # buc_params.set_reference_data(
    #    pdbref="/opt/xtal/ccp4-8.0/lib/data/reference_structures/reference-EMD-4116.pdb",  # noqa E501
    #    mtzref="/opt/xtal/ccp4-8.0/lib/data/reference_structures/reference-EMD-4116.mtz",  # noqa E501
    # )
    # setting X-ray reference data
    # buc_params.set_reference_data(
    #    pdbref="/opt/xtal/ccp4-8.0/lib/data/reference_structures/reference-1tqw.pdb",  # noqa E501
    #    mtzref="/opt/xtal/ccp4-8.0/lib/data/reference_structures/reference-1tqw.mtz",  # noqa E501
    # )
    buc = Buccaneer(buc_params)
    print(f"\nmtzin : {buc_params.mtzin}")
    print(f"pdbin : {buc_params.pdbin}")
    print(f"seqin : {buc_params.ipseq_wrk}\n")
    buc.run()
    print(f"\n### Script ended : {datetime.datetime.now()} ###\n")


if __name__ == "__main__":
    main()
