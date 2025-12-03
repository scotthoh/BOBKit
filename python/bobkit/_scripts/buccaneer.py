"""
Python script to run pybind11 version of Buccaneer
Author: S.W.Hoh, University of York, 2024
"""

import faulthandler
import bobkit.clipper as clipper
import bobkit.buccaneer as buccaneer
import bobkit.util as util
import numpy as np
import sys
import gemmi
from .logger import log2file
from .set_parameters import BuccaneerParams, BucArgParse
import datetime
# from sklearn.cluster import DBSCAN
from os import getcwd, environ
from typing import List
# from functools import wraps
faulthandler.enable()


class Buccaneer:
    """Class for running buccaneer."""

    def __init__(self, buccaneer_args: BuccaneerParams):
        """Initialise

        Args:
            buccaneer_args (BuccaneerParams): Input arguments
        """
        self.args = buccaneer_args
        self.use_ml_seq = False
        if self.args.sequence_method in ["mlinput", "hybrid"]:
            self.use_ml_seq = True
            print("Using ML sequence")
        self.log = buccaneer.Log(buccaneer_args.title)
        self.proteintools = buccaneer.ProteinTools()

    def _debug_print(self, msg, value):
        """Internal debug print out

        Args:
            msg (str): Message to print
            value (Union[str, int, float]): Value to print
        """
        print(f"\n{msg} : {value}\n")

    def run(self):
        """Run setup and Buccaneer"""
        llktgt, llkcls, xwrk, hkls_wrk = self._setup()
        mol_wrk, mol_mr, knownstruc = self.get_input_models_and_knownstructure(hkls_wrk)  # fmt: skip # noqa: E501
        self.run_buccaneer(llktgt, llkcls, xwrk, mol_wrk, mol_mr, knownstruc)

    def _setup(self):
        """Set up all that is needed before running Buccaneer's
        model building steps.

        Returns:
            llktgt (buccaneer.LLK_map_target): Log-likelihood map targets object
            llkcls (buccaneer.LLK_TargetList): Log-likelihood Target list object
            xwrk (clipper.Xmap_float): Xmap object
            hkls_wrk (clipper.HKL_info): Work HKL info object
        """  # noqa: E501
        buccaneer.Ca_prep.set_cpus(self.args.ncpu)
        buccaneer.Ca_find.set_cpus(self.args.ncpu)
        buccaneer.Ca_grow.set_cpus(self.args.ncpu)
        buccaneer.Ca_sequence.set_cpus(self.args.ncpu)
        buccaneer.Ca_sequence.set_semet(self.args.semet)
        self.args.check_steps()
        self.args.set_refcol_fo("FP.F_sigF.F,FP.F_sigF.sigF")
        self.args.set_refcol_hl("FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D")  # noqa: E501
        if self.args.findtype == buccaneer.Ca_find.TYPE.SECSTRUC:
            print("Fast mode selected.\n")
        if self.args.correl:
            print("Correlation mode selected.\n")
        sys.stdout.flush()
        mtz = gemmi.read_mtz_file(self.args.mtzin_ref)
        mtzwrk = gemmi.read_mtz_file(self.args.mtzin)
        # resolution as set by clipper when using clipper::MTZfile
        res_ref = max(0.9999 / np.sqrt(mtz.max_1_d2), self.args.inresol)
        res_wrk = max(0.9999 / np.sqrt(mtzwrk.max_1_d2), self.args.inresol)
        self.resol = clipper.Resolution(max(res_ref, res_wrk))
        if self.args.verbose > 8:
            print("\nResolutions read : ")
            print(
                "Ref: {0:.5f}, Wrk: {1:.5f}, Input: {2:.5f}, "
                "Chosen: {3:.5f}\n".format(
                    res_ref,
                    res_wrk,
                    self.args.inresol,
                    self.resol.limit(),
                )
            )
            sys.stdout.flush()
        # get reference data
        hkls_ref = clipper.HKL_info()
        hkls_ref.init(
            clipper.Spacegroup(mtz.spacegroup.ccp4),
            clipper.Cell(mtz.cell.parameters),
            self.resol,
            True,
        )
        sys.stdout.flush()

        ref_f = clipper.HKL_data_F_sigF_float(hkls_ref)
        ref_hl = clipper.HKL_data_ABCD_float(hkls_ref)
        ref_f.import_from_gemmi(mtz, self.args.ipcol_ref_fo, True)
        ref_hl.import_from_gemmi(mtz, self.args.ipcol_ref_hl, True)
        hkls_wrk = clipper.HKL_info()
        hkls_wrk.init(
            clipper.Spacegroup(mtzwrk.spacegroup.ccp4),
            clipper.Cell(mtzwrk.cell.parameters),
            self.resol,
            True,
        )
        if self.args.verbose > 8:
            self._debug_print("Mtzfile num_reflections", mtzwrk.nreflections)
            self._debug_print(
                "HKL info num_reflections", hkls_wrk.num_reflections()
            )  # noqa: E501
        sys.stdout.flush()
        wrk_f = clipper.HKL_data_F_sigF_float(hkls_wrk)
        wrk_hl = clipper.HKL_data_ABCD_float(hkls_wrk)
        wrk_pw = clipper.HKL_data_Phi_fom_float(hkls_wrk)
        wrk_fp = clipper.HKL_data_F_phi_float(hkls_wrk)
        flag = clipper.HKL_data_Flag(hkls_wrk)
        wrk_f.import_from_gemmi(mtzwrk, self.args.ipcol_fo, True)

        if self.args.ipcol_hl != "NONE":
            wrk_hl.import_from_gemmi(mtzwrk, self.args.ipcol_hl, True)
        if self.args.ipcol_phifom != "NONE":
            wrk_pw.import_from_gemmi(mtzwrk, self.args.ipcol_phifom, True)
        if self.args.ipcol_fc != "NONE":
            wrk_fp.import_from_gemmi(mtzwrk, self.args.ipcol_fc, True)
        if self.args.ipcol_free != "NONE":
            flag.import_from_gemmi(mtzwrk, self.args.ipcol_free, True)
        # aniso correction
        self.u_aniso_correction(wrk_f, wrk_fp, ipcol_fc=self.args.ipcol_fc)
        # apply freer-flag
        wrk_fwrk = self.apply_freerflag(wrk_f, hkls_wrk, flag, self.args.freerindex)  # noqa: E501
        # fill in hl
        if self.args.ipcol_hl == "NONE":
            wrk_hl.compute_from_phi_fom(wrk_pw)
        sim_f, sim_hl = self.simulate_map(hkls_ref, ref_f, ref_hl, wrk_f, wrk_hl)  # noqa: E501
        # calculate targets
        llktgt, llkcls = self.calculate_target_from_ref_data(hkls_ref, sim_f, sim_hl)  # noqa: E501
        xwrk = self.apply_llk_targets(
            llktgt, hkls_wrk, wrk_pw, wrk_hl, wrk_fp, wrk_fwrk
        )
        # initial number of fragments/residues to find
        vol = xwrk.cell.volume / float(xwrk.spacegroup.num_symops)
        nres = int(vol / 320.0)  # 320A^3/residue on average (inc solvent)
        self.args.nfrag = min(self.args.nfrag, int((self.args.nfragr * nres) / 100))  # noqa: E501

        return llktgt, llkcls, xwrk, hkls_wrk

    def get_input_models_and_knownstructure(self, hkls_wrk: clipper.HKL_info):
        """Return input models and knownstructure object

        Args:
            hkls_wrk (clipper.HKL_info): Work HKL info

        Returns:
            tuple: A tuple containing:
                - mol_wrk (clipper.MiniMol): Work MiniMol
                - mol_mr (clipper.MiniMol): Molecular replacement MiniMol
                - knownstruc (buccaneer.KnownStructure): KnownStructure object
        """
        # if self.args.pdbin != "NONE":
        mol_wrk = clipper.MiniMol(hkls_wrk.spacegroup, hkls_wrk.cell)
        mol_mr = clipper.MiniMol(hkls_wrk.spacegroup, hkls_wrk.cell)
        mol_seq = clipper.MiniMol(hkls_wrk.spacegroup, hkls_wrk.cell)
        util.read_structure(mol_wrk, self.args.pdbin)
        util.read_structure(mol_mr, self.args.pdbin_mr)
        util.read_structure(mol_seq, self.args.pdbin_seq)
        if mol_seq.size() > 0:
            buccaneer.Ca_sequence.set_prior_model(mol_seq)
        # check input files match
        if not self.args.find and mol_wrk.is_null():
            ValueError("Missing work model!")
            sys.stdout.flush()
        # store a copy of the input model
        mol_wrk_in = mol_wrk.clone()
        # prepare known structure
        knownstruc = buccaneer.KnownStructure(
            mol_wrk_in,
            self.args.known_ids,
            self.args.nprad,
        )
        if self.args.verbose >= 1:
            msg = ""
            knownstruc.debug(msg)
            print(msg)
        buccaneer.ProteinTools.split_chains_at_gap(mol_wrk)
        return mol_wrk, mol_mr, knownstruc

    @classmethod
    def _print_steps_Casummaries(cls, step, ca_num):
        """Print summaries from each step.

        Args:
            step (str): Step name
            ca_num (int): Number of C-alphas/residues.
        """
        msg = "C-alphas " + step
        print(" {0: <25} : {1: >7}".format(msg, ca_num))
        sys.stdout.flush()

    @classmethod
    def print_steps_Casummaries(cls, step, ca_num):
        """Print summaries from each step.

        Args:
            step (str): Step name
            ca_num (int): Number of C-alphas/residues.
        """
        cls._print_steps_Casummaries(step, ca_num)

    # def run_buccaneer_step(self, func):
    #    @wraps(func)
    #    def wrapper(*args, **kwds):
    #        func(*args, **kwds)
    #        self._print_steps_Casummaries(
    #                "after growing",
    #                len(mol_wrk.select("*/*/CA").atom_list()),
    #            )
    #            sys.stdout.flush()
    #            self.log.log("GROW", mol_wrk, self.args.verbose > 9)

    def u_aniso_correction(
        self,
        wrk_f: clipper.HKL_data_F_sigF_float,
        wrk_fp: clipper.HKL_data_F_phi_float,
        ipcol_fc: str = "NONE",
    ):
        """Perform U-anisotropic corrections

        Args:
            wrk_f (clipper.HKL_data_F_sigF_float, clipper.HKL_data_F_sigF_double): F_sigF hkl data
            wrk_fp (clipper.HKL_data_F_phi_float, clipper.HKL_data_F_phi_double): F_phi hkl data
            ipcol_fc (str, optional): Column name. Defaults to "NONE".
        """  # noqa: E501
        uaniso = clipper.U_aniso_orth(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        sfscl = clipper.SFscale_aniso_float(
            3.0,
            clipper.SFscale_aniso_float.SHARPEN,
        )
        sfscl(wrk_f, -1.0, 12)  # this is the default for sfscal(wrk_f)
        uaniso = sfscl.u_aniso_orth(clipper.SFscale_aniso_float.F)
        # scale map coeffs
        if ipcol_fc != "NONE":
            wrk_fp.compute_scale_u_aniso_fphi(1.0, -uaniso, wrk_fp)
        print(f"Applying anisotropy correction:\n{uaniso.format()}\n")

    def apply_freerflag(
        self,
        wrk_f: clipper.HKL_data_F_sigF_float,
        hkls_wrk: clipper.HKL_info,
        flag: clipper.HKL_data_Flag,
        freerindex: int = 0,
    ):
        """Flag data according to FreeR flags

        Args:
            wrk_f (clipper.HKL_data_F_sigF_float): Work F_sigF HKL data
            hkls_wrk (clipper.HKL_info): Work HKL info
            flag (clipper.HKL_data_Flag): FreeR flags HKL data
            freerindex (int, optional): FreeR index. Defaults to 0.

        Returns:
            clipper.HKL_dat_F_sigF_float: Flagged work F_sigF HKL data
        """
        wrk_fwrk = clipper.HKL_data_F_sigF_float()
        wrk_fwrk.copy_from(wrk_f)
        ih = hkls_wrk.first()
        sys.stdout.flush()
        while not ih.last():
            if flag[ih].flag == freerindex:
                wrk_fwrk[ih] = clipper.HKL_data_F_sigF_float()
            ih.next()
        return wrk_fwrk

    def simulate_map(
        self,
        hkls_ref: clipper.HKL_info,
        ref_f: clipper.HKL_data_F_sigF_float,
        ref_hl: clipper.HKL_data_ABCD_float,
        wrk_f: clipper.HKL_data_F_sigF_float,
        wrk_hl: clipper.HKL_data_ABCD_float,
    ):
        """Initial map simulation

        Args:
            hkls_ref (clipper.HKL_info): Reference HKL info
            ref_f (clipper.HKL_data_F_sigF_float): Reference F_sigF HKL data
            ref_hl (clipper.HKL_data_ABCD_float): Reference Hendrickson-Lattman coefficients HKL data
            wrk_f (clipper.HKL_data_F_sigF_float): Work F_sigF HKL data
            wrk_hl (clipper.HKL_data_ABCD_float): Reference Hendrickson-Lattman coefficients HKL data

        Returns:
            tuple: A tuple containing:
                - clipper.HKL_data_F_sigF_float: Simulated F_sigF HKL data
                - clipper.HKL_data_ABCD_float: Simulated Hendrickson-Lattman coefficents HKL data
        """  # noqa: E501
        sim_f = clipper.HKL_data_F_sigF_float(hkls_ref)
        sim_hl = clipper.HKL_data_ABCD_float(hkls_ref)
        mapsim = buccaneer.MapSimulate(100, 20)
        mapsim(sim_f, sim_hl, ref_f, ref_hl, wrk_f, wrk_hl)

        return sim_f, sim_hl

    def calculate_target_from_ref_data(
        self,
        hkls_ref: clipper.HKL_info,
        sim_f: clipper.HKL_data_F_sigF_float,
        sim_hl: clipper.HKL_data_ABCD_float,
    ):
        """Stage 1. Calculate target from reference data

        Args:
            hkls_ref (clipper.HKL_info): Reference HKL info
            sim_f (clipper.HKL_data_F_sigF_float): Simulated F_sigF HKL data
            sim_hl (clipper.HKL_data_ABCD_float): Simulated Hendrickson-Lattman coefficients HKL data

        Returns:
            tuple: A tuple containing:
                - buccaneer.LLK_map_target: Log-likehood map target
                - buccaneer.LLK_TargetList: List of log-likehood map targets
        """  # noqa: E501
        ref_fp = clipper.HKL_data_F_phi_float(hkls_ref)
        ref_pw = clipper.HKL_data_Phi_fom_float(hkls_ref)
        ref_pw.compute_from_abcd(sim_hl)
        ref_fp.compute_from_fsigf_phifom(sim_f, ref_pw)
        grid = clipper.Grid_sampling(
            hkls_ref.spacegroup, hkls_ref.cell, hkls_ref.resolution
        )
        xref = clipper.Xmap_float(hkls_ref.spacegroup, hkls_ref.cell, grid)
        xref.fft_from(ref_fp)
        # write out work map (optional)
        if self.args.refmapout != "NONE":
            ccp4 = gemmi.Ccp4Map()
            ccp4.grid = gemmi.FloatGrid(
                np.zeros(
                    (
                        xref.grid_sampling.nu,
                        xref.grid_sampling.nv,
                        xref.grid_sampling.nw,
                    ),
                    dtype=np.float32,
                )
            )
            ccp4.grid.unit_cell = clipper.Cell.to_gemmi_cell(hkls_ref.cell)
            ccp4.grid.spacegroup = clipper.Spacegroup.to_gemmi_spacegroup(
                hkls_ref.spacegroup
            )
            ccp4.update_ccp4_header()
            xref.export_to_gemmi(ccp4)
            ccp4.write_ccp4_map(self.args.refmapout)
        sys.stdout.flush()
        # make llk target objects
        llktgt = buccaneer.LLK_map_target()
        llkcls = buccaneer.LLK_TargetList(20)
        # prepare llk targets
        caprep = buccaneer.Ca_prep(
            self.args.main_tgt_rad,
            self.args.side_tgt_rad,
            self.args.rama_fltr,
            self.args.correl,
            self.args.seqnc,
            self.args.verbose > 3,
        )
        mol_ref = clipper.MiniMol()
        gfile_ref = clipper.GEMMIfile()
        gfile_ref.read_file(self.args.pdbin_ref)
        gfile_ref.import_minimol(mol_ref)
        caprep(llktgt, llkcls, mol_ref, xref)
        sys.stdout.flush()
        self.log.log("PREP")

        return llktgt, llkcls

    def apply_llk_targets(
        self,
        llktgt: buccaneer.LLK_map_target,
        hkls_wrk: clipper.HKL_info,
        wrk_pw: clipper.HKL_data_Phi_fom_float,
        wrk_hl: clipper.HKL_data_ABCD_float,
        wrk_fp: clipper.HKL_data_F_phi_float,
        wrk_fwrk: clipper.HKL_data_F_sigF_float,
    ):
        """Stage 2. Apply target to work data

        Args:
            llktgt (buccaneer.LLK_map_target): Log likelihood targets
            hkls_wrk (clipper.HKL_info): Work HKL_info
            wrk_pw (clipper.HKL_data_Phi_fom_float): Work Phi_fom HKL data
            wrk_hl (clipper.HKL_data_ABCD_float): Work Hendrickson Lattmaan HKL data
            wrk_fp (clipper.HKL_data_F_phi_float): Work F_phi HKL data
            wrk_fwrk (clipper.HKL_data_F_sig_float): Work F_sigF HKL_data

        Returns:
            clipper.Xmap: Work Xmap
        """  # noqa: E501
        print("Applying target to work data\n")
        wrk_pw.compute_from_abcd(wrk_hl)
        if self.args.ipcol_fc == "NONE":
            wrk_fp.compute_from_fsigf_phifom(wrk_fwrk, wrk_pw)
        gridxwrk = clipper.Grid_sampling(
            hkls_wrk.spacegroup, hkls_wrk.cell, hkls_wrk.resolution
        )
        xwrk = clipper.Xmap_float(hkls_wrk.spacegroup, hkls_wrk.cell, gridxwrk)
        xwrk.fft_from(wrk_fp)
        # write out work map (optional)
        if self.args.mapout != "NONE":
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
            ccp4.grid.unit_cell = clipper.Cell.to_gemmi_cell(hkls_wrk.cell)
            ccp4.grid.spacegroup = clipper.Spacegroup.to_gemmi_spacegroup(
                hkls_wrk.spacegroup
            )
            ccp4.update_ccp4_header()
            xwrk.export_to_gemmi(ccp4)
            ccp4.write_ccp4_map(self.args.mapout)

        # generate llk distribution of target values for cutoff
        llktgt.prep_llk_distribution(xwrk)
        sys.stdout.flush()
        return xwrk

    def get_work_sequence(self):
        """Return work sequence read from file

        Returns:
            clipper.MMoleculeSequence: Work molecule sequence
        """
        # get work sequence
        seq_wrk = clipper.MMoleculeSequence()
        if self.args.ipseq_wrk != "NONE":
            seqf_wrk = clipper.SEQfile()
            seqf_wrk.read_file(self.args.ipseq_wrk)
            seqf_wrk.import_molecule_sequence(seq_wrk)
        if self.args.seqnc and seq_wrk.is_null():
            ValueError("Missing work sequence!")
        sys.stdout.flush()
        return seq_wrk

    def prepare_ca_find(
        self,
        xwrk: clipper.Xmap_float,
        seq_wrk: clipper.MMoleculeSequence,
        llktgt: buccaneer.LLK_map_target,
        seqlen_multiplier: int = 1,
    ):
        """Prepare Ca find
        Args:
            xwrk (clipper.Xmap_float): Target Xmap
            return_map_index (bool): Return map index if true

        Returns:
            tuple: A tuple containing:
                - cafind (buccaneer.Ca_find): Ca find instance
                - caseq_ml (buccaneer.Ca_sequence_ml): Ca sequence ml instance
                - osaka (buccaneer.Osaka): Osaka instance
        """
        wrkcell = np.array(
                    [
                        xwrk.cell.a,
                        xwrk.cell.b,
                        xwrk.cell.c,
                        xwrk.cell.alpha,
                        xwrk.cell.beta,
                        xwrk.cell.gamma,
                    ]
                )
        if self.args.aa_instance_directory != "NONE" and self.args.use_ml_find:
            print("Using amino acid instance!\n")
            osaka = util.HelperMTStackNetOsaka(
                datapath=self.args.aa_instance_directory,
                workcell=wrkcell,
                ncpu=self.args.ncpu,
            )
            osaka.set_map_parameters(
                self.args.mapin,
                fix_axis_positions=True,
                fix_origin=True,
                shiftback=False,
            )
            seqlen = 0
            for i in range(0, seq_wrk.size()):
                seqlen += len(seq_wrk[i].sequence) * seqlen_multiplier
            aa_instance_coords = osaka.get_map_coords_from_predicted_instance(
                mode="kmeans",
                mapin_path=self.args.mapin,
                fix_origin=True,
                write_npy=True,
                seqlen=seqlen,
                verbose=self.args.verbose,
            )
            inst_mol = clipper.MiniMol(xwrk.spacegroup, xwrk.cell)
            osaka.map_coords_to_ca_atom(aa_instance_coords, inst_mol)
            util.write_structure(
                inst_mol,
                "instance_mol.pdb",
                cif_format=False,
            )
            print(len(aa_instance_coords))
            if len(aa_instance_coords) > self.args.nfrag:
                cafind = buccaneer.Ca_find(
                    len(aa_instance_coords),
                    self.resol.limit(),
                )
            else:
                cafind = buccaneer.Ca_find(self.args.nfrag, self.resol.limit())
            # cafind = buccaneer.Ca_find(args.nfrag, resol.limit())
            print("set starting instance")
            sys.stdout.flush()
            cafind.set_starting_instance_coords(aa_instance_coords, xwrk,
                                                llktgt,
                                                self.args.findtype, True)
            print("after set starting instance")
            sys.stdout.flush()
        else:
            # prepare search target
            cafind = buccaneer.Ca_find(self.args.nfrag, self.resol.limit())
            if (
                self.use_ml_seq
                and not self.args.use_ml_find
                and self.args.aa_instance_directory != "NONE"
            ):
                osaka = util.HelperMTStackNetOsaka(
                    datapath=self.args.aa_instance_directory,
                    workcell=wrkcell,
                    ncpu=self.args.ncpu,
                )
                osaka.set_map_parameters(
                    self.args.mapin,
                    fix_axis_positions=True,
                    fix_origin=True,
                    shiftback=False,
                )
            else:
                osaka = util.HelperMTStackNetOsaka()

        return cafind, osaka

    @staticmethod
    def merge_models(
        mol_wrk: clipper.MiniMol,
        xwrk: clipper.Xmap_float,
        llkcls: buccaneer.LLK_TargetList,
        seq_wrk: clipper.MMoleculeSequence,
        seq_rel: float,
        verbose: int = 0,
        log: str = None,
    ):
        """Merge models

        Args:
            mol_wrk (clipper.MiniMol): Work MiniMol object
            xwrk (clipper.Xmap_float): Work Xmap object
            llkcls (buccaneer.LLK_TargetList): List of log-likelihood targets
            seq_wrk (clipper.MMoleculeSequence): Work molecule sequence
            seq_rel (float): Sequence reliability
            verbose (int, optional): Verbosity. Defaults to 0.
            log (str, optional): Buccaneer log profiler. Defaults to None.
        """
        camerge = buccaneer.Ca_merge(seq_rel)
        Buccaneer.print_steps_Casummaries(
            "before model merge", len(mol_wrk.select("*/*/CA").atom_list())
        )
        camerge(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
        Buccaneer.print_steps_Casummaries(
            "after model merge", len(mol_wrk.select("*/*/CA").atom_list())
        )
        # print(
        #    " C-alphas before model merge: {0}".format(
        #        mol_wrk.select("*/*/CA").atom_list().size()
        #    )
        # )
        # camerge(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
        # print(
        #    " C-alphas after model merge: {0}".format(
        #        mol_wrk.select("*/*/CA").atom_list().size()
        #    )
        # )
        sys.stdout.flush()
        log.log("MRG ", mol_wrk, verbose > 9)
        sys.stdout.flush()

    @staticmethod
    def filter_model(
        mol_wrk: clipper.MiniMol,
        xwrk: clipper.Xmap_float,
        model_filter_sig: float,
        verbose: int = 0,
        log: str = None,
    ):
        """Filter input model

        Args:
            mol_wrk (clipper.MiniMol): Work MiniMol object
            xwrk (clipper.Xmap_float): Work Xmap object
            model_filter_sig (float): Model filter sigma value
            verbose (int, optional): Verbosity. Defaults to 0.
            log (str, optional): Buccaneer log profiler. Defaults to None.
        """
        Buccaneer.print_steps_Casummaries(
            "before model filter",
            len(mol_wrk.select("*/*/CA").atom_list()),
        )
        buccaneer.Ca_filter.filter(mol_wrk, xwrk, model_filter_sig)
        Buccaneer.print_steps_Casummaries(
            "after model filter",
            len(mol_wrk.select("*/*/CA").atom_list()),
        )
        sys.stdout.flush()
        log.log("FLT ", mol_wrk, verbose > 9)
        sys.stdout.flush()

    @staticmethod
    def augment_input_model(
        mol_wrk: clipper.MiniMol,
        mol_mr: clipper.MiniMol,
        mr_filter_sig: float,
        mr_model_filter: bool,
        mr_model_seed: bool,
    ):
        """Augment input model with MR model

        Args:
            mol_wrk (clipper.MiniMol): Work MiniMol object
            mol_mr (clipper.MiniMol): Molecular replacement MiniMol object
            mr_filter_sig (float): Molecular replacement filter sigma cutoff
            mr_model_filter (bool): Filter using molecular replacement model
            mr_model_seed (bool): Use molecular replacement model as seed
        """
        result = buccaneer.Ca_merge.merge_mr(
            mol_wrk,
            mol_mr,
            mr_filter_sig,
            3,
            mr_model_filter,
            mr_model_seed,
        )
        print(
            " MR residues input: {0}, after filter: {1}, "
            "after seeding: {2}".format(result[0], result[1], result[2])
        )
        sys.stdout.flush()

    def run_buccaneer(
        self,
        llktgt: buccaneer.LLK_map_target,
        llkcls: buccaneer.LLK_TargetList,
        xwrk: clipper.Xmap_float,
        mol_wrk: clipper.MiniMol,
        mol_mr: clipper.MiniMol,
        knownstruc: buccaneer.KnownStructure,
    ):
        """Run Buccaneer

        Args:
            llktgt (buccaneer.LLK_map_target): Log-likelihood map target
            llkcls (buccaneer.LLK_TargetList): List of log-likelihood map targets
            xwrk (clipper.Xmap_float): Work Xmap object
            mol_wrk (clipper.MiniMol): Work MiniMol object
            mol_mr (clipper.MiniMol): Molecular replacement MiniMol object
            knownstruc (buccaneer.KnownStructure): KnownStructure object
        """  # noqa: E501
        seq_wrk = self.get_work_sequence()
        mol_wrk_in = mol_wrk.clone()
        cafind, osaka = self.prepare_ca_find(xwrk, seq_wrk, llktgt)
        # merge models
        if self.args.merge:
            buccaneer.Ca_merge.merge(
                mol_wrk,
                xwrk,
                llkcls,
                seq_wrk,
                self.args.seq_rel,
                stdout=sys.stdout,
            )
            sys.stdout.flush()
            self.log.log("MRG ", mol_wrk, self.args.verbose > 9)
            sys.stdout.flush()
        # filter input model
        if self.args.model_filter:
            Buccaneer.filter_model(mol_wrk, xwrk, self.args.model_filter_sig, self.args.verbose, self.log)
            #self._print_steps_Casummaries(
            #    "before model filter",
            #    len(mol_wrk.select("*/*/CA").atom_list()),
            #)
            #buccaneer.Ca_filter.filter(mol_wrk, xwrk, self.args.model_filter_sig)  # fmt: skip # noqa: E501
            #self._print_steps_Casummaries(
            #    "after model filter",
            #    len(mol_wrk.select("*/*/CA").atom_list()),
            #)
            #sys.stdout.flush()
            #self.log.log("FLT ", mol_wrk, self.args.verbose > 9)
            #sys.stdout.flush()
        # augment input model with mr model
        if self.args.mr_model:
            result = buccaneer.Ca_merge.merge_mr(
                mol_wrk,
                mol_mr,
                self.args.mr_filter_sig,
                3,
                self.args.mr_model_filter,
                self.args.mr_model_seed,
            )
            print(
                " MR residues input: {0}, after filter: {1}, "
                "after seeding: {2}".format(result[0], result[1], result[2])
            )
            sys.stdout.flush()

        # trim input model to only protein residues
        buccaneer.ProteinTools.trim_to_protein(mol_wrk)

        # model building loop
        for cyc in range(0, self.args.ncyc):
            print(f"\nCycle: {cyc+1}\n")
            history = ""
            sys.stdout.flush()
            # find C-alphas by slow likelihood/fast secondary structure search
            if self.args.find:
                # self._print_steps_Casummaries(
                #    "before finding",
                #    len(mol_wrk.select("*/*/CA").atom_list()),
                # )
                cafind(
                    mol_wrk,
                    knownstruc,
                    xwrk,
                    llktgt,
                    self.args.findtype,
                    self.args.modelindex,
                )
                self._print_steps_Casummaries(
                    "after finding",
                    len(mol_wrk.select("*/*/CA").atom_list()),
                )
                sys.stdout.flush()
                self.log.log("FIND", mol_wrk, self.args.verbose > 9)
                util.write_structure(mol_wrk, "find.pdb", cif_format=False)
                sys.stdout.flush()
            # grow Ca
            if self.args.grow:
                buccaneer.Ca_grow.grow(mol_wrk, xwrk, llktgt, 25, ncpus=self.args.ncpu, stdout=sys.stdout)
                self.log.log("GROW", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # join Ca
            if self.args.join:
                buccaneer.Ca_join.join(mol_wrk, 2.0, 2.0, stdout=sys.stdout)
                self.log.log("JOIN", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # link Ca
            if self.args.link:
                buccaneer.Ca_link.link(mol_wrk, xwrk, llktgt, 10.0, 24, stdout=sys.stdout)
                self.log.log("LINK", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
                util.write_structure(mol_wrk, "linked.pdb", cif_format=False)
            # sequence
            if self.args.seqnc:
                # caseq_ml = buccaneer.Ca_sequence_ml(aa_pred, corrections)
                if self.use_ml_seq and self.args.aa_instance_directory != "NONE":  # fmt: skip # noqa: E501
                    # print(f"Grid in corrections : {osaka.map_params.grid}")  # fmt: skip # noqa: E501
                    osaka(mol_wrk, correlation_mode=self.args.correl)
                    util.write_out_model_with_seq(
                        mol_wrk, "testwriteout_model_with_seqprob.pdb"
                    )
                    if not osaka.check_is_seqprob_set():
                        print("SEQPROB not set!")
                buccaneer.Ca_sequence.set_use_ml_sequence_probability(
                    self.use_ml_seq, self.args.sequence_method == "hybrid"
                )
                caseq = buccaneer.Ca_sequence(self.args.seq_rel)
                # caseq.set_use_set_use_ml_sequence_probability(True)
                caseq(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                self._print_steps_Casummaries(
                    "sequenced", caseq.num_sequenced()
                )  # noqa: E501
                history = caseq.format()
                # sys.stdout.flush()
                self.log.log("SEQU", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
                # osaka.print_seq_dat(mol_wrk)
            # correct
            if self.args.corct:
                buccaneer.Ca_correct.correct(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk, 12, stdout=sys.stdout)
                # buccaneer.Ca_correct.correct(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk, 12, show_summary=True)
                # cacor = buccaneer.Ca_correct(12)
                # cacor(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                # self._print_steps_Casummaries("corrected", cacor.num_corrected)
                # sys.stdout.flush()
                self.log.log("CORR", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # filter
            if self.args.filtr:
                buccaneer.Ca_filter.filter(mol_wrk, xwrk, 1.0, True, stdout=sys.stdout)
                # buccaneer.Ca_filter.filter(mol_wrk, xwrk, 1.0, True, show_summary=True)
                # cafiltr = buccaneer.Ca_filter(1.0)
                # cafiltr(mol_wrk, xwrk)
                # self._print_steps_Casummaries(
                #    "after filtering",
                #    len(mol_wrk.select("*/*/CA").atom_list()),
                # )
                # sys.stdout.flush()
                self.log.log("FILT", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # ncs build
            if self.args.ncsbd:
                buccaneer.Ca_ncsbuild.ncsbuild(
                    mol_wrk,
                    xwrk,
                    llkcls.get_vector(),
                    seq_wrk,
                    self.args.seq_rel,
                    1.0,
                    12,
                    stdout=sys.stdout,
                )
                # buccaneer.Ca_ncsbuild.ncsbuild(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk, self.args.seq_rel, 1.0, 12, show_summary=True)
                # cancsbuild = buccaneer.Ca_ncsbuild(self.args.seq_rel, 1.0, 12)
                # cancsbuild(mol_wrk, xwrk, llkcls.get_vector(), seq_wrk)
                # self._print_steps_Casummaries(
                #    "after NCS build",
                #    len(mol_wrk.select("*/*/CA").atom_list()),
                # )
                # sys.stdout.flush()
                self.log.log("NCSB", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # prune
            if self.args.prune:
                buccaneer.Ca_prune.prune(mol_wrk, xwrk, 3.0, stdout=sys.stdout)
                # caprune = buccaneer.Ca_prune(3.0)
                # caprune(mol_wrk, xwrk)
                # self._print_steps_Casummaries(
                #    "after pruning",
                #    len(mol_wrk.select("*/*/CA").atom_list()),
                # )
                # sys.stdout.flush()
                self.log.log("PRUN", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()
            # build
            if self.args.build:
                buccaneer.Ca_build.build(mol_wrk, xwrk, self.args.newrestype, False, stdout=sys.stdout)
                # buccaneer.Ca_build.build(mol_wrk, xwrk, self.args.newrestype, False, show_summary=True)
                # cabuild = buccaneer.Ca_build(self.args.newrestype)
                # cabuild(mol_wrk, xwrk)
                # self._print_steps_Casummaries(
                #    "after rebuilding",
                #    len(mol_wrk.select("*/*/CA").atom_list()),
                # )
                # sys.stdout.flush()
                self.log.log("REBU", mol_wrk, self.args.verbose > 9)
                sys.stdout.flush()

            knownstruc.prune(mol_wrk)
            buccaneer.ProteinTools.split_chains_at_gap(mol_wrk)
            buccaneer.ProteinTools.chain_number(mol_wrk)
            buccaneer.ProteinTools.chain_label(mol_wrk, self.args.chainid_2char)  # noqa: E501
            self.log.log("TDYI", mol_wrk, self.args.verbose > 9)
            if self.args.verbose > 7:
                print(history)
            sys.stdout.flush()
            msg = self.log.summary(mol_wrk, mol_mr, seq_wrk)
            print(f"\nInternal cycle {cyc+1:3d}\n {msg}")
            sys.stdout.flush()
            # file output
            if self.args.xmlout != "NONE":
                self.log.xml(self.args.xmlout)
            sys.stdout.flush()
            # output intermediates (optional)
            if self.args.optemp:
                gfile = clipper.GEMMIfile()
                gfile.export_minimol(mol_wrk)
                gfile.write_pdb(self.args.get_outfile_name(cyc + 1))
                if self.args.write_cif:
                    gfile.write_cif(self.args.get_outfile_name(cyc + 1, True))
            # end cycle
        # move model
        if self.args.fixpos:
            buccaneer.ProteinTools.symm_match(mol_wrk, mol_wrk_in)
        # model tidy
        if self.args.tidy:
            mtidy = buccaneer.ModelTidy(
                1.0,
                12,
                self.args.newrestype,
                self.args.verbose > 6,
            )
            chk = mtidy.tidy(mol_wrk, mol_mr, seq_wrk)
            if not chk:
                print("ModelTidy error")
            sys.stdout.flush()
            self.log.log("TIDY", mol_wrk, self.args.verbose > 9)
            sys.stdout.flush()
        # assign default B-factors to missing values
        default_u_iso = buccaneer.ProteinTools.main_chain_u_mean(mol_wrk_in)
        for c in mol_wrk:
            for r in c:
                for a in r:
                    if clipper.Util.is_nan(a.u_iso):
                        a.u_iso = default_u_iso

        if self.args.newresname != "NONE":
            for c in mol_wrk:
                for r in c:
                    if buccaneer.ProteinTools.residue_index_3(r.type) < 0:
                        r.type = self.args.newresname

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
        buccaneer.ProteinTools.chain_label(mol_wrk, self.args.chainid_2char)
        sys.stdout.flush()
        self.log.log("LABE", mol_wrk, self.args.verbose > 9)
        sys.stdout.flush()

        # write answers
        gfile = clipper.GEMMIfile()
        gfile.export_minimol(mol_wrk)
        gfile.write_pdb(self.args.get_outfile_name())
        if self.args.write_cif:
            gfile.write_cif(self.args.get_outfile_name(cifout=True))

        if self.args.verbose > 8:
            self.log.profile(stdout=sys.stdout)
        # end buccaneer


def main(args: List[str] = None, shell: bool = True):
    """Minimal buccaneer run script

    Args:
        args (List[str], optional): List of arguments. Defaults to None.
        shell (bool, optional): Output stdout stderr to terminal/shell. Defaults to True.
    """  # noqa: E501
    cwd = getcwd()
    logger = log2file(f"{cwd}/pybuc.log", "w", shell=shell)
    sys.stdout = sys.stderr = logger
    print(f"\n### Job started : {datetime.datetime.now()} ###\n")
    raw_args = args or sys.argv[1:]
    parser = BucArgParse(sys.argv[0])
    parsed_args = parser.parse_args(raw_args)
    parser.print_command_with_args()
    buc_params = BuccaneerParams(parsed_args)
    if "CLIBD" not in environ:
        venv_path = environ.get("VIRTUAL_ENV")
        environ["CLIB"] = venv_path + "/lib"
        environ["CLIBD"] = venv_path + "/lib/data"
    
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
    print(f"logfile : {cwd}/pybuc.log\n")
    buc.run()
    print(f"\n### Job ended : {datetime.datetime.now()} ###\n")
    logger.close()
    log2file.reset_stdouterr_to_sys()
    # sys.stdout = sys.__stdout__
    # sys.stderr = sys.__stderr__


if __name__ == "__main__":
    main()
