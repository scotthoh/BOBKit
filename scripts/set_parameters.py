# Argument parser for pybuccaneer
# Author: S.W.Hoh
# 2023 -
# York Structural Biology Laboratory
# The University of York
from bobkit.buccaneer import Ca_find, Ca_prep
from bobkit import __version__

from dataclasses import dataclass
from typing import List
import argparse
from typing import Union


class BucArgParse:
    """
    PyBuccaneer command line parser
    """

    def __init__(self, prog=None):
        self.script = prog
        # self.raw_args = args
        # self.args = args

    def _str2bool(self, arg: str):
        if arg.lower() in ("yes", "true", "t", "y", "1", "on"):
            return True
        elif arg.lower() in ("no", "false", "f", "n", "0", "off"):
            return False
        else:
            raise argparse.ArgumentTypeError("Boolean value expected.")

    def parse_args(self, args=None):
        self.raw_args = args
        parser = argparse.ArgumentParser(prog=self.script)
        # parser.prog = self.script
        parser.description = """Buccaneer wrapped in python."""
        parser.add_argument(
            "-v", "--version", action="version", version=__version__
        )  # noqa 501
        self.add_args(parser)
        self.args = parser.parse_args(args)
        return self.args

    def print_args(self):
        """Print input arguments"""
        if self.args is None:
            print("No arguments given.\n")
            return False
        else:
            print("Input arguments: \n")
            for arg in vars(self.args):
                arg_val = getattr(self.args, arg)
                if isinstance(arg_val, list):
                    print(f" {arg}: {arg_val}")
                if isinstance(arg_val, float):
                    if (arg_val > 0.0) and (not None):
                        print(f" {arg}: {arg_val}")
                if isinstance(arg_val, int):
                    if arg_val > 0:
                        print(f" {arg}: {arg_val}")
                if isinstance(arg_val, str):
                    if arg_val is not None:
                        print(f" {arg}: {arg_val}")
            print("\n")
            return True

    def print_command_with_args(self):
        """Print command line arguments as it is from input"""
        s = "Input arguments : \n\n"
        for arg in self.raw_args:
            s += arg + " "
            # print(f"{arg}", end=" ")
        ##s += "\n"
        print(s)

    def add_args(self, parser):
        """Parse input arguments."""
        _GROUP = parser.add_argument_group("required arguments")
        _GROUP.add_argument(
            "--mtzin",
            type=str,
            dest="mtzin",
            help="Input mtz.",
            metavar="X",
            default=None,
            required=True,
        )

        _GROUP.add_argument(
            "--seqin",
            type=str,
            dest="seqin",
            help="Input sequence file.",
            metavar="X",
            default=None,
            required=True,
        )

        _GROUP.add_argument(
            "--resolution",
            type=float,
            dest="resolution",
            help="High resolution limit in Angstroms.",
            metavar="X",
            default=2.0,
            required=True,
        )

        _GROUP = parser.add_argument_group("optional arguments")
        _GROUP.add_argument(
            "--title",
            type=str,
            metavar="X",
            default="Buccaneer job",
            help="Simple job title.",
        )
        _GROUP.add_argument(
            "--mtzin-ref",
            dest="mtzin_ref",
            type=str,
            default="NONE",
            metavar="X",
            help="MTZ file for reference data.",
        )
        _GROUP.add_argument(
            "--pdbin-ref",
            dest="pdbin_ref",
            type=str,
            default="NONE",
            metavar="X",
            help="PDB file for reference data.",
        )
        _GROUP.add_argument(
            "--mapin",
            dest="mapin",
            type=str,
            default="NONE",
            metavar="X",
            help="Input working map for model building.",
        )
        _GROUP.add_argument(
            "--model",
            dest="model_in",
            type=str,
            default="NONE",
            metavar="X",
            help=("A starting model in PDB or mmCIF format."),
        )
        _GROUP.add_argument(
            "--cycles",
            dest="cycles",
            default=3,
            type=int,
            metavar="X",
            help="Number of cycles to run Buccaneer.",
        )

        _GROUP.add_argument(
            "--model-in-mr",
            dest="model_in_mr",
            type=str,
            default="NONE",
            metavar="X",
            help="Input molecular replacement model.",
        )
        _GROUP.add_argument(
            "--output-model",
            dest="model_out",
            type=str,
            default="buccaneer_build",
            metavar="X",
            help="Filename used for output model.",
        )
        _GROUP.add_argument(
            "--output-sim-map",
            dest="sim_mapout",
            type=str,
            default="NONE",
            metavar="X",
            help="Output filename used for simulated map.",
        )
        _GROUP.add_argument(
            "--xmlout",
            dest="xmlout",
            type=str,
            default="summary.xml",
            metavar="X",
            help="Filename for XML output.",
        )
        _GROUP.add_argument(
            "--colin-ref-fo",
            dest="colin_ref_fo",
            type=str,
            default="NONE",
            metavar="X",
            help="Reference column labels for Fo and sigFo.",
        )
        _GROUP.add_argument(
            "--colin-ref-hl",
            dest="colin_ref_hl",
            type=str,
            default="NONE",
            metavar="X",
            help=(
                "Reference column labels for Hendriksson Lattman coefficients."
            ),  # noqa 501
        )
        _GROUP.add_argument(
            "--colin-fo",
            dest="colin_fo",
            type=str,
            default="NONE",
            metavar="X",
            help=(
                "Column labels for Fo,sigFo of working data. "
                "(e.g 'Fout,sigF')"  # noqa 501
            ),
        )
        _GROUP.add_argument(
            "--colin-hl",
            dest="colin_hl",
            type=str,
            default="NONE",
            metavar="X",
            help="Column labels for Hendrikkson Lattman coefficients.",
        )
        _GROUP.add_argument(
            "--colin-phifom",
            dest="colin_phifom",
            type=str,
            default="NONE",
            metavar="X",
            help=(
                "Column labels for phases and figure of merits. "  # noqa 501
                "(e.g. Phi,FOM)"
            ),
        )
        _GROUP.add_argument(
            "--colin-fc",
            dest="colin_fc",
            type=str,
            default="NONE",
            metavar="X",
            help="Column labels for calculated F and sigF.",
        )
        _GROUP.add_argument(
            "--colin-free",
            dest="colin_free",
            type=str,
            default="NONE",
            metavar="X",
            help="Column label for the free-R flag.",
        )
        _GROUP.add_argument(
            "--merge",
            dest="merge",
            # action=argparse.BooleanOptionalAction, py>3.9
            # action="store_true",
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Merge multi-model.",
        )
        _GROUP.add_argument(
            "--find",
            dest="find",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on find step.",
        )
        _GROUP.add_argument(
            "--grow",
            dest="grow",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on grow step.",
        )
        _GROUP.add_argument(
            "--join",
            dest="join",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on join step.",
        )
        _GROUP.add_argument(
            "--link",
            dest="link",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on link step.",
        )
        _GROUP.add_argument(
            "--sequence",
            dest="sequence",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on sequence step.",
        )
        _GROUP.add_argument(
            "--correct",
            dest="correct",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on correct step.",
        )
        _GROUP.add_argument(
            "--filter",
            dest="filter",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on filter step.",
        )
        _GROUP.add_argument(
            "--ncsbuild",
            dest="ncsbuild",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on non crystallographic symmetry building step.",
        )
        _GROUP.add_argument(
            "--prune",
            dest="prune",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on prune step.",
        )
        _GROUP.add_argument(
            "--rebuild",
            dest="rebuild",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on rebuild step.",
        )
        _GROUP.add_argument(
            "--tidy",
            dest="tidy",
            # action=argparse.BooleanOptionalAction,
            type=self._str2bool,
            default="False",
            metavar="true/false",
            help="Turn on tidy step.",
        )
        _GROUP.add_argument(
            "--fast",
            dest="fast",
            # action=argparse.BooleanOptionalAction,
            action="store_true",
            # type=self._str2bool,
            # default="False",
            # metavar="Fast",
            help="Fast mode will be used when finding C-alphas.",
        )
        _GROUP.add_argument(
            "--build-semet",
            dest="build_semet",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Build MSE",
            help="Build Selenmethionine instead of methionine.",
        )
        _GROUP.add_argument(
            "--anisotropy-correction",
            dest="anisotropy_correction",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Anisotropy correction",
            help="Apply anisotropy correction.",
        )
        _GROUP.add_argument(
            "--fix-position",
            dest="fix_position",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Fix position",
            help="Turn on fix position by symmetry match.",
        )
        _GROUP.add_argument(
            "--output-intermediates",
            dest="output_intermediates",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Output intermediates",
            help="Output intermediates model files.",
        )
        _GROUP.add_argument(
            "--threads",
            dest="nthreads",
            type=int,
            metavar="X",
            default=1,
            help=(
                "The number of threads to use in some internal steps of Buccaneer."  # noqa 501
            ),  # noqa 501
        )
        _GROUP.add_argument(
            "--fragments",
            dest="fragments",
            type=int,
            default=500,
            metavar="X",
            help="Set starting number of fragments.",
        )
        _GROUP.add_argument(
            "--fragments-per-100-residues",
            dest="fragments_per_100_residues",
            type=int,
            default=20,
            metavar="X",
            help="Set number of fragments per 100 residues.",
        )
        _GROUP.add_argument(
            "--free-flag",
            dest="free_flag",
            type=int,
            default=0,
            metavar="X",
            help="Free-R flag index.",
        )
        _GROUP.add_argument(
            "--model-index",
            dest="model_index",
            type=int,
            default=False,
            metavar="X",
            help="Model index.",
        )
        _GROUP.add_argument(
            "--random-seed",
            dest="random_seed",
            type=int,
            metavar="X",
            help="Set random seed.",
        )
        _GROUP.add_argument(
            "--ramachandran-filter",
            dest="ramachandran_filter",
            default="all",
            choices=["all", "helix", "strand", "nonhelix"],
            metavar="all | helix | strand | nonhelix",
            help="Specify ramachandran filter type.",
        )
        _GROUP.add_argument(
            "--known-structure",
            dest="known_structure",
            action="append",
            type=str,
            metavar="X",
            help=(
                "This keyword can be used multiple times. "
                "Use this to specify atoms or chains from the input model to "
                "be retained in the output structure. "
                "Syntax: known-structure coodinateID:radius "
                "(E.g. known-structure /A///:2.0) - keep all atoms in chain A "
                "and don't build within 2A."
            ),
        )
        _GROUP.add_argument(
            "--nonprotein-radius",
            dest="nonprotein_radius",
            type=float,
            metavar="X",
            default=-1.0,
            help=(
                "Set non protein radius, so buccaneer will not build within "
                "the radius."
            ),
        )
        _GROUP.add_argument(
            "--main-chain-likelihood-radius",
            dest="main_chain_likelihood_radius",
            type=float,
            default=4.0,
            metavar="X",
            help="Set main chain likelihood radius.",
        )
        _GROUP.add_argument(
            "--side-chain-likelihood-radius",
            dest="side_chain_likelihood_radius",
            type=float,
            default=5.5,
            metavar="X",
            help="Set side chain likelihood radius.",
        )
        _GROUP.add_argument(
            "--sequence-reliability",
            dest="sequence_reliability",
            type=float,
            default=0.95,
            metavar="X",
            help="Set sequence reliability.",
        )
        _GROUP.add_argument(
            "--offset",
            dest="offset",
            type=float,
            default=0.0,
            metavar="X",
            help="Set offset for map density.",
        )
        _GROUP.add_argument(
            "--correlation-mode",
            dest="correlation_mode",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Correlation mode",
            help="Turn on correlation mode.",
        )
        _GROUP.add_argument(
            "--two-char-chain-id",
            dest="two_char_chain_id",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Two character chain id",
            help="Label chain IDs with two character IDs.",
        )
        _GROUP.add_argument(
            "--no-write-cif",
            dest="no_write_cif",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="True",
            # metavar="Write cif",
            help="Prevent writing output model in mmCIF format.",
        )
        # _GROUP.add_argument(
        #    "--no-write-pdb",
        #    dest="no_write_pdb",
        #    action="store_true",
        #    # action=argparse.BooleanOptionalAction,
        #    # type=self._str2bool,
        #    # default="True",
        #    # metavar="Write pdb",
        #    help="Write output model in PDB format.",
        # )
        _GROUP.add_argument(
            "--model-filter",
            dest="model_filter",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="Model filter",
            help="Turn on to filter input model.",
        )
        _GROUP.add_argument(
            "--model-filter-sigma",
            dest="model_filter_sigma",
            type=float,
            default=3.0,
            metavar="X",
            help="Set filter sigma for model.",
        )
        _GROUP.add_argument(
            "--mr-model",
            dest="mr_model",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="MR model",
            help="Turn on MR model.",
        )
        _GROUP.add_argument(
            "--mr-model-filter",
            dest="mr_model_filter",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="MR Model filter",
            help="Turn on to filter input molecular replacement model.",
        )
        _GROUP.add_argument(
            "--mr-model-seed",
            dest="mr_model_seed",
            action="store_true",
            # action=argparse.BooleanOptionalAction,
            # type=self._str2bool,
            # default="False",
            # metavar="MR Model seed",
            help="Turn on molecular replacement model seed.",
        )
        _GROUP.add_argument(
            "--mr-model-filter-sigma",
            dest="mr_model_filter_sigma",
            type=float,
            default=3.0,
            metavar="X",
            help="Set filter sigma for molecular replacement model.",
        )
        _GROUP.add_argument(
            "--verbose",
            dest="verbose",
            type=int,
            default=0,
            metavar="Verbose",
            help="Set verbosity level.",
        )


@dataclass
class BuccaneerParams:
    """Class to hold arguments needed to run Buccaneer"""

    ncpu: int
    ncyc: int  # = 1
    nfrag: int
    nfragr: int
    freerindex: int  # = 0
    write_pdb: bool
    write_cif: bool
    outfile_name: str
    xmlout: str  # = "output.xml"
    mapout: str  # output simulated map, what buccaneer sees
    pdbin: str
    mtzin: str
    mapin: str
    pdbin_ref: str
    mtzin_ref: str
    pdbin_mr: str

    # columns
    title: str
    ipcol_ref_fo: str  # "FP.F_sigF.F,FP.F_sigF.sigF"
    ipcol_ref_hl: str  # "FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D"
    ipcol_fo: str  # = "Fout, SIGF"
    ipcol_phifom: str  # = "Pout, FOM"
    ipcol_free: str  # = "NONE"
    ipcol_hl: str  # = "NONE"
    ipcol_fc: str  # = "NONE"
    ipseq_wrk: str  # = "NONE"
    newresname: str
    newrestype: str  # = "ALA"

    inresol: float  # = 3.2  # 1.85
    main_tgt_rad: float  # = 4.0
    side_tgt_rad: float  # = 5.5
    seq_rel: float  # = 0.95
    moffset: float  # = 0.0
    nprad: float  # = -1.0
    verbose: int  # = 8
    known_ids: List[str]  # = []
    rama_fltr: Ca_prep.Rama_flt
    merge: bool  # = False multimodel merge
    # buccaneer steps
    find: bool  # = True  # False  #
    grow: bool  # = True  # False  #
    join: bool  # = True  # False  #
    link: bool  # = True  # False  #
    seqnc: bool  # = True  # False  #
    corct: bool  # = True  # False  #
    filtr: bool  # = True  # False  #
    ncsbd: bool  # = True  # False  #
    prune: bool  # = True  # False  #
    build: bool  # = True  # False  #
    semet: bool  # further options
    optemp: bool  # false
    # fast mode Ca_find.TYPE.SECSTRUC else Ca_find.TYPE.LIKELIHOOD
    findtype: bool
    correl: bool
    tidy: bool  # = True
    fixpos: bool  # = False
    doanis: bool  # = True
    optemp: bool  # output intermediate models = True
    chainid_2char: bool  # = False
    profiling: bool
    model_filter: bool
    model_filter_sig: float  # 3.0
    mr_model: bool  # False
    mr_model_filter: bool  # = False
    mr_model_seed: bool  # = False
    mr_filter_sig: float  # = 3.0
    modelindex: int  # 0

    def __init__(self, args: argparse.Namespace = None):
        if args is not None:
            self.set_parameters_with_commandline_args(args)
        else:
            self.set_parameters()

    def set_parameters(
        self,
        title: str = "NONE",
        ncyc: int = 0,
        res_in: float = 1.0,
        mtzin: str = "NONE",
        seqin: str = "NONE",
        mapin: str = "NONE",
        pdbin: str = "NONE",
        pdbin_mr: str = "NONE",
        outfile_name: str = "buccaneer_build.pdb",
        xmlout: str = "summary.xml",
        write_pdb: bool = True,
        write_cif: bool = True,
        verbose: int = 0,
    ):
        self.mtzin = mtzin
        self.mapin = mapin
        self.pdbin = pdbin
        self.pdbin_mr = pdbin_mr
        self.verbose = verbose
        self.title = title
        self.ipseq_wrk = seqin
        self.set_reference_data()
        self.set_ncpu()
        self.set_ncyc(ncyc)
        self.set_steps()
        self.write_pdb = write_pdb
        self.write_cif = write_cif
        self.set_outfile(outfile_name)
        self.set_xmlout(xmlout)
        self.set_write_simulated_map("NONE")
        self.set_resolution(res_in)
        self.set_nfrag()
        self.set_nfragr()
        self.set_freerindex()
        self.set_further_options()
        self.set_rama_fltr()
        self.set_newres()
        self.set_refcol_fo()
        self.set_refcol_hl()
        self.set_wrkcol_fc()
        self.set_wrkcol_phifom()
        self.set_wrkcol_fo()
        self.set_wrkcol_free()
        self.set_wrkcol_hl()
        self.set_nprad()
        self.set_target_radius()
        self.set_sequence_reliability()
        self.set_model_filter()
        self.set_mr_model()
        self.modelindex = 0
        self.moffset = 0.0
        self.set_output_intermediate_models()
        self.set_known_structure()
        if verbose >= 5:
            self.profiling = True
        else:
            self.profiling = False

    def set_parameters_with_commandline_args(self, args: BucArgParse):
        self.mtzin = args.mtzin
        self.mapin = args.mapin
        self.pdbin = args.model_in
        self.pdbin_mr = args.model_in_mr
        self.verbose = args.verbose
        self.title = args.title
        self.ipseq_wrk = args.seqin
        self.set_reference_data(args.pdbin_ref, args.mtzin_ref)
        self.set_ncpu(args.nthreads)
        self.set_ncyc(args.cycles)
        self.set_steps(
            find=args.find,
            grow=args.grow,
            join=args.join,
            link=args.link,
            seqnc=args.sequence,
            corct=args.correct,
            filtr=args.filter,
            ncsbd=args.ncsbuild,
            prune=args.prune,
            build=args.rebuild,
        )
        self.write_pdb = True
        self.write_cif = not args.no_write_cif
        self.set_outfile(args.model_out)
        self.set_xmlout(args.xmlout)
        self.set_write_simulated_map(args.sim_mapout)
        self.set_resolution(args.resolution)
        self.set_nfrag(args.fragments)
        self.set_nfragr(args.fragments_per_100_residues)
        self.set_freerindex(args.free_flag)
        self.set_further_options(
            merge=args.merge,
            semet=args.build_semet,
            fast=args.fast,
            correl=args.correlation_mode,
            tidy=args.tidy,
            fixpos=args.fix_position,
            doanis=args.anisotropy_correction,
            optemp=args.output_intermediates,
            chainid_2char=args.two_char_chain_id,
        )
        self.set_rama_fltr_str(args.ramachandran_filter)
        self.set_newres()
        self.set_refcol_fo(args.colin_ref_fo)
        self.set_wrkcol_fc(args.colin_fc)
        self.set_wrkcol_phifom(args.colin_phifom)
        self.set_wrkcol_fo(args.colin_fo)
        self.set_wrkcol_free(args.colin_free)
        self.set_wrkcol_hl(args.colin_hl)
        self.set_nprad(args.nonprotein_radius)
        self.set_target_radius(
            main_tgt_rad=args.main_chain_likelihood_radius,
            side_tgt_rad=args.side_chain_likelihood_radius,
        )
        self.set_sequence_reliability(args.sequence_reliability)
        self.set_model_filter(args.model_filter, args.model_filter_sigma)
        self.set_mr_model(
            mr_model=args.mr_model,
            mr_model_filter=args.mr_model_filter,
            mr_model_seed=args.mr_model_seed,
            mr_filter_sig=args.mr_model_filter_sigma,
        )
        self.modelindex = args.model_index
        self.set_output_intermediate_models(args.output_intermediates)
        self.set_known_structure(args.known_structure)
        self.verbose = args.verbose
        if self.verbose >= 5:
            self.profiling = True
        else:
            self.profiling = False

    def set_outfile(self, outfile_name: str = "buccaneer_build"):
        self.outfile_name = outfile_name

    def set_write_simulated_map(self, mapout: str = "NONE"):
        self.mapout = mapout

    def set_xmlout(self, xmlout: str = "NONE"):
        self.xmlout = xmlout

    def set_steps(
        self,
        find: bool = False,
        grow: bool = False,
        join: bool = False,
        link: bool = False,
        seqnc: bool = False,
        corct: bool = False,
        filtr: bool = False,
        ncsbd: bool = False,
        prune: bool = False,
        build: bool = False,
    ):
        """turn on/off steps in buccaneer"""
        self.find = find
        self.grow = grow
        self.join = join
        self.link = link
        self.seqnc = seqnc
        self.corct = corct
        self.filtr = filtr
        self.ncsbd = ncsbd
        self.prune = prune
        self.build = build

    def check_steps(self):
        """Turn all buccaneer steps on if all are off"""
        if [
            x
            for x in (
                self.find,
                self.grow,
                self.join,
                self.link,
                self.seqnc,
                self.corct,
                self.filtr,
                self.ncsbd,
                self.prune,
                self.build,
                self.tidy,
            )  # noqa 501
            if not x
        ]:
            self.set_steps(
                find=True,
                grow=True,
                join=True,
                link=True,
                seqnc=True,
                corct=True,
                filtr=True,
                ncsbd=True,
                prune=True,
                build=True,
            )
            self.tidy = True

    def set_further_options(
        self,
        merge: bool = False,
        semet: bool = False,
        fast: bool = False,
        correl: bool = False,
        tidy: bool = False,
        fixpos: bool = False,
        doanis: bool = False,
        optemp: bool = False,
        chainid_2char: bool = False,
    ):
        self.merge = merge
        self.semet = semet
        self.set_findtype(fast)
        self.correl = correl
        self.tidy = tidy
        self.fixpos = fixpos
        self.doanis = doanis
        self.optemp = optemp
        self.chainid_2char = chainid_2char

    def set_resolution(self, resol: float = 1.0):
        self.inresol = resol

    def set_rama_fltr(
        self,
        rama_fltr: Ca_prep.Rama_flt = Ca_prep.rama_flt_all,  # noqa 501
    ):
        self.rama_fltr = rama_fltr

    def set_rama_fltr_str(self, rama_fltr: str = "all"):
        self.rama_fltr = {
            "all": Ca_prep.rama_flt_all,
            "helix": Ca_prep.rama_flt_helix,
            "strand": Ca_prep.rama_flt_strand,
            "nonhelix": Ca_prep.rama_flt_nonhelix,
        }.get(rama_fltr, Ca_prep.rama_flt_all)

    def set_ncpu(self, cpu: int = 0):
        self.ncpu = cpu

    def set_nfrag(self, nfrag: int = 500):
        self.nfrag = nfrag

    def set_nfragr(self, nfragr: int = 20):
        self.nfragr = nfragr

    def set_freerindex(self, index: int = 0):
        self.freerindex = index

    def set_newres(self, name: str = "UNK", type: str = "ALA"):
        self.newresname = name
        self.newrestype = type

    def set_wrkcol_fo(self, colin: str = "NONE"):
        self.ipcol_fo = colin

    def set_wrkcol_phifom(self, colin: str = "NONE"):
        self.ipcol_phifom = colin

    def set_wrkcol_free(self, colin: str = "NONE"):
        self.ipcol_free = colin

    def set_wrkcol_hl(self, colin: str = "NONE"):
        self.ipcol_hl = colin

    def set_wrkcol_fc(self, colin: str = "NONE"):
        self.ipcol_fc = colin

    def set_refcol_fo(self, colin: str = "/*/*/FP"):
        self.ipcol_ref_fo = colin

    def set_refcol_hl(self, colin: str = "/*/*/FC"):
        self.ipcol_ref_hl = colin

    def set_ncyc(self, ncyc: int = 3):
        self.ncyc = ncyc

    def set_nprad(self, nprad: float = -1.0):
        self.nprad = nprad

    def set_target_radius(
        self,
        main_tgt_rad: float = 4.0,
        side_tgt_rad: float = 5.5,
    ):
        self.main_tgt_rad = main_tgt_rad
        self.side_tgt_rad = side_tgt_rad

    def set_sequence_reliability(self, seq_rel: float = 0.95):
        self.seq_rel = seq_rel

    def set_findtype(self, fast: bool = False):
        self.findtype = (
            Ca_find.TYPE.SECSTRUC if fast else Ca_find.TYPE.LIKELIHOOD
        )  # noqa 501

    def set_model_filter(self, filter: bool = False, sigma: float = 3.0):
        self.model_filter = filter
        self.model_filter_sig = sigma

    def set_mr_model(
        self,
        mr_model: bool = False,
        mr_model_filter: bool = False,
        mr_model_seed: bool = False,
        mr_filter_sig: float = 3.0,
    ):
        self.mr_model = mr_model
        self.mr_model_filter = mr_model_filter
        self.mr_model_seed = mr_model_seed
        self.mr_filter_sig = mr_filter_sig

    def set_output_intermediate_models(self, optemp: bool = False):
        self.optemp = optemp

    def get_outfile_name(self, cycle: int = -1, cifout: bool = False):
        if cycle != -1:
            cyc = str(cycle)
        else:
            cyc = ""
        if self.outfile_name.rfind(".") == -1:
            out_fname = self.outfile_name
        else:
            out_fname = self.outfile_name[: self.outfile_name.rfind(".")]
        if cifout:
            return out_fname + cyc + ".cif"
        else:
            return out_fname + cyc + ".pdb"

    def set_reference_data(self, pdbref: str = "NONE", mtzref: str = "NONE"):
        self.pdbin_ref = pdbref
        self.mtzin_ref = mtzref

    def set_known_structure(self, ids: Union[List[str], None] = []):
        if isinstance(ids, list):
            if len(ids) != 0:
                self.known_ids = ids
            else:
                self.known_ids = []
        else:
            self.known_ids = []


# if __name__ == "__main__":
#    import sys
#
#    parser = BucArgParse()  # args=args)
#    parser.parse_args(sys.argv[1:])
#    parser.print_args()
#    parser.print_command_with_args()
