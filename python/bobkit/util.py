from __future__ import annotations
from typing import List as _List
from bobkit._util import *
import scipy.spatial as _spatial
import itertools as _IT
import sys

__all__ = [
    "read_structure",
    "write_structure",
]

from bobkit.buccaneer import (
    Ca_build as _Ca_build,
    Ca_group as _Ca_group,
    Ca_sequence as _Ca_sequence,
    ProteinTools as _ProteinTools,
)

from bobkit.clipper import (
    Cell as _Cell,
    Grid as _Grid,
    Grid_sampling as _Grid_sampling,
    Grid_range as _Grid_range,
    Coord_orth as _Coord_orth,
    Coord_grid as _Coord_grid,
    Coord_map as _Coord_map,
    MiniMol as _MiniMol,
    MModel as _MModel,
    MChain as _MChain,
    MResidue as _MRes,
    MAtom as _MAtom,
    Xmap_base as _Xmap_base,
    Xmap_float as _Xmap_float,
    Xmap_int as _Xmap_int,
    Vec3_double as _Vec3_double,
    Property_int as _Property_int,
    Ca_sequence_data_single_atom as _Ca_sequence_data_single_atom,
    Property_sequence_data as _Property_sequence_data,
    Property_sequence_data_single_atom as _Property_sequence_data_single_atom,
    NXmap_base as _NXmap_base,
    NXmap_float as _NXmap_float,
    EDcalc_mask_float as _EDcalc_mask_float,
    RTop_orth as _RTop_orth,
    Mat33_double as _Mat33_double,
    Spacegroup as _Spacegroup,
    Coord_frac as _Coord_frac,
    Resolution as _Resolution,
    NX_operator as _NX_operator
)

import numpy as _np
import gemmi as _gemmi
from sklearn.cluster import HDBSCAN as _HDBSCAN
from sklearn.cluster import KMeans as _KMeans
from sklearn.neighbors import KDTree
from typing import Sequence as _Sequence, Union as _Union, Literal as _Literal
from os import path as _path
from dataclasses import dataclass as _dataclass, field as _field
import re as _re
#from multiprocessing import Process as _Process
#from multiprocessing.managers import BaseManager as _BaseManager
from threading import Thread
from concurrent.futures import ThreadPoolExecutor as _ThreadPoolExecutor
_ProteinTools()

AA_TO_CH_DICT = {
    "ALA": 0,
    "GLY": 1,
    "ILE": 2,
    "LEU": 3,
    "PRO": 4,
    "VAL": 5,
    "PHE": 6,
    "TRP": 7,
    "TYR": 8,
    "ASP": 9,
    "GLU": 10,
    "ARG": 11,
    "HIS": 12,
    "LYS": 13,
    "SER": 14,
    "THR": 5,
    "CYS": 15,
    "MET": 16,
    "ASN": 9,
    "GLN": 10,
}

CH_TO_AA_DICT = {
    0: "ALA",
    1: "GLY",
    2: "ILE",
    3: "LEU",
    4: "PRO",
    5: "VAL/THR",
    6: "PHE",
    7: "TRP",
    8: "TYR",
    9: "ASP/ASN",
    10: "GLU/GLN",
    11: "ARG",
    12: "HIS",
    13: "LYS",
    14: "SER",
    15: "CYS",
    16: "MET",
}

MLI_TO_BUC_DICT = {
    0: 0,           # ALA
    1: 7,           # GLY
    2: 9,           # ILE
    3: 10,          # LEU
    4: 14,          # PRO
    5: [16, 19],    # VAL/THR
    6: 13,          # PHE
    7: 17,          # TRP
    8: 18,          # TYR
    9: [2, 3],      # ASP/ASN
    10: [5, 6],     # GLU/GLN
    11: 1,          # ARG
    12: 8,          # HIS
    13: 11,         # LYS
    14: 15,         # SER
    15: 4,          # CYS
    16: 12,         # MET
}

BUC_TO_MLI_DICT = {
    0: 0,       # ALA
    1: 11,      # ARG
    2: 9,       # ASN
    3: 9,       # ASP
    4: 15,      # CYS
    5: 10,      # GLN
    6: 10,      # GLU
    7: 1,       # GLY
    8: 12,      # HIS
    9: 2,       # ILE
    10: 3,      # LEU
    11: 13,     # LYS
    12: 16,     # MET
    13: 6,      # PHE
    14: 4,      # PRO
    15: 14,     # SER
    16: 5,      # THR
    17: 7,      # TRP
    18: 8,      # TYR
    19: 5,      # VAL
    12: 16,     # MSE
}


@_dataclass
class MapParameters:
    """Class for keeping the map used in ML spacing, origin and correction
    values
    """

    workcell: _np.ndarray = _field(
        default_factory=lambda: _np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])
    )
    workgrid: _np.ndarray = _field(
        default_factory=lambda: _np.ones(3, dtype=_np.int64)
    )
    cell: _np.ndarray = _field(
        default_factory=lambda: _np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])
    )
    grid: _np.ndarray = _field(
        default_factory=lambda: _np.ones(3, dtype=_np.int64)
    )
    grid_asu: _np.ndarray = _field(
        default_factory=lambda: _np.ones(3, dtype=_np.int64)
    )
    origin: _np.ndarray = _field(
        default_factory=lambda: _np.zeros(3, dtype=_np.int64)
    )
    spacing: _np.ndarray = _field(default_factory=lambda: _np.ones(3))
    ncorrect: int = 0
    fix_origin: bool = False
    fix_axis_positions: _np.ndarray = _field(
        default_factory=lambda: _np.array([0, 1, 2], dtype=int)
    )
    # rotation and translation from transformed back to initial model/xtal space
    rotation: _np.ndarray = _field(
        default_factory=lambda: _np.zeros((3, 3))
    )
    translation: _np.ndarray = _field(default_factory=lambda: _np.zeros(3))

    def set_param_values_3(self, param: _np.ndarray, x: int, y: int, z: int):
        param[0] = x
        param[1] = y
        param[2] = z

    def set_grid(self, nu: int, nv: int, nw: int):
        self.grid[0] = nu
        self.grid[1] = nv
        self.grid[2] = nw

    def set_grid_asu(self, nu: int, nv: int, nw: int):
        self.grid_asu[0] = nu
        self.grid_asu[1] = nv
        self.grid_asu[2] = nw

    def set_origin(self, x: int, y: int, z: int):
        self.origin[0] = x
        self.origin[1] = y
        self.origin[2] = z

    def set_spacing(self, x: float, y: float, z: float):
        self.spacing[0] = x
        self.spacing[1] = y
        self.spacing[2] = z

    def set_translation(self, x: float, y: float, z: float):
        self.translation[0] = x
        self.translation[1] = y
        self.translation[2] = z


@_dataclass
class MLOutputFilenames:
    """Data class for quick storing of filenames used in MTStackNet
    """
    predicted_seq: str = "pred.npy"
    predicted_instance: str = "inst_pred.npy"
    density_array: str = "density.npy"
    cutout_input_map: str = "cutout_2mFo-DFc.ccp4"
    cutout_model: str = "cutout.pdb"


class HelperMTStackNetOsaka:
    """Helper class with methods and static methods to process the outputs from
    single shot CNN to predict amino acid instance and sequence from map
    segmentation by Osaka group.
    """
    def __init__(
        self, datapath: str = None,
        filenames: MLOutputFilenames = MLOutputFilenames(),
        workcell: _Union[_Sequence[float], _np.ndarray] = None,
        ncpu: int = 1,
        kmeans_ninit: int = 6,
        verbose: bool = False,
    ):
        """Initialise object

        Args:
            datapath (str, optional): Path to ML output. Defaults to None.
            filenames (MLOutputFilenames, optional): Input and output filenames from MTStackNet
            workcell (Union[Sequence[float], _np.ndarray], optional): Work cell parameters. Defaults to None.
            kmeans_ninit: (int, optional): Number of times K-means algorithm is run with different centroid seeds, 6 seems to be stable. Defaults to 6.
            verbose: (bool, optional): Flag to turn on/off debugging printouts. Defaults to False.

        Raises:
            ValueError: Workcell array is expected to have 6 elements
        """  # noqa: E501
        self.datapath = datapath
        self.map_params = MapParameters()
        self.filenames = filenames
        self.verbose = verbose
        if workcell is not None:
            tmp = _np.asarray(workcell, dtype=float)
            if tmp.shape != (6,):
                raise ValueError("Expected array with 6 elements")
            self.map_params.workcell = tmp
        if self.datapath is not None:
            self.__sequence_array = _np.load(f"{self.datapath}/{self.filenames.predicted_seq}")  # noqa:E501
        else:
            self.__sequence_array = _np.full((1, 1), None, dtype=object)
        self.__ncpu = ncpu
        self.__mol = None
        self.__done = []
        self.xwrk_int = _Xmap_int()
        self.kmeans_ninit = kmeans_ninit
        
    def _debug(self, msg: str, debug: bool):
        """Prints debug message

        Args:
            msg (str): Debug message
        """
        if debug:
            print(f"DEBUG>> {msg}")
            sys.stdout.flush()

    def _debug_list(self, msg: str, data: _Union[_Sequence, _np.ndarray], debug: bool):
        """Prints debug message with data

        Args:
            msg (str): Debug message
            data (Union[Sequence, _np.ndarray]): Data values
        """
        if debug:
            print(f"DEBUG>> {msg} : ", end="")
            for v in data:
                print(f"{v}, ", end="")
            print("\n")
            sys.stdout.flush()

    def __call__(self, mol: _Minimol, correlation_mode: bool = False, single_atom: bool = False):  # , llktargets: LLK_TargetList):  # noqa: E501
        """Apply sequence probabilities from MTStackNet output to model properties

        Args:
            mol (_Minimol): _description_
            xmap (_Xmap): _description_
        """  # noqa: E501
        self.__mol = mol
        # print(f"grid {self.__corrections.grid}")
        # print(f"asu grid {self.__corrections.grid_asu}")
        # print(f"function in grid {grid}")
        grid_samp = _Grid_sampling(self.map_params.grid[0],
                                   self.map_params.grid[1],
                                   self.map_params.grid[2])
        if self.__ncpu > 1:
            #CustomManager.register("MiniMol", _MiniMol)
            #CustomManager.register("list", list)
            #manager = CustomManager()
            #manager.start()
            #shared_mol = manager.MiniMol(mol.spacegroup, mol.cell)
            #shared_mol.model = mol.model.clone()
            #print(f"shared_mol size {self.__mol.size()}")
            #print(f"mol size {mol.size()}")
            threads = []
            #shared_done = manager.list([False] * self.__mol.size())
            self.__done = [False] * self.__mol.size()
            chunk_size = max(1, self.__mol.size() // (self.__ncpu * 4))
            for i in range(0, self.__mol.size(), chunk_size):
                start_chn = i
                end_chn = start_chn + chunk_size
                if end_chn > self.__mol.size():
                    end_chn = self.__mol.size()
                p = Thread(
                    target=self.prepare_scores,
                    args=(start_chn, end_chn, grid_samp, correlation_mode, single_atom), #, shared_mol, shared_done),
                )
                threads.append(p)
                p.start()
            # wait
            for p in threads:
                p.join()
            
            #self.__mol.model = shared_mol.model.clone()
            #mol.model = shared_mol.model.clone()
            #print(f"shared mol size {shared_mol.size()}")
            #print(f"mol size {mol.size()}")
            #print(f"self.__mol size {self.__mol.size()}")
            if self.check_is_seqprob_set():
                print("inside, set")
            else:
                print("inside, not set")
            #print(self.__done)
        else:
            self.__done = [False] * self.__mol.size()
            self.prepare_scores(0, self.__mol.size(), grid_samp, correlation_mode, single_atom),# self.__mol, self.__done)
            for i in self.__done:
                if not i:
                    print("Did not manage to fully assign sequence probability!")  # noqa: E501
        
        for chn in mol:
                for res in chn:
                    if not res.exists_property("SEQPROB"):
                        print("NO SEQPROB")

    #def _run_single_atom_assignment(self, shiftback: bool = False, correlation_mode: bool = False):
    #    grid_samp = _Grid_sampling(self.map_params.grid[0],
    #                               self.map_params.grid[1],
    #                               self.map_params.grid[2])
    #    if self.__ncpu > 1:
    #        processes = []
    #        self.__done = [False] * self.__mol.size()
    #        chunk_size = max(1, self.__mol.size() // (self.__ncpu * 4))
    #        for i in range(0, self.__mol.size(), chunk_size):
    #            start_chn = i
    #            end_chn = start_chn + chunk_size
    #            if end_chn > self.__mol.size():
    #                end_chn = self.__mol.size()
    #            p = _Process(
    #                target=self.prepare_scores,
    #                args=(start_chn, end_chn, grid_samp, shiftback, correlation_mode, True),
    #            )
    #            processes.append(p)
    #            p.start()
    #        # wait
    #        for p in processes:
    #            p.join()
    #        if self.check_is_seqprob_set():
    #            print("inside, set")
    #        else:
    #            print("inside, not set")
    #        print(self.__done)
    #    else:
    #        self.__done = [False] * self.__mol.size()
    #        self.prepare_scores(0, self.__mol.size(), grid_samp, shiftback, correlation_mode, True)
    #        for i in self.__done:
    #            if not i:
    #                print("Did not manage to fully assign sequence probability!")  # noqa: E501

    @staticmethod
    def ml_aa2bucindex(ml_aa_ind: int):
        return MLI_TO_BUC_DICT[ml_aa_ind]

    @staticmethod
    def buc_aa2mlindex(buc_aa_ind: int):
        return BUC_TO_MLI_DICT[buc_aa_ind]

    @staticmethod
    def mlindex2res(ml_aa_ind: int):
        return CH_TO_AA_DICT[ml_aa_ind]

    @staticmethod
    def bucindex2mlres(ml_aa_ind: int):
        return CH_TO_AA_DICT[HelperMTStackNetOsaka.buc_aa2mlindex(ml_aa_ind)]

    @staticmethod
    def get_probability_values(cg: _Coord_grid, seq_arrays: _np.ndarray, correlation_mode: bool = False):
        """Get probability values of all sequence from predicted sequence array.
        For compatibility with Buccaneer's sequencing scores, a list of negative probability values (correlation mode) or
        1.-probability values (non-correlation mode) is returned as ordered in AA_TO_CH_DICT.

        Args:
            cg (_Coord_grid): Coordinate grid position
            seq_arrays (_np.ndarray): Array containing probability values of sequence predicted

        Returns:
            probvals: A list of probability values as ordered in AA_TO_CH_DICT
        """  # noqa: E501
        # pad the 0. scores a bit in correlation mode
        if correlation_mode:
            probvals = -(seq_arrays[:-1, cg.u, cg.v, cg.w] + 0.1)
        else:
            probvals = 1. - seq_arrays[:-1, cg.u, cg.v, cg.w]
        probvals[5] /= 2.0
        probvals[9] /= 2.0
        probvals[10] /= 2.0
        buc_indices = [BUC_TO_MLI_DICT[i] for i in BUC_TO_MLI_DICT.keys()]
        #print('###')
        #print(buc_indices)
        #print(probvals)
        #print(probvals[buc_indices])
        return probvals[buc_indices]

    @staticmethod
    def get_probability_value(ind: int, cg: _Coord_grid, seq_array: _np.ndarray, correlation_mode: bool = False):  # noqa: E501
        """Get probability value from a selected predicted sequence array.
        For compatibility with Buccaneer's sequencing scores, negative probability value (correlation mode) or
        1.-probability value (non-correlation mode) is returned.

        Args:
            ind (int): Index to amino acid according to CH_TO_AA_DICT
            pos (_Coord_orth): Orthogonal coordinates of the atom
            seq_array (_np.ndarray): Array containing probability for the target amino acid
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False
            
        Returns:
            probval: Negative probability value in correlation mode, and 1.-probability for non-correlation mode for compatibility with buccaneer sequencing scores
        """  # noqa: E501
        # print(cg)
        # pad the 0. scores a bit in correlation mode
        if correlation_mode:
            probval = -(seq_array[cg.u, cg.v, cg.w] + 0.1)
        else:
            probval = 1. - seq_array[cg.u, cg.v, cg.w]
        if ind in [5, 9, 10]:
            return probval / 2.0
        else:
            return probval
    
    def get_map_data(self, cg: _Coord_grid):
        return self.xwrk_int.get_data(cg)

    def get_probability_values_coord_orth(self, pos: _Coord_orth, correlation_mode: bool = False, fix_origin=True, shiftback=False):  # noqa: E501
        """Get probability value from a selected predicted sequence array.
        For compatibility with Buccaneer's sequencing scores, negative probability value (correlation mode) or
        1.-probability value (non-correlation mode) is returned.

        Args:
            pos (_Coord_orth): Orthogonal coordinates of the atom
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
            fix_origin (bool, optional): Flag to fix origin. Defaults to True.
            shiftback (bool, optional): Flag to apply translation. Defaults to False.

        Returns:
            probval: Negative probability value in correlation mode, and 1.-probability for non-correlation mode for compatibility with buccaneer sequencing scores
        """  # noqa: E501
        # print(cg)
        grid = _Grid_sampling(self.map_params.grid[0],
                                   self.map_params.grid[1],
                                   self.map_params.grid[2])
        cell = _Cell(self.map_params.cell)
        #if shiftback:
        pos = pos@self.map_params.rotation + _Coord_orth(self.map_params.translation)
        cg = pos.coord_frac(cell).coord_grid(grid)
        #if self.map_params.fix_origin:
        cg = cg - _Coord_grid(self.map_params.origin[0],
                              self.map_params.origin[1],
                              self.map_params.origin[2])
        cg = cg.unit(_Grid_sampling(self.map_params.grid[0],
                                    self.map_params.grid[1],
                                    self.map_params.grid[2]))
        # pad the 0. scores a bit in correlation mode
        if correlation_mode:
            probvals = -(self.__sequence_array[:-1, cg.u, cg.v, cg.w] + 0.1)
        else:
            probvals = 1. - self.__sequence_array[:-1, cg.u, cg.v, cg.w]
        probvals[5] /= 2.0
        probvals[9] /= 2.0
        probvals[10] /= 2.0
        buc_indices = [BUC_TO_MLI_DICT[i] for i in BUC_TO_MLI_DICT.keys()]
        return probvals[buc_indices]

    @staticmethod
    def print_seq_dat(mol: _MiniMol, show_dat: bool = True, show_prob: bool = True):
        """Print SEQDAT or/and SEQPROB values if exists

        Args:
            mol (_MiniMol): MiniMol object
            show_dat (bool, optional): Flag to print SEQDAT. Defaults to True.
            show_prob (bool, optional): Flag to print SEQPROB. Defaults to True.
        """
        for chn in range(0, mol.size()):
            for r in range(0, mol[chn].size()):
                res = mol[chn][r]
                ca = _Ca_group(res)
                if not ca.is_null():
                    if show_dat:
                        if res.exists_property("SEQDAT"):
                            seqdat = res.get_property("SEQDAT").value
                            print(f"{seqdat.ca.coord_ca}")
                            print(f"SD> {seqdat.data}")
                    if show_prob:
                        if res.exists_property("SEQPROB"):
                            seqprob = res.get_property("SEQPROB").value
                            print(f"SP> {seqprob.data}")

    def prepare_score_single_atom(
        self,
        res: _MRes,
        grid: _Grid_sampling,
        #shiftback: bool = False,
        correlation_mode: bool = False,
    ):
        """Prepare and set scores property for each residue

        Args:
            res (_MRes): MResidue object
            grid (_Grid_sampling): Grid sampling
            shiftback (bool, optional): Flag to apply translation. Defaults to False.
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
        """
        cached = False
        ca = res["CA"]
        #cell = _Cell(self.map_params.cell)
        workcell = _Cell(self.map_params.workcell)
        workgrid = _Grid(*self.map_params.workgrid)
        # print(cell.format())
        pos = _Coord_orth.null()
        nxm_grid = _Grid(self.map_params.grid_asu[0], self.map_params.grid_asu[1], self.map_params.grid_asu[2])
        if not ca.pos.is_null():
            if res.exists_property("SEQPROB"):
                seqprob_val = res.get_property("SEQPROB").value
                if ((ca.pos - seqprob_val.pos).lengthsq() < 1.0e-3):
                    cached = True
            pos = ca.pos
            if not cached:
                if res.exists_property("SEQPROB"):
                    res.delete_property("SEQPROB")
                # transform from buccaneer xmap space to grid space in arrays
                #pos = _Coord_orth(pos.array@self.map_params.rotation + self.map_params.translation)
                #cg = pos.coord_frac(cell).coord_grid(grid)
                cg = pos.coord_frac(workcell).coord_grid(workgrid)
                nxm_cg = nxm_grid.deindex(self.xwrk_int.get_data(cg))
                #if self.map_params.fix_origin:
                #    cg = cg - _Coord_grid(self.map_params.origin[0],
                #                          self.map_params.origin[1],
                #                          self.map_params.origin[2])
                #cg = cg.unit(_Grid_sampling(self.map_params.grid_asu[0],
                #                            self.map_params.grid_asu[1],
                #                            self.map_params.grid_asu[2]))
                scores = self.get_probability_values(nxm_cg, self.__sequence_array, correlation_mode)
                seqprob_val = _Ca_sequence_data_single_atom()
                seqprob_val.pos = ca.pos
                seqprob_val.data = scores
                res.set_property("SEQPROB", _Property_sequence_data_single_atom(seqprob_val))

    def prepare_score(
        self,
        res: _MRes,
        grid: _Grid_sampling,
        #shiftback: bool = False,
        correlation_mode: bool = False,
    ):
        """Prepare and set scores property for each residue

        Args:
            res (_MRes): MResidue object
            grid (_Grid_sampling): Grid sampling
            shiftback (bool, optional): Flag to apply translation. Defaults to False.
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
        """
        cached = False
        ca = _Ca_group(res)
        workcell = _Cell(self.map_params.workcell)
        workgrid = _Grid(*self.map_params.workgrid)
        # print(cell.format())
        pos = _Coord_orth.null()
        nxm_grid = _Grid(self.map_params.grid_asu[0], self.map_params.grid_asu[1], self.map_params.grid_asu[2])
        if not ca.is_null():
            if res.exists_property("SEQPROB"):
                seqprob_val = res.get_property("SEQPROB").value
                if (
                    (ca.coord_n - seqprob_val.ca.coord_n).lengthsq() < 1.0e-3
                    and (ca.coord_ca - seqprob_val.ca.coord_ca).lengthsq() < 1.0e-3
                    and (ca.coord_c - seqprob_val.ca.coord_c).lengthsq() < 1.0e-3
                ):
                    cached = True
            pos = ca.coord_ca
            #pos = ca.coord_cb
            #if pos.is_null():
            #    pos = res["CA"].pos
            #else:
            #    a = res.lookup(" CA ")
            #    if a >= 0:
            #        pos = res["CA"].pos
            # print(pos)
            if not cached:
                if res.exists_property("SEQPROB"):
                    res.delete_property("SEQPROB")
                # ntyp = llktargets_size
                # scores = [0.0] * ntyp
                # transform from buccaneer build xmap space to grid array space
                #pos = _Coord_orth(pos.array@self.map_params.rotation + self.map_params.translation)
                #poscb = poscb + _Coord_orth(self.map_params.shiftback)
                # don't have to transform positions as we have xmap to nxmap indices mapping
                cg = pos.coord_frac(workcell).coord_grid(workgrid)
                nxm_cg = nxm_grid.deindex(self.xwrk_int.get_data(cg))
                #nxm_cg = nxm_grid.deindex(self.get_map_data(cg))
                #cg = poscb.coord_frac(workcell).coord_grid(workgrid)
                #nxm_cg = nxm_grid.deindex(self.xwrk_int.get_data(cg))
                #self._debug_list("cg : ", cg, True)
                #self._debug_list("index: ", [self.xwrk_int.get_data(cg)], True)
                #self._debug_list("nxm_cg : ", nxm_cg, True)
                #self._debug_list("check pos/cg", [pos, cg])
                #cgb = pos.coord_frac(cell).coord_grid(grid)
                #if self.map_params.fix_origin:
                # cg = cg - _Coord_grid(self.map_params.origin[0],
                #                       self.map_params.origin[1],
                #                       self.map_params.origin[2])
                #     #cgb = cgb - _Coord_grid(self.map_params.origin[0],
                #     #                      self.map_params.origin[1],
                #     #                      self.map_params.origin[2])
                # cg = cg.unit(_Grid_sampling(self.map_params.grid[0],
                #                             self.map_params.grid[1],
                #                             self.map_params.grid[2]))
                #self._debug_list("  check cg_unit", [cg])
                #cgb = cgb.unit(_Grid_sampling(self.map_params.grid[0],
                #                            self.map_params.grid[1],
                #                            self.map_params.grid[2]))
                scores = self.get_probability_values(nxm_cg, self.__sequence_array, correlation_mode)
                #scores += self.get_probability_values(cgb, self.__sequence_array, correlation_mode)
                #scores /= 2
                # for t in range(0, ntyp):
                #    ml_aa_ind = self.buc_aa2mlindex(t)
                #    scores[t] = self.get_probability_value(
                #        ml_aa_ind,
                #        pos, #ca.coord_ca,
                #        self.__sequence_array[ml_aa_ind],
                #        self.map_params,
                #        cell,
                #        grid,
                #        shiftback,
                #        # self.__fix_axis_positions,
                #    )
                self._debug_list(res, scores, self.verbose)
                seqprob_val = _Ca_sequence.Sequence_data(ca, scores)
                res.set_property("SEQPROB", _Property_sequence_data(seqprob_val))

    def prepare_scores(
        self,
        ichain_start: int,
        ichain_end: int,
        grid: _Grid_sampling,
        #shiftback: bool = False,
        correlation_mode: bool = False,
        single_atom: bool = False,
        #mol: _MiniMol = None,
        #done: list[bool] = [],
    ):
        """Prepare and set scores for specified chains in a structure

        Args:
            ichain_start (int): Starting chain index
            ichain_end (int): Ending chain index
            grid (_Grid_sampling): Grid sampling
            shiftback (bool, optional): Flag to apply translation. Defaults to False.
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
        """
        if not single_atom:
            for chn in range(ichain_start, ichain_end):
                for r in range(0, self.__mol[chn].size()):
                    self.prepare_score(self.__mol[chn][r], grid, correlation_mode)  # noqa: E501
                self.__done[chn] = True
        else:
            for chn in range(ichain_start, ichain_end):
                for r in range(0, self.__mol[chn].size()):
                    self.prepare_score_single_atom(self.__mol[chn][r], grid, correlation_mode)  # noqa: E501
                self.__done[chn] = True

    def map_xmap_to_nxmap_grid(self, vol_points: _np.ndarray, vol_indices: _np.ndarray, rt: _RTop_orth):
        if len(vol_points) != len(vol_indices):
            raise ValueError(f"Error map_xmap_to_nxmap_grid: lengths of vol_points and vol_indices are different!")
        
        if (self.map_params.grid_asu == self.map_params.grid).all():
            nxm_cell = self.map_params.cell
        else:
            celldim = _np.array([self.map_params.spacing[0]*self.map_params.grid_asu[0],
                                 self.map_params.spacing[1]*self.map_params.grid_asu[1],
                                 self.map_params.spacing[2]*self.map_params.grid_asu[2]])
            print(celldim)
            nxm_cell = _Cell([*celldim, 90., 90., 90.])

        cf_points = (vol_points + self.map_params.origin)/self.map_params.grid_asu
        # row vector-matrix product
        co_points = cf_points@nxm_cell.matrix_orth.array
        co_points = (co_points@rt.rot.array) + rt.trn.array
        #cg_points = _np.rint((co_points@nxm_cell.matrix_frac.array) * self.xwrk_int.grid_sampling.array.astype(_np.float64)).astype(int)
        cg_points = _np.rint((_np.matvec(self.xwrk_int.cell.matrix_frac.array, co_points)) * self.xwrk_int.grid_sampling.array.astype(_np.float64)).astype(int)
        for cg, ind in zip(cg_points, vol_indices):
            ix = self.xwrk_int.map_reference_coord(_Coord_grid(*cg))
            self.xwrk_int[ix] = ind

    def set_map_parameters(
        self,
        xwrk: _Xmap_base, #mapin_path: str,
        #fix_axis_positions=False,
        #fix_origin=False,
        initial_model: str = "NONE",
        shiftback=False,
    ):
        """Set map parameters

        Args:
            xwrk (_Xmap_base): Path to ML input map file
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
            fix_origin (bool, optional): Flag to set origin. Defaults to False.
            shiftback (bool, optional): Flag to set translation parameters. Defaults to False.
        """  # noqa: E501
        self.xwrk_int.init(xwrk.spacegroup, xwrk.cell, xwrk.grid_sampling)
        self.map_params.set_param_values_3(self.map_params.workgrid, xwrk.grid_sampling.nu, xwrk.grid_sampling.nv, xwrk.grid_sampling.nw)
        gmap = _gemmi.read_ccp4_header(f"{self.datapath}/{self.filenames.cutout_input_map}")  # noqa: E501
        #if fix_axis_positions: # and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):    # noqa: E501
        # get axis positions, header 1,2,3 and 5,6,7 depends on axis positions
        axispos = _np.array(gmap.axis_positions())
        # fix to 0, 1, 2
        fixed_pos = [_np.where(axispos == 0)[0][0], _np.where(axispos == 1)[0][0], _np.where(axispos == 2)[0][0]]  # noqa: E501
        origin = [gmap.header_i32(5), gmap.header_i32(6), gmap.header_i32(7)]
        self._debug_list(f"origin before", [origin], self.verbose) #{gmap.header_i32(5)}, {gmap.header_i32(6)}, {gmap.header_i32(7)}", self.verbose)   # noqa: E501
        #gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        self.map_params.set_origin(origin[fixed_pos[0]], origin[fixed_pos[1]], origin[fixed_pos[2]])
        if (self.map_params.origin == 0).all():
            self.ncorrect = 0
        else:
            self.map_params.fix_origin = True
        self._debug_list(f"origin after", [self.map_params.origin], self.verbose)# {gmap.header_i32(5)}, {gmap.header_i32(6)}, {gmap.header_i32(7)}", self.verbose)   # noqa: E501
        self.map_params.fix_axis_positions = fixed_pos
        self.map_params.spacing = _np.array(
            [
                gmap.header_float(11) / gmap.header_i32(8),
                gmap.header_float(12) / gmap.header_i32(9),
                gmap.header_float(13) / gmap.header_i32(10),
            ]
        )
        # straight from map header, not repositioned
        grid_asu = _np.array([gmap.header_i32(1), gmap.header_i32(2), gmap.header_i32(3)])  # noqa: E501
        self.map_params.set_grid_asu(grid_asu[fixed_pos[0]], grid_asu[fixed_pos[1]], grid_asu[fixed_pos[2]])
        #grid = _np.array([gmap.header_i32(8), gmap.header_i32(9), gmap.header_i32(10)])  # noqa: E501
        self.map_params.set_grid(gmap.header_i32(8), gmap.header_i32(9), gmap.header_i32(10))  # noqa: E501
        self.map_params.cell = _np.array([gmap.header_float(11), gmap.header_float(12), gmap.header_float(13), gmap.header_float(14), gmap.header_float(15), gmap.header_float(16)])
        # self.corrections.spacing = _np.asarray(spacing)
        #if fix_origin:
        #    self.map_params.origin = _np.array([gmap.header_i32(5), gmap.header_i32(6), gmap.header_i32(7)])  # noqa: E501
        #    self._debug_list("origin", self.map_params.origin, self.verbose)

        if shiftback or initial_model != "NONE":
            init_mol = read_structure(initial_model)  # noqa: F405
            cut_mol = read_structure(f"{self.datapath}/{self.filenames.cutout_model}")  # noqa: F405, E501
            # rtop transforming cut model back to initial model
            rtop = self.get_rt_from_models(cut_mol, init_mol, True)
            if rtop.is_null():
                trn = self.get_translation_from_cutout(set_map_params=True)
                #self.map_params.set_translation(trn[0], trn[1], trn[2])
                rtop = _RTop_orth(_Mat33_double(self.map_params.rotation),
                                  _Vec3_double(self.map_params.translation))
                #self.map_params.shiftback = _np.asarray(tr, dtype=float)

        density = _np.load(f"{self.datapath}/{self.filenames.density_array}")
        if any(p1 != p2 for p1, p2 in zip(fixed_pos, [0, 1, 2])):
            density = _np.swapaxes(density, axispos[0], fixed_pos[0])  # fmt: skip  # noqa: E501
            density = _np.swapaxes(density, axispos[1], fixed_pos[1])  # fmt: skip  # noqa: E501
        xyz = _np.mgrid[
                    0 : density.shape[0],  # fmt: skip  # noqa: E203
                    0 : density.shape[1],  # fmt: skip  # noqa: E203
                    0 : density.shape[2],  # fmt: skip  # noqa: E203
              ].astype(_np.int64)

        vol_points = (xyz[:, density > 0].T).astype(_np.int64)
        vol_indices = _np.ravel_multi_index(xyz, density.shape)[density > 0]
        self.map_xmap_to_nxmap_grid(vol_points, vol_indices, rtop)  # noqa: E501
        #print("check xwrk_int", self.xwrk_int.grid_asu, self.xwrk_int.grid_sampling)
        # need to continue from here fri 19 2025

    def set_map_parameters_from_array(
        self,
        cell: _Union[_Sequence[float], _np.ndarray],
        fix_axis_positions: bool = False,
        fix_origin: bool = False,
    ):
        """Set map parameters from density array

        Args:
            cell (Union[Sequence[float], _np.ndarray]): cell parameters
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
            fix_origin (bool, optional): Flag to set origin. Defaults to False.
        """  # noqa: E501
        density = _np.load(f"{self.datapath}/density.npy")
        if fix_axis_positions:  # will swap the 1st and 3rd axes of the array
            density = _np.swapaxes(density, 0, 2)
        self.map_params.spacing = _np.array(
            [
                cell[0] / density.shape[0],
                cell[1] / density.shape[1],
                cell[2] / density.shape[2],
            ]
        )
        self.map_params.cell = _np.asarray(cell)
        if fix_origin:
            nxs = -density.shape[0] // 2
            nys = -density.shape[1] // 2
            nzs = -density.shape[2] // 2
            self.map_params.origin = _np.array(nxs, nys, nzs)
            self.map_params.ncorrect = 1

    def get_rt_from_models(self, source_mol: _MiniMol, target_mol: _MiniMol, set_map_params: bool = True) -> _RTop_orth:
        """Get rotational translation operator from cutout and original models

        Args:
            source_mol (_MiniMol): Model to be transformed
            target_mol (_MiniMol): Fixed model
            set_map_params (bool, optional): Flag to set rotational and translation arrays in map_parameters field. Defaults to True.

        Returns:
            _RTop_orth: Rotation-translation operator
        """
        # src = atoms to be transformed, tgt = fixed model
        rt = _RTop_orth(source_mol.atom_list(), target_mol.atom_list())
        print(f"rtop :\n{rt.format()}")
        if set_map_params:
            self.map_params.rotation = rt.rot.array
            self.map_params.translation = rt.trn.array
        return rt

    def get_translation_from_cutout(self, set_map_params: bool = True):
        """Get translation coordinates from model cut out used for machine learning

        Returns:
            numpy.ndarray : Array of translation coordinates
        """  # noqa: E501
        xyz = _np.array([0.0, 0.0, 0.0])
        shiftback_set = False
        with open(_path.join(self.datapath, "cutout.pdb"), "r") as fopen:
            for line in fopen:
                if "TRANSLATED BY" in line:
                    pattern = r"TRANSLATED BY \(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)"
                    match = _re.search(pattern, line)
                    if match:
                        xyz[0], xyz[1], xyz[2] = map(float, match.groups())
                        if set_map_params:
                            self.map_params.rotation = _np.eye(3, 3)
                            self.map_params.set_translation(-xyz[0], -xyz[1], -xyz[2])
                        shiftback_set = True
        
        if not shiftback_set:
            print("Error: HelperMTStackNetOsaka.get_translation_from_cutout- Unable to find translation coordinates used in PDB header.\n")
            print("       shiftback not set.")

        return xyz  # _Coord_orth(xyz[0],xyz[1],xyz[2])

    ClusterMode = _Literal["dbscan", "kmeans"]  # noqa: F405

    def get_map_coords_from_predicted_instance(
        self,
        mode: ClusterMode = "kmeans",
        fix_origin: bool = True,
        mapin_path: str = "NONE",
        write_npy: bool = False,
        write_allpoints_npy: bool = False,
        seqlen: int = 0,
        log = None,
        # shiftback=True,
    ):
        """Return coordinates of amino acid instances from machine learning numpy output.

        Args:
            mode (ClusterMode, optional): Clustering mode to use {dbscan, kmeans}. Defaults to dbscan
            fix_origin (bool, optional): Flag to fix non-zero origin. Defaults to True.
            mapin_path (str, optional): Path to map file. Defaults to None.
            write_npy (bool, optional): Flag to write amino acid coordinates as .npy file. Defaults to False.
            seqlen (int, optional): Sequence length, number of cluster centers for kmeans clustering, compulsory is kmeans method is used. Defaults to 0.
            verbose (int, optional): Verbosity. Defaults to 0
        Returns:
            List: list containing the orthogonal coordinates amino acid instances
        """  # noqa: E501
        offsets = _np.load(f"{self.datapath}/{self.filenames.predicted_instance}")# inst_pred.npy")
        density = _np.load(f"{self.datapath}/{self.filenames.density_array}") #density.npy")
        mp = self.map_params
        #if mapin_path != "NONE":
        #    # assuming the same map is used in the NN to predict
        #    # amino acid segmentations and instances output
        #    gmap = _gemmi.read_ccp4_map(mapin_path)
        #    if mp.fix_axis_positions:
        #        axis_pos = _np.array(gmap.axis_positions())
        #        offsets = _np.swapaxes(offsets, axis_pos[0], _np.where(axis_pos == 0)[0][0])  # fmt: skip  # noqa: E501
        #        offsets = _np.swapaxes(offsets, axis_pos[1], _np.where(axis_pos == 1)[0][0])  # fmt: skip  # noqa: E501
        #        density = _np.swapaxes(density, axis_pos[0], _np.where(axis_pos == 0)[0][0])  # fmt: skip  # noqa: E501
        #        density = _np.swapaxes(density, axis_pos[1], _np.where(axis_pos == 1)[0][0])  # fmt: skip  # noqa: E501
        #        gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        #else:
        #    # swap the x and z axes from the ML npy output
        #    if mp.fix_axis_positions:
        #        # offset array is of shape (3,nu,nv,nw)
        #        offsets = _np.swapaxes(offsets, 1, 3)
        #        density = _np.swapaxes(density, 0, 2)
        #if fix_origin:
        #    xyz = _np.mgrid[
        #        mp.origin[0] : density.shape[0] + (mp.origin[0] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
        #        mp.origin[1] : density.shape[1] + (mp.origin[1] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
        #        mp.origin[2] : density.shape[2] + (mp.origin[2] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
        #    ].astype(_np.float64)
        #else:
        xyz = _np.mgrid[
            mp.origin[0] : density.shape[0] + (mp.origin[0] + mp.ncorrect),  # fmt: skip  # noqa: E203
            mp.origin[1] : density.shape[1] + (mp.origin[1] + mp.ncorrect),  # fmt: skip  # noqa: E203
            mp.origin[2] : density.shape[2] + (mp.origin[2] + mp.ncorrect),  # fmt: skip  # noqa: E203
        ].astype(_np.float64)

        self._debug(f"xyz shape : {xyz.shape}", self.verbose)
        self._debug(f"offsets shape : {offsets.shape}", self.verbose)
        self._debug(f"density shape : {density.shape}", self.verbose)
        self._debug_list("spacing", mp.spacing, self.verbose)

        #for i in (0, 1, 2):
        #    offsets[i] = offsets[i] * mp.spacing[i]
        #    xyz[i] = xyz[i] * mp.spacing[i]
        points = (xyz + offsets) * (density > 0)
        points = _np.round(points[:, density > 0].T).astype(int)
        log.log("ROUND")
        if write_allpoints_npy:
            _np.savez_compressed("aa_all_points.npz", pos=points, spacing=mp.spacing)
            #_np.save("aa_all_points.npy", points, allow_pickle=False)
        # if shiftback:
        #    points -= self.corrections.shiftback
        # tr = self.get_translation_from_cutout()
        # points -= tr
        log.log("NPZ")
        aa_instance_coordinates = []
        if mode == "kmeans":
            if seqlen == 0:
                raise ValueError("Please provide total number of sequence 'seqlen=' for KMeans clustering.")  # noqa: E501
            aa_instance_coordinates = self._kmeans_clustering(points, seqlen, write_npy, log)  # noqa: E501
        elif mode == "dbscan": # dbscan
            aa_instance_coordinates = self._dbscan_clustering(points, 10, write_npy)  # noqa: E501
        else:
            tmp = _np.zeros((len(points), 3))
            uniq = _np.unique(
            points, return_inverse=False, return_counts=False, axis=0
            )
            # translate back from nxmap space to xmap space
            uniq = uniq.astype(float) * mp.spacing
            tmp = _np.matvec(mp.rotation, uniq) + mp.translation
            for i in range(0, len(tmp)):
                aa_instance_coordinates.append(_Coord_orth(tmp[i]))
            #for i in range(0, len(uniq)):
            #    cg = _np.array([uniq[i][0], uniq[i][1], uniq[i][2]], dtype=float)
            #    #tmp[i] = co - self.map_params.shiftback
            #    tmp[i] = (cg * self.map_params.spacing)@mp.rotation + mp.translation
            #    aa_instance_coordinates.append(_Coord_orth(tmp[i]))
        
            if write_npy:
                _np.save("aa_inst_kmeans_mean.npy", tmp, allow_pickle=False)
        log.log("CLUS")
        return aa_instance_coordinates

    def get_map_coord_grid_from_predicted_instance(
            self,
            mode: ClusterMode = "dbscan",
            fix_origin: bool = True,
            mapin_path: str = "NONE",
            write_npy: bool = False,
            seqlen: int = 0,
    ):
        """Return coordinate grid of amino acid instances from machine learning numpy output.

        Args:
            mode (ClusterMode, optional): Clustering mode to use {dbscan, kmeans}. Defaults to dbscan
            fix_origin (bool, optional): Flag to fix non-zero origin. Defaults to True.
            mapin_path (str, optional): Path to map file. Defaults to None.
            write_npy (bool, optional): Flag to write amino acid coordinates as .npy file. Defaults to False.
            seqlen (int, optional): Sequence length, number of cluster centers for kmeans clustering, compulsory is kmeans method is used. Defaults to 0.
            verbose (int, optional): Verbosity. Defaults to 0
        Returns:
            List: list containing the coordinate grid for amino acid instances
        """  # noqa: E501
        offsets = _np.load(f"{self.datapath}/inst_pred.npy")
        density = _np.load(f"{self.datapath}/density.npy")
        mp = self.map_params
        if mapin_path != "NONE":
            # assuming the same map is used in the NN to predict
            # amino acid segmentations and instances output
            gmap = _gemmi.read_ccp4_map(mapin_path)
            if mp.fix_axis_positions:
                axis_pos = _np.array(gmap.axis_positions())
                offsets = _np.swapaxes(offsets, axis_pos[0], _np.where(axis_pos == 0)[0][0])  # fmt: skip  # noqa: E501
                offsets = _np.swapaxes(offsets, axis_pos[1], _np.where(axis_pos == 1)[0][0])  # fmt: skip  # noqa: E501
                density = _np.swapaxes(density, axis_pos[0], _np.where(axis_pos == 0)[0][0])  # fmt: skip  # noqa: E501
                density = _np.swapaxes(density, axis_pos[1], _np.where(axis_pos == 1)[0][0])  # fmt: skip  # noqa: E501
                gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        else:
            # swap the x and z axes from the ML npy output
            if mp.fix_axis_positions:
                # offset array is of shape (3,nu,nv,nw)
                offsets = _np.swapaxes(offsets, 1, 3)
                density = _np.swapaxes(density, 0, 2)
        if fix_origin:
            xyz = _np.mgrid[
                mp.origin[0] : density.shape[0] + (mp.origin[0] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
                mp.origin[1] : density.shape[1] + (mp.origin[1] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
                mp.origin[2] : density.shape[2] + (mp.origin[2] + mp.ncorrect),  # fmt: skip  # noqa: E203, E501
            ].astype(_np.float64)
        else:
            xyz = _np.mgrid[
                mp.origin[0] : density.shape[0],  # fmt: skip  # noqa: E203
                mp.origin[1] : density.shape[1],  # fmt: skip  # noqa: E203
                mp.origin[2] : density.shape[2],  # fmt: skip  # noqa: E203
            ].astype(_np.float64)
        
        self._debug(f"xyz shape : {xyz.shape}", self.verbose)
        self._debug(f"offsets shape : {offsets.shape}", self.verbose)
        self._debug(f"density shape : {density.shape}", self.verbose)
        self._debug_list("spacing", mp.spacing, self.verbose)

        points = (xyz + offsets) * (density > 0)
        points = points[:, density > 0].T
        aa_instance_coordinates = []
        if mode == "kmeans":
            if seqlen == 0:
                raise ValueError("Please provide total number of sequence 'seqlen=' for KMeans clustering.")  # noqa: E501
            aa_instance_coordinates = self._kmeans_clustering(points, seqlen, write_npy)  # noqa: E501
        else:  # dbscan
            aa_instance_coordinates = self._dbscan_clustering(points, 10, write_npy)  # noqa: E501

        return aa_instance_coordinates

    def _dbscan_clustering(self,
                           points: _np.ndarray,
                           min_samples: int = 10,
                           write_npy: bool = False,
                           ):
        """Find cluster centers for given points using HDBSCAN method

        Args:
            points (_np.ndarray): Array containing coordinates shaped (N,3)
            min_samples (int, optional): Minimum cluster size . Defaults to 10.
            write_npy (bool, optional): Flag to write output as .npy file. Defaults to False.

        Returns:
            List: A list of coordinates of cluster centers
        """
        # 23 Dec 2025 fix this
        points_int = _np.round(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        mp = self.map_params
        aa_instance_coordinates = []
        counts = count[indices]
        rad = 2.9 * (_np.sqrt(3)/_np.sqrt(_np.sum(mp.spacing**2)))
        db = _HDBSCAN(cluster_selection_epsilon=rad, #_np.sqrt(_np.sum(mp.spacing * mp.spacing)),
                      min_cluster_size=min_samples,
                      algorithm="ball_tree",
                      store_centers="centroid").fit(points)
        cluster = db.centroids_
        aa_instance_coordinates = []
        tmp = _np.zeros((len(cluster), 3))
        rt = _RTop_orth()
        for i in range(0, len(cluster)):
            co = _np.array([cluster[i][0], cluster[i][1], cluster[i][2]], dtype=float)
            #tmp[i] = co - mp.shiftback
            tmp[i] = (co * mp.spacing) - mp.shiftback
            aa_instance_coordinates.append(_Coord_orth(tmp[i]))
        #    eps=_np.sqrt(_np.sum(mp.spacing * mp.spacing)),
        #    min_samples=min_samples,
        #    algorithm="ball_tree",
        #).fit(points)
        
        #labels = db.labels_
        #groups = {}
        #weights = {}
        #for label in _np.unique(labels):
        #    if label != -1:
        #        groups[label] = points[labels == label]
        #        weights[label] = counts[labels == label]
        ## calculate center of the points with weights
        #tmp = _np.zeros((len(groups), 3))
        #ind = 0
        #for label, group in groups.items():
        #    sum_x = 0.0
        #    sum_y = 0.0
        #    sum_z = 0.0
        #    w = 0.0
        #    for p in zip(group, weights[label]):
        #        sum_x += p[0][0] * p[1]
        #        sum_y += p[0][1] * p[1]
        #        sum_z += p[0][2] * p[1]
        #        w += p[1]
        #    tmp[ind] = _np.array([sum_x / w, sum_y / w, sum_z / w]) - mp.shiftback  # noqa: E501
        #    aa_instance_coordinates.append(_Coord_orth(tmp[ind]))
        #    ind += 1
        if write_npy:
            _np.save("aa_inst_dbscan_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates
    
    def _kmeans_clustering(self, 
                           points: _np.ndarray,
                           seqlen: int,
                           write_npy: bool = False,
                           log = None):
        """Find cluster centers for given points using KMeans++ method

        Args:
            points (_np.ndarray): Array containing coordinates grid shaped (N,3)
            seqlen (int): Sequence length, needed for the number of cluster center to expect
            write_npy (bool, optional): Flag to write output as .npy file. Defaults to False.
            
        Returns:
            List: A list of coordinates of cluster centers
        """
        mp = self.map_params
        kmeans = _KMeans(n_clusters=seqlen, init="k-means++", n_init=self.kmeans_ninit, random_state=0).fit(points)  # noqa: E501
        log.log("KMNS")
        centroids = kmeans.cluster_centers_
        labels = kmeans.labels_
        #grid distance= 3.0 Angstrom  / (ratio of grid to real space)
        rad = 2.9 * (_np.sqrt(3)/_np.sqrt(_np.sum(mp.spacing**2)))
        #inv_rot = _np.linalg.inv(mp.rotation)
        centroids, labels = self._merge_centroids_within_radius(points, centroids, labels, rad)
        log.log("MGCT")
        aa_instance_coordinates = []
        tmp = _np.zeros((len(centroids), 3))
        # translate back from nxmap space to xmap space
        for i in range(0, len(centroids)):
            cg = _np.array([centroids[i][0], centroids[i][1], centroids[i][2]], dtype=float)
            #tmp[i] = co - self.map_params.shiftback
            tmp[i] = (cg * self.map_params.spacing)@mp.rotation + mp.translation
            aa_instance_coordinates.append(_Coord_orth(tmp[i]))
        log.log("TRNS")
        if write_npy:
            _np.save("aa_inst_kmeans_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates
    
    def _decide_centroids_to_merge(self, i, centroids, labels, points, radius, tree):
        neighbours, dist = tree.query_radius([centroids[i]], r=radius, return_distance=True)
        neighbours = neighbours[0]
        mask = _np.isin(labels, neighbours)
        merged_points = points[mask]
        # recluster centroids > 2 into just 2
        if len(neighbours) > 2 and len(merged_points) > 15:
            kmeans = _KMeans(n_clusters=2, init="k-means++", n_init=self.kmeans_ninit, random_state=0).fit(merged_points)  # noqa: E501
            return {"i" : i, "neighbours" : neighbours, "type": "kmeans", "centroids": kmeans.cluster_centers_, "labels": kmeans.labels_, "mask": mask}
        # if count of points is <= 15 just get mean value
        elif len(neighbours) > 1:
            mean = _np.mean(merged_points, axis=0)
            return { "i": i, "neighbours": neighbours, "type": "mean", "centroids": mean, "labels": None, "mask": mask}
        else: # keep, no change
            mask = labels == i
            return { "i": i, "neighbours": neighbours, "type": "default", "centroids": centroids[i], "labels": None, "mask": mask}


    def _merge_centroids_within_radius(self,
                                       points: _np.ndarray,
                                       centroids: _np.ndarray,
                                       labels: _np.ndarray,
                                       radius: float = 2.5):
        """Merge centroids that are too close to each other, within given radius

        Args:
            points (_np.ndarray): Array containing coordinates shaped (N,3)
            centroids (_np.ndarray): Array containing centroids
            labels (_np.ndarray): Array containing centroid labels for each points
            radius (float, optional): Radius cutoff for merging centroids. Defaults to 2.5.

        Returns:
            _type_: _description_
        """
        tree = KDTree(centroids)
        visited = set()
        new_centroids = []
        new_labels = _np.zeros_like(labels)
        new_cluster_id = 0
        dbscan_eps = radius * (_np.sqrt(3)/_np.sqrt(_np.sum(self.map_params.spacing**2)))
        # print('dbscan_eps : ', dbscan_eps)

        with _ThreadPoolExecutor(max_workers=self.__ncpu) as executor:
            workers = [
                executor.submit(
                    self._decide_centroids_to_merge, i, centroids, labels, points, radius, tree
                )
                for i in range(0, len(centroids))
            ]
        results = [w.result() for w in workers]
        # process results
        for res in results:
            i = res["i"]
            if i in visited:
                continue
            visited.update(res["neighbours"])
            if res["type"] == "kmeans":
                for c in res["centroids"]:
                    new_centroids.append(c)
                lbls = res["labels"]
                mask = res["mask"]
                for p, lbl in enumerate(lbls):
                    if lbl >= 0:
                        new_labels[mask][p] = new_cluster_id + lbl
                    else:
                        new_labels[mask][p] = lbl

                new_cluster_id += _np.sum(lbls >= 0)
            elif res["type"] == "mean":
                new_centroids.append(res["centroids"])
                new_labels[res["mask"]] = new_cluster_id
                new_cluster_id += 1
            else:
                new_centroids.append(res["centroids"])
                new_labels[res["mask"]] = new_cluster_id
                new_cluster_id += 1
        '''
        for i in range(0, len(centroids)):
            if i in visited:
                continue
            neighbours, dist = tree.query_radius([centroids[i]], r=dbscan_eps, return_distance=True)
            neighbours = neighbours[0]
            dist = dist[0]
            visited.update(neighbours)
            #print(f"neighbours {neighbours}")
            mask = _np.isin(labels, neighbours)
            #print(dist)
            merged_points = points[mask]
            #print(f"labels {labels}")
            if len(neighbours) > 2 and len(merged_points) > 15:
                #print(len(merged_points))
                #print(len(points))
                #print("Shape:", merged_points.shape)
                #print("Dtype:", merged_points.dtype)
                #print("Is object array:", merged_points.dtype == object)
                #print("Contains lists:\n")
                #for i in merged_points:
                #    print(i)
                #_np.save("merged_points.npy", merged_points, allow_pickle=False)
                #sys.stdout.flush()
                # if have a centroid is too close with 2 other centroids, recluster to 2 centroids
                kmeans = _KMeans(n_clusters=2, init="k-means++", n_init=self.kmeans_ninit, random_state=0).fit(merged_points)  # noqa: E501
                #centroids = kmeans.cluster_centers_
                #labels = kmeans.labels_
        
                #db = _HDBSCAN(cluster_selection_epsilon=dbscan_eps, #_np.sqrt(_np.sum(mp.spacing * mp.spacing)),
                #              min_cluster_size=2,
                #              algorithm="kd_tree",
                #              store_centers="centroid",
                #              copy=True).fit(merged_points)
                for c in kmeans.cluster_centers_:
                    new_centroids.append(c)
                for p, lbl in enumerate(kmeans.labels_):
                    if lbl >= 0:
                        new_labels[mask][p] = new_cluster_id + lbl
                    else:
                        new_labels[mask][p] = lbl
                new_cluster_id += len(lbl[lbl >= 0])
                continue
            elif len(neighbours) > 1:
                mean = _np.mean(merged_points, axis=0)
                new_centroids.append(mean)
                new_labels[mask] = new_cluster_id
                new_cluster_id += 1
                continue
            
            mask = labels == i
            new_centroids.append(centroids[i])
            new_labels[mask] = new_cluster_id
            new_cluster_id += 1 '''
        
        return _np.array(new_centroids), new_labels

    @staticmethod
    def group_within_distance(uniq: _List, cutoff: float = 2.0):
        """Returns a list of points grouped within a given distance cutoff

        Args:
            uniq (_List): List of coordinates
            cutoff (float, optional): Cutoff distance. Defaults to 2.0.

        Returns:
            List: A list of coordinates for each group of points
        """
        indices = _np.where(_spatial.distance.cdist(uniq, uniq, 'euclidean') <= cutoff)[0]  # noqa: E501
        indices = _np.unique(indices)
        return indices

    def map_coords_to_ca_atom(
        self,
        aa_instance_coordinates: _List,
        mol: _MiniMol,
        aa_instance_indices: _List = [],
        ncluster: int = 1,
        globularise: bool = False,
    ):
        """Fill a MiniMol object with Ca atoms at position given list of coordinates

        Args:
            aa_instance_coordinates (_List): List of orthogonal coordinates
            mol (_MiniMol): An initialised MiniMol object to be filled
            aa_instance_indices (_List, optional): List of grid indices corresponding to each point. Defaults to [].
            ncluster (int, optional): Number of separated domains expected, KMeans will be run for clustering. Defaults to 1.
        """
        mmodel = _MModel()
        have_index = True if len(aa_instance_indices) > 0 else False
        aa_inst_array = _np.array(aa_instance_coordinates)
        if ncluster > 1:
            km = _KMeans(n_clusters=ncluster, init="k-means++", random_state=0).fit(aa_inst_array)
            labels = km.labels_
            for n in range(0, ncluster):
                d = aa_inst_array[labels == n]
                chn = _MChain()
                count = 1
                for i in range(0, len(d)):
                    res = _MRes()
                    res.set_seqnum(count)
                    res.type = "UNK"
                    atm = _MAtom.null()
                    matom = _MAtom(atm)
                    matom.element = "C"
                    matom.id = "CA"
                    # matom.set_name("CA")
                    matom.b_iso = 40.0
                    matom.occupancy = 1.0
                    # matom.pos = _Coord_orth((co*correction))
                    matom.pos = _Coord_orth(aa_inst_array[i])
                    if have_index:
                        matom.set_property("INDEX", _Property_int(aa_inst_array[i]))  # noqa: E501
                    # print(matom)
                    res.insert(matom, -1)
                    chn.insert(res, -1)
                    count += 1
                mmodel.insert(chn, -1)
        else:
            for i in range(0, len(aa_instance_coordinates)):
                chn = _MChain()
                res = _MRes()
                res.set_seqnum(1)
                res.type = "UNK"
                atm = _MAtom.null()
                matom = _MAtom(atm)
                matom.element = "C"
                matom.id = "CA"
                matom.b_iso = 40.0
                matom.occupancy = 1.0
                matom.pos = _Coord_orth(aa_instance_coordinates[i])
                if have_index:
                    matom.set_property("INDEX", _Property_int(aa_instance_indices[i]))  # noqa: E501
                res.insert(matom, -1)
                chn.insert(res, -1)
                mmodel.insert(chn, -1)
        mol.model = mmodel
        _ProteinTools.chain_label(mol, True)
        if globularise:
            _ProteinTools.globularise(mol, _Coord_frac(0, 0, 0))

    def coord_frac_to_ca_atom(self, aa_instance_frac: _List, mol: _MiniMol):
        """Fill MiniMol object with Ca atoms at position given a list of fractional coordinates

        Args:
            aa_instance_frac (_List): A list of fractional coordinates
            mol (_MiniMol): An initiliased MiniMol object to fill
        """
        mmodel = _MModel()
        for cf in aa_instance_frac:
            chn = _MChain()
            res = _MRes()
            res.set_seqnum(1)
            res.type = "UNK"
            atm = _MAtom.null()
            matom = _MAtom(atm)
            matom.element = "C"
            matom.id = "CA"
            # matom.set_name("CA")
            matom.b_iso = 40.0
            matom.occupancy = 1.0
            matom.pos = cf.coord_orth(mol.cell)  # *correction.spacing))
            # print(matom)
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
            # print("len model ", len(mmodel))
        mol.model = mmodel
        _ProteinTools.chain_label(mol, True)

    def map_grid_to_ca_atom(
        self,
        aa_instance_coordinates: _List,
        mol: _MiniMol,
        cell: _Cell,
        grid: _Grid_sampling,
    ):
        """Fill a MiniMol object with Ca atoms at positions given a list of map grid coordinates

        Args:
            aa_instance_coordinates (_List): A list of map grid coordinates
            mol (_MiniMol): An initialised MiniMol object to be filled
            cell (_Cell): Cell object
            grid (_Grid_sampling): Grid sampling
        """
        mmodel = _MModel()
        for co in aa_instance_coordinates:
            chn = _MChain()
            res = _MRes()
            res.set_seqnum(1)
            res.type = "UNK"
            atm = _MAtom.null()
            matom = _MAtom(atm)
            matom.element = "C"
            matom.id = "CA"
            # matom.set_name("CA")
            matom.b_iso = 40.0
            matom.occupancy = 1.0
            # matom.pos = _Coord_orth((co*correction))
            matom.pos = co.coord_frac(grid).coord_orth(cell)
            # print(matom)
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
            # print("len model ", len(mmodel))
        mol.model = mmodel
        _ProteinTools.chain_label(mol, True)

    def check_is_seqprob_set(self):
        """Check if SEQPROB property is set for all residue

        Returns:
            bool: True if set, False if not
        """
        if self.__mol.size() == 0:
            return False
        else:
            is_set = True
            for chn in self.__mol:
                for res in chn:
                    if not res.exists_property("SEQPROB"):
                        is_set = False
            return is_set

    def sequence_array_to_map(self, gmap):
        """Turn arrays of amino acid probability output from machine learning into corresponding maps

        Args:
            gmap (gemmi.Ccp4Map): An initiliase Gemmi map object with cell and spacegroup
        """
        if (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
            gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        cell = gmap.grid.unit_cell
        spg = gmap.grid.spacegroup
        for i in range(0, self.__sequence_array.shape[0]-1):
            a = self.__sequence_array[i]
            gmap.grid = _gemmi.FloatGrid(_np.array(a, dtype=_np.float32))
            gmap.grid.set_unit_cell(cell)
            gmap.grid.spacegroup = spg
            if i in [5, 9, 10]:
                name = str(CH_TO_AA_DICT[i]).replace("/","")
                gmap.write_ccp4_map(f"channel_{str(name)}_pred.ccp4")
            else:
                gmap.write_ccp4_map(f"channel_{str(CH_TO_AA_DICT[i])}_pred.ccp4")

#@staticmethod
def write_out_model_with_seq(mol: _Minimol = None, outfile: str = "model_with_seqprob.pdb"):
    """Write a model file (PDB format) with residue type deduced from probability values

    Args:
        mol (_Minimol, optional): An initialised MiniMol object. Defaults to None.
        outfile (str, optional): Output file name. Defaults to "NONE".
    """
    for chn in mol:
        for res in chn:
            if res.exists_property("SEQPROB"):
                seqprob = res.get_property("SEQPROB").value
                if all(d == 0.0 for d in seqprob.data):
                    res.type = "UNK"
                else:
                    seqindex = _np.argmin(seqprob.data) # values are -ved
                    res.type = str(_ProteinTools.residue_code_3(seqindex))
    
    write_structure(mol, outfile, False)


def cut_fragment_density(xmap: _Xmap_float, fragment: _Union[_MModel, _MChain, _MRes, _MAtom], radius: float, resolution: float, offset: bool=True, center: bool=True, respad: bool=False, cell_multiplier: float=1.0, return_array: bool=False, verbose: bool=False):
    """Cut out density of given fragment and place in a box

    Args:
        xmap (_Xmap_float): Original density map of model
        fragment (_Union[_MModel, _MChain, _MRes, _MAtom]): Fragment where density is to be cut
        radius (float): Radius around given model to include
        resolution (float): Resolution
        offset (bool, optional): Offset density/fragment to center of mass. Defaults to True.
        center (bool, optional): Offset density/fragment to center of box. Defaults to True.
        respad (bool, optional): Turn on resolution padding. Defaults to False.
        cell_multiplier (float, optional): Cell multiplier. Defaults to 1.0.
        return_array (bool, optional): If true, returns only array contain cut density else Xmap_float instance. Defaults to False.
        verbose: (bool, optional): Verbosity. Defaults to False.

    Returns:
        Xmap_float/Ndarray: If return_array is turned on, a 3D array with cut density values is returned, else Xmap_float instance
        restmp: Shifted or non-shifted fragment
        cellcut: Cut cell
    """  # noqa: E501
    cell = xmap.cell
    grid = xmap.grid_sampling
    # bounds from fragment
    atoms = []
    if not isinstance(fragment, _MAtom):
        atoms = fragment.atom_list()
    else:
        atoms.append(fragment)
    cm0 = _Coord_map(1.0e20,  1.0e20,  1.0e20)
    cm1 = _Coord_map(-1.0e20, -1.0e20, -1.0e20)
    # find boundaries for atom positions
    for i in range(0, len(atoms)):
        cm = xmap.coord_map(atoms[i].pos)
        for j in range(0, 3):
            cm0[j] = min(cm0[j], cm[j])
            cm1[j] = max(cm1[j], cm[j])
    if verbose:
        print(f"Atom positions lower and upper boundaries :\n {cm0}\n {cm1}")
    # grid range
    gr0 = _Grid_range(cell, grid, radius+1.0)
    gr1 = _Grid_range(cm0.floor()+gr0.min, cm1.ceil()+gr0.max)
    if verbose:
        print(f"Grid range :\n {gr0.min}\n {gr0.max}")
        print(f"Grid range :\n {gr1.min}\n {gr1.max}")
        
    mask = _NXmap_float()
    mask.init(cell, grid, gr1)
    edcalc = _EDcalc_mask_float(radius)
    edcalc(mask, atoms)

    co0 = _Coord_orth(1.0e20, 1.e20, 1.e20)
    co1 = _Coord_orth(-1.e20, -1e20, -1e20)
    com = _Coord_orth(0., 0., 0.)
    count = 0.0
    ix = mask.first()
    while not ix.last():
        if mask[ix] > 0.0:
            co = ix.coord_orth()
            for j in range(0, 3):
                co0[j] = min(co0[j], co[j])
                co1[j] = max(co1[j], co[j])
            com += co
            count += 1.0
        ix.next()
    com = (1.0/count) * com
    print(f"Centre of mass : {com.format()}")
    print(f"Masked region bounds :\n{co0.format()}\n{co1.format()}")
    coff = _Coord_orth(0., 0., 0.,)
    restmp = fragment.clone()    
    #co1 = gr1.max.coord_frac(grid).coord_orth(cell)
    if (offset):
        coff = com
        if verbose:
            print(f"CoM Offset : {coff}")
        restmp.transform(_RTop_orth(_Mat33_double.identity(), _Vec3_double(-coff.x,-coff.y,-coff.z)))
    #else:
    #    coff1 = co1
    #    print(f"coff {coff1}")
    #    restmp.transform(clipper.RTop_orth(clipper.Mat33_double.identity(), clipper.Vec3_double(-coff1.x, -coff1.y, -coff1.z)))

    cd = cell_multiplier*(co1-co0)
    cellcut = _Cell([cd[0], cd[1], cd[2], 90., 90., 90.])
    spgrcut = _Spacegroup.p1()
    if verbose:
        print(f"Cut cell : {cellcut}")
    ctr = _Coord_orth(0.,0.,0.)
    if center:
        ctr = _Coord_frac(0.5, 0.5, 0.5).coord_orth(cellcut)
        print(f"Box center offset : {ctr}")
        restmp.transform(_RTop_orth(_Mat33_double.identity(), _Vec3_double(ctr.x,ctr.y,ctr.z)))
    
    # get resol from map
    resol = min( resolution, 2.0/(cell.a_star*float(grid.nu)) )
    resol = min( resolution, 2.0/(cell.b_star*float(grid.nv)) )
    resol = min( resolution, 2.0/(cell.c_star*float(grid.nw)) )
    resol *= 1.5
    if respad:
        resol =  1.0 / ( 1.0/resol + _np.sqrt( pow(cellcut.a(),-2.0) +
                                      pow(cellcut.b(),-2.0) +
                                      pow(cellcut.c(),-2.0) ) )
    rescut = _Resolution(resol)
    print(rescut)
    gridcut = _Grid_sampling(spgrcut, cellcut, rescut)
    xcut = _Xmap_float(spgrcut, cellcut, gridcut)
    xcut.fill_map_with(0.0)
    if verbose:
        print(f"Cut grid : {gridcut}")
    # find bounds of masked region
    cg0 = _Coord_grid(999999,999999,999999)
    cg1 = _Coord_grid(-999999,-999999,-999999)
    ix = mask.first()
    while not ix.last():
        if mask[ix] > 0.0:
            cg = xcut.coord_map(ix.coord_orth()-coff+ctr).coord_grid()
            for j in range(0, 3):
                cg0[j] = min(cg0[j], cg[j]-1)
                cg1[j] = max(cg1[j], cg[j]+1)
        ix.next()

    #spacing = [cell.x/float(grid.nu), cell.y/float(grid.nv), cell.z/float(grid.nw)]
    # fill cut map
    #c0 = clipper.Coord_grid(0,0,0)
    i0 = xcut.map_reference_coord(cg0)
    iu = xcut.map_reference_coord(i0.coord)
    while iu.coord.u <= cg1.u:
        iv = xcut.map_reference_coord(iu.coord)
        while iv.coord.v <= cg1.v:
            iw = xcut.map_reference_coord(iv.coord)
            while iw.coord.w <= cg1.w:
                co = iw.coord_orth()+coff-ctr
                cg = mask.coord_map(co).coord_grid()
                if mask.in_map(cg):
                    if mask.get_data(cg) != 0.0:
                        xcut[iw] = xmap.interp_cubic_orth(co)
                iw.next_w()
            iv.next_v()
        iu.next_u()

    #mapfile = _CCP4MAPfile()
    #fname = f"{fragment.type}_{str(fragment.seqnum)}"
    #mapfile.open_write(f"{fname}.ccp4")
    #mapfile.export_xmap_float( xcut )
    #mapfile.close_write()
    return xcut.array, restmp, cellcut


def interpolate_map_grid(xmap: _Xmap_base, spacing: float, box_size: int, overlap: int, interp_mode: int = 1, whole_unitcell: bool=False, mapout: bool=False):
    """Interpolate map onto new map grid

    Args:
        xmap (_Xmap_float): Xmap object
        spacing (float): Grid spacing
        box_size (int): New grid size in one dimension, this will be use for all three dimensions to make a cubic box
        overlap (int): Number of grid points a two dimension box overlaps
        interp_mode (int): Mode of interpolation. 0=nearest, 1=linear, 3=cubic. Defaults to 1
        whole_unitcell (bool, optional): Set True to interpolate whole unit cell, otherwise only cover the asymmetric unit. Defaults to False.
        mapout (bool, optional): Set True to return Map object, otherwise 3D array. Defaults to False

    Returns:
        Xmap_float or Ndarray: If return_array is turned on, a 3D array is returned, else Xmap_float instance.
        RTop: Rotation-translation operator
    """  # noqa: E501
    cf0 = _Coord_frac(0, 0, 0)
    cf1 = _Coord_frac(1, 1, 1)
    cell = xmap.cell
    grid = xmap.grid_sampling
    if whole_unitcell:
        gr = _Grid_range(grid, cf0, cf1)
    else:
        gr = _Grid_range(xmap.grid_asu.min, xmap.grid_asu.max)

    border = box_size // 2
    co_min = gr.min.coord_frac(grid).coord_orth(cell)
    co_max = gr.max.coord_frac(grid).coord_orth(cell)
    print("before")
    print(co_min.format(), co_max.format())
    gr.add_border(border)
    co_min = gr.min.coord_frac(grid).coord_orth(cell)
    co_max = gr.max.coord_frac(grid).coord_orth(cell)
    print("after add border")
    print(co_min.format(), co_max.format())
    cd = co_max - co_min
    newcell = _Cell([cd[0], cd[1], cd[2], 90., 90., 90.])
    numx = -(-int(newcell.a / spacing) // overlap * overlap)
    numy = -(-int(newcell.b / spacing) // overlap * overlap)
    numz = -(-int(newcell.c / spacing) // overlap * overlap)
    newgrid = _Grid_sampling(numx, numy, numz)
    # RTop to relate grids
    rtop = _RTop_orth(_Mat33_double.identity(), co_min)
    newcell = _Cell([(spacing*numx), (spacing*numy), (spacing*numz), 90., 90., 90.])
    grnew = _Grid_range(_Coord_grid(0, 0, 0), _Coord_grid(numx-1, numy-1, numz-1))
    nxm = _NXmap_float(newcell, newgrid, grnew)
    nxop = _NX_operator(xmap, nxm, rtop)
    nxm.fill_map_with(0.0)
    cg0 = _Coord_grid(0, 0, 0)
    cg1 = _Coord_grid(numx-1, numy-1, numz-1)
    nxop.xmap_interp_all_points(xmap, nxm, cg0, cg1, interp_mode)
    
    if mapout:
        return nxm, rtop
    else:
        return nxm.array, rtop


def reinterpolate_map_grid(work_array: _np.ndarray, rtop: _RTop_orth, spacing: float, input_xmap: _Xmap_base, interp_mode: int = 1, whole_unitcell: bool=False, mapout: bool=False):
    """Reinterpolate results in array to map grid

    Args:
        work_array (_np.ndarray): Results in array
        rtop (_RTop_orth): Rotation-translation operator
        spacing (float): Grid spacing
        input_xmap (_Xmap_base): Initial Xmap object
        interp_mode (int): Mode of interpolation. 0=nearest, 1=linear, 3=cubic. Defaults to 1
        whole_unitcell (bool, optional): Set True to reinterpolate whole unit cell, otherwise only cover the asymmetric unit. Defaults to False.
        mapout (bool, optional): Set True to return Map object, otherwise 3D array. Defaults to False

    Returns:
        Xmap_float or ndarray: If return_array is turned on, a 3D array is returned, else Xmap_float instance.
    """  # noqa: E501
    cf0 = _Coord_frac(0, 0, 0)
    cf1 = _Coord_frac(1, 1, 1)
    input_spg = input_xmap.spacegroup
    input_cell = input_xmap.cell
    input_grid = input_xmap.grid_sampling
    if whole_unitcell:
        gr = _Grid_range(input_grid, cf0, cf1)
    else:
        gr = input_xmap.grid_asu
    
    cellx = work_array.shape[0] * spacing
    celly = work_array.shape[1] * spacing
    cellz = work_array.shape[2] * spacing
    cell = _Cell([cellx, celly, cellz, 90., 90., 90.])
    nxmap = _NXmap_float(work_array, cell)
    xmap = _Xmap_float(input_spg, input_cell, input_grid)
    
    nxop = _NX_operator(xmap, nxmap, rtop)
    cg0 = gr.min #.coord_grid(input_grid)
    cg1 = gr.max #_Coord_frac(1, 1, 1).coord_grid(input_grid)
    nxop.nxmap_interp_all_points(nxmap, xmap, cg0, cg1, 3)
    if mapout:
        return xmap
    else:
        return xmap.array

def get_principal_axes(points: _np.ndarray):
    centered_points = points - _np.mean(points, axis=0)
    _, _, eigenv = _np.linalg.svd(centered_points)
    return eigenv

def reorientate_amino_acid_to_density(aa_atoms, com, xmap):
    # get grid points around com

    axes_aa = get_principal_axes(aa_atoms)
    axes_rho = get_principal_axes(density)
    # 4 possible flips
    flips = [ _np.diag([1,1,1]), # original
             _np.diag([-1,-1,1]), # flip x and y
             _np.diag([-1, 1, -1]), # flip x and z
             _np.diag([1, -1, -1])] # flip y and z
    
    best_score = -_np.inf
    best_coords = None
    com = _np.mean(aa_atoms, axis=0)
    aa_centered = aa_atoms - com
    for flip in flips:
        rot = (axes_rho.T@flip@axes_aa)
        if _np.linalg.det(rot) < 0:
            continue
        #rotate atoms
        rot_atoms = (aa_centered@rot.T) + com
        current_score = score_orientation(xmap, aa_atoms)
        if current_score > best_coords:
            best_score = current_score
            best_coords = trials_coords
    
    return best_scores

    ## calculate rotation matrix
    #rot = axes_rho.T@axes_aa
    ## ensure right handed coordinate system, no reflections
    #axes_aa[-1,:] *= -1
    #rot = axes_rho.T@axes_aa
    ## rotate atoms around their COM
    #com = _np.mean(aa_atoms, axis=0)
    #rotated_atoms = (aa_atoms - com)@rot.T + com
    return rotated_atoms

def score_orientation(res, xmap: _Xmap_float, rtop: _RTop_orth):
    scores = 0.
    pos = _np.matvec(rtop.rot.array, res) + rtop.trn.array
    for p in pos:
        scores += xmap.interp_cubic_orth(_Coord_orth(*p))
    return scores

def grid_points_in_radius(point: _np.ndarray, density: _np.ndarray, radius: float = 6.0):
    r2 = radius * radius
    grmin = _np.floor(point - radius).astype(int)
    grmax = _np.ceil(point + radius).astype(int)
    XYZ = _np.mgrid[
                    grmin[0] : grmax[0],  # fmt: skip  # noqa: E203
                    grmin[1] : grmax[1],  # fmt: skip  # noqa: E203
                    grmin[2] : grmax[2],  # fmt: skip  # noqa: E203
              ].astype(_np.int64)
    mask = ((XYZ[0] - point[0])**2 + (XYZ[1] - point[1])**2 + (XYZ[2] - point[2])**2) <= r2
    points = _np.column_stack((XYZ[0, mask], XYZ[1, mask], XYZ[2, mask]))
    rhomask = density[XYZ[0, mask], XYZ[1, mask], XYZ[2, mask]] > 0.
    return points[rhomask]

def create_temp_residue(n, ca, c, seqnum, type):
    residue = _MRes()
    residue.type = type
    atm = _MAtom.null()
    atom = _MAtom(atm)
    atom.occupancy = 1.0
    atom.b_iso = 40.
    atom.element = "N"
    atom.id = "N"
    atom.pos = n
    residue.insert( atom )

    atom.element = "C"
    atom.id = "CA"
    atom.pos = ca
    residue.insert( atom )

    atom.id = "C"
    atom.pos = c
    residue.insert( atom )

    residue.set_seqnum(seqnum)
    residue.build_sidechain_numbered_rotamer(0)
    co = _Coord_orth(0.,0.,0.)
    no = 0.0
    for a in range(0, residue.size()):
        co += residue[a].pos
        no += 1.0
    co = (1.0/no) * co
    com = residue["CA"].pos
    diff = com - co
    for a in residue:
        a.pos += diff
    
    return residue

def get_principal_axes(points: _np.ndarray, com: _np.ndarray):
    centered_points = points - com
    u, s, vh = _np.linalg.svd(centered_points)
    return vh

def res_to_pos(res: _MRes):
    points = []
    for a in res:
        points.append(a.pos.array)
    return _np.array(points)

def update_residue_pos(res: _MRes, pos: _np.ndarray):
    for a in range(0, res.size()):
        res[a].pos = pos[a]
    
def set_flips():
    rots = []
    rots.append(_RTop_orth(_Mat33_double.identity(), _Vec3_double.zero()))
    rots.append(_RTop_orth(_np.diag([-1,-1,1]), _Vec3_double.zero()))
    rots.append(_RTop_orth(_np.diag([-1,1,-1]), _Vec3_double.zero()))
    rots.append(_RTop_orth(_np.diag([1,-1,-1]), _Vec3_double.zero()))
    return rots


def build_rotated_aa(aa_coms: _np.ndarray, density_file: str, sequence: list, xmap: _Xmap_float, rtop: _RTop_orth): #, threshold: float = 1.0):
    # get points that is within density threshold
    density = _np.load(density_file)
    count = 1
    mol = _MiniMol(xmap.spacegroup, xmap.cell)
    chn = _MChain()
    chn.id = "A"
    for aa_com, seq in zip(aa_coms, sequence):
        #density_points = grid_points_in_radius(aa_com, density)
        nxm_aa_com = aa_com@rtop.rot.inverse + rtop.inverse().trn
        #cagroup = _Ca_group(aa_com * _Ca_group.std_coord_n() , aa_com * _Ca_group.std_coord_ca(), aa_com * _Ca_group.std_coord_c())
        res = create_temp_residue(nxm_aa_com * _Ca_group.std_coord_n() , nxm_aa_com * _Ca_group.std_coord_ca(), nxm_aa_com * _Ca_group.std_coord_c(), count, seq)
        aa_pos = res_to_pos(res)
        axes_aa = get_principal_axes(aa_pos, aa_com)
        axes_rho = get_principal_axes(grid_points_in_radius(aa_com, density), aa_com)
        flips = set_flips()
        best_score = -_np.inf
        best_coords = None
        for flip in flips:
            rot = axes_rho.T@axes_aa
            if _np.linalg.det(rot) < 0:
                continue
            
            rotated_atms = (aa_pos - aa_com)@rot.T + aa_com
            # score
            current_score = score_orientation(rotated_atms, rtop, xmap)
            if current_score > best_score:
                best_score = current_score
                best_coords = rotated_atms

        update_residue_pos(res, best_coords)
        chn.insert(res, -1)

    mol.insert(chn)
    return mol
        


# del _Ca_group, _Cell, _Grid_sampling, _Coord_orth, _Coord_map
# del _MiniMol #_MAtom, _MModel, _MChain, _MRes
del annotations
