from __future__ import annotations
from typing import List as _List
from bobkit.util import *
import scipy.spatial as _spatial
import itertools as _IT

__all__ = [
    "read_structure",
    "write_structure",
    "get_coordinates_from_predicted_instance",
]

from bobkit.buccaneer import (
    Ca_group as _Ca_group,
    Ca_sequence as _Ca_sequence,
    ProteinTools as _ProteinTools,
)

from bobkit.clipper import (
    Cell as _Cell,
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
    Xmap_float as _Xmap_float,
    Vec3_double as _Vec3_double,
    Property_int as _Property_int,
    Property_sequence_data as _Property_sequence_data,
    NXmap_float as _NXmap_float,
    EDcalc_mask_float as _EDcalc_mask_float,
    RTop_orth as _RTop_orth,
    Mat33_double as _Mat33_double,
    Spacegroup as _Spacegroup,
    Coord_frac as _Coord_frac,
    Resolution as _Resolution
)

import numpy as _np
import gemmi as _gemmi
from sklearn.cluster import HDBSCAN as _HDBSCAN
from sklearn.cluster import KMeans as _KMeans
from typing import Sequence as _Sequence, Union as _Union, Literal as _Literal
from os import path as _path
from dataclasses import dataclass as _dataclass
import re as _re
from multiprocessing import Process as _Process

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
    value
    """
    workcell: _np.ndarray = _np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])
    cell: _np.ndarray = _np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])
    grid: _np.ndarray = _np.array([1, 1, 1])
    grid_asu: _np.ndarray = _np.array([1, 1, 1])
    origin: _np.ndarray = _np.array([0, 0, 0])
    spacing: _np.ndarray = _np.array([1.0, 1.0, 1.0])
    ncorrect: int = 0
    fix_origin: bool = False
    fix_axis_positions: bool = False
    shiftback: _np.ndarray = _np.array([0.0, 0.0, 0.0])

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

    def set_shiftback(self, x: float, y: float, z: float):
        self.shiftback[0] = x
        self.shiftback[1] = y
        self.shiftback[2] = z


class HelperMTStackNetOsaka:
    """Helper class with methods and static methods to process the outputs from
    single shot CNN to predict amino acid instance and sequence from map
    segmentation by Osaka group.
    """
    def __init__(
        self, datapath: str = None,
        workcell: _Union[_Sequence[float], _np.ndarray] = None,
        ncpu: int = 1,
    ):
        """Initialise object

        Args:
            datapath (str, optional): Path to ML output. Defaults to None.
            workcell (Union[Sequence[float], _np.ndarray], optional): Work cell parameters. Defaults to None.

        Raises:
            ValueError: Workcell array is expected to have 6 elements
        """  # noqa: E501
        self.datapath = datapath
        self.map_params = MapParameters()
        if workcell is not None:
            tmp = _np.asarray(workcell, dtype=float)
            if tmp.shape != (6,):
                raise ValueError("Expected array with 6 elements")
            self.map_params.workcell = tmp
        if self.datapath is not None:
            self.__sequence_array = _np.load(f"{self.datapath}/pred.npy")
        else:
            self.__sequence_array = _np.full((1, 1), None, dtype=object)
        self.__ncpu = ncpu
        self.__mol = None
        self.__done = []
        
    def _debug(self, msg: str):
        """Prints debug message

        Args:
            msg (str): Debug message
        """
        print(f"DEBUG>> {msg}")

    def _debug_triple(self, msg: str, data: _Union[_Sequence, _np.ndarray]):
        """Prints debug message with data

        Args:
            msg (str): Debug message
            data (Union[Sequence, _np.ndarray]): Data values
        """
        print(f"DEBUG>> {msg} : ", end="")
        for v in data:
            print(f"{v}, ", end="")
        print("\n")

    def __call__(self, mol: _Minimol, shiftback: bool = False, correlation_mode: bool = False):  # , llktargets: LLK_TargetList):  # noqa: E501
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
            processes = []
            self.__done = [False] * self.__mol.size()
            chunk_size = max(1, self.__mol.size() // (self.__ncpu * 4))
        
            for i in range(0, self.__mol.size(), chunk_size):
                start_chn = i
                end_chn = start_chn + chunk_size
                p = _Process(
                    target=self.prepare_scores,
                    args=(start_chn, end_chn, grid_samp, shiftback, correlation_mode),
                )
                processes.append(p)
                p.start()
            # wait
            for p in processes:
                p.join()
        else:
            self.__done = [False] * self.__mol.size()
            self.prepare_scores(0, self.__mol.size(), grid_samp, shiftback, correlation_mode)
            for i in self.__done:
                if not i:
                    print("Did not manage to fully assign sequence probability!")  # noqa: E501

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
    
    def get_probability_values_coord_orth(self, pos: _Coord_orth, correlation_mode: bool = False, fix_origin=True, shiftback=False):  # noqa: E501
        """Get probability value from a selected predicted sequence array.
        For compatibility with Buccaneer's sequencing scores, negative probability value (correlation mode) or
        1.-probability value (non-correlation mode) is returned.

        Args:
            pos (_Coord_orth): Orthogonal coordinates of the atom
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
            fix_origin (bool, optional): Flag to fix origin. Defaults to True.
            shiftback (bool, optional): Flag to shift positions back. Defaults to False.

        Returns:
            probval: Negative probability value in correlation mode, and 1.-probability for non-correlation mode for compatibility with buccaneer sequencing scores
        """  # noqa: E501
        # print(cg)
        grid = _Grid_sampling(self.map_params.grid[0],
                                   self.map_params.grid[1],
                                   self.map_params.grid[2])
        cell = _Cell(self.map_params.cell)
        if shiftback:
            pos = pos + _Coord_orth(self.map_params.shiftback)
        cg = pos.coord_frac(cell).coord_grid(grid)
        if self.map_params.fix_origin:
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


    def prepare_score(
        self,
        res: _MRes,
        grid: _Grid_sampling,
        shiftback: bool = False,
        correlation_mode: bool = False,
    ):
        """Prepare and set scores property for each residue

        Args:
            res (_MRes): MResidue object
            grid (_Grid_sampling): Grid sampling
            shiftback (bool, optional): Flag to shift positions back. Defaults to False.
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
        """
        cached = False
        ca = _Ca_group(res)
        cell = _Cell(self.map_params.cell)
        # print(cell.format())
        pos = _Coord_orth.null()
        if not ca.is_null():
            if res.exists_property("SEQPROB"):
                seqprob_val = res.get_property("SEQPROB").value
                if (
                    (ca.coord_n - seqprob_val.ca.coord_n).lengthsq()
                    and (ca.coord_ca - seqprob_val.ca.coord_ca).lengthsq()
                    and (ca.coord_c - seqprob_val.ca.coord_c).lengthsq()
                ):
                    cached = True
            pos = ca.coord_ca
            poscb = ca.rtop_beta_carbon()
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
                if shiftback:
                    pos = pos + _Coord_orth(self.map_params.shiftback)
                    poscb = poscb + _Coord_orth(self.map_params.shiftback)
                cg = pos.coord_frac(cell).coord_grid(grid)
                cgb = pos.coord_frac(cell).coord_grid(grid)
                if self.map_params.fix_origin:
                    cg = cg - _Coord_grid(self.map_params.origin[0],
                                          self.map_params.origin[1],
                                          self.map_params.origin[2])
                    cgb = cgb - _Coord_grid(self.map_params.origin[0],
                                          self.map_params.origin[1],
                                          self.map_params.origin[2])
                cg = cg.unit(_Grid_sampling(self.map_params.grid[0],
                                            self.map_params.grid[1],
                                            self.map_params.grid[2]))
                cgb = cgb.unit(_Grid_sampling(self.map_params.grid[0],
                                            self.map_params.grid[1],
                                            self.map_params.grid[2]))
                scores = self.get_probability_values(cg, self.__sequence_array, correlation_mode)
                scores += self.get_probability_values(cgb, self.__sequence_array, correlation_mode)
                scores /= 2
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
                # self._debug_triple(res, scores)
                seqprob_val = _Ca_sequence.Sequence_data(ca, scores)
                res.set_property("SEQPROB", _Property_sequence_data(seqprob_val))

    def prepare_scores(
        self,
        ichain_start: int,
        ichain_end: int,
        grid: _Grid_sampling,
        shiftback: bool = False,
        correlation_mode: bool = False,
    ):
        """Prepare and set scores for specified chains in a structure

        Args:
            ichain_start (int): Starting chain index
            ichain_end (int): Ending chain index
            grid (_Grid_sampling): Grid sampling
            shiftback (bool, optional): Flag to shift positions back. Defaults to False.
            correlation_mode (bool, optional): Flag to use correlation mode. Defaults to False.
        """
        for chn in range(ichain_start, ichain_end):
            for r in range(0, self.__mol[chn].size()):
                self.prepare_score(self.__mol[chn][r], grid, shiftback, correlation_mode)  # noqa: E501
            self.__done[chn] = True

    def set_map_parameters(
        self,
        mapin_path: str,
        fix_axis_positions=False,
        fix_origin=False,
        shiftback=False,
    ):
        """Set map parameters

        Args:
            mapin_path (str): Path to map file
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
            fix_origin (bool, optional): Flag to set origin. Defaults to False.
            shiftback (bool, optional): Flag to set shiftback translation parameters. Defaults to False.
        """  # noqa: E501
        gmap = _gemmi.read_ccp4_map(mapin_path)
        if fix_axis_positions and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):  # fmt: skip  # noqa: E501
            self._debug(f"origin before {gmap.header_i32(5)}, {gmap.header_i32(6)}, {gmap.header_i32(7)}")  # fmt: skip  # noqa: E501
            gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
            self._debug(f"origin after {gmap.header_i32(5)}, {gmap.header_i32(6)}, {gmap.header_i32(7)}")  # fmt: skip  # noqa: E501
            self.map_params.fix_axis_positions = fix_axis_positions
        self.map_params.spacing = _np.array(
            [
                gmap.grid.unit_cell.a / gmap.header_i32(8),
                gmap.grid.unit_cell.b / gmap.header_i32(9),
                gmap.grid.unit_cell.c / gmap.header_i32(10),
            ]
        )
        self.map_params.grid_asu = _np.array([gmap.grid.nu, gmap.grid.nv, gmap.grid.nw])  # noqa: E501
        self.map_params.grid = _np.array([gmap.header_i32(8), gmap.header_i32(9), gmap.header_i32(10)])  # noqa: E501
        self.map_params.cell = _np.asarray(gmap.grid.unit_cell.parameters)
        # self.corrections.spacing = _np.asarray(spacing)
        if fix_origin:
            self.map_params.origin = _np.array([gmap.header_i32(5), gmap.header_i32(6), gmap.header_i32(7)])  # noqa: E501
            self._debug_triple("origin", self.map_params.origin)
            self.ncorrect = 0
            self.map_params.fix_origin = fix_origin
        if shiftback:
            tr = self.get_translation_from_cutout()
            self.map_params.shiftback = _np.asarray(tr, dtype=float)

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

    def get_translation_from_cutout(self):
        """Get translation coordinates from model cut out used for machine learning

        Returns:
            numpy.ndarray : Array of translation coordinates
        """  # noqa: E501
        xyz = _np.array([0.0, 0.0, 0.0])
        with open(_path.join(self.datapath, "cutout.pdb"), "r") as fopen:
            for line in fopen:
                if "TRANSLATED BY" in line:
                    pattern = r"TRANSLATED BY \(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)"
                    match = _re.search(pattern, line)
                    if match:
                        xyz[0], xyz[1], xyz[2] = map(float, match.groups())
        return xyz  # _Coord_orth(xyz[0],xyz[1],xyz[2])

    ClusterMode = _Literal["dbscan", "kmeans"]  # noqa: F405

    def get_map_coords_from_predicted_instance(
        self,
        mode: ClusterMode = "dbscan",
        fix_origin: bool = True,
        mapin_path: str = "NONE",
        write_npy: bool = False,
        seqlen: int = 0,
        verbose: int = 0,
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
        if verbose > 5:
            self._debug(f"xyz shape : {xyz.shape}")
            self._debug(f"offsets shape : {offsets.shape}")
            self._debug(f"density shape : {density.shape}")
            self._debug_triple("spacing", mp.spacing)

        for i in (0, 1, 2):
            offsets[i] = offsets[i] * mp.spacing[i]
            xyz[i] = xyz[i] * mp.spacing[i]
        points = (xyz + offsets) * (density > 0)
        points = points[:, density > 0].T
        # if shiftback:
        #    points -= self.corrections.shiftback
        # tr = self.get_translation_from_cutout()
        # points -= tr
        aa_instance_coordinates = []
        if mode == "kmeans":
            if seqlen == 0:
                raise ValueError("Please provide total number of sequence 'seqlen=' for KMeans clustering.")  # noqa: E501
            aa_instance_coordinates = self._kmeans_clustering(points, seqlen, write_npy)  # noqa: E501
        else:  # dbscan
            aa_instance_coordinates = self._dbscan_clustering(points, 10, write_npy)  # noqa: E501

        return aa_instance_coordinates

    def get_map_coord_grid_from_predicted_instance(
            self,
            mode: ClusterMode = "dbscan",
            fix_origin: bool = True,
            mapin_path: str = "NONE",
            write_npy: bool = False,
            seqlen: int = 0,
            verbose: int = 0,
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
        if verbose > 5:
            self._debug(f"xyz shape : {xyz.shape}")
            self._debug(f"offsets shape : {offsets.shape}")
            self._debug(f"density shape : {density.shape}")
            self._debug_triple("spacing", mp.spacing)

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
        points_int = _np.round(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        mp = self.map_params
        aa_instance_coordinates = []
        counts = count[indices]
        db = _HDBSCAN(cluster_selection_epsilon=_np.sqrt(_np.sum(mp.spacing * mp.spacing)),
                      min_cluster_size=min_samples,
                      algorithm="ball_tree",
                      store_centers="centroid").fit(points)
        cluster = db.centroids_
        aa_instance_coordinates = []
        tmp = _np.zeros((len(cluster), 3))
        for i in range(0, len(cluster)):
            co = _np.array([cluster[i][0], cluster[i][1], cluster[i][2]], dtype=float)
            tmp[i] = co - self.map_params.shiftback
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
                           write_npy: bool = False):
        """Find cluster centers for given points using KMeans++ method

        Args:
            points (_np.ndarray): Array containing coordinates shaped (N,3)
            seqlen (int): Sequence length, needed for the number of cluster center to expect
            write_npy (bool, optional): Flag to write output as .npy file. Defaults to False.

        Returns:
            List: A list of coordinates of cluster centers
        """
        kmeans = _KMeans(n_clusters=seqlen, init="k-means++", random_state=0).fit(points)  # noqa: E501
        cluster = kmeans.cluster_centers_

        aa_instance_coordinates = []
        tmp = _np.zeros((len(cluster), 3))
        for i in range(0, len(cluster)):
            co = _np.array([cluster[i][0], cluster[i][1], cluster[i][2]], dtype=float)
            tmp[i] = co - self.map_params.shiftback
            aa_instance_coordinates.append(_Coord_orth(tmp[i]))
        
        if write_npy:
            _np.save("aa_inst_kmeans_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates

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
    ):
        """Fill a MiniMol object with Ca atoms at position given list of coordinates

        Args:
            aa_instance_coordinates (_List): List of orthogonal coordinates
            mol (_MiniMol): An initialised MiniMol object to be filled
            aa_instance_indices (_List, optional): List of grid indices corresponding to each point. Defaults to [].
        """
        mmodel = _MModel()
        have_index = True if len(aa_instance_indices) > 0 else False
        for i in range(0, len(aa_instance_coordinates)):
            chn = _MChain()
            res = _MRes()
            res.set_seqnum(1)
            res.type = "UNK"
            atm = _MAtom.null()
            matom = _MAtom(atm)
            matom.element = "C"
            matom.set_id("CA")
            matom.b_iso = 40.0
            matom.occupancy = 1.0
            matom.pos = _Coord_orth(aa_instance_coordinates[i])
            if have_index:
                matom.set_property("INDEX", _Property_int(aa_instance_indices[i]))  # noqa: E501
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
        mol.set_model(mmodel)
        _ProteinTools.chain_label(mol, True)

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
            matom.set_id("CA")
            # matom.set_name("CA")
            matom.b_iso = 40.0
            matom.occupancy = 1.0
            matom.pos = cf.coord_orth(mol.cell)  # *correction.spacing))
            # print(matom)
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
            # print("len model ", len(mmodel))
        mol.set_model(mmodel)
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
            matom.set_id("CA")
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
        mol.set_model(mmodel)
        _ProteinTools.chain_label(mol, True)

    def check_is_seqprob_set(self):
        """Check if SEQPROB property is set for all residue

        Returns:
            bool: True if set, False if not
        """
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
                seqindex = _np.argmin(seqprob.data) # values are -ved
                res.type = str(_ProteinTools.residue_code_3(seqindex))
    
    write_structure(mol, outfile, False)

   
# del _Ca_group, _Cell, _Grid_sampling, _Coord_orth, _Coord_map
# del _MiniMol #_MAtom, _MModel, _MChain, _MRes
del annotations
