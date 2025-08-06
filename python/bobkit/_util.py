from __future__ import annotations
from typing import List as _List
from bobkit.util import *
import scipy.spatial as _spatial
import itertools as _IT

__all__ = [
    "test_array",
    "read_structure",
    "write_structure",
    "get_coordinates_from_predicted_instance",
]

from bobkit.buccaneer import (
    Ca_group as _Ca_group,
    ProteinTools as _ProteinTools,
)

from bobkit.clipper import (
    Cell as _Cell,
    Grid_sampling as _Grid_sampling,
    Coord_orth as _Coord_orth,
    Coord_map as _Coord_map,
    MiniMol as _MiniMol,
    MModel as _MModel,
    MChain as _MChain,
    MResidue as _MRes,
    MAtom as _MAtom,
    Xmap_float as _Xmap_float,
    Vec3_double as _Vec3_double,
    Property_int as _Property_int,
)

import numpy as _np
import gemmi as _gemmi
from sklearn.cluster import DBSCAN as _DBSCAN
from sklearn.cluster import KMeans as _KMeans
from typing import Sequence, Union, Literal
from os import path as _path
from dataclasses import dataclass as _dataclass
import re

_ProteinTools()

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
        workcell: Union[Sequence[float], _np.ndarray] = None,
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

    def _debug(self, msg: str):
        """Prints debug message

        Args:
            msg (str): Debug message
        """
        print(f"DEBUG>> {msg}")

    def _debug_triple(self, msg: str, data: Union[Sequence, _np.ndarray]):
        """Prints debug message with data

        Args:
            msg (str): Debug message
            data (Union[Sequence, _np.ndarray]): Data values
        """
        print(f"DEBUG>> {msg} : ", end="")
        for v in data:
            print(f"{v}, ", end="")
        print("\n")

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
            self.map_params.shiftback = _np.asarray(tr)


    def set_map_parameters_from_array(
        self,
        cell: Union[Sequence[float], _np.ndarray],
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
        """
        xyz = _np.array([0.0, 0.0, 0.0])
        with open(_path.join(self.datapath, "cutout.pdb"), "r") as fopen:
            for line in fopen:
                if "TRANSLATED BY" in line:
                    pattern = r"TRANSLATED BY \(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)"
                    match = re.search(pattern, line)
                    if match:
                        xyz[0], xyz[1], xyz[2] = map(float, match.groups())
        return xyz  # _Coord_orth(xyz[0],xyz[1],xyz[2])

    ClusterMode = Literal["dbscan", "kmeans"]
    def get_map_coords_from_predicted_instance(
        self,
        mode: ClusterMode = "dbscan",
        fix_origin: bool =True,
        mapin_path: str = "NONE",
        write_npy: bool = False,
        seqlen: int = 0,
        verbose: int = 0,
        # shiftback=True,
    ):
        """Return coordinates of amino acid instances from machine learning numpy output.

        Args:
            cell (clipper.Cell): clipper Cell object
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
            fix_origin (bool, optional): Flag to fix non-zero origin. Defaults to True.
            mapin_path (str, optional): Path to map file. Defaults to None.
            write_npy (bool, optional): Flag to write amino acid coordinates as .npy file. Defaults to False.
            return_map_index (bool, optional): Flag to return map indices instead of orthogonal coordinates. Defaults to True.
            verbose (int, optional): Verbosity. Defaults to 0
        Returns:
            List: list containing the orthogonal coordinates or map indices of amino acid instances
        """  # noqa: E501
        offsets = _np.load(f"{self.datapath}/inst_pred.npy")
        density = _np.load(f"{self.datapath}/density.npy")
        mp = self.map_params
        # nxs = 0
        # nys = 0
        # nzs = 0
        # ncorrect = 0
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
            #if fix_origin:
            #    nxs = gmap.header_i32(5)
            #    nys = gmap.header_i32(6)
            #    nzs = gmap.header_i32(7)
            #    ncorrect = 0
            #spacing = _np.array(
            #    [
            #        self.map_params.cell.a / gmap.header_i32(8),
            #        self.map_params.cell.b / gmap.header_i32(9),
            #        self.map_params.cell.c / gmap.header_i32(10),
            #    ]
            #)
        else:
            if mp.fix_axis_positions:  # will swap the x and z axes from the ML npy output
                # offset array is of shape (3,nu,nv,nw)
                offsets = _np.swapaxes(offsets, 1, 3)
                density = _np.swapaxes(density, 0, 2)
            #if fix_origin:
            #    nxs = -density.shape[0] // 2
            #    nys = -density.shape[1] // 2
            #    nzs = -density.shape[2] // 2
            #    ncorrect = 1
            #spacing = _np.array(
            #    [
            #        cell.a / float(density.shape[0]),
            #        cell.b / float(density.shape[1]),
            #        cell.c / float(density.shape[2]),
            #    ]
            #)
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
                raise ValueError("Please provide total number of sequence 'seqlen=' for KMeans clustering.")
            aa_instance_coordinates = self._kmeans_clustering(points, seqlen, write_npy)
        else:  # dbscan
            aa_instance_coordinates = self._dbscan_clustering(points, 10, write_npy)

        return aa_instance_coordinates

    def _dbscan_clustering(self, points: _np.ndarray, min_samples: int = 10, write_npy: bool = False):
        points_int = _np.round(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        mp = self.map_params
        aa_instance_coordinates = []
        counts = count[indices]
        db = _DBSCAN(
            eps=_np.sqrt(_np.sum(mp.spacing * mp.spacing)),
            min_samples=min_samples,
            algorithm="ball_tree",
        ).fit(points)
        # db = _DBSCAN(eps=_np.max(_np.sqrt(3 * (spacing * spacing))), min_samples=10).fit(
        #    points
        # )
        labels = db.labels_
        groups = {}
        weights = {}
        for label in _np.unique(labels):
            if label != -1:
                groups[label] = points[labels == label]
                weights[label] = counts[labels == label]
        # calculate center of the points with weights
        tmp = _np.zeros((len(groups), 3))
        ind = 0
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
            tmp[ind] = _np.array([sum_x / w, sum_y / w, sum_z / w]) - mp.shiftback
            aa_instance_coordinates.append(_Coord_orth(tmp[ind]))
            ind += 1
        if write_npy:
            _np.save("aa_inst_dbscan_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates
    
    def _kmeans_clustering(self, points: _np.ndarray, seqlen: int, write_npy: bool = False):
        kmeans = _KMeans(n_clusters=seqlen, init="k-means++", random_state=1).fit(points)
        cluster = kmeans.cluster_centers_

        aa_instance_coordinates = []
        tmp = _np.zeros((len(cluster), 3))
        for i in range(0, len(cluster)):
            co = _np.array([cluster[i][0],cluster[i][1],cluster[i][2]])
            tmp[i] = co - self.map_params.shiftback
            aa_instance_coordinates.append(_Coord_orth(tmp[i]))
        
        if write_npy:
            _np.save("aa_inst_kmeans_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates


    @staticmethod
    def iterative_cluster(points: _List): # n: int):
        points_int = _np.round(points).astype(int)
        uniq = _np.unique(points_int, axis=0)

        tmp = _np.zeros((len(uniq),3))
        for u in range(0, len(uniq)):
            pts = points[_np.all(points_int == uniq[u], axis=1)]
            sumpts = _np.sum(pts,axis=0)
            co = [sumpts[0]/len(pts), sumpts[1]/len(pts), sumpts[2]/len(pts)]
            tmp[u] = _np.array(co)
        return tmp

    @staticmethod
    def group_within_distance(uniq: _List, cutoff: float = 2.0):
        indices = _np.where(_spatial.distance.cdist(uniq,uniq,'euclidean') <= cutoff)[0]
        indices = _np.unique(indices)
        return indices

    def map_coords_to_ca_atom(
        self,
        aa_instance_coordinates: _List,
        mol: _MiniMol,
        aa_instance_indices: _List = [],
    ):  # , correction: Corrections):
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
            matom.pos = _Coord_orth(aa_instance_coordinates[i])  # *correction.spacing))
            if have_index:
                matom.set_property("INDEX", _Property_int(aa_instance_indices[i]))
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
        mol.set_model(mmodel)
        _ProteinTools.chain_label(mol, True)

    def coord_frac_to_ca_atom(self, aa_instance_frac: _List, mol: _MiniMol):
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

    


# del _Ca_group, _Cell, _Grid_sampling, _Coord_orth, _Coord_map
# del _MiniMol #_MAtom, _MModel, _MChain, _MRes
del annotations
