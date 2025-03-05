from __future__ import annotations
from typing import List as _List
from bobkit.util import *
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

from os import path as _path
from dataclasses import dataclass as _dataclass
import re


@_dataclass
class Corrections:
    """Class for keeping the spacing, origin and correction value"""

    workcell: _Cell = _Cell()
    cell: _Cell = _Cell()  # _np.ndarray = _np.array([0.,0.,0.])
    grid: _Grid_sampling = None
    grid_asu: _Grid_sampling = None
    origin: _np.ndarray = _np.array([0, 0, 0])
    spacing: _np.ndarray = _np.array([1.0, 1.0, 1.0])
    ncorrect: int = 0
    fix_origin: bool = False
    fix_axis_positions: bool = False
    shiftback: _Coord_orth = _Coord_orth(0.0, 0.0, 0.0)

    def set_grid(self, nu: int, nv: int, nw: int):
        self.grid = _Grid_sampling(nu, nv, nw)

    def set_grid_asu(self, nu: int, nv: int, nw: int):
        self.grid_asu = _Grid_sampling(nu, nv, nw)

    def set_origin(self, x: int, y: int, z: int):
        self.origin[0] = x
        self.origin[1] = y
        self.origin[2] = z

    def set_spacing(self, x: float, y: float, z: float):
        self.spacing[0] = x
        self.spacing[1] = y
        self.spacing[2] = z

    def set_shiftback(self, x: float, y: float, z: float):
        self.shiftback = _Coord_orth(x, y, z)


class Osaka:
    """Class with methods and static methods to process the outputs from
    machine learning to predict amino acid instance and sequence from map segmentation
    by Osaka group.
    """
    def __init__(
        self, datapath: str = None, workcell: type[_Cell | _gemmi.UnitCell] = None
    ):
        if datapath is not None:
            self.datapath = datapath
        else:
            self.datapath = None
        self.corrections = Corrections()
        self.corrections.workcell.init(workcell)

    def grid_spacing_from_map(
        self,
        mapin_path: str,
        fix_axis_positions=False,
        fix_origin=False,
        shiftback=True,
    ):
        """Get grid spacing from map

        Args:
            cell (clipper.Cell): clipper Cell object
            mapin_path (str): Path to map file
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
        """
        gmap = _gemmi.read_ccp4_map(mapin_path)
        if fix_axis_positions and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
            gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
            self.corrections.fix_axis_positions = fix_axis_positions
        spacing = _np.array(
            [
                gmap.grid.unit_cell.a / gmap.header_i32(8),
                gmap.grid.unit_cell.b / gmap.header_i32(9),
                gmap.grid.unit_cell.c / gmap.header_i32(10),
            ]
        )
        self.corrections.set_grid_asu(gmap.grid.nu, gmap.grid.nv, gmap.grid.nw)
        self.corrections.set_grid(
            gmap.header_i32(8), gmap.header_i32(9), gmap.header_i32(10)
        )
        self.corrections.cell.init(gmap.grid.unit_cell)
        self.corrections.set_spacing(spacing[0], spacing[1], spacing[2])
        if fix_origin:
            self.corrections.set_origin(
                gmap.header_i32(5), gmap.header_i32(6), gmap.header_i32(7)
            )
            self.ncorrect = 0
            self.corrections.fix_origin = fix_origin
        if shiftback:
            tr = self.get_translation_from_cutout()
            self.corrections.set_shiftback(tr[0], tr[1], tr[2])
        # return spacing

    def get_spacing_from_array(
        self, cell: _Cell, fix_axis_positions=False, fix_origin=False
    ):
        """Get grid spacing from density array

        Args:
            cell (clipper.Cell): clipper Cell object
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.

        Returns:
            _type_: _description_
        """
        density = _np.load(f"{self.datapath}/density.npy")
        if fix_axis_positions:  # will swap the x and z axes from the ML npy output
            density = _np.swapaxes(density, 0, 2)
        spacing = _np.array(
            [
                cell.a / density.shape[0],
                cell.b / density.shape[1],
                cell.c / density.shape[2],
            ]
        )
        self.corrections.set_cell(cell.a, cell.b, cell.c)
        self.corrections.set_spacing(spacing[0], spacing[1], spacing[2])
        if fix_origin:
            nxs = -density.shape[0] // 2
            nys = -density.shape[1] // 2
            nzs = -density.shape[2] // 2
            self.corrections.set_origin(nxs, nys, nzs)
            self.corrections.ncorrect = 1
        # return spacing

    def get_translation_from_cutout(self):
        xyz = [0.0] * 3
        with open(_path.join(self.datapath, "cutout.pdb"), "r") as fopen:
            for line in fopen:
                if "TRANSLATED BY" in line:
                    pattern = r"TRANSLATED BY \(\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)"
                    match = re.search(pattern, line)
                    if match:
                        xyz[0], xyz[1], xyz[2] = map(float, match.groups())
        return _np.array(xyz)  # _Coord_orth(xyz[0],xyz[1],xyz[2])

    def write_aa_instance_file(
        self,
        cell: _Cell,
        grid_sampling: _Grid_sampling,
        fix_axis_positions=False,
        fix_origin=False,
        mapin_path: str = None,
    ):

        offsets = _np.load(f"{self.datapath}/inst_pred.npy")
        density = _np.load(f"{self.datapath}/density.npy")
        nxs = 0
        nys = 0
        nzs = 0
        if mapin_path != "NONE":
            gmap = _gemmi.read_ccp4_map(mapin_path)
            if fix_axis_positions and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
                axis_pos = _np.array(gmap.axis_positions())
                offsets = _np.swapaxes(
                    offsets, axis_pos[0], _np.where(axis_pos == 0)[0][0]
                )
                offsets = _np.swapaxes(
                    offsets, axis_pos[1], _np.where(axis_pos == 1)[0][0]
                )
                density = _np.swapaxes(
                    density, axis_pos[0], _np.where(axis_pos == 0)[0][0]
                )
                density = _np.swapaxes(
                    density, axis_pos[1], _np.where(axis_pos == 1)[0][0]
                )
                gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
            if fix_origin:
                nxs = gmap.header_i32(5)
                nys = gmap.header_i32(6)
                nzs = gmap.header_i32(7)
                ncorrect = 0
        else:
            if fix_axis_positions:  # will swap the x and z axes from the ML npy output
                offsets = _np.swapaxes(offsets, 1, 3)
                density = _np.swapaxes(density, 0, 2)
            if fix_origin:
                nxs = -density.shape[0] // 2
                nys = -density.shape[1] // 2
                nzs = -density.shape[2] // 2
                ncorrect = 1

        if fix_origin:
            xyz = _np.mgrid[
                nxs : density.shape[0] + (nxs + ncorrect),
                nys : density.shape[1] + (nys + ncorrect),
                nzs : density.shape[2] + (nzs + ncorrect),
            ].astype(_np.float64)
        else:
            xyz = _np.mgrid[
                nxs : density.shape[0], nys : density.shape[1], nzs : density.shape[2]
            ].astype(_np.float64)

        xyz0 = _np.mgrid[
            0 : density.shape[0], 0 : density.shape[1], 0 : density.shape[2]
        ].astype(_np.float64)
        print(f"xyz shape : {xyz.shape}")
        print(f"offsets shape : {offsets.shape}")
        print(f"density shape : {density.shape}")
        print(f"grid sampling : {grid_sampling}")

        points = (xyz + offsets) * (density > 0)
        points0 = (xyz0 + offsets) * (density > 0)
        # grid map equivalent
        points = points[:, density > 0].T
        points0 = points0[:, density > 0].T
        points_int = _np.round(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        aa_instance_coordinates = []
        aa_instance_grid_indices = []
        aa_cg = []
        counts = count[indices]
        # use dbscan to group nearby points together
        db = _DBSCAN(
            eps=1.0,
            min_samples=10,
            algorithm="ball_tree",
        ).fit(points)
        # db = _DBSCAN(eps=_np.max(_np.sqrt(3 * (spacing * spacing))), min_samples=10).fit(
        #    points
        # )
        labels = db.labels_
        groups = {}
        groups0 = {}
        weights = {}
        for label in _np.unique(labels):
            if label != -1:
                groups[label] = points[labels == label]
                groups0[label] = points0[labels == label]
                weights[label] = counts[labels == label]
        tmp = _np.zeros((len(groups), 3))
        count = 0
        # for label, group in groups.items():
        for i in zip(groups.items(), groups0.items()):
            label = i[0][0]
            group = i[0][1]
            label0 = i[1][0]
            group0 = i[1][1]
            sum_x = 0.0
            sum_y = 0.0
            sum_z = 0.0
            w = 0.0
            sum_x0 = 0.0
            sum_y0 = 0.0
            sum_z0 = 0.0
            w0 = 0.0
            for p in zip(group, weights[label]):
                sum_x += p[0][0] * p[1]
                sum_y += p[0][1] * p[1]
                sum_z += p[0][2] * p[1]
                w += p[1]
            for p in zip(group0, weights[label0]):
                sum_x0 += p[0][0] * p[1]
                sum_y0 += p[0][1] * p[1]
                sum_z0 += p[0][2] * p[1]
                w0 += p[1]
            tmp[count] = _np.array([sum_x / w, sum_y / w, sum_z / w])
            cm0 = _Coord_map(_Vec3_double([sum_x0 / w0, sum_y0 / w0, sum_z0 / w0]))
            cm = _Coord_map(_Vec3_double(tmp[count]))
            co = cm.coord_frac(grid_sampling).coord_orth(cell)
            print(f"{cm0[0]:10.6f}, {cm0[1]:10.6f}, {cm0[2]:10.6f}")
            print(
                f"{count:>4}, [ {tmp[count][0]:10.6f}, {tmp[count][1]:10.6f}, {tmp[count][2]:10.6f} ] || {co} || ",
                end="",
            )
            aa_instance_coordinates.append(co)
            # check if points are in density
            cg = _Coord_map(
                _Vec3_double(tmp[count] - _np.array([nxs, nys, nzs]))
            ).coord_grid()
            # cg =cm.coord_grid()
            density_grid = _Grid_sampling(
                density.shape[0], density.shape[1], density.shape[2]
            )
            aa_instance_grid_indices.append(cg.unit(density_grid).index(density_grid))
            aa_cg.append([cg.u, cg.v, cg.w])
            print(f"{cg} | {cg.unit(density_grid)}")
            # if density[cg.u, cg.v, cg.w] <= 0.:
            #    print(f"{cg} is not in density!")
            # else:
            #   print(f"{cg} : {density[cg.u, cg.v, cg.w]}")
            # aa_instance_coordinates.append(_Coord_orth(tmp[count]))
            count += 1

        _np.save("cg_inst.npy", aa_cg, allow_pickle=False)
        # if write_npy:
        #    _np.save("aa_inst_mean.npy", tmp, allow_pickle=False)
        # self.corrections.set_cell(spacing[0], spacing[1], spacing[2])
        # self.corrections.set_origin(nxs, nys, nzs)
        # self.corrections.ncorrect = ncorrect
        return aa_instance_coordinates, aa_instance_grid_indices

    def get_map_coords_from_predicted_instance(
        self,
        cell: _Cell,
        grid_sampling: _Grid_sampling,
        fix_axis_positions=False,
        fix_origin=True,
        mapin_path: str = None,
        write_npy=False,
        return_map_index=False,
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
        """
        offsets = _np.load(f"{self.datapath}/inst_pred.npy")
        density = _np.load(f"{self.datapath}/density.npy")
        nxs = 0
        nys = 0
        nzs = 0
        ncorrect = 0
        if mapin_path != "NONE":
            gmap = _gemmi.read_ccp4_map(mapin_path)
            if fix_axis_positions and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
                axis_pos = _np.array(gmap.axis_positions())
                offsets = _np.swapaxes(offsets, axis_pos[0], _np.where(axis_pos == 0)[0][0])
                offsets = _np.swapaxes(offsets, axis_pos[1], _np.where(axis_pos == 1)[0][0])
                density = _np.swapaxes(density, axis_pos[0], _np.where(axis_pos == 0)[0][0])
                density = _np.swapaxes(density, axis_pos[1], _np.where(axis_pos == 1)[0][0])
                gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
            if fix_origin:
                nxs = gmap.header_i32(5)
                nys = gmap.header_i32(6)
                nzs = gmap.header_i32(7)
                ncorrect = 0
            spacing = _np.array(
                [
                    cell.a / gmap.header_i32(8),
                    cell.b / gmap.header_i32(9),
                    cell.c / gmap.header_i32(10),
                ]
            )
        else:
            if fix_axis_positions:  # will swap the x and z axes from the ML npy output
                offsets = _np.swapaxes(offsets, 1, 3)
                density = _np.swapaxes(density, 0, 2)
            if fix_origin:
                nxs = -density.shape[0] // 2
                nys = -density.shape[1] // 2
                nzs = -density.shape[2] // 2
                ncorrect = 1
            spacing = _np.array(
                [
                    cell.a / float(density.shape[0]),
                    cell.b / float(density.shape[1]),
                    cell.c / float(density.shape[2]),
                ]
            )
        if fix_origin:
            xyz = _np.mgrid[
                nxs : density.shape[0] + (nxs + ncorrect),
                nys : density.shape[1] + (nys + ncorrect),
                nzs : density.shape[2] + (nzs + ncorrect),
            ].astype(_np.float64)
        else:
            xyz = _np.mgrid[
                nxs : density.shape[0], nys : density.shape[1], nzs : density.shape[2]
            ].astype(_np.float64)
        if verbose > 5:
            print(f"xyz shape : {xyz.shape}")
            print(f"offsets shape : {offsets.shape}")
            print(f"density shape : {density.shape}")
            print(f"spacing: {spacing}")

        for i in (0, 1, 2):
            offsets[i] = offsets[i] * spacing[i]
            xyz[i] = xyz[i] * spacing[i]
        points = (xyz + offsets) * (density > 0)
        points = points[:, density > 0].T
        # if shiftback:
        #    points -= self.corrections.shiftback
        # tr = self.get_translation_from_cutout()
        # points -= tr

        points_int = _np.round(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        aa_instance_coordinates = []
        counts = count[indices]
        # use dbscan to group nearby points together
        db = _DBSCAN(
            eps=_np.sqrt(_np.sum(spacing * spacing)),
            min_samples=10,
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

        tmp = _np.zeros((len(groups), 3))
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
            tmp[count] = _np.array([sum_x / w, sum_y / w, sum_z / w])
            # aa_instance_coordinates.append(_Coord_orth(tmp[count]).coord_frac(cell).coord_grid(grid_sampling))
            if return_map_index:
                aa_instance_coordinates.append(
                    (_Coord_orth(tmp[count]) - self.corrections.shiftback)
                    .coord_frac(cell)
                    .coord_grid(grid_sampling)
                )
            else:
                aa_instance_coordinates.append(
                    (_Coord_orth(tmp[count]) - self.corrections.shiftback)
                )
            count += 1
        if write_npy:
            _np.save("aa_inst_mean.npy", tmp, allow_pickle=False)
        # if return_map_index:
        #    db = _DBSCAN(eps=2., min_samples=10, algorithm='ball_tree').fit(points)
        #    if write_npy:
        #        _np.save("aa_inst_mean.npy", uniq, allow_pickle=False)
        #    labels = db.labels_
        #    groups = {}
        #    weights = {}
        #    for label in _np.unique(labels):
        #        if label != -1:
        #            groups[label] = points_int[labels == label]
        #            weights[label] = counts[labels == label]
        #
        #    tmp = _np.zeros((len(groups), 3))
        #    count = 0
        #    for label, group in groups.items():
        #        sum_x = 0
        #        sum_y = 0
        #        sum_z = 0
        #        w = 0
        #        for p in zip(group, weights[label]):
        #            sum_x += p[0][0] * p[1]
        #            sum_y += p[0][1] * p[1]
        #            sum_z += p[0][2] * p[1]
        #            w += p[1]
        #        tmp[count] = _np.array([sum_x / w, sum_y / w, sum_z / w])
        #        #aa_instance_coordinates.append(_Coord_orth(tmp[count]).coord_frac(cell).coord_grid(grid_sampling))
        #        aa_instance_coordinates.append(_Coord_orth(tmp[count]))
        #        count += 1
        #    for point_ind in uniq:
        #        cg = _Coord_orth(point_ind.astype(float)).coord_frac(cell).coord_grid(grid_sampling)
        #        aa_instance_coordinates.append(cg)
        #        #aa_instance_coordinates.append(_Coord_map(point_ind))
        # else:

        self.corrections.set_spacing(spacing[0], spacing[1], spacing[2])
        self.corrections.set_origin(nxs, nys, nzs)
        self.corrections.ncorrect = ncorrect
        return aa_instance_coordinates

    # @staticmethod
    def map_coords_to_ca_atom(
        self,
        aa_instance_coordinates: _List,
        mol: _MiniMol,
        aa_instance_indices: _List = [],
    ):  # , correction: Corrections):
        _ProteinTools()
        # mmodel = mol.model()
        # correction = _np.array([1.,1.,1.])
        # if correction.spacing is not None:
        #    correction = _np.array([spacing[0],spacing[1],spacing[2]])
        mmodel = _MModel()
        have_index = True if len(aa_instance_indices) > 0 else False
        for i in range(0, len(aa_instance_coordinates)):
            # for co in aa_instance_coordinates:
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
            matom.pos = _Coord_orth(aa_instance_coordinates[i])  # *correction.spacing))
            if have_index:
                matom.set_property("INDEX", _Property_int(aa_instance_indices[i]))
            # print(matom)
            res.insert(matom, -1)
            chn.insert(res, -1)
            mmodel.insert(chn, -1)
            # print("len model ", len(mmodel))
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
    ):  # , spacing: _List = None):
        _ProteinTools()
        # mmodel = mol.model()
        # correction = _np.array([1.,1.,1.])
        # if spacing is not None:
        #    correction = _np.array([spacing[0],spacing[1],spacing[2]])
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

    def orth_to_map(self, aa_instance_coordinates: _List, xmap: _Xmap_float):
        """Convert orthogonal coordinates to map coordinates

        Args:
            aa_instance_coordinates (_List): List of orthogonal coordinates
            xmap (_Xmap_float): Xmap object
        """
        cm_list = []
        for co in aa_instance_coordinates:
            cm = xmap.coord_map(co)
            print(f"{co} || {cm}")
            cm_list.append(cm)

        return cm_list


# del _Ca_group, _Cell, _Grid_sampling, _Coord_orth, _Coord_map
# del _MiniMol #_MAtom, _MModel, _MChain, _MRes
del annotations
