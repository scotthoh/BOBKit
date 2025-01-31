from bobkit.util import *
__all__ = [
    "test_array",
    "read_structure",
    "write_structure",
    "get_coordinates_from_predicted_instance",
]

from bobkit.buccaneer import Ca_group as _Ca_group

from bobkit.clipper import (
    Cell as _Cell,
    Grid_sampling as _Grid_sampling,
    Coord_orth as _Coord_orth,
    Coord_map as _Coord_map,
)
import numpy as _np
import gemmi as _gemmi
from sklearn.cluster import DBSCAN as _DBSCAN

from os import path as _path

class Osaka:
    """Class with methods and static methods to process the outputs from
    machine learning to predict amino acid instance and sequence from map segmentation
    by Osaka group.
    """
    def __init__(self, datapath: str = None):
        if datapath is not None:
            self.datapth = _path.abspath(datapath)
        else:
            self.datapath = None
    
    def grid_spacing_from_map(self, cell: _Cell, mapin_path: str, fix_axis_positions=False):
        """Get grid spacing from map

        Args:
            cell (clipper.Cell): clipper Cell object
            mapin_path (str): Path to map file
            fix_axis_positions (bool, optional): Flag to fix axis positions to XYZ. Defaults to False.
        """
        gmap = _gemmi.read_ccp4_map(mapin_path)
        if fix_axis_positions and (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
            gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        spacing = _np.array(
            [
                cell.a / gmap.header_i32(8),
                cell.b / gmap.header_i32(9),
                cell.c / gmap.header_i32(10),
            ]
        )
        return spacing
        
    def get_spacing_from_array(self, cell: _Cell, fix_axis_positions=False):
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
        return spacing

    def get_map_coords_from_predicted_instance(
        self,
        cell: _Cell,
        fix_axis_positions=False,
        fix_origin=True,
        mapin_path: str = None,
        write_npy=False,
        return_map_index=True,
        verbose: int = 0,
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
                    cell.a / density.shape[0],
                    cell.b / density.shape[1],
                    cell.c / density.shape[2],
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
        if not return_map_index:
            for i in (0, 1, 2):
                offsets[i] = offsets[i] * spacing[i]
                xyz[i] = xyz[i] * spacing[i]
        points = (xyz + offsets) * (density > 0)
        points = points[:, density > 0].T
        points_int = _np.floor(points).astype(int)
        uniq, indices, count = _np.unique(
            points_int, return_inverse=True, return_counts=True, axis=0
        )
        aa_instance_coordinates = []
        if return_map_index:
            if write_npy:
                _np.save("aa_inst_mean.npy", uniq, allow_pickle=False)
            for point_ind in uniq:
                aa_instance_coordinates.append(_Coord_map(point_ind))  
        else:
            counts = count[indices]
            # use dbscan to group nearby points together
            db = _DBSCAN(eps=_np.max(_np.sqrt(3 * (spacing * spacing))), min_samples=10).fit(
                points
            )
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
                aa_instance_coordinates.append(_Coord_orth(tmp[count]))
                count += 1
            if write_npy:
                _np.save("aa_inst_mean.npy", tmp, allow_pickle=False)
        return aa_instance_coordinates

del _Ca_group, _Cell, _Grid_sampling, _Coord_orth, _Coord_map
del _np
del _gemmi
del _DBSCAN
del _path