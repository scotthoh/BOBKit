from typing import List, Tuple, Union
import bobkit.clipper as clipper
import bobkit.util as util
import bobkit.buccaneer as buccaneer
import gemmi
import numpy as np
import dataclasses
import functools
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path

@dataclasses.dataclass
class Configuration:
    """Configuration for map grids preparation"""
    n_threads: Union[int, None] = None
    spacing: float = 0.7
    box_size: int = 128  # 32?
    overlap: int = 64  # 16??
    channels: int = 2  # one type of map a.a. backbone
    resolution: float = 2.0
    normalise: bool = True

class PrepareKaggleDataset:
    def __init__(self, datafile_path: str, structure_path: str, config: Configuration = Configuration(), column_names: Union[str, List, None] = None):
        self.input_data = Path(datafile_path)
        self.input_model = Path(structure_path)
        self.config = config
        self.column_names = column_names
        self.xmap = self._read_file_to_xmap()
        self.structure = util.read_structure(structure_path)

    def _read_file_to_xmap(self):
        xmap = clipper.Xmap_float()
        if self.input_data.suffix == ".mtz":
            mtz = gemmi.read_mtz_file(str(self.input_data))
            reso = clipper.Resolution(mtz.resolution_high())
            hklinfo = clipper.HKL_info.from_gemmi_mtz(mtz, generate=False)
            sample_rate = reso.limit() / self.config.spacing
            cell = clipper.Cell(mtz.cell.parameters)
            spg = clipper.Spacegroup(mtz.spacegroup_number)
            f_phi = clipper.HKL_data_F_phi_float(hklinfo)
            f_phi.import_from_gemmi(mtz, ",".join(self.column_names), True)
            grid = clipper.Grid_sampling(spg, cell, reso, sample_rate)
            xmap.init(clipper.Spacegroup(mtz.spacegroup_number), cell, grid)
            xmap.fft_from(f_phi)
        elif self.input_data.suffix == ".cif":
            cif = clipper.CIFfile()
            cif.open_read(str(self.input_data))
            reso = cif.resolution()
            cell = cif.cell
            spg = cif.spacegroup
            hklinfo = clipper.HKL_info()
            cif.import_hkl_info(hklinfo)
            f_phi = clipper.HKL_data_F_phi_float(hklinfo)
            cif.import_hkl_data(f_phi, self.column_names)
            cif.close_read()
            sample_rate = reso.limit() / self.config.spacing
            grid = clipper.Grid_sampling(spg, cell, reso, sample_rate)
            xmap.init(spg, cell, grid)
            xmap.fft_from(f_phi)
        else:
            gmap = gemmi.read_ccp4_map(str(self.input_data))
            xmap.import_from_gemmi(gmap)
        # normalise data
        xmap.normalise()
        return xmap

    def _mtz_to_map(self, mtz_path, normalise: bool = True):
        gemmi.read_mtz_file()

    def precompute_slices(self, grid_shape: np.ndarray, overlap: int = 16) -> List[List[int]]:
        """Precompute indices of slices for inference"""
        slices = []

        for i in range(0, grid_shape[0] - overlap, overlap):
            for j in range(0, grid_shape[1] - overlap, overlap):
                for k in range(0, grid_shape[2] - overlap, overlap):
                    slices.append([i, j, k])

        return slices

    def chop_map_grids(self, array: np.ndarray[np.float32], box_size: int, translation: Tuple[int, int, int]):
        i, j, k = translation
        subarr = array[i : i + box_size, j : j + box_size, k : k + box_size]
        # subarr = subarr[np.newaxis, ..., np.newaxis].astype(np.float16)
        return subarr, translation

    def label_volume(self, array: np.ndarray[np.float32], transform_stucture: clipper.MiniMol, radius: float = 1.5):
        cx = array.shape[0] * self.config.spacing
        cy = array.shape[1] * self.config.spacing
        cz = array.shape[2] * self.config.spacing

        cell = clipper.Cell([cx, cy, cz, 90., 90., 90.])
        zeros = np.zeros(array.shape, dtype=np.float32)
        mask = clipper.NXmap_float(zeros, cell)
        edcalc = clipper.EDcalc_mask_float(radius)
        edcalc(mask, transform_stucture.atom_list())
        array[mask.array > 0.0] = 1.0

    def make_map_grids(self, output_path: Union[Path, str], write_npz: bool=True):
        work_grid, rtop = util.interpolate_map_grid(
            self.xmap,
            self.config.spacing,
            self.config.box_size,
            self.config.overlap,
            3,
            whole_unitcell=False,
            mapout=False,
        )
        transform_st = self.structure.clone()
        transform_st.transform(rtop.inverse())
        trim_st = transform_st.clone()
        trim_st_nooxy = transform_st.clone()
        buccaneer.ProteinTools.remove_sidechain(trim_st)
        buccaneer.ProteinTools.remove_sidechain(trim_st_nooxy, False)
        label_array1 = np.zeros(work_grid.shape, dtype=np.float32)
        label_array2 = np.zeros(work_grid.shape, dtype=np.float32)
        self.label_volume(label_array1, trim_st)
        self.label_volume(label_array2, trim_st_nooxy)
        output_file = ""
        if (write_npz):
            outfile_name = self.input_data.stem + "_interp.npz"
            output_file = Path(output_path) / outfile_name
            np.savez_compressed(output_file, volume=work_grid, label_ncaco=label_array1, label_ncac=label_array2)
        # pytorch channel at the front; Shape (3, n, n, n)
        starting_array = np.stack([work_grid, label_array1, label_array2], axis=0)
        # work_grid_shape = np.array(work_grid.shape)
        # slices = self.precompute_slices(work_grid_shape, overlap=self.config.overlap)
        return starting_array, output_file


def prepare_maps(
    file_list: Union[List[Path], List[str]],
    struture_list: Union[List[Path], List[str]],
    output_path: str,
    amplitude: str = "FWT",
    phase: str = "PHWT",
    config: Configuration = Configuration(),
):
    column_names = {amplitude, phase}
    done = np.array(([False] * len(file_list)))
    count = 0
    for inputdata, inputstructure in zip(file_list, struture_list):
        print(f"{inputdata}, {inputstructure}")
        prepare_kaggle = PrepareKaggleDataset(
            inputdata, inputstructure, config, column_names
        )
        starting_array, output_file = prepare_kaggle.make_map_grids(output_path, write_npz=True)
        done[count] = output_file.exists()
        count += 1
        
    notprocessed = np.array(inputdata)[~done]
    if notprocessed.size != 0:
        print("The following are not processed: ")
        for fname in notprocessed:
            print(fname)



if __name__ == "__main__":
    import sys
    flist = sys.argv[1]
    output_path = sys.argv[2]
    with open(flist, 'r') as fopen:
        file_list = []
        st_list = []
        lines = fopen.readlines()
        for line in lines:
            f, s = line.strip().split(',')
            file_list.append(f.strip())
            st_list.append(s.strip())

    prepare_maps(file_list, st_list, output_path)
    # def prepare_map_slices(self, work_array: np.ndarray, slices: List[List[int]]):
    #    channels = self.config.channels
    #    box_size = self.config.box_size
    #    output_shape = (channels, box_size, box_size, box_size)
    #    process_grid_worker = functools.partial(
    #        self.chop_map_grids,
    #        work_array,
    #        box_size,
    #    )
    #    # get a list of x-ray maps to run separate files for labels and volume
    #    miniter = 1000 if len(slices) > 10000 else 1
    #    max_workers = self.config.n_threads
    #    if max_workers == 1:
    #        grid_samples = list(map(process_grid_worker, slices))
    #    else:
    #        with ThreadPoolExecutor(max_workers=max_workers) as executor:
    #            grid_samples = list(executor.map(process_grid_worker, slices))

    # def cut_amino_acid_backbones(self, radius:)
