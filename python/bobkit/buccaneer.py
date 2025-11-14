#from __future__ import annotations
#from typing import List as _List
from bobkit._buccaneer import *

__all__ = [
    "Ca_build",
    "Ca_chain",
    "Ca_correct",
    "Ca_filter",
    "Ca_find",
    "Ca_group",
    "Ca_grow",
    "Ca_join",
    "Ca_link",
    "Ca_merge",
    "Ca_ncsbuild",
    "Ca_prep",
    "Ca_prune",
    "Ca_sequence",
    "CoordList_5",
    "CoordList_8",
    "Grow_threaded",
    "KnownStructure",
    "LLK_Targetlist",
    "LLK_map_target",
    "Log",
    "MapSimulate",
    "ModelTidy",
    "Optimiser_simplex",
    "Pr_group",
    "Prep_threaded",
    "ProteinLoop",
    "ProteinTools",
    "SSfind",
    "Score_list_RTop_orth",
    "Score_list_String",
    "Score_list_float_array",
    "ScoreResult",
    "Search_threaded",
    "Sequence_score_threaded",
    "Sequence_threaded",
    "Target_fn_order_zero",
    "Target_fn_refine_llk_map_target",
    "Target_fn_refine_n_terminal_build",
    "set_reference",
    #"Ca_sequence_ml",
]

# from _bobkit._buccaneer import (
#    Ca_sequence as _Caseq,
#    Ca_group as _Cagroup,
#    LLK_map_target as _LLKtgt,
#    LLK_TargetList as _LLK_TargetList,
#    ProteinTools as _PT,
# )
#from bobkit.clipper import (
#    Cell as _Cell,
#    Spacegroup as _Spacegroup,
#    Grid_sampling as _Grid_sampling,
#    Coord_orth as _Coord_orth,
#    Coord_grid as _Coord_grid,
#    MiniMol as _Minimol,
#    MChain as _MChain,
#    MResidue as _MRes,
#    Xmap_float as _Xmap,
#    MMoleculeSequence as _MMolSeq,
#    MChain as _MChain,
#    Thread_base as _Thread_base,
#    Property_sequence_data as _Property_sequence_data,
#)

#from bobkit.util import write_structure as _write_structure
#from ._util import MapParameters as _MapParameters
#import numpy as _np
#import gemmi as _gemmi
#from multiprocessing import Process as _Process
'''
ProteinTools()

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
    0: 0,   # ALA
    1: 7,   # GLY
    2: 9,   # ILE
    3: 10,  # LEU
    4: 14,  # PRO
    5: [16, 19], # VAL/THR
    6: 13,  # PHE
    7: 17,  # TRP
    8: 18,  # TYR
    9: [2, 3],   # ASP/ASN
    10: [5, 6],  # GLU/GLN
    11: 1,  # ARG
    12: 8,  # HIS
    13: 11, # LYS
    14: 15, # SER
    15: 4,  # CYS
    16: 12, # MET
}

BUC_TO_MLI_DICT = {
    0: 0,   # ALA
    1: 11,  # ARG
    2: 9,   # ASN
    3: 9,   # ASP
    4: 15,  # CYS
    5: 10,  # GLN
    6: 10,  # GLU
    7: 1,   # GLY
    8: 12,  # HIS
    9: 2,   # ILE
    10: 3,  # LEU
    11: 13, # LYS
    12: 16, # MET
    13: 6,  # PHE
    14: 4,  # PRO
    15: 14, # SER
    16: 5,  # THR
    17: 7,  # TRP
    18: 8,  # TYR
    19: 5,  # VAL
    12: 16, # MSE
}


class Ca_sequence_ml:
    """Class for sequence Ca chains using density and output predictions from machine learning.
    This is based on the 
    """

    def __init__(
        self,
        sequence_array: _np.ndarray = _np.full((1, 1), None, dtype=object),
        corrections: _MapParameters = _MapParameters(),
        reliability: float = 1.0,
        correl: bool = False,  # , fix_axis_positions: bool = False,
    ):
        # for arr in sequence_array:
        #    print(arr.shape)
        self.__sequence_array = (
            sequence_array  # self.process_sequence_array(sequence_array, corrections)
        )
        # for arr in self.__sequence_array:
        #    print(arr.shape)

        # self.__sequence_array = sequence_array
        #self.__reliability = reliability
        self.__ncpu = 1
        self.__corrections = corrections
        self.__mol = None
        self.__done = []
        self.__correl = correl
        ## how to utilise this p1 and sg maps to map the grid coords to the array.
        ##self.__xlookp1 = _Xmap(_Spacegroup.p1(), corrections.workcell, corrections.)
        ##clipper::Xmap<int> xlookp1( clipper::Spacegroup::p1(), cell, grid );
        ##int lresult = 0;
        ##{
        ##  clipper::Xmap<int> xlooksg( spgr, cell, grid );
        ##  for ( MRI ix = xlooksg.first(); !ix.last(); ix.next() )
        ##    xlooksg[ix] = lresult++;
        ##  for ( MRI ix = xlookp1.first(); !ix.last(); ix.next() )
        ##    xlookp1[ix] = xlooksg.get_data(ix.coord());
        ##}
        # self.__fix_axis_positions = fix_axis_positions

    # def __init__(self, sequence_array: _np.ndarray, reliability: float = 0.5):
    #    """Initialise class with numpy array containing probability of amino acid sequence from
    #    machine learning segmentation.
    #
    #    Args:
    #        sequence_array (_np.ndarray): Numpy array of arrays containing probability of amino acid sequence
    #        reliability (float, optional): Reliability of sequence provided. Defaults to 0.5.
    #    """
    #    self.__sequence_array = sequence_array
    #    self.__reliability = reliability
    #    self.__ncpu = 1
    #    self.__semet = False
    @staticmethod
    def process_sequence_array(sequence_array, corrections):
        arraylist = []
        if corrections.fix_axis_positions:
            for arr in sequence_array:
                arr = _np.swapaxes(arr, 0, 2)
                arraylist.append(arr)
                print(arr.shape)
        return arraylist

    @staticmethod
    def pad_sequence_array(sequence_array, corrections):
        arraylist = []
        grid = _np.array(corrections.grid[0], corrections.grid[1],corrections.grid[2])
        arrayshape = sequence_array[0].shape
        if (not _np.array_equal(grid,arrayshape)):
            diff = [grid[i]-arrayshape[i] for i in range(0,3)]
            #newarray = 

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
        return CH_TO_AA_DICT[Ca_sequence_ml.buc_aa2mlindex(ml_aa_ind)]

    @staticmethod
    def get_probability_value(
        ind: int,
        pos: _Coord_orth,
        seq_array: _np.ndarray,
        corrections: _MapParameters,
        cell: _Cell,
        grid: _Grid_sampling,
        correl: bool = False,
        shiftback: bool = False,
        grid_index: int = -1,
    ):
        coords = _Coord_orth(pos)
        if shiftback:
            coords = coords + corrections.shiftback
    
        #if corrections.fix_origin:
        #    coords = coords - _Coord_orth(corrections.origin)
        # if corrections.fix_axis_positions:
        #    coords = _Coord_orth(pos.z,pos.y,pos.x)
        # else:
        #    coords = _Coord_orth(pos.x,pos.y,pos.z)
        # cannot use the line below because the grid
        # cg1 = coords.coord_frac(corrections.cell).coord_grid(grid).unit(grid)
        if grid_index != -1:
            cg = corrections.grid_asu.deindex(grid_index)
        else:
            #cg = _Coord_grid(int(coords.x/corrections.spacing[0]), int(coords.y/corrections.spacing[1]), int(coords.y/corrections.spacing[2]))
            cg = coords.coord_frac(cell).coord_grid(grid)
            #if ind == 0:
            #    print(f"DEBUG : {coords.format()}, {cg.format()}", end=",")
            if corrections.fix_origin:
                cg = cg - _Coord_grid(corrections.origin[0], corrections.origin[1], corrections.origin[2])
            cg = cg.unit(_Grid_sampling(corrections.grid_asu[0], corrections.grid_asu[1], corrections.grid_asu[2]))
            #if ind == 0:
            #    print(f" {cg.format()}, {corrections.origin}" )
        # print(f"{cg}")
        # if corrections.fix_origin:
        #    cg = cg - _Coord_grid(corrections.origin[0] + corrections.ncorrect, corrections.origin[1] + corrections.ncorrect, corrections.origin[2] + corrections.ncorrect)
        # print(f"DEBUG : {pos}, {cg}")
        # mapind = [int(pos.x), int(pos.y), int(pos.z)]
        #try:
        probval = seq_array[cg.u, cg.v, cg.w]
        ##except IndexError:
        ##    cg = cg.unit(corrections.grid_asu)
        ##    probval = seq_array[cg.u, cg.v, cg.w]
        #print(f"{probval}, ", end="")
        # probval = seq_array[cg.u%seq_array.shape[0], cg.v%seq_array.shape[1], cg.w%seq_array.shape[2]]
        # print(f"{coords} || {cg.u%seq_array.shape[0], cg.v%seq_array.shape[1], cg.w%seq_array.shape[2]} || {probval}")
        #if correl:
        if ind in [5, 9, 10]:
            return -probval / 2.0
        else:
            return -probval
        #else:
        #    if probval > 0.0:
        #        if ind in [5, 9, 10]:
        #            llkval = _np.log10(probval / 2.0)
        #        else:
        #            llkval = _np.log10(probval)
        #    elif probval == 0.0:
        #        llkval = 0.0
        #
        #    return llkval

    @staticmethod
    def shift_model(mol: _Minimol, corrections: _MapParameters, grid: _Grid_sampling):
        for chn in range(0, mol.size()):
            for r in range(0, mol[chn].size()):
                for a in range(0, mol[chn][r].size()):
                    coords = mol[chn][r][a].pos
                    cg = coords.coord_frac(corrections.cell).coord_grid(
                        grid
                    )  # corrections.grid)
                    if corrections.fix_origin:
                        cg = cg - _Coord_grid(
                            corrections.origin[0],
                            corrections.origin[1],
                            corrections.origin[2],
                        )
                    coords2 = cg.coord_frac(grid).coord_orth(corrections.cell)
                    mol[chn][r][a].pos = coords2
                    # print(cg)
            #    break
            # break

    def apply_ml_sequence(
        self,
        mol: _Minimol,
        cell: _Cell,
        grid: _Grid_sampling,
        ntype: int = 20,
        shiftback: bool = False,
    ):
        done = [False] * mol.size()
        for chn in range(0, mol.size()):
            for r in range(0, mol[chn].size()):
                scores = [0.0] * ntype
                ind = -1
                if mol[chn][r]["CA"].exists_property("INDEX"):
                    ind = int(mol[chn][r]["CA"].get_property("INDEX").value)
                    # cg= self.__corrections.grid_asu.deindex(ind)
                for t in range(0, ntype):
                    ml_aa_ind = self.buc_aa2mlindex(t)
                    scores[t] = self.get_probability_value(
                        ml_aa_ind,
                        mol[chn][r]["CA"].pos,
                        self.__sequence_array[ml_aa_ind],
                        self.__corrections,
                        cell,
                        grid,
                        self.__correl,
                        shiftback,
                        ind,
                    )
                # if _np.argmax(scores) == 0:
                #    if scores[0] > 0.:
                #        print(f"{cg}")
                #        print(f"DEBUG: ALA, {scores}")
                #    continue
                # if not _np.all(scores):
                #    mol[chn][r].type = "UNK"
                #    print(f"DEBUG: UNK , {scores}")
                #    continue
                print(f"DEBUG: SCORES : {scores}")
                ires = _np.argmax(scores)
                #if ires in [2, 3, 5, 6, 16, 19]:
                #    mol[chn][r].type = "UNK"
                #    print(f"DEBUG: {ires}, {scores}")
                #else:
                #    mol[chn][r].type = ProteinTools.residue_code_3(ires)
                #    print(f"DEBUG: {ires}, {scores}")
            done[chn] = True

        for i in done:
            if not i:
                print("Did not manage to fully assign sequence from ML probability!")
        print("done")

    def apply_ml_sequence_symmcheck(
        self,
        mol: _Minimol,
        grid: _Grid_sampling,
        aa_instance: _List,
        ntype: int = 20,
        shiftback: bool = False,
    ):
        done = [False] * mol.size()
        for chn in range(0, mol.size()):
            for r in range(0, mol[chn].size()):
                scores = [0.0] * ntype
                pos = mol[chn][r]["CA"].pos
                # if shiftback:
                #    pos = pos - self.__corrections.shiftback
                if mol[chn][r].exists_property("INDEX"):
                    ind = int(mol[chn][r].get_property("INDEX").value)
                    pos = aa_instance[ind]
                    ### cotinue 25feb
                for t in range(0, ntype):
                    ml_aa_ind = self.buc_aa2mlindex(t)
                    scores[t] = self.get_probability_value(
                        ml_aa_ind,
                        pos,  # mol[chn][r]['CA'].pos,
                        self.__sequence_array[ml_aa_ind],
                        self.__corrections,
                        grid,
                        self.__correl,
                        shiftback,
                    )
                #if not _np.all(scores):
                #    mol[chn][r].type = "UNK"
                #    print(f"DEBUG: UNK , {scores}")
                #    continue
                #ires = _np.argmax(scores)
                #if ires in [2, 3, 5, 6, 16, 19]:
                #    mol[chn][r].type = "UNK"
                #    print(f"DEBUG: {ires}, {scores}")
                #else:
                #    mol[chn][r].type = ProteinTools.residue_code_3(ires)
                #    print(f"DEBUG: {ires}, {scores}")
            done[chn] = True

        for i in done:
            if not i:
                print("Did not manage to fully assign sequence from ML probability!")

    # def get_probability_value(ind: int, mapind: _List, seq_array: _np.ndarray):
    #    if ind in [5, 9, 10]:
    #        return seq_array[mapind[0], mapind[1], mapind[2]] / 2.0
    #    else:
    #        return seq_array[mapind[0], mapind[1], mapind[2]]

    def __call__(
        self, mol: _Minimol, grid: _np.ndarray, shiftback: bool = False
    ):  # , llktargets: LLK_TargetList):
        """Run sequence

        Args:
            mol (_Minimol): _description_
            xmap (_Xmap): _description_
            llkcls (_LLK_TargetList): _description_
            seq (_MMolSeq): _description_
        """
        self.__mol = mol
        print(f"grid {self.__corrections.grid}")
        print(f"asu grid {self.__corrections.grid_asu}")
        print(f"function in grid {grid}")
        grid_samp = _Grid_sampling(grid[0], grid[1], grid[2])
        if self.__ncpu > 1:
            processes = []
            self.__done = [False] * self.__mol.size()
            chunk_size = self.__mol.size() // self.__ncpu
            if chunk_size == 0:
                chunk_size = 1
                if self.__mol.size() == 1:
                    self.__ncpu = 1

            for i in range(0, self.__ncpu):
                start_chn = i * chunk_size
                end_chn = start_chn + chunk_size
                if i == self.__ncpu - 1:
                    end_chn = mol.size()
                p = _Process(
                    target=self.prepare_scores,
                    args=(start_chn, end_chn, 20, grid_samp, shiftback),
                )
                processes.append(p)
                p.start()
            # wait
            for p in processes:
                p.join()
        else:
            self.__done = [False] * self.__mol.size()
            self.prepare_scores(0, self.__mol.size(), 20, grid_samp, shiftback)
            for i in self.__done:
                if not i:
                    print("Did not manage to fully assign sequence probability!")

        ## extract the necessary bits of likelihood targets
        # llksample = [_LLKtgt.Sampled] * len(llkcls)
        # for t in range(0, len(llkcls)):
        #    llksample[t] = llkcls[t].sampled()

    def prepare_score(
        self,
        res: _MRes,
        llktargets_size: int,
        grid: _Grid_sampling,
        shiftback: bool = False,
    ):  # LLK_TargetList):
        cached = False
        # should this method be called every cycle or just at the start?
        ca = Ca_group(res)
        cell = _Cell(self.__corrections.cell)
        if not ca.is_null():
            if res.exists_property("SEQPROB"):
                seqprob_val = res.get_property("SEQPROB").value
                if (
                    (ca.coord_n - seqprob_val.ca.coord_n).lengthsq()
                    and (ca.coord_ca - seqprob_val.ca.coord_n).lengthsq()
                    and (ca.coord_c - seqprob_val.ca.coord_c).lengthsq()
                ):
                    cached = True
        if not cached:
            if res.exists_property("SEQPROB"):
                res.delete_property("SEQPROB")
            ntyp = llktargets_size
            scores = [0.0] * ntyp
            for t in range(0, ntyp):
                ml_aa_ind = self.buc_aa2mlindex(t)
                scores[t] = self.get_probability_value(
                    ml_aa_ind,
                    ca.coord_ca,
                    self.__sequence_array[ml_aa_ind],
                    self.__corrections,
                    cell,
                    grid,
                    self.__correl,
                    shiftback,
                    # self.__fix_axis_positions,
                )
            #self.print_debug(res, scores)
            seqprob_val = Ca_sequence.Sequence_data(ca, scores)
            res.set_property("SEQPROB", _Property_sequence_data(seqprob_val))

    def prepare_scores(
        self,
        ichain_start: int,
        ichain_end: int,
        llktargets_size: int,
        grid: _Grid_sampling,
        shiftback: bool = False,  # LLK_TargetList
    ):
        for chn in range(ichain_start, ichain_end):
            for r in range(0, self.__mol[chn].size()):
                self.prepare_score(self.__mol[chn][r], llktargets_size, grid, shiftback)
            self.__done[chn] = True

    def write_out_model_with_seq(self, mol: _Minimol = None, outfile: str = "NONE"):
        tmpmol = self.__mol.copy()
        if mol is not None:
            print("use mol!\n")
            tmpmol = mol.copy()

        for chn in tmpmol:
            for res in chn:
                if res.exists_property("SEQPROB"):
                    seqprob = res.get_property("SEQPROB").value
                    seqindex = _np.argmax(seqprob.data)
                    res.type = str(ProteinTools.residue_code_3(seqindex))
        if outfile != "NONE":
            _write_structure(tmpmol, outfile, True)
        else:
            _write_structure(tmpmol, "model_with_seqprob.pdb", True)

    def check_is_seqprob_set(self):
        is_set = True
        for chn in self.__mol:
            for res in chn:
                if not res.exists_property("SEQPROB"):
                    is_set = False
        return is_set

    def print_debug(self, res: _MRes, scores: _List):
        print(f"DEBUG: {res}, {scores}")

    def sequence_array_to_map(self, gmap):
        if (gmap.grid.axis_order != _gemmi.AxisOrder.XYZ):
            gmap.setup(float("nan"), _gemmi.MapSetup.ReorderOnly)
        for i in range(0, self.__sequence_array.shape[0]-1):
            a = self.__sequence_array[i]
            gmap.grid = _gemmi.FloatGrid(_np.array(a, dtype=_np.float32))
            gmap.grid.set_unit_cell(gmap.grid.unit_cell)
            gmap.grid.spacegroup = gmap.grid.spacegroup
            if i in [5, 9, 10]:
                name = str(CH_TO_AA_DICT[i]).replace("/","")
                gmap.write_ccp4_map(f"channel_{str(name)}_pred.ccp4")
            else:
                gmap.write_ccp4_map(f"channel_{str(CH_TO_AA_DICT[i])}_pred.ccp4")
'''
#
## split into separate chains
# _PT.split_chains_at_gap(mol)
#
## score residues
# for chn in mol:
#    _Caseq.prepare_scores(chn, xmap, llksample)
#
## sequence
# seqnc = _Caseq.Sequence_threaded(mol, seq, self.__reliability)
# seqnc(self.__ncpu)
# mol = seqnc.result()
# history = seqnc.history()
#
## break chains where sequence is broken
# _PT.split_chains_at_unk(mol, xmap)
#
## count sequenced residues
# num_seq = 0
# for chn in mol:
#    if mol.size() > 5:
#        for res in chn.size():
#            if res.type != "UNK":
#                num_seq += 1
#
# 16Jan2025
#del annotations
