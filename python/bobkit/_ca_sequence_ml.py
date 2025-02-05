from __future__ import annotations
from typing import List as _List
from bobkit.buccaneer import *

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
    "Ca_sequence_ml",
]

# from _bobkit._buccaneer import (
#    Ca_sequence as _Caseq,
#    Ca_group as _Cagroup,
#    LLK_map_target as _LLKtgt,
#    LLK_TargetList as _LLK_TargetList,
#    ProteinTools as _PT,
# )
from bobkit.clipper import (
    Cell as _Cell,
    Grid_sampling as _Grid_sampling,
    Coord_orth as _Coord_orth,
    MiniMol as _Minimol,
    MChain as _MChain,
    MResidue as _MRes,
    Xmap_float as _Xmap,
    MMoleculeSequence as _MMolSeq,
    MChain as _MChain,
    Thread_base as _Thread_base,
    Property_sequence_data as _Property_sequence_data,
)
import numpy as _np
import gemmi as _gemmi
from multiprocessing import Process as _Process, Manager as _Manager


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
    0: 0,
    1: 7,
    2: 9,
    3: 10,
    4: 14,
    5: [16, 19],
    6: 13,
    7: 17,
    8: 18,
    9: [2, 3],
    10: [5, 6],
    11: 1,
    12: 8,
    13: 11,
    14: 15,
    15: 4,
    16: 12,
}

BUC_TO_MLI_DICT = {
    0: 0,
    1: 11,
    2: 9,
    3: 9,
    4: 15,
    5: 14,
    6: 10,
    7: 1,
    8: 12,
    9: 2,
    10: 3,
    11: 13,
    12: 16,
    13: 6,
    14: 4,
    15: 14,
    16: 5,
    17: 7,
    18: 8,
    19: 5,
    12: 16,
}

ProteinTools()
# const char ProteinTools::rtype3[21][4] =
# {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
#  "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
#  "MSE"};
# const int ProteinTools::tindex[21] =
# {    0,    1,    2,    3,    4,    5,    6,    7,    8,    9,
#     10,   11,   12,   13,   14,   15,   16,   17,   18,   19,
#     12};


class Cathread(_Thread_base):
    """Derive class

    Args:
        _Thread_base (_type_): _description_
    """

    def __init__(
        self,
        data: _List[int] = None,
        chain: _MChain = None,
        xmap: _Xmap = None,
        sequence_array: _np.ndarray = None,
        llksample: _List[LLK_map_target.Sampled] = None,
    ):
        super().__init__()
        self._chain = chain
        self._xmap = xmap
        self._seq_array = sequence_array
        self._llksample = llksample
        self.count = 0
        self.sum = 0
        self.data = data.copy()

    def __call__(self, nthread=0):
        thr = True if nthread > 0 else False
        if thr:
            threads = [Cathread] * (nthread - 1)
            self.run()
            for i in range(0, len(threads)):
                threads[i].run()
            self.join()
            for i in range(0, len(threads)):
                threads[i].join()
            # if self.count >= len(threads):
            #    for i in range(0, len(threads)):
            #        self.merge(threads[i])
            # else:
            #    thr = False
        else:
            for i in range(0, len(self.data)):
                self.sum += self.data[i]

    def sequence_score(self):
        pass

    def result(self):
        print(f"sum = {self.sum}")

    # def merge(self, other: Cathread = None):
    #    for i in range(0, len(self.data)):
    #
    #    pass

    def Run(self):
        self.sum = 0
        while 1:
            self.lock()
            c = self.count + 1
            self.unlock()
            if c >= 100:
                break
            self.sum += self.data[c]

        # print("hello")


class Ca_sequence_ml:
    """Class for sequence Ca chains using density and output predictions from machine learning"""

    def __init__(
        self,
        sequence_array: _np.ndarray,
        apix: _List,
        reliability: float = 1.0,
        correl: bool = False,
    ):
        self.__sequence_array = sequence_array
        self.__reliability = reliability
        self.__ncpu = 1
        self.__apix = apix
        self.__mol = None
        self.__done = []
        self.__correl = correl

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
    def ml_aa2bucindex(ml_aa_ind: int):
        return MLI_TO_BUC_DICT[ml_aa_ind]

    @staticmethod
    def buc_aa2mlindex(buc_aa_ind: int):
        return BUC_TO_MLI_DICT[buc_aa_ind]

    @staticmethod
    def get_probability_value(
        ind: int,
        pos: _Coord_orth,
        seq_array: _np.ndarray,
        corrections: _List,
        correl: bool = False,
    ):
        llkval = -1.0
        x = (
            int(pos.x // corrections[0]) + int(corrections[3] + corrections[6])
        ) % seq_array.shape[0]
        y = (
            int(pos.y // corrections[1]) + int(corrections[4] + corrections[6])
        ) % seq_array.shape[1]
        z = (
            int(pos.z // corrections[2]) + int(corrections[5] + corrections[6])
        ) % seq_array.shape[2]

        # mapind = [int(pos.x), int(pos.y), int(pos.z)]
        probval = seq_array[x, y, z]
        if correl:
            if ind in [5, 9, 10]:
                return probval / 2.0
            else:
                return probval
        else:
            if probval > 0.0:
                if ind in [5, 9, 10]:
                    llkval = _np.log10(seq_array[x, y, z] / 2.0)
                else:
                    llkval = _np.log10(seq_array[x, y, z])
            elif probval == 0.0:
                llkval = -1.0
            return -llkval

    # def get_probability_value(ind: int, mapind: _List, seq_array: _np.ndarray):
    #    if ind in [5, 9, 10]:
    #        return seq_array[mapind[0], mapind[1], mapind[2]] / 2.0
    #    else:
    #        return seq_array[mapind[0], mapind[1], mapind[2]]

    def __call__(self, mol: _Minimol):  # , llktargets: LLK_TargetList):
        """Run sequence

        Args:
            mol (_Minimol): _description_
            xmap (_Xmap): _description_
            llkcls (_LLK_TargetList): _description_
            seq (_MMolSeq): _description_
        """
        self.__mol = mol
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
                p = _Process(target=self.prepare_scores, args=(start_chn, end_chn, 20))
                processes.append(p)
                p.start()
            # wait
            for p in processes:
                p.join()
        else:
            self.__done = [False] * self.__mol.size()
            self.prepare_scores(0, self.__mol.size(), 20)
            for i in self.__done:
                if not i:
                    print("Did not manage to fully assign sequence probability!")

        ## extract the necessary bits of likelihood targets
        # llksample = [_LLKtgt.Sampled] * len(llkcls)
        # for t in range(0, len(llkcls)):
        #    llksample[t] = llkcls[t].sampled()

    def prepare_score(self, res: _MRes, llktargets_size: int):  # LLK_TargetList):
        cached = False
        # should this method be called every cycle or just at the start?
        ca = Ca_group(res)
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
                    self.__apix,
                    self.__correl,
                )
            seqprob_val = Ca_sequence.Sequence_data(ca, scores)
            res.set_property("SEQPROB", _Property_sequence_data(seqprob_val))

    def prepare_scores(
        self, ichain_start: int, ichain_end: int, llktargets_size: int  # LLK_TargetList
    ):
        for chn in range(ichain_start, ichain_end):
            for r in range(0, self.__mol[chn].size()):
                self.prepare_score(self.__mol[chn][r], llktargets_size)
            self.__done[chn] = True

    def check_is_seqprob_set(self):
        is_set = True
        for chn in self.__mol:
            for res in chn:
                if not res.exists_property("SEQPROB"):
                    is_set = False
        return is_set


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
del annotations
