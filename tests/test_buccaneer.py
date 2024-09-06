import os
import urllib.request
import bobkit
import unittest

import sys

# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append("../scripts")
import buccaneer


def TestInsulin():
    mtzin = os.path.join(os.environ["CCP4"], "examples/data/insulin.mtz")
    seqin = "2bn3_seq.fasta"
    url = f"https://www.rcsb.org/fasta/entry/2BN3"
    urllib.request.urlretrieve(url, seqin)
    pdb_ref = os.path.join(
        os.environment["CLIBD"], "reference_structures/reference-1tqw.pdb"
    )
    mtz_ref = os.path.join(
        os.environment["CLIBD"], "reference_structures/reference-1tqw.mtz"
    )
    ipcolin_fsigf = "FP,SIGF"
    ipcolin_freer = "FreeR_flag"
    input_args = ["pybuccaneer"]
    input_args += ["--mtzin", mtzin]
    input_args += ["--seqin", seqin]
    input_args += [
        "--resolution",
    ]
    input_args += ["--cycles", "1"]
    input_args += ["--colin-fo", ipcolin_fsigf]
    input_args += ["--colin-free", ipcolin_freer]
    input_args += ["--pdbin-ref", pdb_ref]
    input_args += ["--mtzin-ref", mtz_ref]
    input_args += ["--resolution", "1.7"]

    buccaneer.main(input_args)


if __name__ == "__main__":
    TestInsulin()
