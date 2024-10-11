import os
import pytest
import sys

top_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
sys.path.append(os.path.join(top_path, "scripts"))
# sys.path.append("../scripts")
import buccaneer


@pytest.fixture
def make_fasta_file():
    sequence_str = (
        ">2BN3_1|Chain A|INSULIN|BOS TAURUS (9913)\n"
        "GIVEQCCTSVCSLYQLENYCN\n"
        ">2BN3_2|Chain B|INSULIN|BOS TAURUS (9913)\n"
        "FVNQHLCGSHLVEALYLVCGERGFFYTPKA\n"
    )
    seqfile_name = "2bn3_seq.fasta"
    with open("2bn3_seq.fasta", "w") as seq_file:
        seq_file.write(sequence_str)
    return seqfile_name


@pytest.fixture
def get_ref1tqw():
    pdb_input = os.path.join(
        top_path, "include/reference_data/reference-1tqw.pdb"
    )  # noqa E501
    mtz_input = os.path.join(
        top_path, "include/reference_data/reference-1tqw.mtz.gz"
    )  # noqa E501
    return [pdb_input, mtz_input]


class TestBuccaneer:
    # def TestInsulin(self, make_fasta_file, get_ref1tqw):
    def test_build_insulin(self, make_fasta_file, get_ref1tqw):
        # mtzin = os.path.join(os.environ["CCP4"], "examples/data/insulin.mtz")
        mtzin = "/opt/xtal/ccp4-8.0/examples/data/insulin.mtz"
        seqin = make_fasta_file
        # url = "https://www.rcsb.org/fasta/entry/2BN3"
        # urllib.request.urlretrieve(url, seqin)
        # pdb_ref = get_ref1tqw[0]
        # mtz_ref = get_ref1tqw[1]
        pdb_ref = os.path.join(
            os.environ["CLIBD"], "reference_structures/reference-1tqw.pdb"
        )
        mtz_ref = os.path.join(
            os.environ["CLIBD"], "reference_structures/reference-1tqw.mtz"
        )
        # ipcolin_fsigf = "FP,SIGF"
        ipcolin_fsigf = "F,SIGF"
        ipcolin_freer = "FreeR_flag"
        # input_args = ["pybuccaneer"]
        # input_args = [""]
        input_args = ["--mtzin", mtzin]
        input_args += ["--seqin", seqin]
        input_args += ["--resolution", "1.7"]
        input_args += ["--cycles", "1"]
        input_args += ["--colin-fo", ipcolin_fsigf]
        input_args += ["--colin-free", ipcolin_freer]
        input_args += ["--pdbin-ref", pdb_ref]
        input_args += ["--mtzin-ref", mtz_ref]

        buccaneer.main(input_args)
        # need some assertion tests/checking
        # remove files written
        files_to_remove = [
            seqin,
            "buccaneer_build.cif",
            "buccaneer_build.pdb",
            "summary.xml",
        ]
        for fname in files_to_remove:
            if os.path.exists(fname):
                os.remove(fname)


# if __name__ == "__main__":
#    TestInsulin()
