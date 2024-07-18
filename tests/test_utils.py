import argparse

import pytest
import numpy as np

import deviaTE.utils


def test_seq2int():
    arr = deviaTE.utils.Seq2Int("ACGT").translate()
    assert len(arr) == 4
    assert np.allclose(arr, np.array([0, 1, 2, 3]))


def test_seq2int_invalid_char():
    arr = deviaTE.utils.Seq2Int("ACGTX").translate()
    assert np.allclose(arr, np.array([0, 1, 2, 3, 40]))


def test_logger():
    deviaTE.utils.init_logger(logfile="test.log", args=argparse.Namespace(test="works"))


@pytest.mark.parametrize("seq,rev", [("ACGT", "ACGT"), ("TTTT", "AAAA"), ("C", "G"), ("X", "X")])
def test_revcomp(seq, rev):
    r = deviaTE.utils.reverse_complement(seq)
    assert r == rev


def test_rawcount(testfq):
    lines = deviaTE.utils.rawcount(testfq)
    assert lines == 19261


def test_rapidgz(testfq_gz):
    lines = deviaTE.utils.rapidgzip_count_lines(testfq_gz)
    assert lines == 19261

