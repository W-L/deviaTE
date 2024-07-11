from pathlib import Path

import pytest

from deviaTE.paf import Paf
from deviaTE.analyse import CoverageConverter
from deviaTE.mapping import Mapper
from deviaTE.utils import readfq, read_fa


@pytest.fixture(scope="module")
def testfq():
    return "../data/jockey_dmel.fastq"


@pytest.fixture(scope="module")
def testfq_gz():
    return "../data/jockey_dmel.fastq.gz"


@pytest.fixture(scope="module")
def testpaf():
    return "jockey_dmel.fastq.paf"


@pytest.fixture(scope="module")
def testref():
    return "../data/transposon_sequence_set_v10.2.fa"


@pytest.fixture(scope="module")
def testgff():
    return "../data/transposon_sequence_set_v10.2.gff"


@pytest.fixture(scope="module")
def testfam():
    return "FBte0000088"


@pytest.fixture(scope="module")
def testresults():
    return "../data/jockey_dmel.fastq.FBte0000088.deviate"


@pytest.fixture(scope="module")
def default_mapping(testfq, testref, testpaf):
    m = Mapper(ref=testref, preset="sr")
    outfile = m.map_file(seq_file=testfq)
    assert outfile == testpaf
    return outfile



@pytest.fixture(scope="module")
def paf_dict(default_mapping, testpaf):
    assert Path(default_mapping).is_file()
    paf_dict = Paf.parse_PAF(paf_file=testpaf)
    return paf_dict


@pytest.fixture(scope="module")
def seqs(testfq):
    seqs = {}
    with open(testfq, 'r') as fq:
        for desc, name, seq, qual in readfq(fq):
            seqs[name] = seq
    return seqs


@pytest.fixture(scope="module")
def quals(testfq):
    quals = {}
    with open(testfq, 'r') as fq:
        for desc, name, seq, qual in readfq(fq):
            quals[name] = qual
    return quals


@pytest.fixture(scope="module")
def increments(paf_dict, seqs, quals):
    cc = CoverageConverter()
    incr = cc.convert_records(paf_dict=paf_dict, seqs=seqs, quals=quals)
    return incr


