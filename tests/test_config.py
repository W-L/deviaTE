import logging
from types import SimpleNamespace

import pytest

import deviaTE.config



testfq = "../data/jockey_dmel.fastq"
testpaf = "jockey_dmel.fastq.paf"
testref = "../data/transposon_sequence_set_v10.2.fa"
testfam = "FBte0000088"


args_default = SimpleNamespace()

args_nonsense_input = SimpleNamespace(
    input="deviate.log"
)

args_non_input = SimpleNamespace(
    input="deviate.lo"
)


args_short = SimpleNamespace(
    input=testfq,
    families=[testfam],
)


args_short_scg = SimpleNamespace(
    input=testfq,
    families=[testfam],
    single_copy_genes=['Dmel_rpl32', 'Dmel_Act5C', 'Dmel_RpII140', 'Dmel_piwi', 'Dmel_p53']
)


args_short_scg_bad = SimpleNamespace(
    input=testfq,
    families=[testfam],
    single_copy_genes=['Dmel_rpl32', 'Dmel_Act5C', 'Dmel_RpII140', 'Dmel_piwi', 'Dmel_p5']
)


args_short_dir = SimpleNamespace(
    input="../data",
    families=[testfam],
)


args_short_dir_badfam = SimpleNamespace(
    input="../data",
    families=['FBte0'],
)

args_no_dir = SimpleNamespace(
    input="__pycache__",
)

args_badlib_notexist = SimpleNamespace(
    input=testfq,
    library="badlib.fa",
)

args_badlib = SimpleNamespace(
    input=testfq,
    library="../data/badlib.fa",
)


#
# args_debug = SimpleNamespace(
#     input="../data/short_reads/jockey_dmel.fastq",
#     library=None,
#     families=['FBte0000088'],
#     annotation=None,
#     min_align_len=1,
#     preset="sr",
#     single_copy_genes=['Dmel_rpl32', 'Dmel_piwi'],
#     rpm=None)


@pytest.mark.xfail(raises=ValueError, strict=True)
def test_no_input():
    _ = deviaTE.config.Config(args_debug=args_default)


@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_nonsense_input():
    _ = deviaTE.config.Config(args_debug=args_nonsense_input)

@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_non_input():
    _ = deviaTE.config.Config(args_debug=args_non_input)


@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_no_dir():
    _ = deviaTE.config.Config(args_debug=args_no_dir)


def test_short_single():
    conf = deviaTE.config.Config(args_debug=args_short)
    for k, v in conf.__dict__.items():
        logging.info(f'{k}')
    assert conf

def test_short_dir():
    conf = deviaTE.config.Config(args_debug=args_short_dir)
    for k, v in conf.__dict__.items():
        logging.info(f'{k}')
    assert conf

@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_badlib_notexist():
    _ = deviaTE.config.Config(args_debug=args_badlib_notexist)

@pytest.mark.xfail(raises=ValueError, strict=True)
def test_badfam():
    _ = deviaTE.config.Config(args_debug=args_short_dir_badfam)

def test_badlib():
    conf = deviaTE.config.Config(args_debug=args_badlib)
    assert conf

def test_scg():
    conf = deviaTE.config.Config(args_debug=args_short_scg)
    assert conf

@pytest.mark.xfail(raises=ValueError, strict=True)
def test_scg_bad():
    _ = deviaTE.config.Config(args_debug=args_short_scg_bad)

