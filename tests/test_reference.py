from types import SimpleNamespace
from pathlib import Path
import sys

import numpy as np
import pytest

import deviaTE.reference
import deviaTE.config
import deviaTE.mapping


testfq = "../data/jockey_dmel.fastq"
testfq_gz = "../data/jockey_dmel.fastq.gz"
testpaf = "../data/jockey_dmel.fastq.paf"
testfam = "FBte0000088"


args_short_scg = SimpleNamespace(
    input=testfq,
    families=[testfam],
    single_copy_genes=['Dmel_rpl32', 'Dmel_Act5C', 'Dmel_RpII140', 'Dmel_piwi', 'Dmel_p53']
)


args_short_rpm = SimpleNamespace(
    input=testfq,
    families=[testfam],
    rpm=True
)

args_short = SimpleNamespace(
    input=testfq,
    families=[testfam],
)


args_short_gz = SimpleNamespace(
    input=testfq_gz,
    families=[testfam],
)


@pytest.fixture
def scg_conf():
    conf = deviaTE.config.Config(args_debug=args_short_scg)
    return conf


@pytest.fixture
def rpm_conf():
    conf = deviaTE.config.Config(args_debug=args_short_rpm)
    return conf


@pytest.fixture
def short_conf():
    conf = deviaTE.config.Config(args_debug=args_short)
    return conf


@pytest.fixture
def short_conf_gz():
    conf = deviaTE.config.Config(args_debug=args_short_gz)
    return conf


@pytest.fixture
def infile_converted_scg(scg_conf):
    infile = deviaTE.reference.InputFile(conf=scg_conf, infile=Path(testfq))
    infile.analyse_coverage()
    return infile


@pytest.fixture
def infile_converted_rpm(rpm_conf):
    infile = deviaTE.reference.InputFile(conf=rpm_conf, infile=Path(testfq))
    infile.analyse_coverage()
    return infile


@pytest.fixture
def infile_converted(short_conf):
    infile = deviaTE.reference.InputFile(conf=short_conf, infile=Path(testfq))
    infile.analyse_coverage()
    return infile


@pytest.fixture
def infile_converted_gz(short_conf_gz):
    infile = deviaTE.reference.InputFile(conf=short_conf_gz, infile=Path(testfq_gz))
    infile.analyse_coverage()
    return infile



def test_scg_normfac(scg_conf, increments):
    sequence_lib = scg_conf.sequences
    scg = deviaTE.reference.single_gene_norm_fac(
        scgs=args_short_scg.single_copy_genes,
        sequence_lib=sequence_lib,
        increments=increments)
    assert np.allclose(scg, 0.0059049079754601224)



def test_inputfile_analyse_coverage(infile_converted_scg):
    infile = infile_converted_scg
    assert len(infile.seqs) == 4815
    assert len(infile.quals) == 4815
    assert isinstance(infile.mapper, deviaTE.mapping.Mapper)
    assert len(infile.incr[testfam]) == 4747
    assert sys.getsizeof(infile.incr[testfam]) == 41880
    assert np.allclose(infile.scg_norm_fac, 0.0059049079754601224)


def test_inputfile_analyse_coverage_gz(infile_converted_gz):
    infile = infile_converted_gz
    assert len(infile.seqs) == 4815
    assert len(infile.quals) == 4815
    assert isinstance(infile.mapper, deviaTE.mapping.Mapper)
    assert len(infile.incr[testfam]) == 4747
    assert sys.getsizeof(infile.incr[testfam]) == 41880


def test_fams_scg(infile_converted_scg):
    infile_converted_scg.analyse_families()
    assert len(infile_converted_scg.results_files) == 1
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate").is_file()
    with open(f"jockey_dmel.fastq.{testfam}.deviate", 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 169.35064935064935)
    # grab estimated insertions
    jj = lines[0]
    ins, inshq = float(jj.split(' ')[2]), float(jj.split(' ')[4])
    assert np.allclose(ins, 33773.88337558856)
    assert np.allclose(inshq, 29360.645935737568)


def test_fams_rpm(infile_converted_rpm):
    infile_converted_rpm.analyse_families()
    assert len(infile_converted_rpm.results_files) == 1
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate").is_file()
    with open(f"jockey_dmel.fastq.{testfam}.deviate", 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 207.68431983385256)


def test_viz(infile_converted):
    infile_converted.analyse_families()
    infile_converted.visualise()
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate.pdf").is_file()


def test_viz_gz(infile_converted_gz):
    infile_converted_gz.analyse_families()
    infile_converted_gz.visualise()
    assert Path(f"jockey_dmel.fastq.{testfam}.gz.deviate.pdf").is_file()





