from pathlib import Path
from types import SimpleNamespace

import pytest

from deviaTE.paf import Paf
from deviaTE.analyse import CoverageConverter
from deviaTE.mapping import Mapper
from deviaTE.config import Config
from deviaTE.reference import InputFile


@pytest.fixture(scope="module")
def testfq():
    return "../data/jockey_dmel.fastq"


@pytest.fixture(scope="module")
def testfq_gz():
    return "../data/jockey_dmel.fastq.gz"


@pytest.fixture(scope="module")
def testfa_nanopore_short():
    return "../data/SRR28208594_1_FBte0000400_270.fa"


@pytest.fixture(scope="module")
def testfq_nanopore_long():
    return "../data/SRR28208594_1_FBte0001206_285.fq"


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
def testfam_short():
    return "FBte0000400"


@pytest.fixture(scope="module")
def testfam_long():
    return "FBte0001206"


@pytest.fixture(scope="module")
def testresults():
    return "../data/jockey_dmel.fastq.FBte0000088.deviate"


##################################################################################


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
def increments(paf_dict):
    cc = CoverageConverter()
    incr = cc.convert_records(paf_dict=paf_dict)
    return incr


##################################################################################

@pytest.fixture(scope="module")
def args_short_scg(testfq, testfam):
    return SimpleNamespace(
        input=testfq,
        families=[testfam],
        single_copy_genes=['Dmel_rpl32', 'Dmel_Act5C', 'Dmel_RpII140', 'Dmel_piwi', 'Dmel_p53']
    )


@pytest.fixture(scope="module")
def args_short_rpm(testfq, testfam):
    return SimpleNamespace(
        input=testfq,
        families=[testfam],
        rpm=True
    )


@pytest.fixture(scope="module")
def args_short(testfq, testfam):
    return SimpleNamespace(
        input=testfq,
        families=[testfam],
    )


@pytest.fixture(scope="module")
def args_short_gz(testfq_gz, testfam):
    return SimpleNamespace(
        input=testfq_gz,
        families=[testfam],
    )


@pytest.fixture(scope="module")
def args_short_tar(testfq, testfam):
    return SimpleNamespace(
        input=testfq,
        families=[testfam],
        tar=True
    )


@pytest.fixture(scope="module")
def args_nanopore_short(testfa_nanopore_short, testfam_short):
    return SimpleNamespace(
        input=testfa_nanopore_short,
        families=[testfam_short],
        preset="map-ont"
    )


@pytest.fixture(scope="module")
def args_nanopore_long(testfq_nanopore_long, testfam_long):
    return SimpleNamespace(
        input=testfq_nanopore_long,
        families=[testfam_long],
        preset="map-ont"
    )


##################################################################################


@pytest.fixture(scope="module")
def scg_conf(args_short_scg):
    conf = Config(args_debug=args_short_scg)
    return conf


@pytest.fixture(scope="module")
def rpm_conf(args_short_rpm):
    conf = Config(args_debug=args_short_rpm)
    return conf


@pytest.fixture(scope="module")
def short_conf(args_short):
    conf = Config(args_debug=args_short)
    return conf


@pytest.fixture(scope="module")
def short_conf_tar(args_short_tar):
    conf = Config(args_debug=args_short_tar)
    return conf


@pytest.fixture(scope="module")
def short_conf_gz(args_short_gz):
    conf = Config(args_debug=args_short_gz)
    return conf


@pytest.fixture(scope="module")
def nanopore_conf_short(args_nanopore_short):
    conf = Config(args_debug=args_nanopore_short)
    return conf


@pytest.fixture(scope="module")
def nanopore_conf_long(args_nanopore_long):
    conf = Config(args_debug=args_nanopore_long)
    return conf


##################################################################################


@pytest.fixture(scope="module")
def infile_converted_scg(scg_conf, testfq):
    infile = InputFile(conf=scg_conf, infile=Path(testfq))
    infile.analyse_coverage()
    return infile


@pytest.fixture(scope="module")
def infile_converted_rpm(rpm_conf, testfq):
    infile = InputFile(conf=rpm_conf, infile=Path(testfq))
    infile.analyse_coverage()
    return infile


@pytest.fixture(scope="module")
def infile_converted_rpm_gz(rpm_conf, testfq_gz):
    infile = InputFile(conf=rpm_conf, infile=Path(testfq_gz))
    infile.analyse_coverage()
    return infile


@pytest.fixture(scope="module")
def infile_converted(short_conf, testfq):
    infile = InputFile(conf=short_conf, infile=Path(testfq))
    infile.analyse_coverage(force_map=True)
    return infile


@pytest.fixture(scope="module")
def infile_converted_gz(short_conf_gz, testfq_gz):
    infile = InputFile(conf=short_conf_gz, infile=Path(testfq_gz))
    infile.analyse_coverage()
    return infile


@pytest.fixture(scope="module")
def infile_converted_tar(short_conf_tar, testfq_gz):
    infile = InputFile(conf=short_conf_tar, infile=Path(testfq_gz))
    infile.analyse_coverage()
    return infile


@pytest.fixture(scope="module")
def infile_converted_nanopore_short(nanopore_conf_short, testfa_nanopore_short):
    infile = InputFile(conf=nanopore_conf_short, infile=Path(testfa_nanopore_short))
    infile.analyse_coverage(force_map=True)
    return infile


@pytest.fixture(scope="module")
def infile_converted_nanopore_long(nanopore_conf_long, testfq_nanopore_long):
    infile = InputFile(conf=nanopore_conf_long, infile=Path(testfq_nanopore_long))
    infile.analyse_coverage(force_map=True)
    return infile


