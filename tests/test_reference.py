from pathlib import Path
import sys
import tarfile

import numpy as np

import deviaTE.reference
import deviaTE.config
import deviaTE.mapping



def test_scg_normfac(scg_conf, increments, args_short_scg):
    sequence_lib = scg_conf.sequences
    scg = deviaTE.reference.single_gene_norm_fac(
        scgs=args_short_scg.single_copy_genes,
        sequence_lib=sequence_lib,
        increments=increments)
    assert np.allclose(scg, 0.0059049079754601224)



def test_inputfile_analyse_coverage(infile_converted_scg, testfam):
    infile = infile_converted_scg
    assert isinstance(infile.mapper, deviaTE.mapping.Mapper)
    assert len(infile.incr[testfam]) == 4747
    assert sys.getsizeof(infile.incr[testfam]) == 41880
    assert np.allclose(infile.scg_norm_fac, 0.0059049079754601224)


def test_inputfile_analyse_coverage_gz(infile_converted_gz, testfam):
    infile = infile_converted_gz
    assert isinstance(infile.mapper, deviaTE.mapping.Mapper)
    assert len(infile.incr[testfam]) == 4747
    assert sys.getsizeof(infile.incr[testfam]) == 41880


def test_fams_scg(infile_converted_scg, testfam):
    infile_converted_scg.analyse_families()
    assert len(infile_converted_scg.results_files) == 1
    result = f"jockey_dmel.fastq.{testfam}.deviate"
    assert Path(result).is_file()
    with open(result, 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 169.35064935064935)
    # grab estimated insertions
    jj = lines[0]
    ins, inshq = float(jj.split(' ')[2]), float(jj.split(' ')[4])
    assert np.allclose(ins, 33773.88337558856)
    assert np.allclose(inshq, 29360.645935737568)
    Path(result).unlink()


def test_fams_rpm(infile_converted_rpm, testfam):
    infile_converted_rpm.analyse_families()
    assert len(infile_converted_rpm.results_files) == 1
    result = f"jockey_dmel.fastq.{testfam}.deviate"
    assert Path(result).is_file()
    with open(result, 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 207.68431983385256)
    Path(result).unlink()


def test_fams_rpm_gz(infile_converted_rpm_gz, testfam):
    infile_converted_rpm_gz.analyse_families()
    assert len(infile_converted_rpm_gz.results_files) == 1
    result = f"jockey_dmel.fastq.gz.{testfam}.deviate"
    assert Path(result).is_file()
    with open(result, 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 207.68431983385256)
    Path(result).unlink()


def test_fams_tar(infile_converted_tar, testfam):
    infile_converted_tar.analyse_families()
    assert len(infile_converted_tar.results_files) == 1
    tar = Path("jockey_dmel.fastq.gz.deviate.tar")
    assert tar.is_file()
    with tarfile.open(str(tar), "r") as tarf:
        members = tarf.getmembers()
        assert len(members) == 1
        assert members[0].name == f'jockey_dmel.fastq.gz.{testfam}.deviate'
        assert members[0].size == 330285
    tar.unlink()



def test_fams_nanopore_short(infile_converted_nanopore_short, testfam_short):
    infile_converted_nanopore_short.analyse_families()
    assert len(infile_converted_nanopore_short.results_files) == 1
    result_file = f"SRR28208594_1_FBte0000400_270.fq.{testfam_short}.deviate"
    assert Path(result_file).is_file()
    with open(result_file, 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    assert np.allclose(test_num, 85)
    Path(result_file).unlink()


def test_fams_nanopore_long(infile_converted_nanopore_long, testfam_long):
    infile_converted_nanopore_long.analyse_families()
    assert len(infile_converted_nanopore_long.results_files) == 1
    result_file = f"SRR28208594_1_FBte0001206_285.fq.{testfam_long}.deviate"
    assert Path(result_file).is_file()
    with open(result_file, 'r') as outf:
        lines = outf.readlines(1000)
    ll = lines[5]
    test_num = float(ll.split(' ')[4])
    print(test_num)
    assert np.allclose(test_num, 0)
    Path(result_file).unlink()


def test_viz(infile_converted, testfam):
    infile_converted.analyse_families()
    infile_converted.visualise()
    res = f"jockey_dmel.fastq.{testfam}.deviate.pdf"
    assert Path(res).is_file()
    Path(res).unlink()


def test_viz_gz(infile_converted_gz, testfam):
    infile_converted_gz.analyse_families()
    infile_converted_gz.visualise()
    res = f"jockey_dmel.fastq.gz.{testfam}.deviate.pdf"
    assert Path(res).is_file()
    Path(res).unlink()


def test_viz_tar(infile_converted_tar, testfam):
    infile_converted_tar.analyse_families()
    infile_converted_tar.visualise()
    assert not Path(f"jockey_dmel.fastq.gz.{testfam}.deviate.pdf").is_file()
    tar = Path("jockey_dmel.fastq.gz.deviate.visualisations.tar")
    assert tar.is_file()
    with tarfile.open(str(tar), "r") as tarf:
        members = tarf.getmembers()
        assert len(members) == 1
        assert members[0].name == f'jockey_dmel.fastq.gz.{testfam}.deviate.pdf'
        assert members[0].size == 139913
    tar.unlink()


def test_viz_nanopore(infile_converted_nanopore_short, testfam_short):
    infile_converted_nanopore_short.analyse_families()
    infile_converted_nanopore_short.visualise()
    res = f"SRR28208594_1_FBte0000400_270.fq.{testfam_short}.deviate.pdf"
    assert Path(res).is_file()
    Path(res).unlink()




