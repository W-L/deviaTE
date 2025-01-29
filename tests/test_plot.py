from pathlib import Path
import os

from deviaTE.plot import visualise


gff_mod = "../data/transposon_sequence_set_v10.2_nojockey.gff"



def test_vis(testresults, testfam):
    visualise(input_file=testresults)
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate.pdf").is_file()
    assert os.path.getsize(f"jockey_dmel.fastq.{testfam}.deviate.pdf") == 138505  # size without annotations


def test_vis_anno(testresults, testfam, testgff):
    visualise(input_file=testresults, annotations=testgff)
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate.pdf").is_file()
    assert os.path.getsize(f"jockey_dmel.fastq.{testfam}.deviate.pdf") == 139729  # with annotations



def test_vis_noanno(testresults, testfam):
    visualise(input_file=testresults, annotations=gff_mod)
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate.pdf").is_file()
    assert os.path.getsize(f"jockey_dmel.fastq.{testfam}.deviate.pdf") == 138505  # size without annotations


def test_vis_anno_notexist(testresults, testfam):
    visualise(input_file=testresults, annotations="invalid")
    assert Path(f"jockey_dmel.fastq.{testfam}.deviate.pdf").is_file()
    assert os.path.getsize(f"jockey_dmel.fastq.{testfam}.deviate.pdf") == 138505  # size without annotations

