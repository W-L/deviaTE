from pathlib import Path

import pytest

from deviaTE.mapping import Mapper




def test_mapping(testfq, testref):
    m = Mapper(ref=testref, preset="sr")
    outfile = m.map_file(seq_file=testfq)
    assert outfile == f"{Path(testfq).name}.paf"
    assert Path(outfile).stat().st_size


def test_mapping_gz(testfq_gz, testref):
    m = Mapper(ref=testref, preset="sr")
    outfile = m.map_file(seq_file=testfq_gz)
    assert outfile == f"{Path(testfq_gz).name}.paf"
    assert Path(outfile).stat().st_size


def test_mapping_mmi(testfq):
    m = Mapper(ref="../data/transposon_sequence_set_v10.2.mmi", preset="sr")
    outfile = m.map_file(seq_file=testfq)
    assert outfile == f"{Path(testfq).name}.paf"
    assert Path(outfile).stat().st_size


def test_mapping_mmi_long(testfq):
    m = Mapper(ref="../data/transposon_sequence_set_v10.2.mmi", preset="map-ont")
    outfile = m.map_file(seq_file=testfq)
    assert outfile == f"{Path(testfq).name}.paf"
    assert Path(outfile).stat().st_size


@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_mapping_notfile():
    m = Mapper(ref="../data/transposon_sequence_set_v10.2.mmi", preset="map-ont")
    _ = m.map_file(seq_file="non-existing.file")


@pytest.mark.xfail(raises=FileNotFoundError, strict=True)
def test_mapping_noref():
    _ = Mapper(ref="non-existing.file", preset="map-ont")


@pytest.mark.xfail(raises=ValueError, strict=True)
def test_mapping_wrongreftype():
    _ = Mapper(ref="../data/transposon_sequence_set_v10.2.gff", preset="map-ont")






