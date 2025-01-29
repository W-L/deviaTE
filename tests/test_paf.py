from io import StringIO
import sys

import pytest

from deviaTE.paf import Paf




def test_paf_init():
    _ = Paf()


def test_parsing(testpaf):
    paf_dict = Paf.parse_PAF(paf_file=testpaf)
    assert paf_dict
    assert len(paf_dict) == 4621
    assert sys.getsizeof(paf_dict) == 103864


def test_parsing_full(testpaf):
    paf_dict = Paf.parse_PAF(paf_file=testpaf)
    paf_dict = Paf.single_rec(paf_dict)
    assert paf_dict
    assert len(paf_dict) == 4621
    assert sys.getsizeof(paf_dict) == 103864


def test_parsing_stringio(testpaf):
    with open(testpaf, 'r') as tp:
        p = ''.join(tp.readlines())
    paf_dict = Paf.parse_PAF(paf_file=StringIO(p))
    assert paf_dict
    assert len(paf_dict) == 4621


def test_parsing_minlen(testpaf):
    paf_dict = Paf.parse_PAF(paf_file=testpaf, min_len=100)
    assert paf_dict
    assert len(paf_dict) == 3275


@pytest.mark.xfail(raises=ValueError, strict=True)
def test_parsing_invalid():
    _ = Paf.parse_PAF(paf_file="invalid_input")





