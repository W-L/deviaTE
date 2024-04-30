import re
import sys

from deviaTE.analyse import CoverageConverter



def test_init_cc():
    cc = CoverageConverter()
    assert cc
    assert isinstance(cc.base2int, dict)
    assert isinstance(cc.cigar_regex, re.Pattern)


def test_convert_cov(testfam, increments):
    incr = increments
    assert incr
    assert len(incr) == 4
    assert len(incr[testfam]) == 4747
    assert sys.getsizeof(incr[testfam]) == 41880





