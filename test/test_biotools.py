import os
import biotools_lib
DataDir = os.path.dirname(os.path.abspath(__file__))


def test_count():
    data = DataDir + '/data/test.fq'
    rc, bc = biotools_lib.readfq(data)
    assert(rc == 748)

test_count()
