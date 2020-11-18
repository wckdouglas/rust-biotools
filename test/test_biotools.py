import os
import biotools_lib
DataDir = os.path.dirname(os.path.abspath(__file__))


def test_count():
    data = DataDir + '/data/test.fq'
    rc, bc = biotools_lib.fq_stat(data)
    assert(rc == 748)
    assert(bc == 47139)


def test_kmer():
    kmer_dict = biotools_lib.kmer_counter('ACTGACTG', 3)
    for k,v in {'GAC': 1, 'CTG': 2, 'ACT': 2, 'TGA': 1}.items():
        assert(kmer_dict[k] == v)


def test_bed():
    data = DataDir + '/data/test.bed'
    with open(data, 'r') as f:
        bedline = next(f)
        r = biotools_lib.BedRecord(bedline.strip())
        assert(r.start == 3192856)
        assert(r.start == "chr1")
        assert(r.coordinate == 'chr1:3192856-3192888')