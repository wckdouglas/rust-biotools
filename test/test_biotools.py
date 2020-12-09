import os
from urllib.request import Request, urlopen
import biotools_lib
import gzip
DataDir = os.path.dirname(os.path.abspath(__file__))


def test_transcriptome():
    refFlat = DataDir + '/data/test.refFlat'
    url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz'
    if not os.path.isfile(refFlat):
        req = Request(url)
        response = urlopen(req)
        with open(refFlat, 'w') as f:
            for line in gzip\
                    .decompress(response.read())\
                    .decode()\
                    .split('\n')[:100]:
                print(line, file = f)


        



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
        assert(r.chrom == "chr1")
        assert(r.coordinate == 'chr1:3192856-3192888')
        assert(r.overlap(3192877, 3192900))
        assert(not r.overlap(3192889, 3192900))

def test_bed12():
    data = DataDir + '/data/test.bed12'
    with open(data, 'r') as f:
        bedline = next(f)
        r = biotools_lib.Bed12Record(bedline.strip())
        assert(r.start == 45691195)
        assert(r.chrom == "chr8")
        assert(r.coordinate == 'chr8:45691195-45809133')
        assert(r.overlap(45691180,45691200))
