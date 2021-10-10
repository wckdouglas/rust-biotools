import os
from urllib.request import Request, urlopen
import biotools_lib
import gzip
import pytest
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


def test_kmer_large_k():
    with pytest.raises(ValueError) as exinfo: 
        biotools_lib.kmer_counter('ACTGACTG', 10)

    assert "k is smaller than sequence length" in str(exinfo)


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
        assert(not r.overlap(45691290,45694369))
        assert(r.exons[0].gstart == 45691195)
        assert(r.exons[-1].gstart == 45691195 + 115101)
        assert(r.exons[0].gend == 45691195 + 94)
        assert(r.exons[-1].gend == (45691195 + 115101 + 2837))
        assert(r.exons[-1].gend == r.end)

        bedline = next(f)
        r = biotools_lib.Bed12Record(bedline.strip())
        assert(r.end == 74484164)
        assert(r.chrom == "chrX")
        assert(r.coordinate == 'chrX:74440264-74484164')
        assert(r.exons[0].gstart == 43858 + 74440264)
        assert(r.exons[-1].gstart == 0 + 74440264 )
        
        assert(r.exons[0].gend == 74440264 + 43858 + 42)
        assert(r.exons[-1].gend == 74440264 + 0 + 1073)
        assert(r.exons[0].gend == r.end)


def test_block():
    data = DataDir + '/data/test.bed12'
    with open(data, 'r') as f:
        bedline = next(f)
        r = biotools_lib.Bed12Record(bedline.strip())
        block_starts, block_ends = r.blocks(10,280)

        assert(len(block_starts) == 3)
        assert(block_starts[0] == r.start + 10 - 1)
        assert(block_starts[1] == r.exons[1].gstart)
        assert(block_starts[2] == r.exons[2].gstart)
        assert(block_ends[0] == r.exons[0].gend)
        assert(block_ends[1] == r.exons[1].gend)
        assert(block_ends[2] == r.exons[2].gstart + 38  - 1)

        bedline = next(f)
        r = biotools_lib.Bed12Record(bedline.strip())
        block_starts, block_ends = r.blocks(10,280)

        assert(len(block_starts) == 4)
        assert(block_ends[-1] == r.end - 10 + 1)
        assert(block_starts[-1] == r.exons[0].gstart )
        assert(block_ends[0] == r.exons[3].gend)
        assert(block_starts[0] == r.exons[3].gend - 73 + 1)
