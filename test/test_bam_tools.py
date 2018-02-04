from sequencing_tools.bam_tools import *
from pysam import AlignedRead



def artificial_read(flag=163):
    aln = AlignedRead()
    aln.qname = 'try_this'
    aln.flag = flag
    aln.cigarstring = '16M'
    aln.seq = 'A' * 16
    aln.pos = 1000000 if flag == 163 else 1000020
    aln.reference_id = 1
    return aln


def test_split_cigar():
    assert(split_cigar('63M4S') == [[63,4],['M','S']])


def test_read_ends():
    aln = artificial_read()
    assert(read_ends(aln) == (1000000,1000015))


def test_fragment_ends():
    read1 = artificial_read(flag = 83)
    read2 = artificial_read(flag = 163)
    assert(fragment_ends(read1, read2) == (1000000, 1000036))


def test_concordant_pairs():
    read1 = artificial_read(flag = 83)
    read2 = artificial_read(flag = 163)
    assert(concordant_pairs(read1, read2))

    read1 = artificial_read(flag = 81)
    read2 = artificial_read(flag = 163)
    assert(concordant_pairs(read1, read2)==False)
    

def test_concordant_alignment():
    read1 = artificial_read(flag = 81)
    read2 = artificial_read(flag = 163)
    assert(concordant_alignment(read1)==False)
    assert(concordant_alignment(read2))


def test_make_regions():
    assert(list(make_regions(100,51)) == [(0, 51), (51, 100)])


def test_make_cigar_seq():
    seq = ''
    for c in make_cigar_seq(['3','5','1','3'],['S','M','I','M']):
        seq += c
    assert(seq == 'SSSMMMMMIMMM')


def test_cigar_to_str():
    cs = '3S5M1I3M'
    assert(cigar_to_str(cs) == 'SSSMMMMMIMMM')


def test_get_strand():
    assert(get_strand(artificial_read(flag = 83)) == '-')
    assert(get_strand(artificial_read(flag = 163)) == '-')
    assert(get_strand(artificial_read(flag = 99)) == '+')
    assert(get_strand(artificial_read(flag = 147)) == '+')


def test_check_concordant():
    assert(check_concordant(artificial_read(flag=83), artificial_read(flag=163)))
    assert(check_concordant(artificial_read(flag=99), artificial_read(flag=147)))


def test_check_primary():
    assert(check_primary(artificial_read(flag=83), artificial_read(flag=163)))
    assert(check_primary(artificial_read(flag=339), artificial_read(flag=163))==False)