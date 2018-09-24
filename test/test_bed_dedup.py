
from __future__ import print_function
from sequencing_tools.bam_tools import bed_dedup 
from sequencing_tools.bam_tools import poisson_umi_tools
from itertools import groupby



def test_fragment_group():
    expected = ['chr1\t10000\t20000\tACG_name\t10000\t+',
                'chr1\t10000\t20000\tGGG_name\t10000\t+',
                'chr1\t10000\t20000\tACT_name\t10000\t+',
                'chr1\t10000\t20000\tGGG_name\t10000\t+',
                'chr1\t10000\t20000\tACT_name\t10000\t+',
                 'chr1\t10000\t20000\tACT_name\t10000\t+']


    member_count = 0
    for coordinates, fragments in groupby(expected, bed_dedup.fragment_coordinates):
        barcode_group = bed_dedup.fragment_group(coordinates, 
                                fragments, 
                                umi_delim = '_', 
                                umi_f = 0, 
                                ct = -1)
        barcode_group.resolve_fragments()
        max_member_count = barcode_group.demultiplexing_barcodes(threshold = 1)
        uc, mc = barcode_group.print_demultiplexed()
    assert(max_member_count == 4)
    assert(mc == len(expected))
    assert(uc == 2)
    barcode_group.unique_barcodes == ['ACG_4_members', 'GGG_2_members']


    max_member_count = barcode_group.demultiplexing_barcodes(threshold = 0)
    assert(max_member_count == 3)
    assert(barcode_group.unique_barcodes == ['ACG_1_members', 'GGG_2_members', 'ACT_3_members'])


def test_umi():
    import os
    import sys
    import string

    template = 'chr1\t100\t200\t{umi}_1_member\t0\t+'
    test_f = []
    for s in string.ascii_letters:
        test_f += [template.format(umi = (s * 3))]

    outfile = open(os.devnull, 'w')
    ic, oc = poisson_umi_tools.parse_dedup(test_f, outfile, umi_nt = 3, read_prefix=None)
    assert(ic == 107 and oc == 52)


