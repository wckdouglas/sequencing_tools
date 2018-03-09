from __future__ import print_function
from sequencing_tools.bam_tools import bed_dedup 
from sequencing_tools.bam_tools import poisson_umi_tools


def test_hamming_barcode():
    assert(bed_dedup.barcode_distance(('ACTG','ACCC')) == 2)


def test_fragment_group():
    fg = bed_dedup.fragment_group('chr1','10000', '20000', '-','ACT','')
    fg.add_member('ACT','')
    fg.add_member('ACG','')
    assert(not fg.check_fragment('chr1','10000','20000','+'))

    fg.demultiplexing_barcodes(0)
    expected = ['chr1\t10000\t20000\tACG_1_members\t10000\t-',
                 'chr1\t10000\t20000\tACT_2_members\t10000\t-']
    expected.sort()
    out = list(fg.output_bed_line())
    out.sort()
    assert(out == expected)


    fg.demultiplexing_barcodes(1)
    expected = ['chr1\t10000\t20000\tACT_3_members\t10000\t-', 
            'chr1\t10000\t20000\tACG_3_members\t10000\t-']
    out = list(fg.output_bed_line())[0]
    assert(out in expected)



    

def test_umi():
    import os
    import sys
    import string
    f = open(os.devnull, 'w')
    test_fragment = poisson_umi_tools.fragment('chr1','1000','2000','-', 'AAA')
    for s in string.ascii_letters:
        test_fragment.add_fragment(s * 3)
    assert(test_fragment.output_fragments(umi_nt = 3, out_file = f) == 107)
    assert(test_fragment.output_fragments(umi_nt = 5, out_file = f) == 53)


