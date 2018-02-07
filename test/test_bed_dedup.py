from sequencing_tools.bam_tools.bed_dedup import *

def test_hamming_barcode():
    assert(hamming_barcode(('ACTG','ACCC')) == 2)


def test_fragment_group():
    fg = fragment_group('chr1','10000', '20000', '-','ACT')
    fg.add_member('ACT')
    fg.add_member('ACG')
    assert(fg.check_fragment('chr1','10000','20000','+') == False)

    fg.demultiplexing_barcodes(0)
    expected = ['chr1\t10000\t20000\tACG_1_members\t10000\t-',
                 'chr1\t10000\t20000\tACT_2_members\t10000\t-']
    assert(list(fg.output_bed_line()) == expected)

    
    