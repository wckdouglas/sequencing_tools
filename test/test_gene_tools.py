from sequencing_tools.gene_tools import Bed12Record


def test_bed12():
    line = 'chr1\t11868\t14409\tENST00000456328.2_1\t3\t+\t14409\t14409\t11868,12612,13220,\t12227,12721,14409,\n'
    br = Bed12Record(line)
    first_intron = next(br.get_introns())
    assert(first_intron == 'chr1\t12227\t12612\tENST00000456328.2_1\t1\t+')