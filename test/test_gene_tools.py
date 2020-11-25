from sequencing_tools.gene_tools import Bed12Record


def test_bed12():
    line = 'chr1\t14403\t29570\tENST00000488147.1\t0\t-\t14403\t14403\t0\t11\t98,34,152,159,198,136,137,147,99,154,37,\t0,601,1392,2203,2454,2829,3202,3511,3864,10334,15130,'
    br = Bed12Record(line)
    first_intron = next(br.get_introns())
    assert(first_intron == 'chr1\t14501\t15004\tENST00000488147.1\t11\t-')