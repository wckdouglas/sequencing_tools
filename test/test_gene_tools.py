from sequencing_tools.gene_tools import Bed12Record


def test_bed12():
    line = 'chr1\t14403\t29570\tENST00000488147.1\t0\t-\t14403\t14403\t0\t11\t98,34,152,159,198,136,137,147,99,154,37,\t0,601,1392,2203,2454,2829,3202,3511,3864,10334,15130,'
    br = Bed12Record(line)
    first_intron = next(br.get_introns())
    assert(first_intron == 'chr1\t15167\t10334\tENST00000488147.1\t11\t-')
    assert(br.genomic_position(300) == 17094)

def test_GTF():
    line = 'chr1\tensGene\texon\t14404\t14501\t.\t-\t.\tgene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "1"; exon_id "ENST00000488147.1"; gene_name "ENSG00000227232";'
    gtf = GTFRecord(line)

    assert(gtf.info['gene_id'] == 'ENSG00000227232')
    assert(gtf.info['transcript_id'] == 'ENST00000488147')
    assert(gtf.info['exon_number'] == '1')
    assert(gtf.info['exon_id'] == 'ENST00000488147.1')
    assert(gtf.info['gene_name'] == 'ENSG00000227232')
    assert(gtf.strand == '-')
    assert(gtf.feature_type == "exon")

