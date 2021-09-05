import pytest

from sequencing_tools.gene_tools import Bed12Record, GTFRecord


def gtfline():
    line = 'chr1\tensGene\texon\t14404\t14501\t.\t-\t.\tgene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "1"; exon_id "ENST00000488147.1"; gene_name "ENSG00000227232";'
    gtf = GTFRecord(line)
    return gtf


def test_bed12():
    line = "chr1\t14403\t29570\tENST00000488147.1\t0\t-\t14403\t14403\t0\t11\t98,34,152,159,198,136,137,147,99,154,37,\t0,601,1392,2203,2454,2829,3202,3511,3864,10334,15130,"
    br = Bed12Record(line)
    first_intron = next(br.get_introns())
    assert first_intron == "chr1\t24892\t29533\tENST00000488147.1\t1\t-"
    assert br.genomic_position(300) == 17904


@pytest.mark.parametrize(
    "attr, expected_result",
    [
        ("gene_id", "ENSG00000227232"),
        ("transcript_id", "ENST00000488147"),
        ("exon_number", "1"),
        ("gene_name", "ENSG00000227232"),
    ],
)
def test_GTF_info(attr, expected_result):
    gtf = gtfline()
    assert gtf.info[attr] == expected_result


def test_GTF():
    gtf = gtfline()
    assert gtf.strand == "-"
    assert gtf.feature_type == "exon"
