from itertools import groupby

import pytest
from sequencing_tools.bam_tools.bed_dedup import fragment_coordinates, fragment_group


@pytest.mark.parametrize(
    "threshold,max_member_count,unique_barcodes,unique_count,molecular_count",
    [
        (1, 4, ["ACT_4_members", "ACG_4_members", "GGG_2_members"], 2, 6),
        (0, 3, ["ACG_1_members", "GGG_2_members", "ACT_3_members"], 3, 6),
    ],
)
def test_fragment_group(
    threshold, max_member_count, unique_barcodes, unique_count, molecular_count
):
    expected = [
        "chr1\t10000\t20000\tACG_name\t10000\t+",
        "chr1\t10000\t20000\tGGG_name\t10000\t+",
        "chr1\t10000\t20000\tACT_name\t10000\t+",
        "chr1\t10000\t20000\tGGG_name\t10000\t+",
        "chr1\t10000\t20000\tACT_name\t10000\t+",
        "chr1\t10000\t20000\tACT_name\t10000\t+",
    ]

    for coordinates, fragments in groupby(expected, fragment_coordinates):
        barcode_group = fragment_group(
            coordinates, fragments, umi_delim="_", umi_f=0, ct=-1
        )
        barcode_group.resolve_fragments()
        max_member_count = barcode_group.demultiplexing_barcodes(threshold=threshold)
        uc, mc = barcode_group.print_demultiplexed()
    assert uc == unique_count
    assert mc == molecular_count
    assert max_member_count == max_member_count
    assert set(barcode_group.unique_barcodes).issubset(unique_barcodes)
