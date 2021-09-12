import os
from pathlib import Path

import pysam
import pytest
from pysam import AlignedRead

from sequencing_tools.bam_tools import (
    split_cigar,
    read_ends,
    fragment_ends,
    concordant_pairs,
    concordant_alignment,
    make_regions,
    make_cigar_seq,
    cigar_to_str,
    get_strand,
    check_concordant,
    check_primary,
    paired_bam
)


def artificial_read(flag=163):
    aln = AlignedRead()
    aln.qname = "try_this"
    aln.flag = flag
    aln.cigarstring = "16M"
    aln.seq = "A" * 16
    aln.pos = 1000000 if flag == 163 else 1000020
    aln.reference_id = 1
    return aln


@pytest.fixture(scope="module")
def read1():
    return artificial_read(flag=83)


@pytest.fixture(scope="module")
def read2():
    return artificial_read(flag=163)


@pytest.mark.parametrize(
    "cigar, result",
    [
        ("63M4S", [[63, 4], ["M", "S"]]),
        ("5M4D6I10M", [[5, 4, 6, 10], ["M", "D", "I", "M"]]),
    ],
)
def test_split_cigar(cigar, result):
    assert split_cigar(cigar) == result


def test_read_ends(read2):
    aln = artificial_read()
    assert read_ends(aln) == (1000000, 1000015)


def test_fragment_ends(read1, read2):
    assert fragment_ends(read1, read2) == (1000000, 1000036)


def test_concordant_pairs(read1, read2):
    assert concordant_pairs(read1, read2) == True


def test_non_concordant_pairs(read2):
    read1 = artificial_read(flag=81)
    assert concordant_pairs(read1, read2) == False


def test_concordant_alignment(read2):
    read1 = artificial_read(flag=81)
    assert concordant_alignment(read1) == False
    assert concordant_alignment(read2) == True


@pytest.mark.parametrize(
    "end, interval, output",
    [(100, 51, [(0, 51), (51, 100)]), (10, 5, [(0, 5), (5, 10)])],
)
def test_make_regions(end, interval, output):
    assert list(make_regions(end, interval)) == output


@pytest.mark.parametrize(
    "cigar_nums, cigar_ops, result",
    [
        (["3", "5", "1", "3"], ["S", "M", "I", "M"], "SSSMMMMMIMMM"),
        (["10", "5", "2", "5"], ["M", "D", "I", "M"], "MMMMMMMMMMIIMMMMM"),
    ],
)
def test_make_cigar_seq(cigar_nums, cigar_ops, result):
    seq = ""
    for c in make_cigar_seq(cigar_nums, cigar_ops):
        seq += c
    assert seq == result


@pytest.mark.parametrize(
    "cigar, cigar_seq",
    [
        ("3S5M1I3M", "SSSMMMMMIMMM"),
        ("5M4D6I10M", "MMMMMIIIIIIMMMMMMMMMM"),
    ],
)
def test_cigar_to_str(cigar, cigar_seq):
    result = cigar_to_str(cigar)
    assert result == cigar_seq


@pytest.mark.parametrize(
    "flag, expected_result", [(83, "-"), (163, "-"), (99, "+"), (147, "+")]
)
def test_get_strand(flag, expected_result):
    assert get_strand(artificial_read(flag=flag)) == expected_result


@pytest.mark.parametrize(
    "r1_flag, r2_flag",
    [
        (83, 163),
        (99, 147),
    ],
)
def test_check_concordant(r1_flag, r2_flag):
    assert check_concordant(
        artificial_read(flag=r1_flag), artificial_read(flag=r2_flag)
    )


def test_check_primary(read1, read2):
    assert check_primary(read1, read2) == True


def test_check_primary_false(read2):
    assert check_primary(artificial_read(flag=339), read2) == False


def test_paired_bam():
    test_data_path = Path(__file__).parents[2] / "data"
    bam = test_data_path / "clipped.bam"
    bam = pysam.Samfile(bam)
    pairs = 0
    for read1, read2 in paired_bam(bam):
        pairs += 1
    assert pairs == 218
