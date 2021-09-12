import os
import string

from sequencing_tools.bam_tools.poisson_umi_tools import parse_dedup


def test_umi():
    template = "chr1\t100\t200\t{umi}_1_member\t0\t+"
    test_f = []
    for s in string.ascii_letters:
        test_f += [template.format(umi=(s * 3))]

    outfile = open(os.devnull, "w")
    ic, oc = parse_dedup(test_f, outfile, umi_nt=3, read_prefix=None)
    assert ic == 107 and oc == 52
