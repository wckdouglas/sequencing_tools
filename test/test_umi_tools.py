#!/usr/bin/env python

from sequencing_tools.bam_tools import umi_network
import time
from collections import Counter


def test_umi_dedup_directional():
    umis = ["ACTGA"] * 10 + ["ACTGG"] * 2 + ["ACTCC"] * 1 + ["ACCCC"] * 17

    umi_counter = Counter(umis)
    out, mc = umi_network.demultiplex_directional(umi_counter)
    assert sum(int(s.split("_")[1]) for s in out) == len(umis)
    umi_set = set(s.split("_")[0] for s in out)
    assert umi_set == {"ACCCC", "ACTGA"}
    assert mc == 18

    umis = ["ACTG", "ACCG"]
    umi_counter = Counter(umis)
    out, mc = umi_network.demultiplex_directional(umi_counter)
    assert sum(int(s.split("_")[1]) for s in out) == len(umis)
    assert mc == 2
