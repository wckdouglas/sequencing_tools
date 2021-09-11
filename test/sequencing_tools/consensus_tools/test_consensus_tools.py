import os

import numpy as np
import pysam

from scipy.special import logsumexp
from sequencing_tools.bam_tools.pileup_errors import analyze_region
from sequencing_tools.consensus_tools import (ConsensusAlignments,
                                              ErrorCorrection,
                                              calculatePosterior)

test_data_path = os.path.dirname(os.path.realpath(__file__)) + "/data"


def test_alignment():
    bam = test_data_path + "/MT_TF.bam"
    ca = ConsensusAlignments(bam)
    con = ca.consensus_seq("MT-TF", 0, 75)
    seq = "NGTNTATGTTGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACCCCN"
    assert con == seq


def test_stranded_base_qual():
    bam = test_data_path + "/MT_TF.bam"
    bam = pysam.Samfile(bam)
    aln_count, base_dict = analyze_region(bam, "MT-TF", 30, False, False, 0, 75)
    assert aln_count == 11
    assert base_dict[9]["+"]["T"] == 8

    # base qual cutoff
    coverages = [
        0,
        1,
        1,
        0,
        1,
        1,
        1,
        1,
        1,
        8,
        9,
        8,
        8,
        8,
        6,
        8,
        5,
        9,
        8,
        8,
        8,
        8,
        7,
        5,
        4,
        8,
        8,
        5,
        7,
        7,
        4,
        9,
        4,
        8,
        6,
        7,
        5,
        8,
        8,
        7,
        8,
        9,
        9,
        10,
        8,
        10,
        7,
        0,
        10,
        6,
        9,
        7,
        8,
        9,
        9,
        1,
        8,
        3,
        8,
        7,
        1,
        9,
        9,
        9,
        9,
        9,
        10,
        9,
        6,
        2,
        7,
        1,
        1,
        1,
        0,
    ]

    inferred_coverages = np.array([sum(base_dict[i]["+"].values()) for i in range(75)])
    assert np.isclose(np.array(coverages), inferred_coverages).all()


def test_stranded_no_base_qual():
    bam = test_data_path + "/MT_TF.bam"
    bam = pysam.Samfile(bam)
    aln_count, base_dict = analyze_region(bam, "MT-TF", 0, False, False, 0, 75)
    assert aln_count == 11
    assert base_dict[3]["+"]["T"] == 1

    # no base qual cutoff
    coverages = [
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        8,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        9,
        7,
        9,
        9,
        10,
        11,
        11,
        11,
        11,
        11,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        9,
        9,
        1,
        1,
        1,
        0,
    ]
    inferred_coverages = np.array([sum(base_dict[i]["+"].values()) for i in range(75)])
    assert np.isclose(np.array(coverages), inferred_coverages).all()


def test_error_correction():
    correction_mode = ErrorCorrection(mode="prob")
    column_bases = np.array(list("CCAAAAAGT"))
    column_qualities_str = np.array(list(")))A--A)A"))
    column_qualities = np.array(list(map(ord, column_qualities_str))) - 33
    possible_bases = np.unique(column_bases)
    log_posteriors = [
        calculatePosterior(column_bases, column_qualities, guess_base)
        for guess_base in possible_bases
    ]

    pos_res = [
        -17.593092907715462,
        -39.35055643122346,
        -42.118680221372976,
        -36.420550581755265,
    ]
    assert np.all(np.isclose(log_posteriors, pos_res))

    log_posterior_sum = logsumexp(pos_res)
    loglik = np.array(log_posteriors) - log_posterior_sum

    base, loglikelihood = correction_mode.__posteriorConcensus__(
        (column_bases, column_qualities_str, "")
    )
    assert np.isclose(loglikelihood, np.exp(loglik).max())
    assert base == "A"
