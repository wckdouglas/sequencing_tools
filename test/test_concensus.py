import os
import numpy as np
import pysam
from sequencing_tools.consensus_tools import ConsensusAlignments
from sequencing_tools.bam_tools.pileup_errors import analyze_region
test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'


def test_alignment():
    bam = test_data_path + '/MT_TF.bam'
    ca = ConsensusAlignments(bam)
    con = ca.consensus_seq('MT-TF',0, 75)
    seq = 'NGTNTATGTTGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACCCCN'
    assert(con == seq)


def test_stranded_base_qual():
    bam = test_data_path + '/MT_TF.bam'
    bam = pysam.Samfile(bam)
    aln_count, base_dict = analyze_region(bam, 'MT-TF', 30, False, False, 0, 75)
    assert(aln_count == 11)
    assert(base_dict[9]['+']['T'] == 8)

    # base qual cutoff
    coverages = [ 0,  1,  1,  0,  1,  1,  1,  1,  1,  8,  9,  8,  8,  8,  6,  8,  5,
        9,  8,  8,  8,  8,  7,  5,  4,  8,  8,  5,  7,  7,  4,  9,  4,  8,
        6,  7,  5,  8,  8,  7,  8,  9,  9, 10,  8, 10,  7,  0, 10,  6,  9,
        7,  8,  9,  9,  1,  8,  3,  8,  7,  1,  9,  9,  9,  9,  9, 10,  9,
        6,  2,  7,  1,  1,  1,  0]

    inferred_coverages = np.array([sum(base_dict[i]['+'].values()) for i in range(75)])
    assert(np.isclose(np.array(coverages),
                    inferred_coverages).all())


def test_stranded_no_base_qual():    
    bam = test_data_path + '/MT_TF.bam'
    bam = pysam.Samfile(bam)
    aln_count, base_dict = analyze_region(bam, 'MT-TF', 0, False, False, 0, 75)
    assert(aln_count == 11)
    assert(base_dict[3]['+']['T'] == 1)

    # no base qual cutoff
    coverages = [ 0,  1,  1,  1,  1,  1,  1,  1,  1,  8,  9,  9,  9,  9,  9,  9,  9,
        9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
        9,  9,  7,  9,  9, 10, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10,
       10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
       10,  9,  9,  1,  1,  1,  0]
    inferred_coverages = np.array([sum(base_dict[i]['+'].values()) for i in range(75)])
    assert(np.isclose(np.array(coverages),
                    inferred_coverages).all())

   

