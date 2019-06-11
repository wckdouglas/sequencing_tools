import os
from sequencing_tools.consensus_tools import ConsensusAlignments
test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'


def test_alignment():
    bam = test_data_path + '/MT_TF.sorted.bam'
    ca = ConsensusAlignments(bam)
    con = ca.consensus_seq('MT-TF',0, 75)
    seq = 'NGTNTATGTTGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACCCCN'
    assert(con == seq)