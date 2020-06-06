
from __future__ import print_function
import os
from sequencing_tools.fasta_tools import readfa, MultiAlignments
from sequencing_tools.fastq_tools import readfq
import numpy as np
test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'


def test_readfa():
    fa_file = test_data_path + '/corrected.qual.fa'
    fq_iter = readfq(open(fa_file.replace('.fa','.fq')))
    with open(fa_file, 'w') as fa:
        for record in fq_iter:
            print('>{id}\n{seq}'.format(id = record.id, seq=record.seq), file = fa)

    fa_iter = readfa(open(fa_file))
    for seqid, seq in fa_iter:
        pass

    assert(seq == record.seq)

def test_multialign():
    ma = MutliAlignments(test_data_path + '/multi.fa')
    consensus_seq, score = ma.concensus()
    assert(consensus_seq == 'GGGGAATTAGCTCAAGCGGTAAAACGCTTGCTTAGCATGCAAGAGGTAGTGGGATCGATG')
    assert( min(max(s) for s in score) == 0.5 )
    assert( min(min(s) for s in score) == 1/6 )
    assert( min(min(s) for s in score) == 1/6 )


    ma.PairMatrix()
    pairwise_dist = [0, 0, 2, 4, 6, 1, 0, 0, 2, 4, 6, 1, 2, 
                    2, 0, 2, 4, 3, 4, 4, 2, 0, 2, 5, 6, 6, 4, 
                    2, 0, 7, 1, 1, 3, 5, 7, 0]
    assert( np.all(ma.pairwise.distance == pairwise_dist ) )
