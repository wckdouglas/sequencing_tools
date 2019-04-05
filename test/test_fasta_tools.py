from sequencing_tools.fasta_tools import readfa
from sequencing_tools.fastq_tools import readfq
import os
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




