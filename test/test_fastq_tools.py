from sequencing_tools.fastq_tools import *
from sequencing_tools.fastq_tools.function_clip import clip_read
from sequencing_tools.fastq_tools.pe_align import ConsensusBuilder
import six
import numpy.random as random


test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'

def test_complement():
    assert(complement('ACTGNactgn') == 'TGACNtgacn')


def test_reverse_complement():
    assert(reverse_complement('ACTGNactgn') == 'ncagtNCAGT')


def test_fastq():
    fq_iter = readfq(open(test_data_path + '/corrected.qual.fq'))
    record = six.next(fq_iter)
    assert(record.id=="AAAAAA_1_member/1")
    assert(record.seq=='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    assert(record.qual=="EAEE6666EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA")

    record.subseq(0,20)
    assert(record.seq == 'AAAAAAAAAAAAAAAAAAAA')
    assert(record.qual=='EAEE6666EEEEEEEEEEEE')

def test_kmer():
    test_seq = 'ACTGACT'
    assert(list(extract_kmer(test_seq, 2)) == ['AC','CT','TG','GA','AC', 'CT'])

    d = kmer_bag(test_seq, k_start = 1, k_end = 4)
    assert(d == {'A': 2,
             'AC': 2,
             'ACT': 2,
             'ACTG': 1,
             'C': 2,
             'CT': 2,
             'CTG': 1,
             'CTGA': 1,
             'G': 1,
             'GA': 1,
             'GAC': 1,
             'GACT': 1,
             'T': 2,
             'TG': 1,
             'TGA': 1,
             'TGAC': 1})

def test_onehot_default():
    dna_encoder = onehot_sequence_encoder(bases = 'ACTGN')
    test_seq = 'ACTGACTGNACTGNACTGN'

    onehot = np.array([
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.]
    ])
    assert(dna_encoder.base_encoder == {'A': 0, 'C': 1, 'G': 2, 'N': 3, 'T': 4})
    assert(np.array_equal(dna_encoder.transform(test_seq),onehot))
    assert(dna_encoder.decode(onehot) == test_seq)


def test_onehot_fit_transform():
    dna_encoder = onehot_sequence_encoder()
    test_seq = 'ACTGACTGNACTGNACTGN'

    onehot = np.array([
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.]
    ])

    out = dna_encoder.fit_transform(test_seq)
    assert(dna_encoder.base_encoder ==  {'A': 0, 'C': 1, 'G': 2, 'N': 3, 'T': 4})
    assert(np.array_equal(out, onehot))    
    assert(dna_encoder.decode(onehot) == test_seq)


def test_onehot_fit():
    dna_encoder = onehot_sequence_encoder()
    test_seq = 'ACTGACTGACTGACTG'
    dna_encoder.fit('ACTG')
    onehot = np.array([
       [ 1.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  0.]])
    
    assert(dna_encoder.base_encoder ==  {'A': 0, 'C': 1, 'G': 2, 'T': 3})
    assert(np.array_equal(dna_encoder.transform(test_seq), onehot))
    assert(dna_encoder.decode(onehot) == test_seq)


def test_insert_trimmer():

    constant = ''.join(random.choice(list('ACTG'), size = 70))
    seq1 = constant + 'TATATATA' 
    seq2 = reverse_complement(seq1) + 'ACTGACTG'
    qual1 = len(seq1) * 'A'
    qual2 = len(seq2) * 'A'
    assert(seq2 != reverse_complement(seq1))
    seq1, seq2, qual1, qual2 = clip_read().__insert_trimmer__(seq1, seq2, qual1, qual2)

    #assert(seq1 == constant)
    assert(seq1 == reverse_complement(seq2))
    assert(len(seq1) == len(seq2) == len(qual1) == len(qual2))



    constant = ''.join(random.choice(list('ACTG'), size = 10))
    seq1 = constant + 'TATATATA' 
    seq2 = reverse_complement(constant) + 'ACTGACTG'
    qual1 = len(seq1) * 'A'
    qual2 = len(seq2) * 'A'
    assert(seq2 != reverse_complement(seq1))
    new_seq1, new_seq2, qual1, qual2 = clip_read().__insert_trimmer__(seq1, seq2, qual1, qual2)

    #assert(seq1 == constant)
    assert(seq2 == new_seq2)
    assert(seq1 == new_seq1)


def test_umi_trimmer():
    clipping = clip_read(umi_bases=6)

    read1 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG',
        'ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAAAGATCGGAAGAGCACACGTCTGAA',
        'A6AA6//EA/AEE/AEAE///EEE/AAA/6/EE/A//EEE/A//EE/AE//E/////AEE///////////////',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG    Read1',
        )

    read2 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG',
        'TTTCCCGCACCGATAGCAGCCGCTCTGGTCTCCCATCCCGGGCAATTGTGTGATCGTCGGACTGTAGAACTCTGA',
        'AAA//E/E/E/A/A///EEE/E///<//EE/EE6/</EEA/A</EE<AAE/E/</<<AAE/E<///EE/EA///A',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG    Read2',
        )
                    
    ret_code, outread_1, outread_2 = clipping.trim_reads(read1, read2)
    assert(ret_code == 1)
    name1, seq1, _, qual1 = outread_1.strip().split('\n')
    name2, seq2, _, qual2 = outread_2.strip().split('\n')
    assert(name1.split('/')[0] == name2.split('/')[0])
    assert(seq2 == 'TTTCCCGCACCGATAGCAGCCGCTCTGGTCTCCCATCCCGGGCAA')
    assert(seq1 == 'TTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAA')
    assert(name1.split('_')[0].replace('@','') == read1.seq[:6])


def test_pe_align():
    consensus_builder = ConsensusBuilder(error_toleration = 0.1,
                                        min_len = 15, report_all=False)
    read1 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG',
        'ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAAAGATCGGAAGAGCACACGTCTGAA',
        'A6AA6//EA/AEE/AEAE///EEE/AAA/6/EE/A//EEE/A//EE/AE//E/////AEE///////////////',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG    Read1',
        )

    read2 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG',
        'TTTCCCGCACCGATAGCAGCCGCTCTGGTCTCCCATCCCGGGCAATTGTGTGATCGTCGGACTGTAGAACTCTGA',
        'AAA//E/E/E/A/A///EEE/E///<//EE/EE6/</EEA/A</EE<AAE/E/</<<AAE/E<///EE/EA///A',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG    Read2',
        )

    out = consensus_builder.run(read1, read2)
    id, seq, _, qual = out.strip().split('\n')

    res_seq ='ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAA'
    assert(res_seq == seq)



def test_pe_align():
    consensus_builder = ConsensusBuilder(error_toleration = 0.1,
                                        min_len = 15, report_all=True)
    read1 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG',
        'ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAAAGATCGGAAGAGCACACGTCTGAA',
        'A6AA6//EA/AEE/AEAE///EEE/AAA/6/EE/A//EEE/A//EE/AE//E/////AEE///////////////',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG    Read1',
        )

    read2 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG',
        'TTTCCCGCACCGATAGCAGCCGCTCTGGTCTCCCATCCCGGGCAATTGTGTGATCGTCGGACTGTAGAACTCTGA',
        'AAA//E/E/E/A/A///EEE/E///<//EE/EE6/</EEA/A</EE<AAE/E/</<<AAE/E<///EE/EA///A',
        'NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG    Read2',
        )

    out = consensus_builder.run(read1, read2)
    id, seq, _, qual = out.strip().split('\n')

    res_seq ='TCAGAGTTCTACAGTCCGACGATC'\
        'ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAA'\
        'AGATCGGAAGAGCACACGTCTGAA'
    assert(res_seq == seq)


