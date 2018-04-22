from sequencing_tools.fastq_tools import *
import six


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

    d = kmer_bag(test_seq, k_start = 1, k_end = 5)
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
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 0.,  0.,  0.,  0.,  1.],
       [ 1.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.],
       [ 0.,  0.,  0.,  0.,  1.]
    ])
    assert(dna_encoder.base_encoder == {'A': 0, 'C': 1, 'G': 3, 'N': 4, 'T': 2})
    assert(np.array_equal(dna_encoder.transform(test_seq),onehot))


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

