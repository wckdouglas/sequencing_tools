import os
import re
import string
import six
from collections import defaultdict

# define fastq record type
cdef class fastqRecord:
    def __init__(self, str id, str seq, str qual, str description):
        self.id = id
        self.seq = seq
        self.qual = qual
        self.description = description

    def subseq(self, int start, int end):
        self.seq = self.seq[start:end]
        self.qual = self.qual[start:end]

    def __str__(self):
        return '@{id}\n{seq}\n+\n{qual}'.format(id = self.id, 
                                                seq = self.seq, 
                                                qual = self.qual)

    def to_string(self):
        return self.__str__()


def readfq(fp): # this is a generator function
    '''
    A fastq iterator
        https://github.com/lh3/readfq/blob/master/readfq.py


    usage: readfq(fp)
    ==============================
    Parameter:

    fp: file handle of a fastq file

    return:

    fastqRecord object
        name: sequence id
        seq: sequence
        qual: quality
    ===============================
    '''
    cdef:
        str l, name, seq, description

    space_strip = re.compile('\s+')
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, description, seqs, last = last[1:].partition(" ")[0], last[1:].split('\n')[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield fastqRecord(name, ''.join(seqs), None, description) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield fastqRecord(name, seq, ''.join(seqs), description) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None, None # yield a fasta record instead
                break



# python 3 compatibility
try:
    complement_seq = string.maketrans('ACTGNactgn','TGACNtgacn')
except AttributeError:
    complement_seq = str.maketrans('ACTGNactgn','TGACNtgacn')

def complement(seq):
    """
    Find complement a sequence.

    ============================
    parameter:
    
    string seq: sequence to be complemented
    

    return:

    complemented sequence
    """
    return seq.translate(complement_seq)


def reverse_complement(seq):
    """
    Reverse complement a sequence.

    ============================
    parameter:
    
    string seq: sequence to be reverse complemented
    

    return:

    reverse complemented sequence
    """
    return complement(seq)[::-1]


def read_interleaved(infile):   
    '''
    A interleaved fastq iterator

    usage: parse_interleaved(fp)
    ==============================
    Parameter:

    fp: file handle of a fastq file

    return:

    R1: fastqRecord object
        name: sequence id
        seq: sequence
        qual: quality

    R2: fastqRecord object
        name: sequence id
        seq: sequence
        qual: quality

    ===============================
    '''

    cdef:
        str r1_id, r2_id
        fastqRecord R1, R2
        
    fastq_file = readfq(infile)

    try:
        while True:
            R1 = six.next(fastq_file)
            R2 = six.next(fastq_file)
    
            r1_id, r2_id = R1.id.split('/')[0], R2.id.split('/')[0]
            assert r1_id == r2_id, 'Not interleaved'
            yield R1, R2

    except StopIteration:
        pass

from collections import defaultdict
def extract_kmer(str sequence, int k):
    '''
    output kmer
    '''
    cdef int i

    for i in range(len(sequence) - k + 1):
        yield(sequence[i:i+k])



def kmer_bag(str sequence, k_start = 1, k_end = 4):
    '''
    K-mer bag method for feature extraction

    ---
    We used k = 1,2,3,4,5, resulting in 41 + 42 + 43 + 44 + 45 = 1364 features. 
    For example, the se- quence ACTGG would produce a length-1364 feature vector
    where the entries corresponding to the k-mers A, C, T, AC, CT, TG, GG, ACT, 
    CTG, TGG, ACTG, CTGG, and ACTGG would each equal 1, 
    and the entry corresponding to G would equal 2. 
    All other entries equal 0.
    (from Zhang and Kamath. Learning the Language of the Genome using RNNs)


    usage:
    kmer_bag(str sequence, k_start = 1, k_end = 4)
    '''

    cdef: 
        int k
        str kmer

    bag = defaultdict(int)
    assert k_start >= 1 and k_end <= len(sequence) and k_start < k_end, \
            'Bad k_range being used!!'

    i = 0
    for k in range(k_start, k_end + 1):
        for kmer in extract_kmer(sequence, k):
            bag[kmer] += 1
    
    return bag

import numpy as np
class onehot_sequence_encoder:
    '''
    A onehot encoder for DNA sequence

    usage:
        dna_encoder = onehot_sequence_encoder()
        onehot_encoded_matrix = dna_encoder.fit_transform(sequence)
        dna_encoder.base_encoder  # check which base each column represents
    '''

    def __init__(self, bases = 'ACTGN'):
        self.fit(sequence = bases)
    

    def fit(self, sequence='ACTGN'):
        self.bases = list(set(sequence))
        self.bases.sort()
        self.base_encoder = {b:i for i, b in enumerate(self.bases)}
        self.base_decoder = {i:b for b, i in self.base_encoder.items()}
        self.acceptable_nuc = set(self.bases)
        self.column_number = len(self.acceptable_nuc)

    def transform(self, sequence):   
        '''
        One hot sequence encoder

        Parameter:
            sequence: a string of sequence, only accept 'ACTGN'

        return:
            onehot encoded array:  len(sequence)-by-distinct(base) matrix
                                    columns represent each base
                                    rows represent each position along the sequence

        '''
        cdef:
            int pos
            str base

        encoded_mat = np.zeros((len(sequence),self.column_number))

        assert set(sequence).issubset(self.acceptable_nuc),  \
            'Sequence contain bases other than "ACTGN:"' + sequence

        for pos, base in enumerate(sequence):
            encoded_mat[pos, self.base_encoder[base]] = 1
        return encoded_mat

    def fit_transform(self, sequence):   
        '''
        One hot sequence encoder

        Parameter:
            sequence: a string of sequence, only accept 'ACTGN'

        return:
            onehot encoded array:  len(sequence)-by-distinct(base) matrix
                                    columns represent each base
                                    rows represent each position along the sequence

        '''
        cdef:
            int pos
            str base

        self.fit(sequence)
        return self.transform(sequence)
 
    def decode(self, encoded_mat):
        '''
        Parameter:
            onehot encoded array: len(sequence)-by-distinct(base) matrix
                                    columns represent each base
                                    rows represent each position along the sequence

        Return:
            sequence: a string of sequence, only accept 'ACTGN'
        '''

        decoded = np.matmul(encoded_mat, np.arange(len(self.bases)))
        return ''.join([self.base_decoder[i] for i in decoded])
