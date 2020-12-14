import os
import re
import string
import six
import numpy as np
from collections import defaultdict
from ..utils import SeqUtilsError

# define fastq record type
cdef class fastqRecord:
    """
    Fastq record

    Args:
        id (str): fastq record id
        seq (str): fastq record seq
        qual (str): fastq record quality string
    """
    def __init__(self, str id, str seq, str qual, str description):
        self.id = id #: sequence id
        self.seq = seq # sequence
        self.qual = qual # quality score string
        self.description = description # description for the fastq record

    @property
    def length(self):
        return self.__len__()
    
    def __len__(self):
        return len(self.seq)

    def subseq(self, int start, int end):
        '''
        Trim off sequences outside of start:end

        Args:
            start (int):  start position on the sequence
            end (str):  end position on the sequence
        '''
        if start < 0 or start > self.__len__():
            raise SeqUtilsError('Start position must be positive and smaller than sequence length: %i' %self.__len__())

        if end < start or end > self.__len__():
            raise SeqUtilsError('End position must be greater than start position and smaller than sequence length: %i' %self.__len__())

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
    A fastq iterator from `Heng li <https://github.com/lh3/readfq/blob/master/readfq.py>`_

    Args:
        fp: file handle of a fastq file

    Returns:
        Generator(FastqRecord): :class:`sequencing_tools.fastq_tools._fastq_tools.fastqRecord`
    
    Usage::

        with open('test.fq') as fastq:
            for record in readfq(fastq):
                print(record.id, record.seq, record.qual, record.description)
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
    complement_seq = string.maketrans('ACTGNactgnRYSWKMDHVryswkmdhv','TGACNtgacnBDHKMWSRYbdhkmwsry')
except AttributeError:
    complement_seq = str.maketrans('ACTGNactgnRYSWKMDHVryswkmdhv','TGACNtgacnBDHKMWSRYbdhkmwsry')

def complement(seq):
    """
    Find complement a sequence.

    Args:
        seq (str): sequence to be complemented
    

    Returns:
        str: complemented sequence
    """
    return seq.translate(complement_seq)


def reverse_complement(seq):
    """
    Reverse complement a sequence.

    Args:
        seq (str): sequence to be reverse complemented
    

    Returns:
        str: reverse complemented sequence
    """
    return complement(seq)[::-1]


def read_interleaved(infile):   
    '''
    A interleaved fastq iterator

    Usage::

        with open('test.fq') as fastq:
            for read1, read2 in parse_interleaved(fastq):
                print(read1.id, read1.seq, read1.qual, read1.description)
                print(read2.id, read2.seq, read2.qual, read2.description)

    Args:
        fp: file handle of a fastq file

    Returns:
        Generator(tuple):  (fastqRecord for R1, fastqRecord for R2) :class:`sequencing_tools.fastq_tools._fastq_tools.fastqRecord` 
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

def extract_kmer(str sequence, int k):
    '''
    Args:
        sequence (str): input sequence for extracting kmer
        k (int): kmer size
    
    Returns:
        generator(str): kmer from the sequence
    '''
    cdef int i

    for i in range(len(sequence) - k + 1):
        yield(sequence[i:i+k])



def kmer_bag(str sequence, k_start = 1, k_end = 4):
    '''
    K-mer bag method for feature extraction

    We used k = 1,2,3,4,5, resulting in 41 + 42 + 43 + 44 + 45 = 1364 features. 
    For example, the sequence ACTGG would produce a length-1364 feature vector
    where the entries corresponding to the k-mers A, C, T, AC, CT, TG, GG, ACT, 
    CTG, TGG, ACTG, CTGG, and ACTGG would each equal 1, 
    and the entry corresponding to G would equal 2. 
    All other entries equal 0.
    (from Zhang and Kamath. Learning the Language of the Genome using RNNs)


    Example::

        sequence = 'ACTGG'
        kmer_bag(sequence, k_start = 1, k_end = 4) 
        #defaultdict(int,
        #    {'A': 1,
        #     'C': 1,
        #     'T': 1,
        #     'G': 2,
        #     'AC': 1,
        #     'CT': 1,
        #     'TG': 1,
        #     'GG': 1,
        #     'ACT': 1,
        #     'CTG': 1,
        #     'TGG': 1,
        #     'ACTG': 1,
        #     'CTGG': 1})
    
    Args:
        sequence (str): the sequence for feature extraction
        k_start (int): k for the smallest kmer
        k_end (int): k for the largest kmer

    Returns:
        dict: key is kmer, values: count
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

class onehot_sequence_encoder:
    """
    A onehot encoder for DNA sequence

    usage::

        sequence = 'AAACTTTG'
        dna_encoder = onehot_sequence_encoder()
        onehot_encoded_matrix = dna_encoder.fit_transform(sequence)
        onehot_encoded_matrix                                                                                                                                                           
        #array([[1., 0., 0., 0.],
        #    [1., 0., 0., 0.],
        #    [1., 0., 0., 0.],
        #    [0., 1., 0., 0.],
        #    [0., 0., 0., 1.],
        #    [0., 0., 0., 1.],
        #    [0., 0., 0., 1.],
        #    [0., 0., 1., 0.]])
        dna_encoder.base_encoder  # check which base each column represents
        # {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    """

    def __init__(self, bases = 'ACTGN'):
        self.fit(sequence = bases)
    

    def fit(self, sequence='ACTGN'):
        """
        Args:
            sequence (str): the bases to be extract
        """
        self.bases = list(set(sequence))
        self.bases.sort()
        self.base_encoder = {b:i for i, b in enumerate(self.bases)}
        self.base_decoder = {i:b for b, i in self.base_encoder.items()}
        self.acceptable_nuc = set(self.bases)
        self.column_number = len(self.acceptable_nuc)

    def transform(self, sequence):   
        '''
        One hot sequence encoder

        Args:
            sequence (str): a string of sequence, only accept the bases suppled in `fit`, with length n

        Returns:
            np.ndarray(n,m): onehot encoded array, len(sequence)-by-distinct(base) matrix, columns represent each base, rows represent each position along the sequence
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

        Args:
            sequence (str): a string of sequence, only accept the bases suppled in `fit`, with length n, with m unique bases

        Returns:
            np.ndarray(n,m): onehot encoded array, len(sequence)-by-distinct(base) matrix, columns represent each base, rows represent each position along the sequence
        '''
        cdef:
            int pos
            str base

        self.fit(sequence)
        return self.transform(sequence)
 
    def decode(self, encoded_mat):
        '''
        Args:
            onehot encoded array: len(sequence)-by-distinct(base) matrix, columns represent each base, rows represent each position along the sequence

        Returns:
            str: Sequence, the decoded sequence
        

        Usage::

            mat = np.array([[1., 0., 0., 0.], 
                            [1., 0., 0., 0.], 
                            [1., 0., 0., 0.], 
                            [0., 1., 0., 0.], 
                            [0., 0., 0., 1.], 
                            [0., 0., 0., 1.], 
                            [0., 0., 0., 1.], 
                            [0., 0., 1., 0.]]) 
            dna_encoder.decode(mat)   #'AAACTTTG'
            
        '''

        decoded = np.matmul(encoded_mat, np.arange(len(self.bases)))
        return ''.join([self.base_decoder[i] for i in decoded])
