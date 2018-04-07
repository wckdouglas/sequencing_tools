import os
import string
import six
from collections import defaultdict

# define fastq record type
cdef class fastqRecord:
    def __init__(self, str id, str seq, str qual):
        self.id = id
        self.seq = seq
        self.qual = qual

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
        str l, name, seq

    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield fastqRecord(name, ''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield fastqRecord(name, seq, ''.join(seqs)) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break



def gzopen(filename, read_flag = 'rb'):
    if 'r' in read_flag:
        return os.popen('zcat '+ filename)
    elif 'w' in read_flag:
        return open(filename, read_flag)

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



def kmer_bag(str sequence, k_range = 5):
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
    '''

    cdef: 
        int k
        str kmer

    bag = defaultdict(int)
    assert k_range > 1 and k_range <= len(sequence), \
            'Bad k_range being used!!'

    i = 0
    for k in range(1, k_range):
        for kmer in extract_kmer(sequence, k):
            bag[kmer] += 1
    
    return bag
