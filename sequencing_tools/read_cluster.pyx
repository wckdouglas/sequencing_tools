from __future__ import division, print_function
import numpy as np
import pysam
from builtins import zip, map
import string
from libc.math cimport log10, exp
import sys
from cpython cimport bool
from scipy.misc import logsumexp

np_len = np.vectorize(len,otypes=[np.int32])
np_ord = np.vectorize(ord, otypes=[np.int16])

cdef:
    double min_q = 33.0
    double max_q = 73.0
    double max_prob = 0.999999


complement =  string.maketrans('ACTGNactgn','TGACNTGACN')
def reverse_complement(seq):
    '''
    Generate reverse complement as name suggested....-,-
    '''
    return seq.translate(complement)[::-1]


def fix_strand(str seq, str qual, bool strand):
    if strand:
        seq = reverse_complement(seq)
        qual = qual[::-1]
    return seq, qual


def mode(any_list):
    any_list = list(any_list)
    count_tuple = map(lambda x: (any_list.count(x), x), set(any_list))
    return max(count_tuple)[1]


cdef double clip(double q, double floor, double ceiling):
    '''
    mimic np.clip
    '''
    cdef:
        double res
    if q < floor:
        res = floor
    elif q > ceiling:
        res = ceiling
    else:
        res = q
    return res


cdef double error_prob_to_phred(double error_prob):
    '''
    probability to Phred score
    '''
    return -10 * log10(1- error_prob)


cdef str phred_to_string(double phred):
    '''
    convert quality score to string
    '''
    cdef:
        double adjusted_phred = phred + 33
        int integer_phred
        str string_phred

    adjusted_phred = clip(adjusted_phred, min_q, max_q)
    integer_phred = int(adjusted_phred)
    string_phred =chr(integer_phred)
    return string_phred


cdef str prob_to_qual_string(posterior):
    '''
    Input a list of probabilities and output quality string
    '''
    cdef:
        double q
        str qual_char

    # convert posterior to quality score
    q = clip(posterior, 0, max_prob)
    q = error_prob_to_phred(q)
    qual_char = phred_to_string(q)
    return qual_char


def qual_to_prob(base_qual):
    '''
    Given a q list,
    return a list of prob
    '''
    cdef:
        double q

    return [10**(-q/10) for q in base_qual]


def cumulative_product_qual(qs, hit=True):
    '''
    cythonize numpy.prod
    '''
    if hit:
        result = np.log10(sum(1-q for q in qs))
    return result

cdef double cumulative_product(qs):
    '''
    cythonize numpy.prod
    '''
    cdef:
        double result = 1.0
        double q
    for q in qs:
        result = q * result
    return result


cdef double calculatePosterior(column_bases, column_qualities, guess_base):
    '''
    From each column of the sequence alignemnt (base position),
    extract the probabilty of being the guess_base
    '''
    cdef:
        double q
        double hit, missed, log_posterior

    correct_base = column_bases==guess_base
    qual_missed = column_qualities[~correct_base]
    qual_hit = column_qualities[correct_base]
    hit = sum(log10(1 - q) for q in qual_to_prob(qual_hit))
    missed = sum(log10(q/3.0) for q in qual_to_prob(qual_missed))
    log_posterior = missed + hit
    return log_posterior

def calculate_concensus_base(arg):
    """Given a list of sequences,
        a list of quality line and
        a position,
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    cdef:
        int number_possible_bases
        double log_posterior
        double total_posterior
        double posterior_correct_probability
    column_bases, in_column_qualities = arg
    column_qualities = np_ord(in_column_qualities) - 33
    possible_bases = np.unique(column_bases)
    number_possible_bases = len(possible_bases)
    if number_possible_bases == 1:
        concensus_base = possible_bases[0]
        posterior_correct_probability = max_prob
    else:
        log_posteriors = [calculatePosterior(column_bases, column_qualities, guess_base) for guess_base in possible_bases]
        total_posterior = logsumexp(log_posteriors)
        likelihoods = [log_posterior - total_posterior for log_posterior in log_posteriors]
        arg_max_likelihood = np.argmax(likelihoods)
        concensus_base = possible_bases[arg_max_likelihood]
        posterior_correct_probability = exp(likelihoods[arg_max_likelihood])
    return concensus_base, posterior_correct_probability


def vote_concensus_base(arg):
    """Given a list of sequences,
        a list of quality line and
        a position,
    return the maximum likelihood base at the given position,
        along with the mean quality of these concensus bases.
    """
    cdef:
        int number_possible_bases
        double log_posterior
        double total_posterior
        double posterior_correct_probability

    column_bases, in_column_qualities = arg
    column_qualities = np_ord(in_column_qualities) - 33
    possible_bases, possible_counts = np.unique(column_bases, return_counts=True)

    best_vote = possible_counts.max()
    best_fraction = best_vote/possible_counts.sum()
    if best_fraction < 0.9:
        concensus_base = 'N'
        posterior_correct_probability = 0

    else:
        best_index = possible_counts==possible_counts.max()
        concensus_base = possible_bases[best_index][0]
        best_quals = column_qualities[column_bases==concensus_base].sum()
        posterior_correct_probability = 1 - 10**(-best_quals/10)

    return concensus_base, posterior_correct_probability


def concensus_sequence(conserved, aln_table):
    """given a list of sequences, a list of quality and sequence length.
        assertion: all seq in seqlist should have same length (see function: selectSeqLength)
    return a consensus sequence and the mean quality line (see function: calculateConcensusBase)
    """
    cdef:
        str sequence = ''
        str quality = ''
        double posterior_correct_prob

    in_seq_list = aln_table[:,0]
    in_qual_list = aln_table[:,1]
    if len(in_seq_list) == 1:
        sequence = str(in_seq_list[0])
        quality = str(in_qual_list[0])
    else:
        len_array = np_len(in_seq_list)
        seq_len = min(len_array)
        len_filter = (len_array ==  seq_len)
        in_seq_list = in_seq_list[len_filter]
        in_qual_list = in_qual_list[len_filter]
        seq_list = np.array(list(map(list, in_seq_list)))
        qual_list = np.array(list(map(list, in_qual_list)))
        iter_list = ((seq_list[:,pos], qual_list[:,pos]) for pos in xrange(seq_len))
        if conserved:
            for base, posterior_correct_prob in map(vote_concensus_base, iter_list):
                sequence += base
                quality += prob_to_qual_string(posterior_correct_prob)
        else:
            for base, posterior_correct_prob in map(calculate_concensus_base, iter_list):
                sequence += base
                quality += prob_to_qual_string(posterior_correct_prob)
    return sequence, quality

cdef class readGroup:
    '''
    Read group object
    '''
    def __init__(self,aln, tag):
        #assert self.barcode == '', 'Cluster already initialzed with %s' %(self.barcode)

        #self.barcode = aln.query_name.split('_')[0]
        self.barcode = aln.get_tag(tag)
        self.R1 = []
        self.R2 = []
        self.R1_flag = []
        self.R2_flag = []
        self.R1_position = []
        self.R2_position = []
        self.R1_chrom = []
        self.R2_chrom = []
        self.concensus_read1 = []
        self.concensus_read2 = []
        self.member_count_list = []
        self.concensus_flag1 = []
        self.concensus_flag2 = []
        self.fastq_record = ''

        self.put_alignment(aln)


    def put_alignment(self, aln):
        '''
            add alignment to read group
        '''

        if aln.is_read1:
            self.R1.append([aln.query_sequence, aln.qual])
            self.R1_flag.append(aln.is_reverse)
            self.R1_position.append(aln.pos)
            self.R1_chrom.append(aln.reference_id)

        elif aln.is_read2:
            self.R2.append([aln.query_sequence, aln.qual])
            self.R2_flag.append(aln.is_reverse)
            self.R2_position.append(aln.pos)
            self.R2_chrom.append(aln.reference_id)


    def cluster(self, conserved):
        '''
            from read group, generate concensus sequence, quality
        '''
        assert self.R2 and self.R1, (self.R1, self.R2, self.barcode)
        assert self.R2.shape == self.R1.shape, 'Unequal R1 list vs R2'

        iterator = set(zip(self.R1_chrom, self.R2_chrom,
                           self.R1_position, self.R2_position,
                           self.R1_flag, self.R2_flag))

        R1_array = np.array(self.R1)
        R2_array = np.array(self.R2)
        R1_chrom_array = np.array(self.R1_chrom)
        R2_chrom_array = np.array(self.R2_chrom)
        R1_flag_array = np.array(self.R1_flag)
        R2_flag_array = np.array(self.R2_flag)
        R1_pos_array = np.array(self.R1_position)
        R2_pos_array = np.array(self.R2_position)

        for _chrom1, _chrom2, _pos1, _pos2, _R1_flag, _R2_flag in iterator:
            chrom_is_right = (R1_chrom_array == _chrom1) & (R2_chrom_array == _chrom2)
            flag_is_right = (R1_flag_array == _R1_flag) & (R2_flag_array == _R2_flag)
            pos_is_right = (R1_pos_array == _pos1) & (R2_pos_array == _pos2)

            cluster = chrom_is_right & flag_is_right & pos_is_right
            R1_filtered, R2_filtered = R1_array[cluster,:], R2_array[cluster,:]
            self.concensus_read1.append(concensus_sequence(conserved, R1_filtered))
            self.concensus_read2.append(concensus_sequence(conserved, R2_filtered))
            self.member_count_list.append(len(R1_filtered[:,1]))
            self.concensus_flag1.append(_R1_flag)
            self.concensus_flag2.append(_R2_flag)

    def generate_fastq_record(self):
        '''
            from concensus sequence generate fastq record
        '''
        cdef:
            str r1_seq, r1_qual, r2_seq, r2_qual
            int member_count
            bool strand1, strand2

        assert self.concensus_read1, 'No concensus read generated'
        iterable = zip(self.concensus_read1, self.concensus_read2,
                        self.member_count_list,
                        self.concensus_flag1, self.concensus_flag2)
        for (r1_seq, r1_qual), (r2_seq, r2_qual), member_count, strand1, strand2 in iterable:
            r1_seq, r1_qual = fix_strand(r1_seq, r1_qual, strand1)
            r2_seq, r2_qual = fix_strand(r2_seq, r2_qual, strand2)
            fastq_record_1 = '@%s_%i_member/1\n%s\n+\n%s\n' %(self.barcode, member_count, r1_seq, r1_qual)
            fastq_record_2 = '@%s_%i_member/2\n%s\n+\n%s\n' %(self.barcode, member_count, r2_seq, r2_qual)
            self.fastq_record += fastq_record_1 + fastq_record_2
