from __future__ import division, print_function
import numpy as np
import pysam
from builtins import zip, map
from libc.math cimport log10, exp, log
import sys
from cpython cimport bool
from scipy.special import logsumexp
from ..fastq_tools import reverse_complement
from six.moves import xrange
from functools import partial


cdef:
    double MIN_Q = 33.0
    double MAX_Q = 73.0
    double MAX_PROB = 0.99999999999
    double MIN_PROB = 1-MAX_PROB


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


cdef str phred_to_string(double phred):
    '''
    convert quality score to string
    '''
    cdef:
        double adjusted_phred = phred + 33
        int integer_phred
        str string_phred

    adjusted_phred = clip(adjusted_phred, MIN_Q, MAX_Q)
    integer_phred = int(adjusted_phred)
    string_phred =chr(integer_phred)
    return string_phred


cdef double error_prob_to_phred(double error_prob):
    '''
    probability to Phred score
    '''
    return -10 * log10(1- error_prob)


cpdef str prob_to_qual_string(double posterior):
    '''
    Input a list of probabilities and output quality string
    '''
    cdef:
        double q
        str qual_char

    # convert posterior to quality score
    q = clip(posterior, 0, MAX_PROB)
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

    for q in base_qual:
        yield 10**(-q/10)  


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


def calculatePosterior(column_bases, column_qualities, guess_base):
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
    hit = sum(log(1 - q) for q in qual_to_prob(qual_hit))
    missed = sum(log(q/3.0) for q in qual_to_prob(qual_missed))
    log_posterior = missed + hit
    return log_posterior


class ErrorCorrection():
    def __init__(self, mode = 'prob', threshold=0.8):
        '''
        a module for error correction in fastq records, included two modes:
            1. prob: using base quality as probability of errors as prior to calculate posterior quality, see:
                https://github.com/fulcrumgenomics/fgbio/wiki/Calling-Consensus-Reads
            2. vote: using a voting scheme as SafeSeq to generate consensus base, see:
                https://www.pnas.org/content/108/23/9530


        params:
            mode:  prob or vote
            threshold: only consider for "vote" mode, as a cutoff for returning a "N" if not enough fraction of bases agree


        example usage:
        
        ec = ErrorCorrection(mode='prob')
        ec.Correct(['AAACA','AAAAA','AAACA','AAACA'],
                    ['IIIII','IIIAI','FFFFF', 'FFFFF'])
        '''
        self.np_len = np.vectorize(len,otypes=[np.int32])
        self.np_ord = np.vectorize(ord, otypes=[np.int16])
        self.threshold = threshold

        assert(mode in {'prob','vote'})
        self.correction_mode = mode

        if self.correction_mode == 'prob':
            self.correction_function = partial(self.__posteriorConcensus__)
        else:
            self.correction_function = partial(self.__voteConcensus__)


    def __posteriorConcensus__(self, arg):
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

        column_bases, in_column_qualities, _ = arg
        column_qualities = self.np_ord(in_column_qualities) - 33
        possible_bases = np.unique(column_bases)
        number_possible_bases = len(possible_bases)
        if number_possible_bases == 1:
            concensus_base = possible_bases[0]
            posterior_correct_probability = MAX_PROB
        else:
            log_posteriors = [calculatePosterior(column_bases, column_qualities, guess_base) for guess_base in possible_bases]
            total_posterior = logsumexp(log_posteriors)
            likelihoods = [log_posterior - total_posterior for log_posterior in log_posteriors]
            arg_max_likelihood = np.argmax(likelihoods)
            concensus_base = possible_bases[arg_max_likelihood]
            posterior_correct_probability = exp(likelihoods[arg_max_likelihood])
        return concensus_base, posterior_correct_probability


    def __voteConcensus__(self, arg):
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
            double fraction_threshold
            int best_vote
            double best_fraction

        column_bases, in_column_qualities, fraction_threshold = arg
        column_qualities = self.np_ord(in_column_qualities) - 33
        possible_bases, possible_counts = np.unique(column_bases, return_counts=True)

        best_vote = possible_counts.max()
        best_fraction = best_vote/possible_counts.sum()
        if best_fraction < fraction_threshold:
            concensus_base = 'N'
            posterior_correct_probability = 0

        else:
            best_index = possible_counts==possible_counts.max()
            concensus_base = possible_bases[best_index][0]
            best_quals = column_qualities[column_bases==concensus_base].sum()
            posterior_correct_probability = min(MAX_PROB, 1 - 10**(-best_quals/10))

        return concensus_base, posterior_correct_probability


    def Correct(self, seq_list, qual_list):
        """given a list of sequences, a list of quality and sequence length.
            return a consensus sequence and a recalibrated quality

        params: 
            seq_list: list of sequences with the same length
            qual_list: list of the corresponding quality strings, must be in the same length as seq_list

        return: tuple
            consensus seq, 
            consensus qual

        """
        cdef:
            str sequence = ''
            str quality = ''
            double posterior_correct_prob
            int num_seq

        num_seq = len(seq_list)
        assert num_seq == len(qual_list), 'Number of sequences not match number of quality sequence'
        if len(seq_list) == 1:
            '''
            return the only sequence if only 1 input
            '''
            sequence = str(seq_list[0])
            quality = str(qual_list[0])
        else:
            seq_len = min(len(s) for s in seq_list)
            seq_list = map(lambda x: x[:seq_len], seq_list)
            qual_list = map(lambda x: x[:seq_len], qual_list)
            seq_list = np.array(list(map(list, seq_list)))
            qual_list = np.array(list(map(list, qual_list)))
            iter_list = ((seq_list[:,pos], qual_list[:,pos], self.threshold) for pos in xrange(seq_len))
            for base, posterior_correct_prob in map(self.correction_function, iter_list):
                sequence += base
                quality += prob_to_qual_string(posterior_correct_prob)
        return sequence, quality



class ConsensusAlignments():
    '''
    generate a consensus with bam alignments by counting aligned base at each position
    in cases where two bases had the same count, report 'N' 


    usage:

        ca = ConsensusAlignments(bam_file)
        ca.concensus_seq('chrM',1000,1100)

    args:
        bam: bam file path
        min_cov: minimum coverage to consider
    '''

    def __init__(self, bam, min_cov = 10):
        self.bam = pysam.Samfile(bam)
        assert self.bam.has_index(), 'BAM input should be indexed'
        self.decoder = {0:'A',1:'C',2:'G',3:'T',       # for decoding onehot encoded sequence matrix
                        4:'N',5:'N',6:'N',7:'N',8:'N'} # account for bases with same count

    def consensus_seq(self, chrom, start, end):
        #four array.arrays of the same length in order A C G T (onehot encode?)
        arr = self.bam.count_coverage(chrom, start, end) 
        arr = np.array(arr)
        arr = np.where(arr==arr.max(axis=0),1,0) # onehot encoded
        seq = np.matmul(np.arange(4),arr)
        seq = [self.decoder[i] for i in seq]
        seq = ''.join(seq)

        return seq