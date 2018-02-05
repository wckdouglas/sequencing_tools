from __future__ import print_function
from itertools import product
from builtins import zip, map
from multiprocessing import Pool
from functools import partial
import numpy as np
import gzip
cimport numpy as np
import sys
from cpython cimport bool
import io
import os
from sequencing_tools.fastq_tools import readfq, gzopen, reverse_complement
from sequencing_tools.fastq_tools._fastq_tools cimport fastqRecord
from sequencing_tools.fastq_tools.cutadapt_align import locate


cpdef int levenshtein_distance(str s1, str s2):
    '''
    Calculating Levenshtein distance from two strings
    algorithm from: http://rosettacode.org/wiki/Levenshtein_distance#Python

    usage: levenshtein_distance(string1, string2)
    ==============================
    Parameter:

    string1
    string2

    return:
    edit distance: the edit distance between two string

    '''

    cdef:
        list distance, new_distances
        int index1, index2
        str char1, char2
        int len_s1 = len(s1)
        int len_s2 = len(s2)


    if len_s1 > len_s2:
        s1, s2 = s2, s1
    distances = range(len_s1 + 1)
    for index2, char2 in enumerate(s2):
        new_distances = [index2+1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                new_distances.append(distances[index1])
            else:
                new_distances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             new_distances[-1])))
        distances = new_distances
    return distances[-1]


cpdef int hamming_distance(str s1, str s2):
    '''
    Calculating hamming distance from two strings

    usage: hamming_distance(string1, string2)
    ==============================
    Parameter:

    string1
    string2

    has to be same length

    return:
    edit distance: the edit distance between two string
    ===============================
    '''

    cdef:
        str i, j
        int hamming = 0
    assert len(s1) == len(s2), 'Wrong barcode extraction'

    for i, j in zip(s1, s2):
        if i != j:
            hamming += 1

    return hamming

def trim_other_read(sequence, qual, barcode, adapter):

    cdef:
        int seq_start, seq_end, clip_start, clip_end, matched, error

    clip_seq = reverse_complement(barcode) + adapter
    located = locate(sequence ,clip_seq, 0.2)
    if located:
        seq_start, seq_end, clip_start, clip_end, matched, error = located
        if clip_start == 0:
            sequence, qual = sequence[:seq_start], qual[:seq_start]
    return sequence, qual


def clip_read1(barcode_cut_off, constant, constant_no_evaluation, prefix_split,
            idx_base, usable_seq, int hamming_threshold, str adapter,
            fastqRecord read1, fastqRecord read2):
    """
    from each read1, clipped barcode and constant regions and make it as ID
    only when 1. barcode quality > cutoff and
              2. constant region matched
    """

    cdef:
        float barcode_mean_qual
        bool no_N_barcode, hiQ_barcode, accurate_constant
        str prefix, seq_record, seq_left, qual_left, seq_right, qual_right


    seq_name = read1.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read1.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read1.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(list(map(ord, barcode_qual))) - 33
    prefix = barcode[:prefix_split]

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual > barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hamming_distance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_left = read1.seq[usable_seq:]
        qual_left = read1.qual[usable_seq:]
        seq_right, qual_right = trim_other_read(read2.seq, read2.qual, barcode, adapter)
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, seq_left, qual_left) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, seq_right, qual_right)
        return 1, prefix, seq_record
    else:
        return 0, None, None


def clip_read2(barcode_cut_off, constant, constant_no_evaluation, prefix_split,
            idx_base, usable_seq, int hamming_threshold, str adapter,
            fastqRecord read1, fastqRecord read2):
    """
    from each read1, clipped barcode and constant regions and make it as ID
    only when 1. barcode quality > cutoff and
              2. constant region matched
    """

    cdef:
        float barcode_mean_qual
        bool no_N_barcode, hiQ_barcode, accurate_constant
        str prefix, seq_record, seq_left, qual_left, seq_right, qual_right


    seq_name = read2.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read2.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read2.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(list(map(ord, barcode_qual))) - 33
    prefix = barcode[:prefix_split]

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual > barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hamming_distance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_right = read2.seq[usable_seq:]
        qual_right = read2.qual[usable_seq:]
        seq_left, qual_left = trim_other_read(read1.seq, read1.qual, barcode, adapter)
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, seq_left, qual_left) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, seq_right, qual_right)
        return 1, prefix, seq_record
    else:
        return 0, None, None


def open_files(output_prefix, prefix_split):
    '''
    making prefixed fastq files
    '''
    barcode_prefix = product(['A','C','T','G'],repeat=prefix_split)
    barcode_prefix = [''.join(prefix) for prefix in barcode_prefix]
    if prefix_split > 0:
        file_dict = {prefix: open('%s_%s.fq' %(output_prefix, prefix), 'w') for prefix in barcode_prefix}
    else:
        file_dict = {prefix: open('%s.fq' %(output_prefix), 'w') for prefix in barcode_prefix}
    return file_dict


def clip_funtion(read, barcode_cut_off, constant,
                constant_no_evaluation, prefix_split, idx_base,
                usable_seq, hamming_threshold):
    if read == 'read1':
        adapter = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
        clipping = partial(clip_read1, barcode_cut_off, constant,
                        constant_no_evaluation, prefix_split, idx_base,
                        usable_seq, hamming_threshold, adapter)

    elif read == 'read2':
        adapter = 'AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
        clipping = partial(clip_read2, barcode_cut_off, constant,
                        constant_no_evaluation, prefix_split, idx_base,
                        usable_seq, hamming_threshold, adapter)

    return clipping




def run_pairs(outputprefix, inFastq1, inFastq2, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname,
            prefix_split, read):
    '''
    Feed in paried end fastq file
    loop over pairs and generate interleaved fastq
    '''
    cdef:
        int usable_seq
        int hamming_threshold
        str out_R1, out_R2
        int out, count, out_count = 0
        fastqRecord read1, read2
        str prefix, seq_record

    constant_length = len(constant)
    constant_no_evaluation = constant_length < 1
    hamming_threshold = allow_mismatch if not constant_no_evaluation else 0
    constant_no_evaluation = constant_length < 2
    usable_seq = idx_base if constant_no_evaluation else idx_base + constant_length


    file_dict = open_files(outputprefix, prefix_split)
    read_flag1 = 'rb' if inFastq1.endswith('.gz') else 'r' 
    read_flag2 = 'rb' if inFastq2.endswith('.gz') else 'r' 
    with gzopen(inFastq1, read_flag = read_flag1) as in1, gzopen(inFastq2, read_flag = read_flag2) as in2:
        clipping = clip_funtion(read, barcode_cut_off, constant,
                        constant_no_evaluation, prefix_split, idx_base,
                        usable_seq, hamming_threshold)

        for count, (read1, read2) in enumerate(zip(readfq(in1), readfq(in2))):
            out, prefix, seq_record = clipping(read1, read2)
            out_count += out
            if out == 1:
                print(seq_record, file = file_dict[prefix])
            if count % 10000000 == 0 and count != 0:
                print('[%s] Parsed %i records'%(programname, count), file = sys.stderr)
    [f.close() for f in file_dict.values()]

    print('[%s] Parsed:           %i sequences' %(programname, count), file = sys.stderr)
    print('[%s] Output:           %i sequences' %(programname, out_count), file = sys.stderr)
    return 0


def run_pairs_stdout(inFastq1, inFastq2, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname,
            prefix_split, read):
    '''
    Feed in paried end fastq file
    loop over pairs and generate interleaved fastq
    '''
    cdef:
        int usable_seq
        int hamming_threshold
        str out_R1, out_R2
        int out, count, out_count = 0
        fastqRecord read1, read2
        str prefix, seq_record

    constant_length = len(constant)
    constant_no_evaluation = constant_length < 1
    hamming_threshold = allow_mismatch if not constant_no_evaluation else 0
    usable_seq = idx_base if constant_no_evaluation else idx_base + constant_length


    read_flag1 = 'rb' if inFastq1.endswith('.gz') else 'r' 
    read_flag2 = 'rb' if inFastq2.endswith('.gz') else 'r' 
    with gzopen(inFastq1, read_flag = read_flag1) as in1, gzopen(inFastq2, read_flag = read_flag2) as in2:
        clipping = clip_funtion(read, barcode_cut_off, constant,
                constant_no_evaluation, prefix_split, idx_base,
                usable_seq, hamming_threshold)

        for count, (read1, read2) in enumerate(zip(readfq(in1), readfq(in2))):
            out, prefix, seq_record = clipping(read1, read2)
            out_count += out
            if out == 1:
                print(seq_record, file = sys.stdout)
            if count % 10000000 == 0 and count != 0:
                print('[%s] Parsed %i records'%(programname, count), file = sys.stderr)

    print('[%s] Parsed:           %i sequences' %(programname, count), file = sys.stderr)
    print('[%s] Output:           %i sequences' %(programname, out_count), file = sys.stderr)
    return 0
