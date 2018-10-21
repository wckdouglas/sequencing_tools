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
from sequencing_tools.fastq_tools._fastq_tools cimport fastqRecord
from sequencing_tools.fastq_tools import readfq, reverse_complement
from sequencing_tools.fastq_tools.cutadapt_align import locate
from sequencing_tools.io_tools import xopen
from sequencing_tools.stats_tools import hamming_distance



def insert_trimmer(str seq1, str seq2, str qual1, str qual2):
    '''
    aligned two reads, remove extensions
    only return overlapping regions
    '''
    cdef:
        int r1_start, r1_end, r2_start, r2_end, match
        double err

    seq2 = reverse_complement(seq2)
    qual2 = qual2[::-1]
    aligned = locate(seq1, seq2, 0.1)

    if aligned:
        r1_start, r1_end, r2_start, r2_end, match, err = aligned 
        if match >= 15:

            # read1
            seq1 = seq1[r1_start:r1_end]
            qual1 = qual1[r1_start:r1_end]

            # read2
            seq2 = seq2[r2_start:r2_end]
            qual2 = qual2[r2_start:r2_end]

    return seq1, reverse_complement(seq2), qual1, qual2[::-1]

def trim_other_read(sequence, qual, barcode, adapter):
    '''
    append UMI (barcode) onto adapter,
    and trim the read without UMI
    '''

    cdef:
        int seq_start, seq_end, clip_start, clip_end, matched, error

    clip_seq = reverse_complement(barcode) + adapter
    located = locate(sequence, clip_seq, 0.15)
    if located:
        seq_start, seq_end, clip_start, clip_end, matched, error = located
        if clip_start == 0 and matched >= 4:
            sequence, qual = sequence[:seq_start], qual[:seq_start]
    return sequence, qual


def clip_read1(barcode_cut_off, constant, constant_no_evaluation,
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
        str seq_record, seq_left, qual_left, seq_right, qual_right


    seq_name = read1.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read1.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read1.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(list(map(ord, barcode_qual))) - 33

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual >= barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hamming_distance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_left = read1.seq[usable_seq:]
        qual_left = read1.qual[usable_seq:]
        seq_right, qual_right = trim_other_read(read2.seq, read2.qual, barcode, adapter)
        seq_left, seq_right, qual_left, qual_right = insert_trimmer(seq_left, seq_right, qual_left, qual_right)
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, seq_left, qual_left) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, seq_right, qual_right)
        return 1, seq_record
    else:
        return 0, None


def clip_read2(barcode_cut_off, constant, constant_no_evaluation,
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
        str seq_record, seq_left, qual_left, seq_right, qual_right


    seq_name = read2.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read2.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read2.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(list(map(ord, barcode_qual))) - 33

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual >= barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hamming_distance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_right = read2.seq[usable_seq:]
        qual_right = read2.qual[usable_seq:]
        seq_left, qual_left = trim_other_read(read1.seq, read1.qual, barcode, adapter)
        seq_left, seq_right, qual_left, qual_right = insert_trimmer(seq_left, seq_right, qual_left, qual_right)
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, seq_left, qual_left) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, seq_right, qual_right)
        return 1, seq_record
    else:
        return 0, None


def clip_funtion(read, barcode_cut_off, constant,
                constant_no_evaluation, idx_base,
                usable_seq, hamming_threshold):
    if read == 'read1':
        adapter = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
        clipping = partial(clip_read1, barcode_cut_off, constant,
                        constant_no_evaluation, idx_base,
                        usable_seq, hamming_threshold, adapter)

    elif read == 'read2':
        adapter = 'AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
        clipping = partial(clip_read2, barcode_cut_off, constant,
                        constant_no_evaluation, idx_base,
                        usable_seq, hamming_threshold, adapter)

    return clipping


def clip_pairs(inFastq1, inFastq2, out_file, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname, read):
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
        str seq_record

    constant_length = len(constant)
    constant_no_evaluation = constant_length < 1
    hamming_threshold = allow_mismatch if not constant_no_evaluation else 0
    usable_seq = idx_base if constant_no_evaluation else idx_base + constant_length

    out_handle = sys.stdout if out_file in ['/dev/stdout','-'] else xopen(out_file, 'w')
    with xopen(inFastq1, mode = 'r') as in1, xopen(inFastq2, mode = 'r') as in2:
        clipping = clip_funtion(read, barcode_cut_off, constant,
                constant_no_evaluation, idx_base,
                usable_seq, hamming_threshold)

        for count, (read1, read2) in enumerate(zip(readfq(in1), readfq(in2))):
            out, seq_record = clipping(read1, read2)
            out_count += out
            if out == 1:
                print(seq_record, file = out_handle)
            if count % 10000000 == 0 and count != 0:
                print('[%s] Parsed %i records'%(programname, count), file = sys.stderr)

    print('[%s] Parsed:           %i sequences' %(programname, count), file = sys.stderr)
    print('[%s] Output:           %i sequences' %(programname, out_count), file = sys.stderr)
    return 0
