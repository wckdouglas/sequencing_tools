from __future__ import print_function, division
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
from libc.math cimport fmin
from ._fastq_tools cimport fastqRecord
from ._fastq_tools import readfq, reverse_complement
from .cutadapt_align import locate
from ..io_tools import xopen
from ..stats_tools import hamming_distance
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


class adapters():
    R1 = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
    R2 = 'AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

class clip_read:
    '''
    implement read trimmer
    '''
    def __init__(self, 
                barcode_cut_off = 20, 
                constant = '', 
                constant_no_evaluation = True,
                umi_bases = 6, 
                usable_seq = 6, 
                hamming_threshold = 6, 
                adapter = adapters.R1, 
                min_length = 12):

        self.barcode_cut_off = barcode_cut_off
        self.constant = constant 
        self.constant_no_evaluation = constant_no_evaluation
        self.umi_bases = umi_bases 
        self.usable_seq = usable_seq
        self.hamming_threshold = hamming_threshold
        self.adapter = adapter 
        self.min_length = min_length
        self.template = '@{UMI}_{READNAME}/1\n{SEQ}\n+\n{QUAL}'

    def trim_reads(self, umi_read, opposite_read):
        """
        from each read1, clipped barcode and constant regions and make it as ID
        only when 1. barcode quality > cutoff and
                2. constant region matched
        """

        ret_code = 0
        umi_out_read = ''
        opposite_out_read = ''

        seq_name = umi_read.id.split(' ')[0]
        umi =  umi_read.seq[:self.umi_bases]

        umi_qual = umi_read.qual[:self.umi_bases]
        constant_region = umi_read.seq[self.umi_bases:self.usable_seq]
        barcode_mean_qual = sum(map(ord, umi_qual))/self.umi_bases - 33

        no_N_barcode = 'N' not in umi
        hiQ_barcode = barcode_mean_qual >= self.barcode_cut_off
        hamming_qual_pass = hamming_distance(self.constant, constant_region) <= self.hamming_threshold
        accurate_constant = self.constant_no_evaluation or hamming_qual_pass

        if no_N_barcode and hiQ_barcode and accurate_constant:


            umi_seq = umi_read.seq[self.usable_seq:]
            umi_qual = umi_read.qual[self.usable_seq:]

            opposite_seq, opposite_qual = self.__trim_other_read__(opposite_read.seq, opposite_read.qual, umi)
            umi_seq, opposite_seq, umi_qual, opposite_qual = self.__insert_trimmer__(umi_seq, opposite_seq, 
                                                                            umi_qual, opposite_qual)



            if fmin(len(umi_seq), len(opposite_seq)) >= self.min_length:
                umi_out_read = self.template.format(UMI = umi,
                                                    READNAME = seq_name,
                                                    SEQ = umi_seq,
                                                    QUAL = umi_qual)
                opposite_out_read = self.template.format(UMI = umi,
                                                    READNAME = seq_name,
                                                    SEQ = opposite_seq,
                                                    QUAL = opposite_qual)
                ret_code = 1
        return ret_code, umi_out_read, opposite_out_read


    def __insert_trimmer__(self, str seq1, str seq2, str qual1, str qual2):
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
            if match >= 15 and r1_start == 0 and r2_end == len(seq2):

                # read1
                seq1 = seq1[:r1_end]
                qual1 = qual1[r1_start:r1_end]

                # read2
                seq2 = seq2[r2_start:r2_end]
                qual2 = qual2[r2_start:r2_end]

        return seq1, reverse_complement(seq2), qual1, qual2[::-1]


    def __trim_other_read__(self, sequence, qual, umi):
        '''
        append UMI (barcode) onto adapter,
        and trim the read without UMI
        '''

        clip_seq = reverse_complement(umi) + self.adapter
        located = locate(sequence, clip_seq, 0.15)
        if located:
            seq_start, seq_end, clip_start, clip_end, matched, error = located
            if clip_start == 0 and matched > self.umi_bases:
                sequence, qual = sequence[:seq_start], qual[:seq_start]
        return sequence, qual


def clip_pairs(inFastq1, inFastq2, out_file, umi_bases,
            barcode_cut_off, constant, allow_mismatch, programname, read, min_length):
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
    usable_seq = umi_bases if constant_no_evaluation else umi_bases + constant_length

    out_handle = sys.stdout if out_file in ['/dev/stdout','-'] else xopen(out_file, 'w')
    with xopen(inFastq1, mode = 'r') as in1, xopen(inFastq2, mode = 'r') as in2:
        if read == 'read1':
            adapter = adapters().R1
            iterable = zip(readfq(in1), readfq(in2))

        else:
            adapter = adapters().R2
            iterable = zip(readfq(in2), readfq(in1))

        clipping = clip_read(barcode_cut_off, constant,
                    constant_no_evaluation, umi_bases,
                    usable_seq, hamming_threshold, adapter, 
                    min_length)

        for count, (umi_read, opposite_read) in enumerate(iterable):
            out, umi_read, opposite_read = clipping.trim_reads(umi_read, opposite_read)
            out_count += out
            if out == 1:
                if read == "read1":
                    print(umi_read +'\n' + opposite_read, file = out_handle)
                if read == "read2":
                    print(opposite_read +'\n' + umi_read, file = out_handle)

            if count % 10000000 == 0 and count != 0:
                logger.info('[%s] Parsed %i records'%(programname, count))

    logger.info('[%s] Parsed:           %i sequences' %(programname, count))
    logger.info('[%s] Output:           %i sequences' %(programname, out_count))
    return 0
