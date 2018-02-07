#!/usr/bin/env python

from __future__ import print_function
from sequencing_tools.fastq_tools import read_interleaved, reverse_complement
from sequencing_tools.fastq_tools._fastq_tools cimport fastqRecord
from sequencing_tools.fastq_tools.cutadapt_align import locate
from sequencing_tools.io_tools import xopen
import sys
from builtins import zip
from functools import partial
from cpython cimport bool


cdef calibrate_qual(str b1, str b2, str q1, str q2):
    '''
    https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md

    1. if same base, append base and sum of qual
    2. else use higher qual base
    3. if same, append N and qual of 0 
    '''
    cdef:
        int qual
        str base 
        int qq1, qq2

    qq1, qq2 = ord(q1) - 33 , ord(q2) - 33
    if b1==b2:
        qual = qq1 + qq2
        qual = 40 if qual > 40 else qual
        base = b1
    elif qq1 > qq2:
        base = b1 
        qual = qq1
    elif qq2 > qq1:
        base = b2 
        qual = qq2
    else:
        base = 'N'
        qual = 0
    return base, chr(qual + 33)

        

cdef correct_error(str r1_seq, str r1_qual, str r2_seq, str r2_qual):
    '''
    
    loop over overlapping region, 
    if same base, append base and sum of qual
    else append N and qual of 0 

    ==============================
    Parameter:
        r1_seq: overlaping region on seq1
        r1_qual: qual of r1_seq
        r2_seq: overlaping region on seq2
        r2_qual: qual of r2_Seq

    return:
        seq: concensus seq
        qual: qual for concensus seq
    ===============================
    '''

    cdef:
        str seq = ''
        str qual = ''
        str b1, b2, q1, q2

    iterator = zip(r1_seq, r1_qual, r2_seq, r2_qual)
    for b1, q1, b2, q2 in iterator:
        b, q = calibrate_qual(b1, b2, q1, q2)
        seq += b
        qual += q
    return seq, qual


cdef make_concensus(float error_toleration, int min_len, 
                bool report_all,
                fastqRecord R1, fastqRecord R2):
    '''
    reverse complement read2 sequence and find matching position on read1
    return concensus sequence
    '''
    
    cdef:
        str seq, qual, r2_seq, r1_id, r2_id
        str out_line = None
        str left_add_seq = ''
        str right_add_seq =''
        str left_add_qual = ''
        str right_add_qual = ''

    
    r1_id = R1.id.split('/')[0]

    r2_seq = reverse_complement(R2.seq)
    r2_qual = R2.qual[::-1]
    aligned = locate(R1.seq, r2_seq, error_toleration)
    if aligned:
        r1_start, r1_end, r2_start, r2_end, match, err = aligned 
        if match >= min_len:
#                print(aligned, file=sys.stdout)
#                print(R1.seq[r1_start:r1_end], file=sys.stdout)
#                print(r2_seq[r2_start:r2_end], file=sys.stdout)
            seq, qual = correct_error(R1.seq[r1_start:r1_end], 
                                    R1.qual[r1_start:r1_end], 
                                    r2_seq[r2_start:r2_end],
                                    r2_qual[r2_start:r2_end])
            if report_all:
                if r1_start != 0:
                    left_add_seq = R1.seq[:r1_start] 
                    left_add_qual = R1.qual[:r1_start]

                if r2_end != len(r2_seq):
                    right_add_seq = r2_seq[r2_end:]
                    right_add_qual = r2_qual[r2_end:]
                seq = left_add_seq + seq + right_add_seq
                qual = left_add_qual + qual + right_add_qual

            out_line = '@%s\n%s\n+\n%s' %(r1_id,seq, qual)
    return out_line


def merge_interleaved(infile, outfile_handle, min_len, error_toleration, report_all):
    cdef:
        fastqRecord R1, R2
        int record_count = 0
        int out_count = 0

    infile_handle = sys.stdin if infile == '-' or infile == '/dev/stdin' else xopen(infile,mode = 'r')

    concensus_builder = partial(make_concensus, error_toleration, min_len, report_all)

    for R1, R2 in read_interleaved(infile_handle):
        out_line = concensus_builder(R1, R2)
        if out_line:
            out_count += 1 
            print(out_line, file=outfile_handle)
        record_count += 1

    print('Parsed %i records' %(record_count), file=sys.stderr)
    print('Merged %i records' %(out_count), file=sys.stderr)
