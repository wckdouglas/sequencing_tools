#!/usr/bin/env python

from __future__ import print_function
from tgirt_seq_tools.fastq_tools import readfq, reverse_complement
from tgirt_seq_tools.fastq_tools cimport fastqRecord
import sys
from tgirt_seq_tools.cutadapt_align import locate
from itertools import izip


cdef calibrate_qual(str b1, str b2, str q1, str q2):
    '''
    https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md

    if same base, append base and sum of qual
    else append N and qual of 0 
    '''
    cdef:
        int qual
        str base 

    if b1==b2:
        qual = ord(q1) + ord(q2) - 66
        qual = 40 if qual > 40 else qual
        base = b1
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

    iterator = izip(r1_seq, r1_qual, r2_seq, r2_qual)
    for b1, q1, b2, q2 in iterator:
        b, q = calibrate_qual(b1, b2, q1, q2)
        seq += b
        qual += q
    return seq, qual

cdef make_concensus(fastqRecord R1, fastqRecord R2, 
            float error_toleration, int min_len):
    
    cdef:
        str seq, qual, r2_seq, r1_id, r2_id
        str out_line = None
    
    r1_id, r2_id = R1.id.split('/')[0], R2.id.split('/')[0]
    assert r1_id == r2_id, 'Not interleaved'

    r2_seq = reverse_complement(R2.seq)
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
                                    R2.qual[::-1][r2_start:r2_end])
            out_line = '@%s\n%s\n+\n%s' %(r1_id,seq, qual)
    return out_line


def merge_interleaved(infile, outfile_handle, min_len, error_toleration):
    cdef:
        fastqRecord R1, R2
        int record_count = 0
        int out_count = 0

    infile_handle = sys.stdin if infile == '-' or infile == '/dev/stdin' else open(infile,'r')

    fastq_file = readfq(infile_handle)
    try:
        while True:
            R1 = fastq_file.next()
            R2 = fastq_file.next()
            out_line = make_concensus(R1, R2, error_toleration, min_len)
            if out_line:
                out_count += 1 
                print(out_line, file=outfile_handle)
            record_count += 1

    except StopIteration:
        print('Parsed %i records' %(record_count), file=sys.stderr)
        print('Merged %i records' %(out_count), file=sys.stderr)

