#!/usr/bin/env python

from __future__ import print_function
from tgirt_seq_tools.fastq_tools import readfq, reverse_complement
from tgirt_seq_tools.fastq_tools cimport fastqRecord
import sys
from tgirt_seq_tools.cutadapt_align import locate
from itertools import izip


cpdef calibrate_qual(str b1, str b2, str q1, str q2):
    '''
    https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md
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


def merge_interleaved(infile, outfile, min_len, error):
    cdef:
        fastqRecord R1, R2
        int record_count = 0
        int out_count = 0
        str seq, qual, r2_seq, r1_id, r2_id

    infile_handle = sys.stdin if infile == '-' or infile == '/dev/stdin' else open(infile,'r')
    outfile_handle = sys.stdout if outfile == '-' or outfile == '/dev/stdin' else open(outfile,'w')

    fastq_file = readfq(infile_handle)
    try:
        while True:
            R1 = fastq_file.next()
            R2 = fastq_file.next()
            r1_id, r2_id = R1.id.split('/')[0], R2.id.split('/')[0]
            assert r1_id == r2_id, 'Not interleaved'

            r2_seq = reverse_complement(R2.seq)
            aligned = locate(R1.seq, r2_seq, error)
            if aligned:
                r1_start, r1_end, r2_start, r2_end, match, error = aligned 
                if match > min_len:
    #                print(aligned, file=sys.stdout)
    #                print(R1.seq[r1_start:r1_end], file=sys.stdout)
    #                print(r2_seq[r2_start:r2_end], file=sys.stdout)
                    seq, qual = correct_error(R1.seq[r1_start:r1_end], 
                                            R1.qual[r1_start:r1_end], 
                                            r2_seq[r2_start:r2_end],
                                            R2.qual[::-1][r2_start:r2_end])
                    print('@%s\n%s\n+\n%s' %(r1_id,seq, qual), file=outfile_handle)
                    out_count += 1 
            record_count += 1

    except StopIteration:
        print('Parsed %i records' %(record_count), file=sys.stderr)
        print('Merged %i records' %(out_count), file=sys.stderr)

