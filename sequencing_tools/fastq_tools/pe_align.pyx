from __future__ import print_function
import sys
from builtins import zip
from functools import partial
from cpython cimport bool
from ._fastq_tools import read_interleaved, reverse_complement
from ._fastq_tools cimport fastqRecord
from .cutadapt_align import locate
from ..io_tools import xopen
from ..consensus_tools import ErrorCorrection
import numpy as np
cdef:
    double EPSILON = 0.999999


class ConsensusBuilder:
    def __init__(self, error_toleration = 0.1, min_len = 15, 
                report_all=False, conserved = False):
        self.error_toleration = error_toleration
        self.min_len = min_len
        self.report_all = report_all
    
        mode = 'vote' if conserved else 'prob'
        self.correction_module = ErrorCorrection(mode = mode)
                
    
    def run(self, fastqRecord R1, fastqRecord R2):
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
            bool no_indel

        
        r1_id = R1.id.split('/')[0]
        r2_seq = reverse_complement(R2.seq)
        r2_qual = R2.qual[::-1]
        aligned = locate(R1.seq, r2_seq, self.error_toleration)
        if aligned:
            r1_start, r1_end, r2_start, r2_end, match, err = aligned 
            no_indel =(r1_end - r1_start) == (r2_end - r2_start)
            if match >= self.min_len and no_indel:
    #                print(aligned, file=sys.stdout)
    #                print(R1.seq[r1_start:r1_end], file=sys.stdout)
    #                print(r2_seq[r2_start:r2_end], file=sys.stdout)
                seq, qual = self.correction_module.Correct(
                                        [R1.seq[r1_start:r1_end], 
                                         r2_seq[r2_start:r2_end]],
                                        [R1.qual[r1_start:r1_end],
                                         r2_qual[r2_start:r2_end]])
                if self.report_all:
                    if r2_end != len(R2.seq):
                        '''
                        R1:       --------------------->|======|
                        R2:           <-------------------------
                        '''
                        right_add_seq = r2_seq[r2_end:] 
                        right_add_qual = r2_qual[r2_end:]

                    if r1_start !=0:
                        '''
                        R1:       ----------------------->
                        R2:       |=====|<-------------------------
                        '''
                        left_add_seq = R1.seq[:r1_start]
                        left_add_qual = R1.qual[:r1_start]
                    
                    if r2_start != 0:
                        '''
                        R1:      |=========|----------------------->
                        R2:      <-------------------------
                        '''
                        left_add_seq = r2_seq[:r2_start]
                        left_add_qual = r2_qual[:r2_start]

                    if r1_end != len(R1.seq):
                        '''
                        R1:                 ----------------------->
                        R2:      <-------------------------|=======|
                        '''
                        right_add_seq = R1.seq[r1_end:] 
                        right_add_qual = R1.qual[r1_end:]

                seq = left_add_seq + seq + right_add_seq
                qual = left_add_qual + qual + right_add_qual

                out_line = '@%s\n%s\n+\n%s' %(r1_id,seq, qual)
        return out_line


def merge_interleaved(infile, outfile_handle, min_len, error_toleration, report_all, conserved=False):
    cdef:
        fastqRecord R1, R2
        int record_count = 0
        int out_count = 0
        str out_line

    infile_handle = sys.stdin if infile == '-' or infile == '/dev/stdin' else xopen(infile,mode = 'r')

    concensus_builder = ConsensusBuilder(error_toleration, min_len, report_all, conserved)

    for R1, R2 in read_interleaved(infile_handle):
        out_line = concensus_builder.run(R1, R2)
        if out_line:
            out_count += 1 
            print(out_line, file=outfile_handle)
        record_count += 1

    print('Parsed %i records' %(record_count), file=sys.stderr)
    print('Merged %i records' %(out_count), file=sys.stderr)
