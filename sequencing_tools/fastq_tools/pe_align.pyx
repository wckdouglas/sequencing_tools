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
import logging 
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('PE align')
cdef:
    double EPSILON = 0.999999


class ConsensusBuilder:
    """
    Merging the two reads from paired end sequencing, and build consensus sequence for the overlapped part, paired end sequencing of a fragments:

    Illustration::

        Reads:
                       | adapter|                            |adapter |
                    R1:         |------------------>
          DNA Fragment: ========|----------------------------|========
                    R2:                     <----------------|
            report all:         |---------------------------->
        not report all:                     |------>
    

    Args:
        error_toleration (float): maximum mismatch fraction at the overlap region, reads are merged if they have mismatch fraction lower than this
        min_len (int): minumum overlapping bases at the overlapping region, if lower than this, reads won't merge
        report_all (boolean): True if the output merged sequence contains non-overlapping bases

    Example::

        consensus_builder = ConsensusBuilder(error_toleration = 0.1,
                                    min_len = 15, report_all=False)
        read1 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG',
            'ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAAAGATCGGAAGAGCACACGTCTGAA',
            'A6AA6//EA/AEE/AEAE///EEE/AAA/6/EE/A//EEE/A//EE/AE//E/////AEE///////////////',
            'NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG    Read1',
            )

        read2 = fastqRecord('NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG',
            'TTTCCCGCACCGATAGCAGCCGCTCTGGTCTCCCATCCCGGGCAATTGTGTGATCGTCGGACTGTAGAACTCTGA',
            'AAA//E/E/E/A/A///EEE/E///<//EE/EE6/</EEA/A</EE<AAE/E/</<<AAE/E<///EE/EA///A',
            'NB501060:148:HNFYCBGX5:1:11101:10036:1116 2:N:0:GAGTGG    Read2',
            )

        out = consensus_builder.run(read1, read2)
        print(out)                                                                                                                                                                     
        # @NB501060:148:HNFYCBGX5:1:11101:10036:1116 1:N:0:GAGTGG
        # ACACAATTGCCCGGGATGGGAGACCAGAGCGGCTGCTATCGGTGCGGGAAA
        # +
        # IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    """
    def __init__(self, error_toleration = 0.1, min_len = 15, 
                report_all=False, conserved = False, highlight=False):

        self.error_toleration = error_toleration #: maximum fraction of error in the overlapped region 
        self.min_len = min_len #: minimum overlapping positions to be considered as overlapped
        self.report_all = report_all #: reporting the non-overlapping bases as well?
    
        mode = 'vote' if conserved else 'prob'
        self.correction_module = ErrorCorrection(mode = mode)
        self.highlight = self.report_all and highlight
                
    
    def run(self, fastqRecord R1, fastqRecord R2):
        """
        reverse complement read2 sequence and find matching position on read1

        Args:
            R1: :class:`sequencing_tools.fastq_tools._fastq_tools.fastqRecord`
            R2: :class:`sequencing_tools.fastq_tools._fastq_tools.fastqRecord`

        Returns: 
            str: fastq line for the merged concensus sequence [+ non-overlapping sequence]
        """
        
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
                left_adapter = False
                right_adapter = False
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
                        left_adapter = True

                    if r1_end != len(R1.seq):
                        '''
                        R1:                 ----------------------->
                        R2:      <-------------------------|=======|
                        '''
                        right_add_seq = R1.seq[r1_end:] 
                        right_add_qual = R1.qual[r1_end:]
                        right_adapter = True

                    if self.highlight:
                        if left_adapter:
                            left_add_seq = self.__highlight__(left_add_seq)
                        if right_adapter:
                            right_add_seq = self.__highlight__(right_add_seq)
                        if not left_adapter and not right_adapter:
                            return ''

                    seq = left_add_seq + seq + right_add_seq
                    qual = left_add_qual + qual + right_add_qual

                out_line = '@%s\n%s\n+\n%s' %(r1_id,seq, qual) 
        return out_line
    
    def __highlight__(self, string):
        return '\x1b[6;30;42m' + string + '\x1b[0m'



def merge_interleaved(infile, outfile_handle, min_len, error_toleration, report_all, conserved=False, highlight=False):
    cdef:
        fastqRecord R1, R2
        int record_count = 0
        int out_count = 0
        str out_line

    infile_handle = sys.stdin if infile == '-' or infile == '/dev/stdin' else xopen(infile,mode = 'r')

    concensus_builder = ConsensusBuilder(error_toleration, min_len, report_all, conserved, highlight)

    for R1, R2 in read_interleaved(infile_handle):
        out_line = concensus_builder.run(R1, R2)
        if out_line:
            out_count += 1 
            print(out_line, file=outfile_handle)
        record_count += 1

    logger.info('Parsed %i records' %(record_count))
    logger.info('Merged %i records' %(out_count))
