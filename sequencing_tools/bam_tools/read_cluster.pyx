from __future__ import division, print_function

import logging
import os
import sys
from builtins import map, zip

import numpy as np
import pysam
from scipy.special import logsumexp
from six.moves import xrange

from libc.math cimport log10, exp, log
from cpython cimport bool

from sequencing_tools.consensus_tools import ErrorCorrection
from sequencing_tools.fastq_tools import reverse_complement
from sequencing_tools.utils import SeqUtilsError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))




def fix_strand(str seq, str qual, bool strand):
    if strand:
        seq = reverse_complement(seq)
        qual = qual[::-1]
    return seq, qual

class readGroup:
    '''
    Read group object
    '''
    def __init__(self,aln, tag, conserved):
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
    
        correction_mode = 'prob' if not conserved else 'vote'
        self.correction_module = ErrorCorrection(mode = correction_mode,
                                               threshold=0.8)

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


    def cluster(self):
        '''
            from read group, generate concensus sequence, quality
        '''

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

        if not (self.R2 and self.R1):
            raise SeqUtilsError('No input: {} {} {}'.format(self.R1, self.R2, self.barcode))

        if  R2_array.shape != R1_array.shape:
            raise SeqUtilsError('Unequal R1 list vs R2')

        for _chrom1, _chrom2, _pos1, _pos2, _R1_flag, _R2_flag in iterator:
            chrom_is_right = (R1_chrom_array == _chrom1) & (R2_chrom_array == _chrom2)
            flag_is_right = (R1_flag_array == _R1_flag) & (R2_flag_array == _R2_flag)
            pos_is_right = (R1_pos_array == _pos1) & (R2_pos_array == _pos2)

            cluster = chrom_is_right & flag_is_right & pos_is_right
            R1_filtered, R2_filtered = R1_array[cluster,:], R2_array[cluster,:]
            self.concensus_read1.append(self.correction_module.Correct(R1_filtered[:,0], R1_filtered[:,1]))
            self.concensus_read2.append(self.correction_module.Correct(R2_filtered[:,0], R2_filtered[:,1]))
            self.member_count_list.append(R1_filtered.shape[0])
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

        if not self.concensus_read1:
            raise SeqUtilsError('No concensus read generated')
        iterable = zip(self.concensus_read1, self.concensus_read2,
                        self.member_count_list,
                        self.concensus_flag1, self.concensus_flag2)
        for (r1_seq, r1_qual), (r2_seq, r2_qual), member_count, strand1, strand2 in iterable:
            r1_seq, r1_qual = fix_strand(r1_seq, r1_qual, strand1)
            r2_seq, r2_qual = fix_strand(r2_seq, r2_qual, strand2)
            fastq_record_1 = '@%s_%i_member/1\n%s\n+\n%s\n' %(self.barcode, member_count, r1_seq, r1_qual)
            fastq_record_2 = '@%s_%i_member/2\n%s\n+\n%s\n' %(self.barcode, member_count, r2_seq, r2_qual)
            self.fastq_record += fastq_record_1 + fastq_record_2
