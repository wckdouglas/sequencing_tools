from __future__ import print_function
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import pysam
import re
import sys
import numpy as np
from numpy cimport ndarray
from cpython cimport bool
from sequencing_tools.bam_tools import concordant_alignment, split_cigar, concordant_pairs


cdef class fragment_pairs:
    '''
    A data structure for pair of reads
    '''
    def __init__(self, AlignedSegment read1, AlignedSegment read2):
        self.read1 = read1
        self.read2 = read2
        self.pass_clip_check = None 
        self.flag_qualify_ok = None
        self.has_soft_clip = None
        self.read1_has_S = None
        self.read2_has_S = None

    def check_flags(self):
        self.flag_qualify_ok = concordant_pairs(self.read1, self.read2)

    def check_soft_clips(self):
        self.read1_has_S = 'S' in self.read1.cigarstring 
        self.read2_has_S = 'S' in self.read2.cigarstring
        self.has_soft_clip = self.read1_has_S or self.read2_has_S
    
    def check_pair_clips(self, s_et, b_et):
        read1_check = True
        read2_check = True
        if self.read1_has_S:
            read1_check = check_aln(self.read1, s_et, b_et)
        if self.read2_has_S:
            read2_check = check_aln(self.read2, s_et, b_et)
        self.pass_clip_check = read1_check and read2_check
        

    def output_aln(self, AlignmentFile out_bam_handle):
        out_bam_handle.write(self.read1)
        out_bam_handle.write(self.read2)


cpdef bool check_aln(AlignedSegment aln, float single_end_thresh,
            float both_end_thresh):
    '''
    Compare soft clip length and threshold and return boolean if alignment has softclipped < threshold

    input: aln, single_end_threshold, both end threshold
    return: boolean if alignment soft clip is below threshold
    '''
    cdef:
        int total_clipped, seq_len
        int max_single_clipped
        float b_et, s_et

    seq_len = len(aln.query_sequence)
    s_et = (single_end_thresh * seq_len)
    b_et = (both_end_thresh * seq_len)
    cigar_array = split_cigar(aln.cigarstring)
    all_soft_clipped = [n for n, c in zip(*cigar_array) if c =='S']
    total_clipped = sum(all_soft_clipped)
    max_single_clipped =  max(all_soft_clipped)
    single_pass = abs(max_single_clipped) <  s_et
    both_pass =  abs(total_clipped) < b_et
    return single_pass and both_pass and not aln.is_unmapped

def filter_bam_single_end(in_bam, out_bam, single_end_thresh,
                    both_end_thresh, inverse):
    '''
    This function filter softclipped sequence by tresholds regarding to the sequence length
    '''
    cdef:
        AlignmentFile inbam, outbam
        AlignedSegment aln
        ndarray cigar_array, all_soft_clipped
        int aln_count
        int output_count
        bool flag_qualify_ok, soft_clipped, clipped_size_right
        bool inverse_ok, non_inverse_ok

    with pysam.Samfile(in_bam,'rb') as inbam:
        with pysam.Samfile(out_bam,'wb',template = inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                flag_qualify_ok = concordant_alignment(aln)
                if not aln.is_unmapped and flag_qualify_ok:
                    soft_clipped = 'S' in aln.cigarstring
                    if not soft_clipped and not inverse:
                        outbam.write(aln)
                        output_count += 1

                    elif soft_clipped:
                        clipped_size_right = check_aln(aln, single_end_thresh, both_end_thresh)

                        inverse_ok = (not clipped_size_right and inverse)
                        non_inverse_ok = (clipped_size_right and not inverse)
                        if  inverse_ok or non_inverse_ok:
                            outbam.write(aln)
                            output_count += 1

                if aln_count % 1000000 == 0 and aln_count != 0:
                    print('Parsed %i alignments' %(aln_count), file = sys.stdout)
    return output_count

def filter_bam_pair_end(in_bam, out_bam, single_end_thresh,
                    both_end_thresh, inverse):
    '''
    This function filter softclipped sequence by tresholds regarding to the sequence length
    '''
    cdef:
        AlignmentFile inbam, outbam
        AlignedSegment read1, read2
        ndarray cigar_array, all_soft_clipped
        int pair_count = 0
        int output_count
        bool flag_qualify_ok, soft_clipped, clipped_size_right
        bool inverse_ok, non_inverse_ok
        fragment_pairs pairs

    with pysam.Samfile(in_bam,'rb') as inbam:
        with pysam.Samfile(out_bam,'wb',template = inbam) as outbam:
            while True:
                try:
                    read1 = inbam.next()
                    read2 = inbam.next()
                    pairs = fragment_pairs(read1, read2)
                    pairs.check_flags()
                    assert read1.query_name == read2.query_name, 'Wrong pairs: %s, %s' %(read1.query_name, read2.query_name)
                    pair_count += 1

                    if pairs.flag_qualify_ok:
                        pairs.check_soft_clips()
                        if not pairs.has_soft_clip and not inverse:
                            pairs.output_aln(outbam)
                            output_count += 1

                        elif pairs.has_soft_clip:
                            pairs.check_pair_clips(single_end_thresh, both_end_thresh)


                            inverse_ok = (not pairs.pass_clip_check  and inverse)
                            non_inverse_ok = (pairs.pass_clip_check and not inverse)
                            if  inverse_ok or non_inverse_ok:
                                pairs.output_aln(outbam)
                                output_count += 1

                    if pair_count % 1000000 == 0 and pair_count != 0:
                        print('Parsed %i alignments' %(pair_count), file = sys.stdout)
                except StopIteration:
                    break
    return output_count
