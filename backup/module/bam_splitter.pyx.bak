from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import pysam
import re
import numpy as np
from numpy cimport ndarray
from cpython cimport bool
from tgirt_seq_tools.bam_tools import concordant_alignment, split_cigar

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
        float bet, set

    seq_len = len(aln.query_sequence)
    set = (single_end_thresh * seq_len)
    bet = (both_end_thresh * seq_len)
    cigar_array = split_cigar(aln.cigarstring)
    all_soft_clipped = [n for n, c in zip(*cigar_array) if c =='S']
    total_clipped = sum(all_soft_clipped)
    max_single_clipped =  max(all_soft_clipped)
    single_pass = abs(max_single_clipped) <  set
    both_pass =  abs(total_clipped) < bet
    return single_pass and both_pass and not aln.is_unmapped

def filter_bam(in_bam, out_bam, single_end_thresh,
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
                    print 'Parsed %i alignments' %(aln_count)
    return output_count
