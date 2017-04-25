import pysam
from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
import re
import numpy as np
from numpy cimport ndarray
from cpython cimport bool


cdef bool qualify_aln(AlignedSegment aln):
    '''
    Check if alignment is properly mapped
    '''
    cdef:
        bool qualify
        int flag

    flag = aln.flag
    qualify = (flag == 99) or (flag == 147) or (flag == 163) or (flag == 83)
    return qualify


cpdef int split_bam(str in_bam, str outputprefix):

    cdef:
        AlignmentFile inbam, out_bam_1, out_bam_2
        AlignedSegment aln
        int count_R1 = 0, count_R2 = 0

    r1_out_bam_file = outputprefix + '_R1.bam'
    r2_out_bam_file = outputprefix + '_R2.bam'
    with pysam.AlignmentFile(in_bam,'rb') as inbam:
        with pysam.Samfile(r1_out_bam_file,'wb',template = inbam) as out_bam_1,\
                pysam.Samfile(r2_out_bam_file,'wb',template = inbam) as out_bam_2:
            for aln in inbam:
                if not aln.is_unmapped:
                    if aln.is_read1:
                        out_bam_1.write(aln)
                        count_R1 += 1
                    if aln.is_read2:
                        out_bam_2.write(aln)
                        count_R2 += 1
    print 'Done splitting'
    print 'Read 1: %i' %(count_R1)
    print 'Read 2: %i' %(count_R2)
    return 0

cdef ndarray split_cigar(cigar_string):
    '''
    split cigar string to numpy array
    return:
        [ [list of numbers],
          [list of cigar operators correspongs to the numbers] ]
    '''
    cdef:
        ndarray cigar_array

    cigar_numbers = re.findall(r'[0-9]+', cigar_string)
    cigar_operator = re.findall(r'[A-Z]', cigar_string)
    cigar_array = np.array([cigar_numbers, cigar_operator])
    return cigar_array

cdef bool check_aln(AlignedSegment aln, float single_end_thresh,
            float both_end_thresh):
    '''
    Compare soft clip length and threshold and write
    '''
    cdef:
        int total_clipped, seq_len
        ndarray cigar_array, all_soft_clipped
        int max_single_clipped
        float bet, set

    seq_len = len(aln.query_sequence)
    set = (single_end_thresh * seq_len)
    bet = (both_end_thresh * seq_len)
    cigar_array = split_cigar(str(aln.cigarstring))
    all_soft_clipped = np.array(cigar_array[0][cigar_array[1]=='S'],dtype=np.int16)
    total_clipped = all_soft_clipped.sum()
    max_single_clipped =  all_soft_clipped.max()
    single_pass = abs(max_single_clipped) <  set
    both_pass =  abs(total_clipped) < bet
    return single_pass and both_pass and not aln.is_unmapped

cpdef int filter_bam(str in_bam, str out_bam, float single_end_thresh,
                    float both_end_thresh, bool inverse):
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

    with pysam.AlignmentFile(in_bam,'rb') as inbam:
        with pysam.Samfile(out_bam,'wb',template = inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                flag_qualify_ok = qualify_aln(aln)
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
