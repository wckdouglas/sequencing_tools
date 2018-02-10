from __future__ import print_function
import sys
import re
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool
from functools import partial
from builtins import zip
from sequencing_tools.bam_tools import split_cigar
from sequencing_tools.bam_tools.fragment_pairs import concordant_pairs, is_split_pair

cdef int softClipSize(AlignedSegment aln):
    # compare softclip size, output maximum softclipped base on either side
    cdef:
        str cigar = str(aln.cigarstring)
        str cigar_char
        int cigar_num

    iterator = zip(*split_cigar(cigar))
    clipped = [cigar_num for cigar_num, cigar_char in iterator if cigar_char=='S']
    return max(clipped) if clipped else 1


cdef bool bowtie2_is_unique(AlignedSegment read1, AlignedSegment read2):
    '''
    bowtie2 output mapq as 255 for uniquely mapped
    '''
    return read1.mapq == 255 and read2.mapq == 255


cpdef bool hisat2_is_unique(AlignedSegment read1, AlignedSegment read2):
    '''
    hisat2 output NH==1 for uniquely mapped
    '''
    return read1.get_tag('NH') == 1 and read2.get_tag('NH') == 1


cpdef bool both_read_good_clip(AlignedSegment read1, AlignedSegment read2):
    '''
    check if both reads has clipped base < 20
    '''
    cdef:
        int read1_clipped = softClipSize(read1)
        int read2_clipped = softClipSize(read2)

    return read1_clipped < 20 and read2_clipped < 20


cpdef int split_bam_pair(AlignmentFile bam, AlignmentFile uniquebam, AlignmentFile multibam, aligner):

    cdef:
        int pair_count = 0
        int uniq_written = 0
        int multi_written = 0
        AlignedSegment read1, read2

    is_unique = partial(bowtie2_is_unique) if aligner == 'bowtie2' else partial(hisat2_is_unique)
    while True:
        try:
            read1 = next(bam)
            read2 = next(bam)
            pair_count += 1
            assert read1.query_name == read2.query_name, 'Not paired end'
            if not read1.is_unmapped and not read2.is_unmapped:
                if both_read_good_clip(read1, read2):
                    if is_unique(read1, read2):
                        uniquebam.write(read1)
                        uniquebam.write(read2)
                        uniq_written += 1
                    else:
                        multibam.write(read1)
                        multibam.write(read2)
                        multi_written += 1
            if pair_count % 5000000 == 0:
                print('Parsed %i alignment pairs' %pair_count, file = sys.stderr)
        except StopIteration:
            break
    print('Written %i to uniq bam, %i to multi bam' %(uniq_written, multi_written), file = sys.stderr)
    return 0

cpdef int split_bam_single(AlignmentFile bam, AlignmentFile uniquebam, AlignmentFile multibam, aligner):

    cdef:
        int count = 0
        int uniq_written = 0
        int multi_written = 0
        AlignedSegment read

    is_unique = partial(bowtie2_is_unique) if aligner == 'bowtie2' else partial(hisat2_is_unique)
    for count, read in enumerate(bam):
        if not read.is_unmapped:
            if (aligner == 'bowtie2' and read.mapq==255) or (aligner=='hisat2' and read.get_tag('NH') == 1):
                uniquebam.write(read)
                uniq_written += 1
            else:
                multibam.write(read)
                multi_written += 1
        if count % 5000000 == 0:
            print('Parsed %i alignment pairs' %count, file = sys.stderr)
    print('Written %i to uniq bam, %i to multi bam' %(uniq_written, multi_written), file = sys.stderr)
    return 0

cpdef int split_N_bam(AlignmentFile inbam,
                    AlignmentFile out_split_bam,
                    AlignmentFile out_fragment_bam):
    cdef:
        int split_out = 0
        int out = 0
        AlignedSegment read1, read2

    while True:
        try:
            read1 = next(inbam)
            read2 = next(inbam)
            if concordant_pairs(read1, read2):
                if is_split_pair(read1, read2):
                    out_split_bam.write(read1)
                    out_split_bam.write(read2)
                    split_out += 2
                else:
                    out_fragment_bam.write(read1)
                    out_fragment_bam.write(read2)
                    out += 2
        except StopIteration:
            break
    print('Written %i unsplit and %i split alignments' %(out, split_out), file = sys.stderr)
    return 0
 
