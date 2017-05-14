import re
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool
from functools import partial
from itertools import izip

cdef int softClipSize(AlignedSegment aln):
    # compare softclip size, output maximum softclipped base on either side
    cdef:
        str cigar = str(aln.cigarstring)
        str cigar_char
        str cigar_num

    cigar_nums = re.findall('[0-9]+',cigar)
    cigar_str = re.findall('[A-Z]',cigar)
    clipped = [int(cigar_num) for cigar_num, cigar_char in izip(cigar_nums, cigar_str) if cigar_char == 'S']
    return max(clipped) if clipped else 0


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



cpdef int split_bam(AlignmentFile bam, AlignmentFile uniquebam, AlignmentFile multibam, aligner):

    cdef:
        int pair_count = 0
        AlignedSegment read1, read2

    is_unique = partial(bowtie2_is_unique) if aligner == 'bowtie2' else partial(hisat2_is_unique)
    while True:
        try:
            read1 = bam.next()
            read2 = bam.next()
            pair_count += 1
            assert read1.qname == read2.qname, 'Not paired end'
            if not read1.is_unmapped and not read2.is_unmapped:
                if both_read_good_clip(read1, read2):
                    if is_unique(read1, read2):
                        uniquebam.write(read1)
                        uniquebam.write(read2)
                    else:
                        multibam.write(read1)
                        multibam.write(read2)
            if pair_count % 5000000 == 0:
                print 'Parsed %i alignment pairs' %pair_count
        except StopIteration:
            break
    return 0
