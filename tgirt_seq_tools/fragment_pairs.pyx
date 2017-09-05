#!/usr/bin/env python

from __future__ import print_function
import pysam
from operator import itemgetter
from cpython cimport bool
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import sys
from tgirt_seq_tools.bam_tools import concordant_pairs, read_ends, fragment_ends

cdef class read_fragment:
    cdef:
        long start, end

    def __init__(self, frag_start, frag_end):
        self.end = frag_end
        self.start = frag_start

    def fragment_ends(self):
        return self.start, self.end

    def start_site(self):
        return self.start

    def end_site(self):
        return self.end

def bam_to_bed(bam_file, out_file, int min_size, int max_size):
    '''
    Read two alignments at a time,
    assume they are pairs,
    make paired- fragments as bed line
    '''
    cdef:
        AlignmentFile in_bam
        AlignedSegment read_1, read_2
        str chrom, strand
        int fragment_size
        str line
        long start, end

    pair_count = 0
    with pysam.Samfile(bam_file,'rb') as in_bam:
        while True:
            try:
                read_1 = in_bam.next()
                read_2 = in_bam.next()
                assert read_1.query_name == read_2.query_name, 'Paired not stored together: %s, %s'  %(read_1.query_name , read_2.query_name)
                if concordant_pairs(read_1, read_2) and not read_1.is_duplicate and not read_2.is_duplicate:
                    chrom = read_1.reference_name
                    strand = '-' if read_1.is_reverse else '+'
                    start, end = fragment_ends(read_1, read_2)
                    fragment_size = end - start
                    if min_size < fragment_size < max_size:
                        line = '%s\t%i\t%i\t%s\t%i\t%s' %(chrom, start, end,
                                                        read_1.query_name,
                                                        fragment_size,strand)
                        pair_count += 1
                        print(line, file=out_file)
            except StopIteration:
                break
    sys.stderr.write('Witten %i fragments\n' %(pair_count))
    return 0


### bedpe_to_bed utils
cpdef str filterBed(str bedline, int min_length, int max_length):

    cdef:
        str chrom1, chrom2
        str start1, start2
        str end1, end2
        str strand1, strand2
        str name
        long start, end
        int length
        bool strandeness_correct, length_correct, chrom_correct
        str alignment = ''

    fields = bedline.lstrip().split('\t')
    chrom1,start1,end1 = fields[:3]
    chrom2, start2, end2, strand1, strand2 = itemgetter(3,4,5,8,9)(fields)
    name = bedline[6]
    start = min(map(long,[start1,start2]))
    end = max(map(long,[end1,end2]))
    length = end - start

    strandeness_correct = strand1 != strand2
    length_correct = min_length < length < max_length
    chrom_correct = chrom1 == chrom2

    if strandeness_correct and length_correct and chrom_correct:
        alignment = '\t'.join([chrom1, str(start), str(end), name, str(length), strand1])
    return alignment


cpdef int process_bedpe(bed_iterator, int min_length, int max_length, out_handle):
    cdef:
        str frag
        str aln

    for frag in bed_iterator:
        aln = filterBed(frag, min_length, max_length)
        if aln:
            print(aln, file = out_handle)
    return 0

cpdef bool is_split_pair(AlignedSegment read1, AlignedSegment read2):
    return 'N' in read1.cigarstring or 'N' in read2.cigarstring
