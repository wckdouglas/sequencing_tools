#!/usr/bin/env python

from __future__ import print_function
import pysam
from operator import itemgetter
from cpython cimport bool
from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
import sys

cdef bool qualify_pairs(AlignedSegment read1, AlignedSegment read2):
    cdef:
        bool reverse_fragment = read1.flag == 83 and read2.flag == 163
        bool forward_fragment = read1.flag == 99 and read2.flag == 147
    return reverse_fragment or forward_fragment

def read_ends(AlignedSegment read):
    positions = read.get_reference_positions()
    start, end = itemgetter(0,-1)(positions)
    return start, end

def fragment_ends(AlignedSegment read1, AlignedSegment read2):
    cdef:
        long start1, end1, start2, end2
        long start, end

    start1, end1 = read_ends(read1)
    start2, end2 = read_ends(read2)
    start = min(start1, start2)
    end = max(end1, end2)
    return start, end


def bam_to_bed(bam_file, out_file):
    cdef:
        AlignmentFile in_bam
        AlignedSegment read_1, read_2
        str chrom, strand
        int fragment_size
        str line
        long start, end

    pair_count = 0
    with pysam.Samfile(bam_file) as in_bam:
        while True:
            try:
                read_1 = in_bam.next()
                read_2 = in_bam.next()
                assert read_1.query_name == read_2.query_name, 'Paired not stored together'
                if qualify_pairs(read_1, read_2):
                    chrom = read_1.reference_name
                    strand = '-' if read_1.is_reverse else '+'
                    start, end = fragment_ends(read_1, read_2)
                    fragment_size = end - start
                    if 10 < fragment_size < 10000:
                        line = '%s\t%i\t%i\t%s\t%i\t%s' %(chrom, start, end,
                                                        read_1.query_name,
                                                        fragment_size,strand)
                        pair_count += 1
                        print(line, file=out_file)
            except StopIteration:
                break
    sys.stderr.write('Witten %i fragments\n' %(pair_count))
    return 0
