#!/usr/bin/env python

from __future__ import print_function
import pysam
from operator import itemgetter
from cpython cimport bool
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import sys
from sequencing_tools.bam_tools import concordant_pairs, read_ends, fragment_ends, check_concordant

class read_paired_fragment:
    def __init__(self, read1, read2, tag, max_size, min_size):
        self.tag = tag
        self.max_size = max_size
        self.min_size = min_size
        self.read_1 = None
        self.read_2 = None

        # output elements
        self.bed_line = None
        self.chrom = None
        self.start = None
        self.end = None
        self.fragment_size = None
        self.strand = None
        self.pass_filter = False


        if read1.is_read1:
            self.read_1 = read1
            self.read_2 = read2
        elif read2.is_read1:
            self.read_1 = read2
            self.read_2 = read1

        if self.tag:
            self.rt1 = self.read_1.get_tag(self.tag)
            self.rt2 = self.read_2.get_tag(self.tag)

        if concordant_pairs(self.read_1, self.read_2) and not (self.read_1.is_duplicate or self.read_2.is_duplicate):
            self.chrom = self.read_1.reference_name
            self.strand = '-' if self.read_1.is_reverse else '+'
            self.start, self.end = fragment_ends(self.read_1, self.read_2)
            self.fragment_size = self.end - self.start
            self.pass_filter =  self.min_size <= self.fragment_size <= self.max_size

    def generate_fragment(self):
        if self.pass_filter:
            self.bed_line = '{chrom}\t{start}\t{end}\t{read_name}\t{fragment_size}\t{strand}'\
                .format(chrom = self.chrom,
                        start = self.start, 
                        end = self.end, 
                        read_name = self.read_1.query_name,
                        fragment_size = self.fragment_size,
                        strand = self.strand)
            if self.tag: 
                self.bed_line = self.bed_line + '\t' + self.rt1
        return self.bed_line


class read_fragment:

    def __init__(self, read, tag, max_size, min_size):
        self.tag = tag
        self.read = read
        self.max_size = max_size
        self.min_size = min_size

        if self.tag:
            self.rt = self.read.get_tag(self.tag)
        
        self.chrom = self.read.reference_name
        self.strand = '-' if self.read.is_reverse else '+'
        self.start, self.end = self.read.reference_start, self.read.reference_end
        self.fragment_size = self.end - self.start

    def generate_fragment(self):
        if self.min_size < self.fragment_size < self.max_size:
            self.bed_line = '{chrom}\t{start}\t{end}\t{read_name}\t{fragment_size}\t{strand}'\
                .format(chrom = self.chrom,
                        start = self.start, 
                        end = self.end, 
                        read_name = self.read_1.query_name,
                        fragment_size = self.fragment_size,
                        strand = self.strand)
            if self.tag: 
                self.bed_line = self.bed_line + '\t' + self.rt1
        return self.bed_line


def bam_to_bed(bam_file, out_file, int min_size, int max_size, tag, output_all):
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
        int pair_count, single_count

    pair_count = 0
    single_count = 0
    with pysam.Samfile(bam_file,'rb') as in_bam:
        while True:
            try:
                read_1 = in_bam.next()
                read_2 = in_bam.next()
                if check_concordant(read_1, read_2):#, 'Paired not stored together: %s, %s'  %(read_1.query_name , read_2.query_name)
                    pair_fragment = read_paired_fragment(read_1, read_2, tag, max_size, min_size)
                    line = pair_fragment.generate_fragment()
                    if line:
                        pair_count += 1
                        print(line, file=out_file)
                else:
                    while not check_concordant(read_1, read_2):
                        if not read_1.is_secondary and not read_2.is_secondary:
                            if output_all:
                                if read_2.is_supplementary:
                                    fragment = read_fragment(read_2, tag, max_size, min_size)
                                    read_2 = in_bam.next()

                                else:
                                    fragment = read_fragment(read_1, tag, max_size, min_size)
                                    read_1 = read_2
                                    read_2 = in_bam.next()

                                line = fragment.generate_fragment()
                                single_count += 1
                                print(line, file=out_file)
                            else:
                                read_1 = read_2
                                read_2 = in_bam.next()
                        elif read_1.is_secondary:
                            read_1 = read_2
                            read_2 = in_bam.next()

                        elif read_2.is_secondary:
                            read_2 = in_bam.next()

                    pair_fragment = read_paired_fragment(read_1, read_2, tag, max_size, min_size)
                    line = pair_fragment.generate_fragment()
                    if line:
                        pair_count += 1
                        print(line, file=out_file)
            except StopIteration:
                break
    print('Witten %i pair fragments and %i multiple jump fragments' %(pair_count, single_count), file = sys.stderr)
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
