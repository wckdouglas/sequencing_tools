from __future__ import print_function

import logging
import os
import sys
from operator import itemgetter

import pysam
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool

from sequencing_tools.utils import SeqUtilsError
from sequencing_tools.bam_tools._bam_tools import (check_concordant, check_primary, concordant_pairs,
                         fragment_ends, read_ends)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))



cdef class read_fragment:
    cdef:
        str tag
        int min_size
        long max_size
        str chrom, strand
        long start, end
        int fragment_size
        AlignedSegment read_1
        str bed_line
        str rt1
        bool pass_filter
        str cigar_field
        str prefix

    def __init__(self, read, tag, max_size, min_size, prefix=None):
        self.tag = tag
        self.read_1 = read
        self.max_size = max_size
        self.min_size = min_size
        self.cigar_field = self.read_1.cigarstring
        self.prefix = prefix

        if self.tag:
            self.rt1 = str(self.read_1.get_tag(self.tag))
        
        self.chrom = self.read_1.reference_name
        self.strand = '-' if self.read_1.is_reverse else '+'
        self.start, self.end = self.read_1.reference_start, self.read_1.reference_end
        self.fragment_size = self.end - self.start
        self.pass_filter =  self.min_size <= self.fragment_size <= self.max_size

    def generate_fragment(self, cigar=False):
        read_name = self.read_1.query_name 
        read_name = read_name if not self.prefix else self.prefix + ':' + read_name
        if self.pass_filter:
            self.bed_line = '{chrom}\t{start}\t{end}\t{read_name}\t{fragment_size}\t{strand}'\
                .format(chrom = self.chrom,
                        start = self.start, 
                        end = self.end, 
                        read_name = read_name,
                        fragment_size = self.fragment_size,
                        strand = self.strand)
            if self.tag: 
                self.bed_line = self.bed_line + '\t' + self.rt1
            if cigar:
                self.bed_line = self.bed_line + '\t' + self.cigar_field
        return self.bed_line



cdef class read_paired_fragment(read_fragment):
    cdef: 
        AlignedSegment read_2
        str rt2

    def __init__(self, read1, read2, tag=None, min_size=0, max_size=1000000, prefix=None):
        self.tag = tag
        self.min_size = min_size
        self.max_size = max_size

        # output elements
        self.bed_line = None
        self.pass_filter = False


        if read1.is_read1:
            self.read_1 = read1
            self.read_2 = read2
        elif read2.is_read1:
            self.read_1 = read2
            self.read_2 = read1

        if self.tag:
            self.rt1 = str(self.read_1.get_tag(self.tag))
            self.rt2 = str(self.read_2.get_tag(self.tag))
            try:
                self.rt1 = str(max(float(self.rt1), float(self.rt2)))
            except ValueError:
                if self.rt1 != self.rt2:
                    raise SeqUtilsError('Wrong tag %s and %s' %(self.rt1, self.rt2))

        #if concordant_pairs(self.read_1, self.read_2) and not (self.read_1.is_duplicate or self.read_2.is_duplicate):
        self.chrom = self.read_1.reference_name
        self.strand = '-' if self.read_1.is_reverse else '+'
        self.start, self.end = fragment_ends(self.read_1, self.read_2)
        self.fragment_size = self.end - self.start
        self.pass_filter =  self.min_size <= self.fragment_size <= self.max_size
        self.cigar_field = self.read_1.cigarstring + ':' + self.read_2.cigarstring


def pair_end_iterator(in_bam):
    '''
    read bam file two alignments by two alignments
    '''
    cdef:
        AlignedSegment read_1
        AlignedSegment read_2

    while True:
        try:
            read_1 = next(in_bam)
            read_2 = next(in_bam)
            if check_concordant(read_1, read_2):
                yield read_1, read_2
            else:
                raise SeqUtilsError('Pairs not sorted together, try samtools sort -n and/or samtools view -F2048 -F256')

        except StopIteration:
            break

def fragment_iterator(in_bam):
    '''
    read bam file two alignments by two alignments and return fragments
    '''
    cdef:
        AlignedSegment read_1
        AlignedSegment read_2

    for read_1, read_2 in pair_end_iterator(in_bam):
            yield read_paired_fragment(read_1, read_2)


def bam_to_bed(bam_file, out_file, int min_size, int max_size, 
                str tag, bool output_all, bool only_primary, bool add_cigar,
               prefix = None):
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
                read_1 = next(in_bam)
                read_2 = next(in_bam)
                concordant = check_concordant(read_1, read_2)
                primary_pair = check_primary(read_1, read_2)
                if concordant and ((primary_pair and only_primary) or not only_primary):#, 'Paired not stored together: %s, %s'  %(read_1.query_name , read_2.query_name)
                    pair_fragment = read_paired_fragment(read_1, read_2, tag = tag, 
                                                        max_size = max_size, min_size = min_size,
                                                        prefix = prefix)
                    line = pair_fragment.generate_fragment(cigar = add_cigar)
                    if line:
                        pair_count += 1
                        print(line, file=out_file)
#                else:
#                    raise SeqUtilsError("There are disconcordant pairs: %s and %s" %(read_1.query_name, read_2.query_name))
#                else:
#                    while not check_concordant(read_1, read_2):
#                        if not read_1.is_secondary and not read_2.is_secondary:
#                            if output_all:
#                                if read_2.is_supplementary:
#                                    fragment = read_fragment(read_2, tag, max_size, min_size)
#                                    read_2 = next(in_bam)
#
#                                else:
#                                    fragment = read_fragment(read_1, tag, max_size, min_size)
#                                    read_1 = read_2
#                                    read_2 = next(in_bam)
#
#                                line = fragment.generate_fragment(cigar = add_cigar)
#                                single_count += 1
#                                print(line, file=out_file)
#                            else:
#                                read_1 = read_2
#                                read_2 = next(in_bam)
#                        elif read_1.is_secondary:
#                            read_1 = read_2
#                            read_2 = next(in_bam)
#
#                        elif read_2.is_secondary:
#                            read_2 = next(in_bam)
#
#                    pair_fragment = read_paired_fragment(read_1, read_2, tag = tag, 
#                                                        max_size = max_size, min_size = min_size)
#                    line = pair_fragment.generate_fragment(cigar = add_cigar)
#                    if line:
#                        pair_count += 1
#                        print(line, file=out_file)
            except StopIteration:
                break
    logger.info('Witten %i pair fragments and %i multiple jump fragments' %(pair_count, single_count))
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


def mate_info(in_bam, out_bam):
    '''
    Read two alignments at a time,
    assume they are pairs,
    add mate cigar and mate mapq
    '''
    cdef:
        AlignmentFile in_bam_handle, out_bam_handle
        AlignedSegment read_1, read_2
        long pair_count

    pair_count = 0
    with pysam.Samfile(in_bam,'rb') as in_bam_handle:
        with pysam.Samfile(out_bam,'w', template = in_bam_handle) as out_bam_handle:
            while True:
                try:
                    read_1 = next(in_bam_handle)
                    read_2 = next(in_bam_handle)
                    concordant = check_concordant(read_1, read_2)
                    primary_pair = check_primary(read_1, read_2)
                    if concordant:#, 'Paired not stored together: %s, %s'  %(read_1.query_name , read_2.query_name)

                        # set mate cigar
                        read_1.set_tag('MC', read_2.cigarstring)
                        read_2.set_tag('MC', read_1.cigarstring)

                        # set mate mapq
                        read_1.set_tag('MQ', read_2.mapping_quality)
                        read_2.set_tag('MQ', read_1.mapping_quality)

                        out_bam_handle.write(read_1)
                        out_bam_handle.write(read_2)

                        pair_count += 1
                except StopIteration:
                    break
    logger.info('Witten %i pair fragments' %(pair_count))
    return 0 
