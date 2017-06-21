import re
import pysam
from itertools import izip
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool

numbers = re.compile(r'[0-9]+')
strings = re.compile(r'[MIS]')
cpdef str cigar_to_str(str cigar_string):
    '''
    cigar string to string
    '''
    cigar_numbers = numbers.findall(cigar_string)
    cigar_operator = strings.findall(cigar_string)
    cigar_str = ''.join(make_cigar_seq(cigar_numbers, cigar_operator))
    return cigar_str

cpdef str get_strand(AlignedSegment aln):
    cdef:
        str strand = ''
        bool read1_rvs
        bool read2_rvs

    read1_rvs = (aln.is_read1) and (aln.is_reverse)
    read2_rvs = (aln.is_read2) and (not aln.is_reverse)
    if read1_rvs or read2_rvs:
        strand = '-'
    else:
        strand = '+'
    return strand

def remove_insert(sequence, qual_seq, cigar):
    cdef:
        str base, op
        int qual

    for base, qual, op in zip(sequence, qual_seq, cigar):
        if op != 'I':
            yield base, qual

def extract_bases(base_dict, pos):
    cdef:
        str strand, base

    base_counts = []
    coverage = 0
    for strand in ['+','-']:
        for base in 'ACGT':
            bcount = base_dict[pos][strand][base]
            coverage += bcount
            base_counts.append(str(bcount))
    return coverage,'\t'.join(base_counts)

def make_cigar_seq(cigar_numbers, cigar_operator):
    cdef:
        str num, op

    for num, op in zip(cigar_numbers, cigar_operator):
        if op != 'S':
            yield int(num)*op


def analyze_region(bam, chromosome, qual_threshold, crop, base_dict, start, end):
    cdef:
        int aln_count
        AlignedSegment aln
        int pos, qual, seq_len
        str base
        int i, crop_end

    for aln_count, aln in enumerate(bam.fetch(chromosome, start, end)):
        strand = get_strand(aln)
        if not aln.is_unmapped and strand:
            positions = aln.get_reference_positions()
            sequence = aln.query_alignment_sequence
            cigar_str = cigar_to_str(aln.cigarstring)
            qual_seq = aln.query_alignment_qualities
            adjusted_sequence = remove_insert(sequence, qual_seq, cigar_str)
            crop_end = len(positions) - crop
            for i, (pos, (base, qual)) in enumerate(izip(positions, adjusted_sequence)):
                if crop_end >= i >= crop and qual >= qual_threshold:
                    base_dict[pos][strand][base] += 1
    return aln_count, base_dict
