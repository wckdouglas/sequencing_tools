import re
import pysam
from itertools import izip
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool

numbers = re.compile(r'[0-9]+')
strings = re.compile(r'[MIS]')


def make_cigar_seq(cigar_numbers, cigar_operator):
    '''
    generator: convert number and operator into a sequence of aligned status of bases

    usage: make_cigar_seq(cigar_numbers, cigar_operator)
    return cigar_base

    ==================================
    parameter:

    cigar_numbers: list of numbers
    cigar_operator: list of single character

    return:
    cigar_base: sequence of cigar base

    example:
    for c in make_cigar_seq('3S5M1I3M'):
        print c

    SSS
    MMMMM
    I
    MMM
    ==================================
    '''
    cdef:
        str num, op

    for num, op in zip(cigar_numbers, cigar_operator):
        if op != 'S':
            yield int(num)*op

cpdef str cigar_to_str(str cigar_string):
    '''
    cigar string to string, only extract cigar op == M, I or S

    usage: cigar_to_str(cigar_string)
    return: cigar_seq

    ==================================
    parameter:

    cigar_string: standard cigar string from BAM file

    return:

    cigar_seq: same length string as sequence, with every character matched the aligned status of the base

    example:
    cigar_seq = cigar_to_str('3S5M1I3M')
    print cigar_seq
    'SSSMMMMMIMMM'

    ==================================
    '''
    cigar_numbers = numbers.findall(cigar_string)
    cigar_operator = strings.findall(cigar_string)
    cigar_str = ''.join(make_cigar_seq(cigar_numbers, cigar_operator))
    return cigar_str

cpdef str get_strand(AlignedSegment aln):
    '''
    get strand of the paired fragment

    usage: get_strand(alignment)
    return: strand

    ==================================
    parameter:

    alignment: an aligned segment in bam (see: pysam)

    return:

    strand: "+" if it is reverse read2 or normal read1
            "-" if it is reverse read1 or normal read2
    ==================================
    '''
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
    '''
    iterator remove insertion base from aligned sequence

    usage: remove_insert(sequence, quality_string, cigar_seq)
    return: base, base_quality

    ==================================
    parameter:

    sequence: DNA sequence from BAM
    quality_string: qual string from BAM
    cigar_seq: cigar seq from cigar_to_str

    yield:

    base:        base that are not annotated as insertion
    base_qual:   BAQ associated with the base
    ==================================
    '''
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
