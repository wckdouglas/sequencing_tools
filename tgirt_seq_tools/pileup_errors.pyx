import re
import pysam
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment    

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
    cdef str strand = ''
    if aln.is_reverse:
        strand = '-'
    else:
        strand = '+'

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
        for base in 'ACTG':
            bcount = base_dict[pos][strand][base]
            coverage += bcount
            base_counts.append(str(bcount))
    return sum(base_counts),'\t'.join(base_counts)

def make_cigar_seq(cigar_numbers, cigar_operator):
    cdef:
        str num, op
        
    for num, op in zip(cigar_numbers, cigar_operator):
        if op != 'S':
            yield int(num)*op
