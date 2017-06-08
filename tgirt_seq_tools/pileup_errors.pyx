import re
import pysam

numbers = re.compile(r'[0-9]+')
strings = re.compile(r'[MIS]')
def cigar_to_str(cigar_string):
    '''
    cigar string to string
    '''
    cigar_numbers = numbers.findall(cigar_string)
    cigar_operator = strings.findall(cigar_string)
    cigar_str = ''.join(make_cigar_seq(cigar_numbers, cigar_operator))
    return cigar_str

def get_strand(aln):
    if aln.is_reverse:
        return '-'
    else:
        return '+'

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
    for strand in ['+','-']:
        for base in 'ACTG':
            base_counts.append(str(base_dict[pos][strand][base]))
    return '\t'.join(base_counts)

def make_cigar_seq(cigar_numbers, cigar_operator):
    cdef:
        str num, op
        
    for num, op in zip(cigar_numbers, cigar_operator):
        if op != 'S':
            yield int(num)*op
