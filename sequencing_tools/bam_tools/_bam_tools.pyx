import re
from builtins import map, range, zip
from operator import itemgetter
from libcpp cimport bool

import pysam
import six
from pysam.libcalignmentfilecimportAlignmentFile cimport AlignedSegment

from ..utils import SeqUtilsError

numbers = re.compile(r'[0-9]+')
strings = re.compile(r'[MIS]')
all_strings = re.compile(r'[A-Z]')
def split_cigar(cigar_string):
    '''
    Splitting cigar string to numpy array

    Args:
        cigar_string: a cigar string e.g. '63M1S4D10M'

    Returns:
        (list,list): ([list of cigar numbers], [list of cigar operators correspongs to the numbers])

    Example::

        $ split_cigar('63M')
        [63], [M]
    '''

    cigar_numbers = list(map(int, numbers.findall(cigar_string)))
    cigar_operator = all_strings.findall(cigar_string)
    cigar_array = [cigar_numbers, cigar_operator]
    return cigar_array


def read_ends(AlignedSegment read):
    '''
    Get read end positions, output start and end position of a read

    Args:
        AlignedSegment : a pysam alignment

    Returns:
        tuple: (start, end), leftmost and rightmost position of the read
    '''
    positions = read.get_reference_positions()
    start, end = itemgetter(0,-1)(positions)
    return start, end


def fragment_ends(AlignedSegment read1, AlignedSegment read2):
    '''
    get start and end position of a pair of reads

    Example::

        fragment_ends(read1, read2)

    Args:
        read1: a pysam alignment
        read2: a pysam alignment

    Returns:
        tuple: (start, end), leftmost and rightmost position of the fragment
    '''
    cdef:
        long start1, end1, start2, end2
        long start, end

    start1, end1 = read_ends(read1)
    start2, end2 = read_ends(read2)
    start = min(start1, start2)
    end = max(end1, end2) + 1
    return start, end


cpdef bool concordant_pairs(AlignedSegment read1, AlignedSegment read2):
    '''
    Check if pair is properly mapped flag == 99, 147, 163, 83

    Example::

        usage: concordant_pairs(read1, read2)

    Args:
        read1: read1 of the pair
        read2: read2 of the pair

    Returns:
        boolean: True if 83 and 163 or 99 and 147 for their flags otherwise False
    '''
    cdef:
        bool reverse_fragment = read1.flag == 83 and read2.flag == 163
        bool forward_fragment = read1.flag == 99 and read2.flag == 147
    return reverse_fragment or forward_fragment

cpdef bool concordant_alignment(AlignedSegment aln):
    '''
    Check if alignment is properly mapped flag == 99, 147, 163, 83

    Args:
        alignment: pysam alignment segment

    Returns:
        boolean: True if 83 or 163 or 99 or 147 for their flags otherwise False
    '''
    cdef:
        bool qualify
        int flag

    flag = aln.flag
    qualify = (flag == 99) or (flag == 147) or (flag == 163) or (flag == 83)
    return qualify

def make_regions(chromosome_length, how_many_bases_to_look_at):
    '''
    generator: segment chromosome in to regions
    usage: make_regions(chromosome_length, how_many_bases_every_time)
    return cigar_base

    Args:
        chromosome_length: last base you want to look at
        how_many_bases_every_time: segment size

    Returns:
        tuple: (start, end): start, end of segment

    Example::

        for start, end in make_regions(100,10):
            print start, end

        0 10
        10 20
        20 30
        30 40
        40 50
        50 60
        60 70
        70 80
        80 90
        90 100
    '''
    cdef:
        int start = 0
        int end = start + how_many_bases_to_look_at

    while end < chromosome_length:
        yield (start, end)
        start = end
        end = end + how_many_bases_to_look_at
    yield (start, chromosome_length)

def make_cigar_seq(cigar_numbers, cigar_operator):
    '''
    Generator convert number and operator into a sequence of aligned status of bases, see split_cigar

    Args:
        cigar_numbers: list of numbers in string format
        cigar_operator: list of single character

    Returns:
        Generator: cigar_base, sequence of cigar base

    Example::

        $ for c in make_cigar_seq(['3','5','1','3'],['S','M','I','M']):
        $     print c
        SSS
        MMMMM
        I
        MMM
    '''
    cdef:
        str num, op

    for num, op in zip(cigar_numbers, cigar_operator):
        if strings.search(op):
            yield int(num)*op

cpdef str cigar_to_str(str cigar_string):
    '''
    Cigarstring to string, only extract cigar op == M, I or S, Deletions (D) are skipepd to return same length cigar_seq as the query_sequence

    Args:
        cigar_string: standard cigar string from BAM file

    Returns:
        str: cigar_seq - same length string as sequence, with every character matched the aligned status of the base

    Example::

        $ cigar_seq = cigar_to_str('3S5M1I3M')
        $ print cigar_seq
        'SSSMMMMMIMMM'

    '''
    cigar_numbers = numbers.findall(cigar_string)
    cigar_operator = all_strings.findall(cigar_string)
    cigar_str = ''.join(make_cigar_seq(cigar_numbers, cigar_operator))
    return cigar_str

cpdef str get_strand(AlignedSegment aln, str direction = 'fr'):
    '''
    Get strand of the paired fragment

    Args:
        alignment: an pysam aligned segment

    Returns:
        str: Strand -  "+" if it is reverse read2 or forward read1 "-" if it is reverse read1 or forward read2
    '''
    cdef:
        str strand = ''
        bool read1_rvs
        bool read2_rvs

    if direction not in ['fr','rf', 'U']:
        raise SeqUtilsError('Wrong direction: only accept either fr or rf')
    read1_rvs = (aln.is_read1) and (aln.is_reverse)
    read2_rvs = (aln.is_read2) and (not aln.is_reverse)

    if direction == 'fr':
        if read1_rvs or read2_rvs:
            strand = '-'
        else:
            strand = '+'

    elif direction == 'rf':
        if read1_rvs or read2_rvs:
            strand = '+'
        else:
            strand = '-'
    elif direction == 'U': #single-end
        strand = '-' if aln.is_reverse else '+'

    return strand

def remove_insert(sequence, qual_seq, cigar):
    '''
    Iterator remove insertion base from aligned sequence

    usage: remove_insert(sequence, quality_string, cigar_seq)
    return: base, base_quality

    Args:
        sequence: DNA sequence from BAM
        quality_string: qual string from BAM
        cigar_seq: cigar seq from cigar_to_str

    Returns:
        Generator(tuple): (base, base_qual) base that are not annotated as insertion, BAQ associated with the base
    '''
    cdef:
        str base, op
        int qual

    for base, qual, op in zip(sequence, qual_seq, cigar):
        if op != 'I':
            yield base, qual

cpdef check_concordant(AlignedSegment read_1, AlignedSegment read_2):
    '''
    Check if read pairs are concordant

    Args:
        read1: first pysam alignment 
        read2: second pysam alignment 

    Returns:
        boolean: True if they have same read ID, ref ID, are read 1 and read2, and is opposite strand
    '''
    cdef:
        bool same_name = read_1.query_name == read_2.query_name
        bool same_ref = read_1.reference_id == read_2.reference_id 
        bool opposite_read = read_1.is_read1 != read_2.is_read1
        bool directional_pair = read_1.is_reverse != read_2.is_reverse
    return same_name and same_ref and opposite_read and check_concordant


cpdef check_primary(AlignedSegment read_1, AlignedSegment read_2):
    '''
    Check if read pairs are both primary

    Args:
        read1: first pysam alignment 
        read2: second pysam alignment 

    Returns:
        boolean: True if both of them are primary alignments
    '''
    return not read_1.is_secondary and not read_2.is_secondary
 
def paired_bam(AlignmentFile bam_handle):
    '''
    iterator for paired end bam file

    Args:
        bam_handle: bam file handle opened by pysam.Samfile

    Returns:
        Generator(tuple): (read1, read2)

    Example::

        with pysam.Samfile('path/to/bam/file') as bam:
            for read1, read2 in paired_bam(bam):
                print(read1.query_name)
                print(read2.query_name)
    '''
    cdef:
        AlignedSegment read1
        AlignedSegment read2

    try:
        while True:
            read1 = six.next(bam_handle)
            read2 = six.next(bam_handle)

            if not (read1.is_read1 and read2.is_read2):
                raise SeqUtilsError('First read is not read 1 or Second read is not read 2')
            if read1.query_name.split('/')[0] != read2.query_name.split('/')[0]:
                raise SeqUtilsError("Read names doesn't match! %s and %s" %(read1.query_name, read2.query_name))
                
            yield read1, read2

    except StopIteration:
        pass
