from __future__ import print_function
import sys
from cpython cimport bool
from operator import itemgetter

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
    chrom2, start2, end2, strand1, strand2 = itemgetter(3,4,5,8,9,10)(fields)
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


cpdef int processFile(bed_iterator, int min_length, int max_length, out_handle):
    cdef:
        str frag
        str aln

    for frag in bed_iterator:
        aln = filterBed(frag, min_length, max_length)
        if aln:
            print(aln, file = out_handle)
    return 0
