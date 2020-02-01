from __future__ import print_function, division
from libc.math cimport log, ceil, round
from libc.stdint cimport uint32_t
from builtins import range
import fileinput
from operator import itemgetter
import sys
import cython
from itertools import groupby
from .bed_dedup import fragment_coordinates


'''
implemented the function from:
Counting individual DNA molecules by the stochastic attachment of diverse labels.
https://www.ncbi.nlm.nih.gov/pubmed/21562209
'''

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef long correct_umi_count(int number_of_unique_umi, int umi_nt = 6):
    ## applying equation
    cdef:
        double m, modelled
        long diversity
        double result
        
    diversity = 4 ** umi_nt 
    m = number_of_unique_umi/diversity if number_of_unique_umi < diversity else (diversity-1)/diversity
    modelled = - diversity * log(1 - m)
    result = round(modelled)

    return max(number_of_unique_umi, long(result))




@cython.boundscheck(False)
@cython.wraparound(False)
cpdef parse_dedup(infile, outfile, int umi_nt, str read_prefix):
    cdef:
        str line, chrom, start, end, umi, strand
        uint32_t fragment_count
        int in_count = 0
        long out_count = 0
        set umis
        long true_fragments
        str outline_template
        str outline
        int distinct_umi 

    outline_template = '{chrom}\t{start}\t{end}\t{read_prefix}:UMI_{distinct_umi}_{frag_num}\t0\t{strand}' 

    for coordinates, lines in groupby(infile, fragment_coordinates):
        chrom, start, end, strand = coordinates
        umis = set()

        # count umis
        for line in lines:
            umi = line.split('\t')[3].split('_')[0]
            umis.add(umi)
            in_count += 1

        # compute true fragment counts
        distinct_umi = len(umis)
        true_fragments = correct_umi_count(distinct_umi, umi_nt = umi_nt)

        # output fragments
        for fragment_count in range(true_fragments):
            outline = outline_template\
                .format(chrom = chrom,
                        start = start,
                        end = end, 
                        read_prefix = read_prefix, 
                        distinct_umi = distinct_umi,
                        frag_num = fragment_count + 1,
                        strand = strand)
            print(outline, file = outfile)
            out_count += 1

    print('[Poisson UMI adjustment] Written %i from %i fragments' %(out_count, in_count), file = sys.stderr)
    return out_count, in_count
