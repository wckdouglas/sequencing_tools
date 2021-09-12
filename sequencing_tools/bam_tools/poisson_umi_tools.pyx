from __future__ import division, print_function

import fileinput
import logging
import os
import sys
from builtins import range
from itertools import groupby
from operator import itemgetter

import cython

from libc.math cimport log, ceil, round
from libc.stdint cimport uint32_t

from sequencing_tools.bam_tools.bed_dedup import fragment_coordinates

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


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
    """
    parsing dedup-ed file and write a poisson-adjusted bed file output, 
    implemented the function from: Counting individual DNA molecules by the stochastic attachment of diverse labels.
    https://www.ncbi.nlm.nih.gov/pubmed/21562209

    Basically taking in how many unique UMI being seen at one position, and using the complexity (number of UMI bases) to check for
    if UMI complexity is lower than what we actually need. Because UMI count would be inflated.

    An example input file looks like:

        chr1	10000	20000	ACG_1_members	10000	+
        chr1	10000	20000	GGG_2_members	10000	+
        chr1	10000	20000	ACT_3_members	10000	+
    
    The output file looks like:

        '{chrom}\t{start}\t{end}\t{read_prefix}:UMI_{distinct_umi}_{frag_num}\t0\t{strand}'

    :param IO[AnyStr] infile: input file handle 
    :param IO[AnyStr] outfile: output file handle 
    :param int umi_nt: number of nucleotides in the UMI
    :param str read_prefix: prefix of the read name

    """
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

    logger.info('[Poisson UMI adjustment] Written %i from %i fragments' %(out_count, in_count))
    return out_count, in_count
