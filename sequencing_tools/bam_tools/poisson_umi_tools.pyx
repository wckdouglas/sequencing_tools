from __future__ import print_function, division
from libc.math cimport log, ceil, round
from libc.stdint cimport uint32_t
from builtins import range
import fileinput
from operator import itemgetter
import sys
import cython


'''
implemented the function from:
Counting individual DNA molecules by the stochastic attachment of diverse labels.
https://www.ncbi.nlm.nih.gov/pubmed/21562209
'''

cdef class fragment:
    cdef: 
        str chrom, start, end, strand
        set umi

    def __init__(self, str chrom, str start, str end, str strand, str umi):
        '''
        Initialize fragment group
        '''
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.umi = set()

        self.umi.add(umi)
    
    def check_fragment(self, str chrom, str start, str end, str strand):
        '''
        check if fragment has same genomic coordination as the fragment group
        '''
        return self.chrom == chrom and \
                self.start == start and \
                self.end == end and \
                self.strand == strand


    def add_fragment(self, str umi):
        '''
        adding unique UMI
        '''
        self.umi.add(umi)

    
    def output_fragments(self, int umi_nt = 6, out_file = sys.stdout):
        '''
        print X fragments according to the predicted number of template
        '''
        cdef:
            str line
            long theoretical 
            uint32_t frag_count
            uint32_t distinct_umi

        distinct_umi = len(self.umi)
        theoretical = correct_umi_count(distinct_umi, umi_nt)
        for frag_count in range(theoretical):
            line = '{chrom}\t{start}\t{end}\tUMI_{distinct_umi}_{frag_num}\t0\t{strand}' \
                    .format(chrom = self.chrom,
                            start = self.start,
                            end = self.end, 
                            distinct_umi = distinct_umi,
                            frag_num = frag_count + 1,
                            strand = self.strand)
            print(line, file = out_file)
        return theoretical 


@cython.boundscheck(False)
@cython.wraparound(False)
cdef long correct_umi_count(int number_of_unique_umi, int umi_nt = 6):
    ## applying equation
    cdef:
        double m, modelled
        long diversity
        double result
        
    diversity = 4 ** umi_nt 
    m = number_of_unique_umi/diversity
    modelled = - diversity * log(1 - m)
    result = round(modelled)

    return max(number_of_unique_umi, long(result))



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int parse_dedup(infile, outfile, int umi_nt):
    cdef:
        str line, chrom, start, end, umi, strand
        fragment _fragment
        uint32_t in_count
        long out_count = 0

    _fragment = fragment('chrom','start','end','strand','umi')
    for in_count, line in enumerate(infile):
        fields = line.strip().split('\t') 
        chrom, start, end, umi, strand = itemgetter(0,1,2,3,5)(fields)
        umi = umi.split('_')[0]
        if _fragment.check_fragment(chrom, start, end, strand):
            _fragment.add_fragment(umi)

        else:
            if _fragment.chrom != "chrom":
                out_count += _fragment.output_fragments(umi_nt = umi_nt, out_file = outfile)

            _fragment = fragment(chrom, start, end, strand, umi)

    out_count += _fragment.output_fragments(umi_nt = umi_nt, out_file = outfile)
    print('[Poisson UMI adjustment] Written %i from %i fragments' %(out_count, in_count), file = sys.stderr)
    return 0
