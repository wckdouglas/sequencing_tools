from __future__ import print_function

from builtins import map, zip
from collections import defaultdict

from sequencing_tools.bam_tools import *

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool


def extract_bases(base_dict, pos):
    cdef:
        str strand, base
        int coverage 
        dict coverage_dict = {}

    base_counts = []
    for strand in ['+','-']:
        coverage = 0
        for base in 'ACGT':
            bcount = base_dict[pos][strand][base]
            coverage += bcount
            base_counts.append(str(bcount))
        coverage_dict[strand] = coverage
    return coverage_dict,'\t'.join(base_counts)

def test_direction(bam):
    aln = next(bam)
    return 'U' if aln.flag in [0,16] else 'fr'

INDEL = re.compile('I\D')
def analyze_region(bam, chromosome, qual_threshold, crop, no_indel, start, end):
    '''
    input:
        bam: pysam SamFile object
        chromosome: str
        qual_threshold: int
        crop: boolean, whether sequence head or tail is accounted for
        no_indel: boolean, whether alignment iwth indels are accounted for
        start: int
        end: int


    return:
        aln_count: number of alignments
        base_dict: multilevel dictionary: refpos/strand/base -> base count

    '''


    cdef:
        int aln_count = 0
        AlignedSegment aln
        long pos
        int qual, seq_len
        str base
        int i, crop_end
        int aln_pos
        long ref_pos
        bool no_indel_condition, with_indel_condition

    direction = test_direction(bam)
    base_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int))) #ref pos//strand//base
    for aln_count, aln in enumerate(bam.fetch(chromosome, start, end)):
        strand = get_strand(aln, direction = direction)
        if not aln.is_unmapped and \
                strand and \
                not aln.is_duplicate and \
                not aln.is_supplementary:
            with_indel = INDEL.search(str(aln.cigarstring))
            no_indel_condition = (not with_indel and no_indel)
            with_indel_condition = (not no_indel)
            if (with_indel_condition or no_indel_condition):

                ### using pysam count_coverage approach ###
                crop_end = aln.query_length - crop
                for aln_pos, ref_pos in aln.get_aligned_pairs(matches_only = True):
                    qual = aln.query_qualities[aln_pos]
                    base = aln.query_sequence[aln_pos]
                    if ref_pos and aln_pos and crop_end >= aln_pos >= crop and qual >= qual_threshold:
                        base_dict[ref_pos][strand][base] += 1

    return aln_count + 1, base_dict
