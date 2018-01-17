from __future__ import print_function
from builtins import zip, map
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool
from sequencing_tools.bam_tools import *


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



INDEL = re.compile('I\D')
def analyze_region(bam, chromosome, qual_threshold, crop, no_indel, base_dict, start, end):
    cdef:
        int aln_count = 0
        AlignedSegment aln
        int pos, qual, seq_len
        str base
        int i, crop_end
        int aln_pos
        long ref_pos
        bool no_indel_condition, with_indel_condition

    for aln_count, aln in enumerate(bam.fetch(chromosome, start, end)):
        strand = get_strand(aln)
        if not aln.is_unmapped and strand:
            with_indel = INDEL.search(str(aln.cigarstring))
            no_indel_condition = (not with_indel and no_indel)
            with_indel_condition = (not no_indel)
            if (with_indel_condition or no_indel_condition):

                ### using pysam count_coverage approach ###
                crop_end = aln.query_length - crop
                for aln_pos, ref_pos in aln.get_aligned_pairs(True):
                    qual = aln.query_qualities[aln_pos]
                    base = aln.query_sequence[aln_pos]
                    if crop_end >= aln_pos >= crop and qual >= qual_threshold:
                        base_dict[ref_pos][strand][base] += 1

                '''
                # old algorithm #
                positions = aln.get_reference_positions()
                sequence = aln.query_alignment_sequence
                cigar_str = cigar_to_str(aln.cigarstring).replace('S','')
                qual_seq = aln.query_alignment_qualities
                adjusted_sequence = remove_insert(sequence, qual_seq, cigar_str)
                crop_end = len(positions) - crop
                for i, (pos, (base, qual)) in enumerate(zip(positions, adjusted_sequence)):
                    if crop_end >= i >= crop and qual >= qual_threshold:
                        base_dict[pos][strand][base] += 1
                ''''
    return aln_count, base_dict
