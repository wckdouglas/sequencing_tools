from itertools import izip
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool
from tgirt_seq_tools.bam_tools import *


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
            cigar_str = cigar_to_str(aln.cigarstring).replace('S','')
            qual_seq = aln.query_alignment_qualities
            adjusted_sequence = remove_insert(sequence, qual_seq, cigar_str)
            crop_end = len(positions) - crop
            for i, (pos, (base, qual)) in enumerate(izip(positions, adjusted_sequence)):
                if crop_end >= i >= crop and qual >= qual_threshold:
                    base_dict[pos][strand][base] += 1
    return aln_count, base_dict
