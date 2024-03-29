from __future__ import print_function

import logging
import os
import sys

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

from sequencing_tools.bam_tools.fragment_pairs import concordant_pairs, is_split_pair

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))

cpdef int parse_bam(AlignmentFile inbam,
                    AlignmentFile out_split_bam,
                    AlignmentFile out_fragment_bam):
    cdef:
        int split_out = 0
        int out = 0
        AlignedSegment read1, read2

    while True:
        try:
            read1 = next(inbam)
            read2 = next(inbam)
            if concordant_pairs(read1, read2):
                if is_split_pair(read1, read2):
                    out_split_bam.write(read1)
                    out_split_bam.write(read2)
                    split_out += 2
                else:
                    out_fragment_bam.write(read1)
                    out_fragment_bam.write(read2)
                    out += 2
        except StopIteration:
            break
    logger.info('Written %i unsplit and %i split alignments' %(out, split_out))
    return 0
