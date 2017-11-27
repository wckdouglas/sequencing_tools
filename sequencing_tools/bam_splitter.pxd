from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport  bool

cdef class fragment_pairs:
    cdef readonly:
        AlignedSegment read1
        AlignedSegment read2
        bool pass_clip_check
        bool flag_qualify_ok
        bool has_soft_clip
        bool read1_has_S
        bool read2_has_S


