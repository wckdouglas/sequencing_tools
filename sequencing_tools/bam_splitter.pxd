from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport  bool

cdef class fragment_pairs:
    cdef readonly:
        AlignedSegment read1
        AlignedSegment read2
        bool pass_clip_check, flag_qualify_ok, has_soft_clip, read1_has_S, read2_has_S


