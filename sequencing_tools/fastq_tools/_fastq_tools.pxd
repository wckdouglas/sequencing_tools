#define fastq record type
cdef class fastqRecord:
    cdef readonly:
        str id
        str seq
        str qual
        str description
