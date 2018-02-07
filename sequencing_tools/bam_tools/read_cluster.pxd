cdef class readGroup:
    cdef public:
        str barcode
        list R1
        list R2
        list R1_flag
        list R2_flag
        list R1_position
        list R2_position
        list R1_chrom
        list R2_chrom
        list concensus_read1
        list concensus_read2
        list member_count_list
        list concensus_flag1
        list concensus_flag2
        str fastq_record
#
