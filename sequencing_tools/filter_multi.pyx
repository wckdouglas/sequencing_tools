import pysam
from cpython cimport bool
import sys
from pysam.libcalignmentfile cimport AlignedSegment
import numpy as np
from itertools import izip
import re


class read_pairs:
    def __init__(self):
        self.read1 = []
        self.read2 = []
        self.out_read1 = ''
        self.out_read2 = ''
        self.read_id = ''

    def initiate_group(self, alignment):
        self.read_id = alignment.query_name
        self.put_in_group(alignment)

    def put_in_group(self, alignment):
        if alignment.is_read2:
            self.read2.append(alignment)
        else:
            self.read1.append(alignment)


    def generate_filter_alingments(self):
        cdef:
            AlignedSegment r1, r2

        chroms = []
        isizes = []
        read1_group = []
        read2_group = []
        for r1, r2 in izip(self.read1, self.read2):
            if r1.reference_name == r2.reference_name and abs(r1.isize) == abs(r2.isize):
                read1_group.append(r1)
                read2_group.append(r2)
                isizes.append(abs(r1.isize))
                chroms.append(r1.reference_name)


        chroms = np.array(chroms)
        isizes = np.array(isizes)

        #start filtering
        size_bool = (isizes == np.min(isizes))
        ribo_bool = is_ribo_chrom(chroms)
        regular_chrom_bool = is_regular_chrom(chroms)

        read1 , read2 = map(np.array, [read1_group, read2_group])

        if len(isizes[size_bool]) == 1:
            self.out_read1 =  read1[size_bool][0]
            self.out_read2 = read2[size_bool][0]
        elif len(chroms[ribo_bool]) == 1:
            self.out_read1 =  read1[ribo_bool][0]
            self.out_read2 = read2[ribo_bool][0]
        elif len(chroms[regular_chrom_bool]) > 0:
            self.out_read1 =  read1[regular_chrom_bool][0]
            self.out_read2 = read2[regular_chrom_bool][0]
        else:
            self.out_read1 =  read1[0]
            self.out_read2 = read2[0]

    def output_read(self):
        read1_aln, read2_aln = fix_flag(self.out_read1, self.out_read2)
        return read1_aln, read2_aln

cigar_num = re.compile('\d+')
cigar_str = re.compile('[A-Z]')
cdef mapped_length(str cigarstring):
    '''
    for cigar string, find sum of number of mapped base
    '''
    cdef:
        str cnum
        int mapped = 0
        str cstr
    

    for cnum, cstr in izip(cigar_num.findall(cigarstring), cigar_str.findall(cigarstring)):
        if cstr == 'M':
            mapped += int(cnum)
    return mapped

class single_read:
    def __init__(self):
        self.reads = []
        self.out_read = None
        self.read_id = None

    def initiate_group(self, alignment):
        self.read_id = alignment.query_name
        self.put_in_group(alignment)

    def put_in_group(self, alignment):
        self.reads.append(alignment)

    def generate_filter_alingments(self):
        cdef:
            AlignedSegment read

        chroms = []
        isizes = []
        read_group = []
        for read in self.reads:
            isizes.append(mapped_length(read.cigarstring))
            chroms.append(read.reference_name)

        chroms = np.array(chroms)
        isizes = np.array(isizes)

        #start filtering
        size_bool = (isizes == np.min(isizes))
        ribo_bool = is_ribo_chrom(chroms)
        regular_chrom_bool = is_regular_chrom(chroms)

        reads = np.array(self.reads)

        if len(isizes[size_bool]) == 1:
            self.out_read =  reads[size_bool][0]
        elif len(chroms[ribo_bool]) == 1:
            self.out_read =  reads[ribo_bool][0]
        elif len(chroms[regular_chrom_bool]) > 0:
            self.out_read =  reads[regular_chrom_bool][0]
        else:
            self.out_read = reads[0]

    def output_read(self):
        read_aln = fix_single(self.out_read)
        return read_aln


regular_chroms = range(1,23)
regular_chroms.extend(list('XY'))
regular_chroms.append('MT')
regular_chroms = map(str, regular_chroms)

def is_regular_chrom(chroms):
    return np.in1d(chroms, regular_chroms)

def is_ribo_chrom(chroms):
    return np.array([True if 'gi' in chrom else False for chrom in chroms])

def fix_single(read):
    read.flag = 16 if read.is_reverse else 0
    read.is_supplementary = False
    read.is_secondary = False
    return read

def fix_flag(read1, read2):
    read1_flag, read2_flag = 99, 147
    if read1.is_reverse:
        read1_flag = 83
        read2_flag = 163

    read1.flag = read1_flag
    read2.flag = read2_flag
    read1.is_supplementary = False
    read1.is_secondary = False
    read2.is_supplementary = False
    read2.is_secondary = False
    return read1, read2

def process_pair_bam(in_bam, out_bam, bam_in_bool, bam_out_bool):

    cdef:
        AlignedSegment alignment, read1_aln, read2_aln
        int read_count
        int out_read_count = 0
        int read_group_count = 0

    read_flag = 'rb' if bam_in_bool else 'r'
    write_flag = 'wb' if bam_out_bool else 'w'
    sys.stderr.write('Start processing bam file: %s\n' %(in_bam))
    sys.stderr.write('Writing to: %s\n' %(out_bam))
    with pysam.Samfile(in_bam, read_flag) as in_sam:
        with pysam.Samfile(out_bam, write_flag, template = in_sam) as out_sam:
            for read_count, alignment in enumerate(in_sam):
                if read_group_count == 0:
                    # initial group for first alignment
                    read_group = read_pairs()
                    read_group.initiate_group(alignment)
                    read_group_count = 1
                else:
                    if alignment.query_name != read_group.read_id:
                        read_group.generate_filter_alingments()
                        read1_aln, read2_aln = read_group.output_read()
                        out_sam.write(read1_aln)
                        out_sam.write(read2_aln)
                        out_read_count += 1
                        read_group = read_pairs()
                        read_group.initiate_group(alignment)
                    else:
                        read_group.put_in_group(alignment)
            #After all reading whole file, clean out memory and output last alignment group
            read_group.generate_filter_alingments()
            read1_aln, read2_aln = read_group.output_read()

            out_sam.write(read1_aln)
            out_sam.write(read2_aln)
            out_read_count += 1
    sys.stderr.write('Writting %i read pairs\n' %(out_read_count))
    return 0

def process_single_bam(in_bam, out_bam, bam_in_bool, bam_out_bool):

    cdef:
        AlignedSegment alignment, read_aln
        int read_count
        int out_read_count = 0
        int read_group_count = 0

    read_flag = 'rb' if bam_in_bool else 'r'
    write_flag = 'wb' if bam_out_bool else 'w'
    sys.stderr.write('Start processing bam file: %s\n' %(in_bam))
    sys.stderr.write('Writing to: %s\n' %(out_bam))
    with pysam.Samfile(in_bam, read_flag) as in_sam:
        with pysam.Samfile(out_bam, write_flag, template = in_sam) as out_sam:
            for read_count, alignment in enumerate(in_sam):
                if read_group_count == 0:
                    # initial group for first alignment
                    read_group = single_read()
                    read_group.initiate_group(alignment)
                    read_group_count = 1
                else:
                    if alignment.query_name != read_group.read_id:
                        read_group.generate_filter_alingments()
                        read_aln = read_group.output_read()
                        out_sam.write(read_aln)
                        out_read_count += 1
                        read_group = single_read()
                        read_group.initiate_group(alignment)
                    else:
                        read_group.put_in_group(alignment)
            #After all reading whole file, clean out memory and output last alignment group
            read_group.generate_filter_alingments()
            read_aln = read_group.output_read()

            out_sam.write(read_aln)
            out_read_count += 1
    sys.stderr.write('Writting %i reads\n' %(out_read_count))
    return 0