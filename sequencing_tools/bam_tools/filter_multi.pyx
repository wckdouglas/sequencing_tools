from __future__ import print_function
import pysam
from cpython cimport bool
import sys
from pysam.libcalignmentfile cimport AlignedSegment
from pysam.libctabix cimport TabixFile
import numpy as np
from builtins import zip, range
import re
from libc.stdlib cimport rand
from itertools import groupby

regular_chroms = list(range(1,23))
regular_chroms.extend(list('XY'))
regular_chroms = list(map(str, regular_chroms))
regular_chroms.extend(list(map(lambda x: 'chr'+x, regular_chroms)))


cpdef int fast_random_number(int scale):
    return rand() % scale


class read_pairs:
#    cdef:
#        list read1, read2
#        AlignedSegment out_read1, out_read2
#        str read_id
#        bool is_group
#        object reads

    def __init__(self, read_id, reads, gene_tab):
        self.read1 = []
        self.read2 = []
        self.out_read1 = None
        self.out_read2 = None
        self.read_id = read_id
        self.is_group = None
        self.gene_tabix = gene_tab

        for read in reads:
            if read.is_read1:
                self.read1.append(read)
            elif read.is_read2:
                self.read2.append(read)


    def generate_filtered_alingments(self):
        cdef:
            AlignedSegment r1, r2, r
            int scale

        read1_group = []
        read2_group = []
        scores = []
        for r1, r2 in zip(self.read1, self.read2):
            if r1.reference_name == r2.reference_name and abs(r1.isize) == abs(r2.isize):
                read1_group.append(r1)
                read2_group.append(r2)
                score = r1.get_tag('AS') + r2.get_tag('AS')
                scores.append(score)

        read1, read2, scores = map(np.array, [read1_group, read2_group, scores])

        #start filtering
        '''
        1. best score?
        2. ribo or mt?
        3. shortest isize
        4. regular chrom
        5. random picked

        if it is a condition that will definitely yield an output 
        (smallese isize, best scores vs mt/ribo chrom), 
        use read1 as filtered, otherwise use new_read1
        '''
        best_score_bool = (scores == scores.max())
        read1, read2 = read1[best_score_bool], read2[best_score_bool]

        if len(read1) == 1:
            self.out_read1 = read1[0]
            self.out_read2 = read2[0]
        
        else:
            # ribo mt reads?
            ribo_bool = is_ribo_chrom([r.reference_name for r in read1])
            new_read1 = read1[ribo_bool]
            new_read2 = read2[ribo_bool]

            if len(new_read1) == 1:
                self.out_read1 = new_read1[0]
                self.out_read2 = new_read2[0]

            elif len(new_read1) > 1:
                scale = len(new_read1)
                selected = fast_random_number(scale)
                self.out_read1 =  new_read1[selected]
                self.out_read2 = new_read2[selected]


            else:
                #is on gene?
                if self.gene_tabix:
                    gene_bool = np.array(list(map(self.is_gene, read1, read2)))
                    if any(gene_bool):
                        read1 = read1[gene_bool]
                        read2 = read2[gene_bool]

                # isize ?
                isizes = np.array([abs(r1.isize) for r1 in read1])
                size_bool = (isizes == np.min(isizes))
                read1, read2 = read1[size_bool], read2[size_bool]

                if len(read1) == 1:
                    self.out_read1 = read1[0]
                    self.out_read2 = read2[0]

                else:
                    # if no rRNA, look at regular chromosome fragments
                    regular_chrom_bool = is_regular_chrom([r.reference_name for r in read1])
                    new_read1, new_read2 = read1[regular_chrom_bool], read2[regular_chrom_bool]

                    if len(new_read1) == 1:
                        self.out_read1 = new_read1[0]
                        self.out_read2 = new_read2[0]

                    elif len(new_read1) > 1:
                        # randomly select one if there are multiple regular chrom fragments
                        scale = len(new_read1)
                        selected = fast_random_number(scale)
                        self.out_read1 =  new_read1[selected]
                        self.out_read2 =  new_read2[selected]
                    
                    else:
                        # randomly pick one if all are same
                        scale = len(read1)
                        selected = fast_random_number(scale)
                        self.out_read1 = read1[selected]
                        self.out_read2 = read2[selected]


    def output_read(self):
        read1_aln, read2_aln = fix_flag(self.out_read1, self.out_read2)
        return read1_aln, read2_aln


    def is_gene(self, read1, read2):
        frag_start = min(read1.pos, read2.pos)
        frag_end = max(read1.pos + read1.pos + read1.alen,
                        read2.pos + read2.pos + read2.alen)
        genes = self.gene_tabix.fetch(read1.reference_name, frag_start, frag_end)
        return True if list(genes) else False
        



cdef int mapped_length(AlignedSegment read):
    '''
    for cigar string, find sum of number of mapped base
    '''
    cdef:
        str cnum
        int mapped = 0
        str cstr
        int nt, op
    
    for op, nt in read.cigartuples:
        if op == 0: #BAM CMATCH
            mapped += nt

    return mapped

class single_read:
    def __init__(self, read_id, reads):
        self.reads = list(reads)
        self.out_read = None
        self.read_id = read_id


    def generate_filtered_alingments(self):
        cdef:
            AlignedSegment read

        chroms = []
        isizes = []
        read_group = []
        for read in self.reads:
            isizes.append(mapped_length(read))
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

def is_regular_chrom(chroms):
    return np.in1d(chroms, regular_chroms)

ribo_mt = re.compile('^gi|^chrM|^MT')
def is_ribo_chrom(chroms):
    return np.array([True if ribo_mt.search(chrom) else False for chrom in chroms])

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

def process_pair_bam(in_bam, out_bam, bam_in_bool, bam_out_bool, gene_file):

    cdef:
        AlignedSegment read1_aln, read2_aln
        int in_read_count
        int out_read_count = 0

    read_flag = 'rb' if bam_in_bool else 'r'
    write_flag = 'wb' if bam_out_bool else 'w'
    print('Start processing bam file: %s' %(in_bam), file = sys.stderr)
    print('Writing to: %s' %(out_bam), file = sys.stderr)
    gene_tabix = pysam.Tabixfile(gene_file) if gene_file else None
    with pysam.Samfile(in_bam, read_flag) as in_sam:
        with pysam.Samfile(out_bam, write_flag, template = in_sam) as out_sam:
            for in_read_count, (read_id, alignments) in enumerate(groupby(in_sam, lambda aln: aln.query_name)):
                # initial group for first alignment
                read_group = read_pairs(read_id, alignments, gene_tabix)
                read_group.generate_filtered_alingments()
                read1_aln, read2_aln = read_group.output_read()
                out_sam.write(read1_aln)
                out_sam.write(read2_aln)
                out_read_count += 1
    print('Writting %i read pairs from %i groups' %(out_read_count, in_read_count), file = sys.stderr)
    return 0

def process_single_bam(in_bam, out_bam, bam_in_bool, bam_out_bool, gene_file):

    cdef:
        AlignedSegment alignment, read_aln
        int in_read_count
        int out_read_count = 0

    read_flag = 'rb' if bam_in_bool else 'r'
    write_flag = 'wb' if bam_out_bool else 'w'
    print('Start processing bam file: %s' %(in_bam), file = sys.stderr)
    print('Writing to: %s' %(out_bam), file = sys.stderr)
    read_group = single_read()

    with pysam.Samfile(in_bam, read_flag) as in_sam:
        with pysam.Samfile(out_bam, write_flag, template = in_sam) as out_sam:
            for in_read_count, (read_id, alignments) in enumerate(in_sam, lambda aln: aln.query_name):
                read_group = single_read(read_id, alignments)
                read_group.generate_filtered_alingments()
                read_aln = read_group.output_read()
                out_sam.write(read_aln)
                out_read_count += 1
    print('Writting %i reads from %i groups' %(out_read_count, in_read_count), file = sys.stderr)
    return 0
