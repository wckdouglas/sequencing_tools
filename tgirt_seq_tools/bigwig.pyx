from cpython cimport array
import pyBigWig as pbw
import numpy as np
import array
import sys


class chrom_depth:

    def __init__(self):
        self.values = []
        self.chrom = ''
        self.bw = None
        self.filename = ''
        self.positions = []


    def initiate_bigwig(self, prefix, chrom, genome):
        self.filename = prefix + '.' + chrom + '.bigWig'
        self.bw = pbw.open(self.filename, 'w')
        self.chrom = chrom
        chrom_size = genome[chrom]
        self.bw.addHeader([(chrom, chrom_size)])
        self.values=array.array('i',np.zeros(chrom_size, dtype='int'))

    def add_value(self, int position, int value):
        self.values[position] = value

    def write_bw(self):
        self.positions = range(len(self.values))
        self.values = map(float, self.values)
        self.bw.addEntries(self.chrom, self.positions, values=self.values, span=1)
        self.bw.close()
        print 'Written %s' %self.filename


def parse_depth_bed(bed_file, genome_file, output_prefix):
    cdef:
        str init_chrom = ''
        str line
        str chrom, position, value

    genome = {}
    for line in open(genome_file,'r'):
        fields = line.strip().split()
        genome[fields[0]] = long(fields[1])

    file_handle = open(bed_file,'r') if bed_file not in ['-','dev/stdin'] else sys.stdin
    for line in file_handle:
        fields = line.rstrip().split('\t')
        chrom, position, value = fields
        if chrom != init_chrom:
            if init_chrom:
                depth.write_bw()
            depth = chrom_depth()
            chrom_size = genome[chrom]
            depth.initiate_bigwig(output_prefix, chrom, genome)
            init_chrom = chrom
            depth.add_value(long(position),int(value))
        else:
            depth.add_value(long(position),int(value))
    depth.write_bw()
    return 0
