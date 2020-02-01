#!/usr/bin/env python

from __future__ import print_function
import pysam
import numpy as np
from functools import partial
from collections import defaultdict
import re
import os
import sys
import string
import argparse
from sequencing_tools.bam_tools.pileup_errors import extract_bases, analyze_region, make_regions
from sequencing_tools.io_tools import xopen
from operator import itemgetter
import six
long = six.integer_types[-1]

def getopt():
    parser = argparse.ArgumentParser(description='Pileup whole genome, only output bases where coverage > 0')
    parser.add_argument('-i', '--bam',required=True, help='Input bam file (indexed)')
    parser.add_argument('-f','--fasta', required=True, help='reference fasta file')
    parser.add_argument('-b','--bases', default = 100000,
                    type=int, help='number of bases to look at every iteration (default: 100000)')
    parser.add_argument('-q','--qual', default=30,
                    type=int, help='base quality to filter (defulat: 30)')

    parser.add_argument('-c','--crop', default=0,
                    type=int, help='Crop how many bases from ends (defulat: 0)')
    parser.add_argument('-r','--bed', default='', help='bed file for regions (default: whole genome)')
    parser.add_argument('--no_indel', action='store_true', help='Not considering alignments with Indel')
    parser.add_argument('--min_coverage', default = 0, type = int,
                        help='Minimum coverage to output')
    parser.add_argument('--concensus_fasta', help = 'Generate concensus fasta (only work if bed file is provided)')
    args = parser.parse_args()
    return args

def bed_generator(bed_file):
    bed = sys.stdin if bed_file == '-' else xopen(bed_file)
    for line in bed:
        fields = line.split('\t')
        chrom, start, end = itemgetter(0,1,2)(fields)
        yield chrom, long(start), long(end)
    
def base_dict_to_fa(base_dict):
    fa_dict = defaultdict(str)
    for pos, pos_dict in six.iteritems(base_dict):
        for strand, strand_dict in six.iteritems(pos_dict):
            base = np.array(list(strand_dict.keys()))
            base_count = np.array(list(strand_dict.values()))
            concensus_base = base[base_count.argmax()] 
            fa_dict[strand] += concensus_base
    return fa_dict


def output_table(fa, chromosome, base_dict, start, end, min_cov):
    '''
    output table: chrom, pos, base, A+, C+, G+, T+, A-, C-, G-, T-
    '''
    for i, base in enumerate(fa.fetch(reference = chromosome,
                                      start=start,
                                      end=end)):
        pos = i + start
        coverage_dict, base_count_string = extract_bases(base_dict, pos)
        if coverage_dict['+'] > min_cov or coverage_dict['-'] > min_cov:
            outline =  '%s\t%i\t%s\t%s' %(chromosome, pos, base, base_count_string)
            print(outline, file = sys.stdout)
    return 0

def analyze_chromosome(chromosome, in_bam, fa, bases_region, qual_threshold, crop, no_indel, min_cov):
    chrom_length = fa.get_reference_length(chromosome)
    get_error = partial(analyze_region, in_bam, chromosome, qual_threshold, crop, no_indel)
    output = partial(output_table, fa, chromosome)
    region_generator = make_regions(chrom_length, bases_region)
    for i, (start, end) in enumerate(region_generator):
        aln_count, base_dict = get_error(start, end)
        out = output(base_dict, start, end, min_cov)
        if i % 10 == 0:
            print('Written %s:%i-%i with %i alignments' %(chromosome, start, end, aln_count), file=sys.stderr)

def analyze_bam(in_bam, fa, bases_region, qual_threshold, crop, 
                bed_file, use_bed, no_indel, min_cov, concensus_fasta):
    chromosomes = fa.references
    chromosomes.sort()
    header = 'chrom\tpos\tbase\t'
    header = header + 'A+\tC+\tG+\tT+\tA-\tC-\tG-\tT-'
    print(header, file=sys.stdout)
    if use_bed:
        con_fa = xopen(concensus_fasta,'w') if concensus_fasta else None
        for chrom, start, end in bed_generator(bed_file):
            aln_count, base_dict = analyze_region(in_bam, chrom, qual_threshold, crop, no_indel, start, end)
            out = output_table(fa, chrom, base_dict, start, end, min_cov)
            if con_fa:
                fa_dict = base_dict_to_fa(base_dict)
                for strand in ['+','-']:
                    seq = fa_dict[strand]
                    seq = seq[::-1] if strand == '-' else seq
                    seq_record = '>{chrom}:{start}-{end}({strand})\n{seq}' \
                            .format(chrom = chrom, 
                                    start = start, 
                                    end = end,
                                    strand = strand,
                                    seq = seq)
                    print(seq_record, file = con_fa)
        if con_fa:
            con_fa.close()
    else:
        for chromosome in chromosomes:
            analyze_chromosome(chromosome, in_bam, fa, bases_region, qual_threshold, crop, no_indel, min_cov)


def main():
    args = getopt()
    bam_file = args.bam
    ref_fasta = args.fasta
    bases_region = args.bases
    qual_threshold = args.qual
    crop = args.crop
    use_bed = True if args.bed != '' else False
    no_indel = args.no_indel
    min_cov = args.min_coverage
    concensus_fasta = args.concensus_fasta

    with pysam.Samfile(bam_file, 'rb') as in_bam, \
            pysam.FastaFile(ref_fasta) as fa:
        analyze_bam(in_bam, fa, bases_region, qual_threshold, crop, args.bed, use_bed, no_indel, min_cov, concensus_fasta)


if __name__ == '__main__':
    main()
