#!/usr/bin/env python

from __future__ import print_function
import pysam
import numpy as np
from pyfaidx import Fasta
from functools import partial
from collections import defaultdict
from itertools import izip
import re
import os
import sys
import string
import argparse
from tgirt_seq_tools.pileup_errors import extract_bases, analyze_region

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
    args = parser.parse_args()
    return args

def make_cigar_seq(cigar_numbers, cigar_operator):
    for num, op in zip(cigar_numbers, cigar_operator):
        if op != 'S':
            yield int(num)*op

def make_regions(chromosome_length, how_many_bases_to_look_at):
    start = 0
    end = start + how_many_bases_to_look_at
    while end < chromosome_length:
        yield (start, end)
        start = end
        end = end + how_many_bases_to_look_at
    yield (start, chromosome_length)


def output_table(fa, chromosome, base_dict, start, end):
    for i, base in enumerate(fa.get_seq(chromosome,start+1,end+1)):
        pos = i + start
        coverage, base_count_string = extract_bases(base_dict, pos)
        if coverage > 0:
            outline =  '%s\t%i\t%s\t%s' %(chromosome, pos, base, base_count_string)
            print(outline+'\n', file = sys.stdout)
    return 0

def analyze_chromosome(chromosome, in_bam, fa, bases_region, qual_threshold, crop):
    chrom_length = len(fa[chromosome])
    get_error = partial(analyze_region, in_bam, chromosome, qual_threshold, crop)
    output = partial(output_table, fa, chromosome)
    region_generator = make_regions(chrom_length, bases_region)
    for i, (start, end) in enumerate(region_generator):
        base_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        aln_count, base_dict = get_error(base_dict, start, end)
        out = output(base_dict, start, end)
        if i % 10 == 0:
            print('Written %s:%i-%i with %i alignments' %(chromosome, start, end, aln_count), file=sys.stderr)

def analyze_bam(in_bam, fa, bases_region, qual_threshold, crop):
    chromosomes = fa.keys()
    header = 'chrom\tpos\tbase\t'
    header = header + 'A+\tC+\tG+\tT+\tA-\tC-\tG-\tT-'
    print(header, file=sys.stdout)
    for chromosome in chromosomes:
        analyze_chromosome(chromosome, in_bam, fa, bases_region, qual_threshold, crop)

def main():
    args = getopt()
    bam_file = args.bam
    ref_fasta = args.fasta
    bases_region = args.bases
    qual_threshold = args.qual
    crops = args.crop
    with pysam.Samfile(bam_file, 'rb') as in_bam, \
            Fasta(ref_fasta) as fa:
        analyze_bam(in_bam, fa, bases_region, qual_threshold, crop)


if __name__ == '__main__':
    main()
