#!/usr/bin/env python

import pysam
from ..bam_tools.split_bam_tools import split_N_bam

def getopt(subparsers):
    parser = subparsers.add_parser(name = 'spliceBam',
                                description = 'Splitting paired-end bam into two bam files: '+\
                                                    '1.  $OUT_PREFIX.fragment.bam (no N in cigar string on neither read1 nor read2) and ' +\
                                                    '2.  $OUT_PREFIX.split.bam (with N in cigar string on neither read1 nor read2)')
    parser.add_argument('-i', '--in_bam', required=True,
                        help = 'BAM file name, or stdin (-) ** name sorted' )
    parser.add_argument('-o','--out_prefix', required=True, help = 'bam file output prefix')


def run(args):
    out_prefix = args.out_prefix
    split_bam = out_prefix + '.split.bam'
    other_bam = out_prefix + '.fragment.bam'
    with pysam.Samfile(args.in_bam, 'rb') as inbam:
        with pysam.Samfile(split_bam,'wb',template = inbam) as out_split_bam, \
                pysam.Samfile(other_bam, 'wb', template = inbam) as out_fragment_bam:
            split_N_bam(inbam, out_split_bam, out_fragment_bam)