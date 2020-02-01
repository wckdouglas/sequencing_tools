#!/usr/bin/env python

import sys
from ..bam_tools.umi_tag import filter_umi

def getopt(subparsers):
    parser = subparsers.add_parser(name = 'filterUMI',description = 'Filter out alignments with UMI having nucleotide runs > 3nt')
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file')
    parser.add_argument('-o','--outbam', default = '-', help = 'output interleaved fastq file (defulat: - )')
    parser.add_argument('-c','--consecutive_bases',default=3, help = 'Threshold for nucleotide runs (defualt: 3, only output alignments with UMI without AAAA,TTTT,CCCT,GGGG)')
    parser.add_argument('-t','--tag',default=None, help = 'BAM tag storing UMI')

def run(args):
    filter_umi(args.inbam, args.outbam, args.consecutive_bases, args.tag)