#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from sequencing_tools.umi_tag import filter_umi

def getopt():
    parser = argparse.ArgumentParser(description = 'Filter out alignments with UMI having nucleotide runs > 3nt')
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file')
    parser.add_argument('-o','--outbam', default = '-', help = 'output interleaved fastq file (defulat: - )')
    parser.add_argument('-c','--consecutive_bases',default=3, help = 'Threshold for nucleotide runs (defualt: 3, only output alignments with UMI without AAAA,TTTT,CCCT,GGGG)')
    parser.add_argument('-t','--tag',default=None, help = 'BAM tag storing UMI')
    args = parser.parse_args()
    return args

def main():
    args = getopt()
    filter_umi(args.inbam, args.outbam, args.consecutive_bases, args.tag)


if __name__ == '__main__':
    main()