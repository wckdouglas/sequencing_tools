#!/usr/bin/env python


import numpy as np
import argparse
import pysam
from sequencing_tools.bam_splitter import filter_bam
import gzip
import sys

def getopt():
    parser = argparse.ArgumentParser(description = 'Filter alignments from sam file by softclip ratio')
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file (can be stdin, use -)')
    parser.add_argument('-o','--outbam', default = '-', help = 'output bam file (defulat: - )')
    parser.add_argument('-s','--single_end',default = 0.2, type=float,
                        help ='Maximum ratio of the whole alignment being clipped in one end (default: 0.2)')
    parser.add_argument('-b','--both_end',default = 0.5, type=float,
                        help ='Maximum ratio of the whole alignment being clipped in sum(each end) (default : 0.5)')
    parser.add_argument('-v','--inverse', action = 'store_true',
                        help ='Only output alignment with clipped base > threshold (like grep -v)')
    args = parser.parse_args()
    return args

def main():
    args = getopt()
    in_bam = args.inbam
    out_bam = args.outbam
    single_end_thresh = args.single_end
    both_end_thresh = args.both_end
    print('Filtering %s' %in_bam, file = sys.stdout)
    output_count = filter_bam(in_bam, out_bam, single_end_thresh, both_end_thresh, args.inverse)
    print('Written %i alignments to %s' %(output_count, out_bam), file = sys.stdout)

if __name__ == '__main__':
    main()
