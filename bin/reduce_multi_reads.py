#!/usr/bin/env python

import pysam
import numpy as np
import sys
import argparse
from sequencing_tools.bam_tools.filter_multi import process_pair_bam, process_single_bam
from functools import partial

def getopt():
    parser = argparse.ArgumentParser(description='Process multiply mapped bam and get either shortest insert size or ribosomal reads')
    parser.add_argument('-i','--infile', required=True, help = 'Input bam/sam file (use < - > for stdin)')
    parser.add_argument('-o','--outfile', required=True, help = 'Output bam/sam file (use < - > for stdout)' )
    parser.add_argument('-b','--bam_in', action='store_true', help = 'Provide this flag if bam instead of sam is used as input' )
    parser.add_argument('-z','--bam_out', action='store_true', help = 'Provide this flag if bam is needed for output')
    parser.add_argument('-s','--single_end', action='store_true', help = 'For single end bam files')
    return parser.parse_args()

def main():
    args = getopt()
    in_bam = args.infile
    out_bam = args.outfile
    bam_in_bool = args.bam_in
    bam_out_bool = args.bam_out
    processBam = partial(process_pair_bam) if not args.single_end else partial(process_single_bam)
    processBam(in_bam, out_bam, bam_in_bool, bam_out_bool)
    return 0

if __name__ == '__main__':
    main()
