#!/usr/bin/env python

import pysam
import numpy as np
import sys
import argparse
import pyximport
import pysam
pyximport.install(setup_args = {'include_dirs':[np.get_include() ]+pysam.get_include()}) #pysam=0.9.0
from tgirt_seq_tools.filter_multi import processBam

def getopt():
    parser = argparse.ArgumentParser(description='Process multiply mapped bam and get either shortest insert size or ribosomal reads')
    parser.add_argument('-i','--infile', required=True, help = 'Input bam/sam file (use < - > for stdin)')
    parser.add_argument('-o','--outfile', required=True, help = 'Output bam/sam file (use < - > for stdout)' )
    parser.add_argument('-b','--bam_in', action='store_true', help = 'Provide this flag if bam instead of sam is used as input' )
    parser.add_argument('-z','--bam_out', action='store_true', help = 'Provide this flag if bam is needed for output')
    return parser.parse_args()

def main():
    args = getopt()
    in_bam = args.infile
    out_bam = args.outfile
    bam_in_bool = args.bam_in
    bam_out_bool = args.bam_out
    processBam(in_bam, out_bam, bam_in_bool, bam_out_bool)
    return 0

if __name__ == '__main__':
    main()
