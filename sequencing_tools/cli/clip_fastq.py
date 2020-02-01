#!/usr/bin/env python

from __future__ import print_function
import argparse
import numpy as np
import time
import sys
import re
from sequencing_tools.fastq_tools.function_clip import clip_pairs


def getOptions():
    '''
    reading input
    '''
    descriptions = 'Clip the barcode sequence and attached to the front of seq id'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-o', '--out_file', default='-',
        help='Interleaved Fastq files (default: -)')
    parser.add_argument('-1', '--fastq1', required=True,
        help='Paired end Fastq file 1 with four line/record')
    parser.add_argument('-2', '--fastq2',required=True,
        help='Paired end Fastq file 2 with four line/record')
    parser.add_argument("-x", "--idxBase", default='XXXXXXXXXXXXX',
        help="how many base in 5' end as index? (default: XXXXXXXXXXXXX) "+\
            "X as umi bases, can also add constant regions at the back, such as XXXXCATGC, if CATGC is the constant region")
    parser.add_argument('-q', '--barcodeCutOff', type=int, default=20,
        help="Average base calling quality for barcode sequence (default=20)")
    parser.add_argument("-a", "--mismatch", type=int,default=1,
        help="Allow how many mismatch in constant region (deflaut: 1)")
    parser.add_argument("-m", "--min_length", type=int,default=15,
        help="Minimum length to report (deflaut: 15)")
#    parser.add_argument("-s", "--prefix_split", type=int,default=0, choices = range(5),
#        help="Using how many bases on the barcode to split the fastq? A choice of 3 will generate 4^3 = 64 files (deflaut: 4)")
    parser.add_argument("-r", "--read", default='read1',choices = ['read1','read2'],
        help="Which read is the UMI on?")
    args = parser.parse_args()
    return args


def main(args):
    """
    main function:
        controlling work flow
        1. generate read clusters by reading from fq1 and fq2
        2. obtain concensus sequence from read clusters
        3. writing concensus sequence to files
    """
    start = time.time()
    out_file = args.out_file
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idx_base = args.idxBase
    barcode_cut_off = args.barcodeCutOff
    constant = args.idxBase.lstrip('X')
    idx_base = len(re.findall('X+', idx_base)[0])
    allow_mismatch = args.mismatch
    UMI_side = args.read

    #print out parameters
    programname = sys.argv[0]
    print('[%s] [Parameters] ' %(programname), file = sys.stderr)
    print('[%s] indexed bases:                     %i' %(programname, idx_base), file = sys.stderr)
    print('[%s] min mean barcode quality:          %i' %(programname, barcode_cut_off), file = sys.stderr)
    print('[%s] output file:                       %s' %(programname, out_file), file = sys.stderr)
    print('[%s] using constant regions:            %s' %(programname, constant), file = sys.stderr)
    print('[%s] allowed mismatches:                %i' %(programname, allow_mismatch), file = sys.stderr)
    print('[%s] Using UMI side:                    %s' %(programname, UMI_side), file = sys.stderr)
    print('[%s] Using min-length:                  %i' %(programname, args.min_length), file = sys.stderr)

    # divide reads into subclusters
    clip_pairs(inFastq1, inFastq2, out_file, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname, UMI_side, args.min_length)
    print('[%s] time lapsed:      %2.3f min' %(programname, np.true_divide(time.time()-start,60)), file = sys.stderr)
    return 0

if __name__ == '__main__':
    args = getOptions()
    main(args)
