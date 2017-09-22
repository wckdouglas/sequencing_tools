#!/usr/bin/env python

import argparse
import numpy as np
import time
import sys
import re
from sys import stderr
from sequencing_tools.function_clip import run_pairs, run_pairs_stdout


def getOptions():
    '''
    reading input
    '''
    descriptions = 'Clip the barcode sequence and attached to the front of seq id'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-o', '--outputprefix', default='-',
        help='Interleaved Fastq files (default: -)')
    parser.add_argument('-1', '--fastq1', required=True,
        help='Paired end Fastq file 1 with four line/record')
    parser.add_argument('-2', '--fastq2',required=True,
        help='Paired end Fastq file 2 with four line/record')
    parser.add_argument("-x", "--idxBase", default='XXXXXXXXXXXXX',
        help="how many base in 5' end as index? (default: XXXXXXXXXXXXX) "+\
            "X as umi bases can add constant regions add back, such as XXXXCATGC")
    parser.add_argument('-q', '--barcodeCutOff', type=int, default=20,
        help="Average base calling quality for barcode sequence (default=20)")
    parser.add_argument("-a", "--mismatch", type=int,default=1,
        help="Allow how many mismatch in constant region (deflaut: 1)")
    parser.add_argument("-s", "--prefix_split", type=int,default=0, choices = range(5),
        help="Using how many bases on the barcode to split the fastq? A choice of 3 will generate 4^3 = 64 files (deflaut: 4)")
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
    outputprefix = args.outputprefix
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idx_base = args.idxBase
    barcode_cut_off = args.barcodeCutOff
    constant = args.idxBase.lstrip('X')
    idx_base = len(re.findall('X+', idx_base)[0])
    allow_mismatch = args.mismatch
    prefix_split = args.prefix_split
    UMI_side = args.read

    #print out parameters
    programname = sys.argv[0]
    stderr.write('[%s] [Parameters] \n' %(programname))
    stderr.write('[%s] indexed bases:                     %i\n' %(programname, idx_base))
    stderr.write('[%s] min mean barcode quality:          %i\n' %(programname, barcode_cut_off))
    stderr.write('[%s] outputPrefix:                      %s\n' %(programname, outputprefix))
    stderr.write('[%s] using constant regions:            %s\n' %(programname, constant))
    stderr.write('[%s] allowed mismatches:                %i\n' %(programname, allow_mismatch))
    stderr.write('[%s] Using prefix bases:                %i\n' %(programname, prefix_split))
    stderr.write('[%s] Using UMI side:                    %s\n' %(programname, UMI_side))

    # divide reads into subclusters
    if outputprefix != '-':
        run_pairs(outputprefix, inFastq1, inFastq2, idx_base,
                barcode_cut_off, constant, allow_mismatch, programname,
                prefix_split, UMI_side)
    else:
        stderr.write('[%s] Using STDOUT, Will not split prefix!!\n' %(programname))
        run_pairs_stdout(inFastq1, inFastq2, idx_base,
                barcode_cut_off, constant, allow_mismatch, programname,
                prefix_split, UMI_side)
    stderr.write('[%s] time lapsed:      %2.3f min\n' %(programname, np.true_divide(time.time()-start,60)))
    return 0

if __name__ == '__main__':
    args = getOptions()
    main(args)
