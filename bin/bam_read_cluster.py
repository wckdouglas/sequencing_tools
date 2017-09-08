#!/usr/bin/env python

import pysam
import argparse
from tools_umi.bam_cluster import cluster_bam
import time
import sys
import glob
from functools import partial

def getopt():
    parser = argparse.ArgumentParser(description = 'Make clustre from BAM file that contains alignments with barcode on their ID')
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file')
    parser.add_argument('-o','--outfq', default = '-', help = 'output interleaved fastq file (defulat: - )')
    parser.add_argument('-c','--conserved', action='store_true', help = 'Use of a voting algorithm')
    parser.add_argument('-t','--tag',default='MI', help = 'UMI tag on bam alignment')
    args = parser.parse_args()
    return args

def main():
    args = getopt()
    in_bam = args.inbam
    out_fastq = args.outfq
    start = time.time()

    sys.stderr.write('Demultiplexing: %s\n' %in_bam)
    out_fastq = '/dev/stdout' if out_fastq == '-' else out_fastq
    with pysam.Samfile(in_bam, 'rb') as inbam:
        with open(out_fastq, 'w') as outfastq:
            out_count = cluster_bam(args.tag, args.conserved, inbam, outfastq)
    sys.stderr.write('Finished clustering: output %i clusters in %.3f min\n' %(out_count, (time.time() - start) /60))


if __name__ == '__main__':
    main()
