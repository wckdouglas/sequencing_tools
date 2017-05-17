#!/usr/bin/env python

import pysam
import argparse
import pyximport
pyximport.install()
from tgirt_seq_tools.split_N_bam import parse_bam

def getopt():
    parser = argparse.ArgumentParser(description = 'Splitting paired-end bam into two bam files:\n'+\
                                                    '1.  $OUT_PREFIX.fragment.bam:    no N in cigar string on neither read1 nor read2\n' +\
                                                    '1.  $OUT_PREFIX.split.bam:    with N in cigar string on neither read1 nor read2')
    parser.add_argument('-i', '--in_bam', required=True,
                        help = 'BAM file name, or stdin (-) ** name sorted' )
    parser.add_argument('-o','--out_prefix', required=True, help = 'bam file output prefix')
    return parser.parse_args()


def main():
    args = getopt()

    out_prefix = args.out_prefix
    split_bam = out_prefix + '.split.bam'
    other_bam = out_prefix + '.fragment.bam'
    with pysam.Samfile(args.in_bam, 'rb') as inbam:
        with pysam.Samfile(split_bam,'wb',template = inbam) as out_split_bam, \
                pysam.Samfile(other_bam, 'wb', template = inbam) as out_fragment_bam:
            parse_bam(inbam, out_split_bam, out_fragment_bam)
    return 0


if __name__ == '__main__':
    main()
