#!/usr/bin/env python

from tgirt_seq_tools.umi_tag import add_umi_tag
import argparse

def getopt():
    parser = argparse.ArgumentParser(description = 'Putting UMI as RX tag in bam, enables picard Markduplicates with BARCODE_TAG=RX')
    parser.add_argument('-i', '--in_bam', required=True,
                        help = 'BAM file name, or stdin (-) ** name sorted' )
    parser.add_argument('-o','--out_bam', default='-', help = 'BAM file output (default: - )')
    return parser.parse_args()

def main():
    args = getopt()
    add_umi_tag(args.in_bam, args.out_bam)
    return 0

if __name__ == '__main__':
    main()
