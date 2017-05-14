#!/usr/bin/env python

import pyximport
pyximport.install()
from tgirt_seq_tools.bigwig import parse_depth_bed
import argparse


def getopt():
    parser = argparse.ArgumentParser(description = 'Splitting genomecov bed into chromosomes and write to bigwig')
    parser.add_argument('-i', '--infile', default='-',
                        help = 'genomecov -d file, or stdin (default: <->)' )
    parser.add_argument('-o','--outprefix', required=True,
                        help='output prefix ($OUTPUTPREFIX.$CHROM.bigWig)')
    parser.add_argument('-g','--genome', required=True,
                        help = 'genome file as genomeCov needed')
    return parser.parse_args()


def main():
    args = getopt()
    parse_depth_bed(args.infile, args.genome, args.outprefix)

if __name__ == '__main__':
    main()
