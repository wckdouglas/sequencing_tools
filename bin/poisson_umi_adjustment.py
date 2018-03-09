#!/usr/bin/env python

import sys
from sequencing_tools.bam_tools import poisson_umi_tools 



def getopt():
    parser = argparse.ArgumentParser(description = 'Making paired-end bam into bed file for every fragment')
    parser.add_argument('-i', '--in_bed', required=True,
                        help = 'BED file name, or stdin (-) ** name sorted' )
    parser.add_argument('-o','--out_bed', default='-', help = 'BED file output (default: - )')
    parser.add_argument('--umi', type=int, help = 'Number of nucleotide as umi (default: 6)', default=6)

    return parser.parse_args()

def main():
    args = getopt()
    in_bed = open(sys.stdin) if args.in_bed in ['-', '/dev/stdin'] else open(args.in_bed,'r')
    out_bed = open(sys.stdout) if args.out_bed in ['-','/dev/stdout'] else open(args.out_bed, 'w')
    nt = args.umi
    poisson_umi_tools.parse_dedup(in_bed, out_bed, nt)

if __name__ == '__main__':
    main()
