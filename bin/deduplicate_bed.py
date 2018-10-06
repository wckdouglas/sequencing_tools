#!/usr/bin/env python

from __future__ import print_function
from sequencing_tools.bam_tools.bed_dedup import  dedup_bed
import sys
import argparse


def get_opt():
    help_message = 'Demultiplexing UMI bed file with the follwing columns:\n'+\
            '1. chrom name\n2. start\n3. end\n4. {$UMI}_{$READ_ID}\n5. score\n6.strand\n'+\
            'The program internally used hamming distance matrix of the '+\
            'barcodes to generate connected network and identified UMI clusters'
    parser = argparse.ArgumentParser(description=help_message)
    parser.add_argument('-i','--infile',
                        help ='input bedfile, can be "-" or "/dev/stdin" for stdin (default: -).' +\
                                'Format should be sorted BED: ' + \
                                '\n1. chrom\n2. start\n3. end\n4. {barcode}_{id}\n6. strand',
                        required=True)
    parser.add_argument('-o','--outfile',
                        help ='output bedfile, can be "-" or "/dev/stdout" for stdin (default: -).', default= '-')
    parser.add_argument('-t','--threshold',
                        help='How many error between barcodes can be tolerated? (default = 1)',
                        type=int, default = 1)
    parser.add_argument('-d','--delim',
                        help='Deliminator separating read id and bc (default: _ )',
                        default='_')
    parser.add_argument('-f',
                        help='after splitting read name using {delim}, which fragmnet is UMI? can use -1 as last piece (default: 0)',
                        default=0, type=int)
    parser.add_argument('--ct',
                        help='cigar tag field (0-based)? (default: None)',
                        default=-1, type=int)
    args = parser.parse_args()
    return args


def main():
    args = get_opt()
    in_filename = args.infile
    out_filename = args.outfile
    threshold = args.threshold

    ## input
    if in_filename != '-' and in_filename != '/dev/stdin':
        in_file_handle = open(in_filename, 'r')
        print('[Deduplicate BED] Using file %s...' %in_filename, file=sys.stderr)
    else:
        in_file_handle = sys.stdin
        print('[Deduplicate BED] Using stdin...', file=sys.stderr)

    ## input
    if out_filename != '-' and out_filename != '/dev/stdin':
        out_file_handle = open(out_filename, 'r')
        print('[Deduplicate BED] Writing file %s...' %out_filename, file=sys.stderr)
    else:
        out_file_handle = sys.stdout
        print('[Deduplicate BED] Writing to stdout...', file=sys.stderr)

    dedup_bed(in_file_handle, out_file_handle, threshold, args.delim, args.f, args.ct)
    return 0


if __name__ == '__main__':
    main()
