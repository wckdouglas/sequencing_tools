#!/usr/bin/env python

from tgirt_seq_tools.pe_align import merge_interleaved
import argparse

def getOptions():
    '''
    reading input
    '''
    descriptions = 'Merging interleaved, paired-end fastq file and output overlapped '+ \
        'regions only with error correction using cutadapt module to find overlapping regions'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-i', '--infile', default='-',
        help='Interleaved Fastq files (default: -)')
    parser.add_argument('-o', '--outfile', default='-',
            help='Merged fastq file (default: -)')
    parser.add_argument('-m', '--min_len', default=18, type=int,
            help='Minimum length of sequence to output (default: 18)')
    parser.add_argument('-e', '--error', default=0.1, type=float,
            help='Maximum error rate of alignment (default: 0.1)')
    return parser.parse_args()

def main():
    args = getOptions()
    merge_interleaved(args.infile, args.outfile, 
            args.min_len, args.error)


if __name__ == '__main__':
    main()
