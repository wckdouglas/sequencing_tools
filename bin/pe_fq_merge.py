#!/usr/bin/env python

from tgirt_seq_tools.pe_align import merge_interleaved
import argparse
import sys

def getOptions():
    '''
    reading input
    '''
    descriptions = 'Merging interleaved, paired-end fastq file and output overlapped '+ \
        'regions only with error correction using cutadapt module to find overlapping regions'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-i', '--interleaved', default='-',
        help='Interleaved Fastq files (default: -)')
    parser.add_argument('-o', '--outfile', default='-',
            help='Merged fastq file (default: -)')
    parser.add_argument('-m', '--min_len', default=18, type=int,
            help='Minimum length of sequence to output (default: 18)')
    parser.add_argument('-e', '--error', default=0.1, type=float,
            help='Maximum error rate of alignment (default: 0.1)')
    parser.add_argument('-a','--all',action='store_true', 
        help='Output all bases (default: only overlapping regions)')
    return parser.parse_args()

def main():
    args = getOptions()
    outfile=args.outfile
    outfile_handle = sys.stdout if outfile == '-' or outfile == '/dev/stdin' else open(outfile,'w')
    merge_interleaved(args.interleaved, outfile_handle, 
            args.min_len, args.error, args.all)


if __name__ == '__main__':
    main()
