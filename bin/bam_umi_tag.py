#!/usr/bin/env python

from sequencing_tools.bam_tools.umi_tag import add_umi_tag
import argparse

def getopt():
    parser = argparse.ArgumentParser(description = 'Putting UMI as RX tag in bam, enables picard Markduplicates with BARCODE_TAG=RX')
    parser.add_argument('-i', '--in_bam', required=True,
                        help = 'BAM file name, or stdin (-)' )
    parser.add_argument('-o','--out_bam', default='-', help = 'BAM file output (default: - )')
    parser.add_argument('-t','--tag', default='RX', help = 'Tag id (default: RX )')
    parser.add_argument('-d','--delim', default='_', help = 'Deliminator separating read id and bc (default: _ )')
    parser.add_argument('-f','--fragment', default='0', 
                        help = 'after splitting read name using {delim}, which' +\
                            'fragment is UMI? can use -1 as last piece (default: 0)')

    return parser.parse_args()

def main():
    args = getopt()
    add_umi_tag(args.in_bam, args.out_bam, args.tag, args.delim, int(args.fragment))
    return 0

if __name__ == '__main__':
    main()
