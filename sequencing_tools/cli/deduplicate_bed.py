#!/usr/bin/env python

import argparse
import logging
import os
import sys

from ..bam_tools.bed_dedup import dedup_bed

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


def getopt(subparser):
    help_message = (
        "Demultiplexing UMI bed file with the follwing columns:\n"
        + "1. chrom name\n2. start\n3. end\n4. {$UMI}_{$READ_ID}\n5. score\n6.strand\n"
        + "The program internally used hamming distance matrix of the "
        + "barcodes to generate connected network and identified UMI clusters"
    )
    parser = subparser.add_parser(name="dedup", description=help_message)
    parser.add_argument(
        "-i",
        "--infile",
        help='input bedfile, can be "-" or "/dev/stdin" for stdin (default: -).'
        + "Format should be sorted BED: "
        + "\n1. chrom\n2. start\n3. end\n4. {barcode}_{id}\n6. strand",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help='output bedfile, can be "-" or "/dev/stdout" for stdin (default: -).',
        default="-",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="How many error between barcodes can be tolerated? (default = 1)",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-d",
        "--delim",
        help="Deliminator separating read id and bc (default: _ )",
        default="_",
    )
    parser.add_argument(
        "-f",
        help="after splitting read name using {delim}, which fragmnet is UMI? can use -1 as last piece (default: 0)",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--ct", help="cigar tag field (0-based)? (default: None)", default=-1, type=int
    )


def run(args):
    in_filename = args.infile
    out_filename = args.outfile
    threshold = args.threshold

    ## input
    if in_filename != "-" and in_filename != "/dev/stdin":
        in_file_handle = open(in_filename, "r")
        logger.info("Using file %s..." % in_filename)
    else:
        in_file_handle = sys.stdin
        logger.info("Using stdin...")

    ## input
    if out_filename != "-" and out_filename != "/dev/stdin":
        out_file_handle = open(out_filename, "w")
        logger.info("Writing file %s..." % out_filename)
    else:
        out_file_handle = sys.stdout
        logger.info("Writing to stdout...")

    dedup_bed(in_file_handle, out_file_handle, threshold, args.delim, args.f, args.ct)


if __name__ == "__main__":
    main()
