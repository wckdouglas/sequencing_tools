#!/usr/bin/env python

from ..fastq_tools.pe_align import merge_interleaved
from ..io_tools import xopen
import sys


def getopt(subparsers):
    """
    reading input
    """
    descriptions = (
        "Merging interleaved, paired-end fastq file and output overlapped "
        + "regions only with error correction using cutadapt module to find overlapping regions"
    )
    parser = subparsers.add_parser("mergepe", description=descriptions)
    parser.add_argument(
        "-i", "--interleaved", default="-", help="Interleaved Fastq files (default: -)"
    )
    parser.add_argument(
        "-o", "--outfile", default="-", help="Merged fastq file (default: -)"
    )
    parser.add_argument(
        "-m",
        "--min_len",
        default=18,
        type=int,
        help="Minimum length of sequence to output (default: 18)",
    )
    parser.add_argument(
        "-e",
        "--error",
        default=0.1,
        type=float,
        help="Maximum error rate of alignment (default: 0.1)",
    )
    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        help="Output all bases (default: only overlapping regions)",
    )
    parser.add_argument(
        "--highlight",
        action="store_true",
        help="Highlight the non overlapping base (adapter sequences; only useful when --all is used)",
    )
    parser.add_argument(
        "-c",
        "--conserved",
        action="store_true",
        help="Use of a voting algorithm, " "otherwise use posterior error from qualit",
    )


def run(args):
    outfile = args.outfile
    outfile_handle = (
        sys.stdout
        if outfile == "-" or outfile == "/dev/stdin"
        else xopen(outfile, mode="w")
    )
    merge_interleaved(
        args.interleaved,
        outfile_handle,
        args.min_len,
        args.error,
        args.all,
        args.conserved,
        args.highlight,
    )
