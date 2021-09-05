#!/usr/bin/env python

import argparse
import sys

from ..bam_tools.fragment_pairs import process_bedpe


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="bedpe2bed", description="To convert bedpe file to bed file"
    )
    parser.add_argument(
        "-i",
        default="-",
        help="Input bedpe file, must have no unmapped bedpe record (default: -)",
    )
    parser.add_argument("-o", default="-", help="Output bed file (default: -)")
    parser.add_argument(
        "--min", default=10, type=int, help="Minimum length of fragment (default = 10)"
    )
    parser.add_argument(
        "--max",
        default=10000,
        type=int,
        help="Maximum length of fragment (default = 10000)",
    )


def run(args):
    bed_file = args.i
    out_file = args.o
    min_length = args.min
    max_length = args.max
    bed_iterator = (
        open(bed_file, "r")
        if bed_file != "/dev/stdin" and bed_file != "-"
        else sys.stdin
    )
    out_handle = (
        open(out_file, "w")
        if out_file != "-" and out_file != "/dev/stdout"
        else sys.stdout
    )
    done = process_bedpe(bed_iterator, min_length, max_length, out_handle)
