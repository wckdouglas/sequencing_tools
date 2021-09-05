#!/usr/bin/env python

import sys

from ..bam_tools import poisson_umi_tools


def getopt(subparsers):
    parser = subparsers.add_parser(
        "calcUMI",
        description="Adjusting fragment count for UMI saturation, "
        + "see paper: Counting individual DNA molecules by the stochastic attachment of diverse labels. ",
    )
    parser.add_argument(
        "-i",
        "--in_bed",
        required=True,
        help="BED file name, or stdin (-) ** name sorted",
    )
    parser.add_argument(
        "-o", "--out_bed", default="-", help="BED file output (default: - )"
    )
    parser.add_argument(
        "--umi", type=int, help="Number of nucleotide as umi (default: 6)", default=6
    )
    parser.add_argument(
        "--prefix",
        help="Prefix adding to fragment name (default: Sample)",
        default="Sample",
    )


def run(args):
    in_bed = sys.stdin if args.in_bed in ["-", "/dev/stdin"] else open(args.in_bed, "r")
    out_bed = (
        sys.stdout if args.out_bed in ["-", "/dev/stdout"] else open(args.out_bed, "w")
    )
    nt = args.umi
    read_prefix = args.prefix
    res = poisson_umi_tools.parse_dedup(in_bed, out_bed, nt, read_prefix)
