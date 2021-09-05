#!/usr/bin/env python

import argparse

from ..bam_tools.bigwig import parse_depth_bed


def getopt(subparser):
    parser = subparser.add_parser(
        name="cov2bw",
        description="Splitting genomecov bed into chromosomes and write to bigwig",
    )
    parser.add_argument(
        "-i", "--infile", default="-", help="genomecov -d file, or stdin (default: <->)"
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        required=True,
        help="output prefix ($OUTPUTPREFIX.$CHROM.bigWig)",
    )
    parser.add_argument(
        "-g", "--genome", required=True, help="genome file as genomeCov needed"
    )


def run(args):
    parse_depth_bed(args.infile, args.genome, args.outprefix)
