#!/usr/bin/env python

from __future__ import print_function

import logging
import os
import sys

from sequencing_tools.fastq_tools import read_interleaved
from sequencing_tools.io_tools import xopen

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="deinterleave", description="Deinterleaving a interleaved fastq file"
    )
    parser.add_argument("-i", "--infq", default="-", help="interleaved fq (defualt: -)")
    parser.add_argument("-1", "--read1", required=True, help="read 1 fastq output")
    parser.add_argument("-2", "--read2", help="read 2 fastq output")
    parser.add_argument(
        "-m",
        "--min_length",
        help="Seq length cutoff (default: 10)",
        default=10,
        type=int,
    )


def run(args):
    fastq_in = args.infq
    fastq_forward = args.read1
    fastq_reverse = args.read2
    r1_count, r2_count = 0, 0
    logger.info("Reading from %s " % fastq_in)
    logger.info("Writing to read1: %s" % (fastq_forward))
    logger.info("Writing to read2: %s" % (fastq_reverse))

    with xopen(fastq_forward, mode="w") as r1, xopen(fastq_reverse, mode="w") as r2:
        seq_id_1 = ""
        in_file = xopen(fastq_in, mode="r") if fastq_in != "-" else sys.stdin
        for R1, R2 in read_interleaved(in_file):

            if len(R1.seq) >= args.min_length and len(R2.seq) >= args.min_length:

                print("@%s\n%s\n+\n%s" % (R1.id, R1.seq, R1.qual), file=r1)
                r1_count += 1
                print("@%s\n%s\n+\n%s" % (R2.id, R2.seq, R2.qual), file=r2)
                r2_count += 1

    assert r1_count == r2_count, "Not equal reads!!!! %i != %i" % (r1_count, r2_count)
    logger.info("Written %i records" % (r1_count))
