#!/usr/bin/env python

import pysam
import time
import sys
import glob
import os
from functools import partial
import logging
from ..bam_tools.bam_cluster import cluster_bam

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="demux",
        description="Make clustre from BAM file that contains alignments with barcode on their ID",
    )
    parser.add_argument("-i", "--inbam", required=True, help="input bam file")
    parser.add_argument(
        "-o", "--outfq", default="-", help="output interleaved fastq file (defulat: - )"
    )
    parser.add_argument(
        "-c",
        "--conserved",
        action="store_true",
        help="Use of a voting algorithm to correct read base errors",
    )
    parser.add_argument("-t", "--tag", default="MI", help="UMI tag on bam alignment")


def run(args):
    in_bam = args.inbam
    out_fastq = args.outfq
    start = time.time()

    logger.info("Demultiplexing: %s" % in_bam)
    out_fastq = "/dev/stdout" if out_fastq == "-" else out_fastq
    with pysam.Samfile(in_bam, "rb") as inbam:
        with open(out_fastq, "w") as outfastq:
            out_count = cluster_bam(args.tag, args.conserved, inbam, outfastq)
    logger.info(
        "Finished clustering: output %i clusters in %.3f min"
        % (out_count, (time.time() - start) / 60)
    )
