#!/usrbin/env python

import pysam
import sys
import re
from ..bam_tools.split_bam_tools import split_bam_pair, split_bam_single


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="splitBam",
        description="Splitting bam file to uniquely mapped and multiple mapped (only support bowtie2 and hisat2)",
    )
    parser.add_argument(
        "-a",
        "--aligner",
        choices=["bowtie2", "hisat2"],
        default="bowtie2",
        help="PRogram generating the bam file (default: bowtie2)",
    )
    parser.add_argument(
        "-i", "--inBam", default="-", help="bam file name, or stdin (default: stdin)"
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        required=True,
        help="output prefix ($OUTPUTPREFIX.unique.bam, $OUTPUTPREFIX.multi.bam)",
    )
    parser.add_argument(
        "-s", "--single_end", action="store_true", help="is single end bam"
    )


def run(args):
    outprefix = args.outprefix
    uniqueBam = outprefix + ".unique.bam"
    multiBam = outprefix + ".multi.bam"
    with pysam.Samfile(args.inBam, "rb") as inbam:
        with pysam.Samfile(uniqueBam, "wb", template=inbam) as uniquebam, pysam.Samfile(
            multiBam, "wb", template=inbam
        ) as multibam:
            if args.single_end:
                split_bam_single(inbam, uniquebam, multibam, args.aligner)
            else:
                split_bam_pair(inbam, uniquebam, multibam, args.aligner)
