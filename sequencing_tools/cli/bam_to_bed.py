import logging
import os
import sys

from sequencing_tools.bam_tools.fragment_pairs import bam_to_bed
from sequencing_tools.utils import SeqUtilsError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="bam2bed", help="Making paired-end bam into bed file for every fragment"
    )
    parser.add_argument(
        "-i",
        "--in_bam",
        required=True,
        help="BAM file name, or stdin (-) ** name sorted",
    )
    parser.add_argument(
        "-o", "--out_bed", default="-", help="BED file output (default: - )"
    )
    parser.add_argument(
        "-m",
        "--min_size",
        default=10,
        type=int,
        help="minimum fragment size to report (default: 10)",
    )
    parser.add_argument(
        "-M",
        "--max_size",
        default=10000,
        type=int,
        help="minimum fragment size to report (default: 10000)",
    )
    parser.add_argument("-t", "--tag", default=None, help="tag to extract")
    parser.add_argument(
        "-a", "--all", action="store_true", help="supplementary as single fragment"
    )
    parser.add_argument("-p", "--primary", action="store_true", help="Only primary")
    parser.add_argument(
        "-c",
        "--add_cigar",
        action="store_true",
        help="Add cigar string on the last field",
    )
    parser.add_argument("--prefix", default=None, help="Prefix for read fragment name")


def run(args):
    in_bam = args.in_bam
    out_file = sys.stdout if args.out_bed == "-" else open(args.out_bed, "w")
    tag = args.tag
    if args.max_size <= args.min_size:
        raise SeqUtilsError("!!!!! Min fragment size > Max fragment size")
    bam_to_bed(
        in_bam,
        out_file,
        args.min_size,
        args.max_size,
        tag,
        args.all,
        args.primary,
        args.add_cigar,
        args.prefix,
    )
    return 0
