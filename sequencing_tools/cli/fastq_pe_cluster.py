import logging
import os
import time
from collections import defaultdict
from sys import stderr

import numpy as np

from sequencing_tools.fastq_tools.cluster_reads import (
    Clustering,
    dictToJson,
    plotBCdistribution,
    recordsToDict,
)
from sequencing_tools.utils import SeqUtilsError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))


def getopt(subparser):
    """
    reading input
    """
    descriptions = (
        "Clustering fastq reads to fasta reads with the first $idx_base bases as cDNA-synthesis barcode. "
        + "Concensus bases are called only when the fraction of reads that contain the concensus base exceed some threshold. "
        + "Quality scores are generated by the average score for the bases that matched concensus base. "
    )
    parser = subparser.add_parser(name="correctFQ", description=descriptions)
    parser.add_argument(
        "-o",
        "--outputprefix",
        required=True,
        help="Paired end Fastq files with R1_001.fastq.gz as suffix for read1, and R2_001.fastq.gz as suffix for read2",
    )
    parser.add_argument(
        "-1",
        "--fastq1",
        required=True,
        help="Paired end Fastq file 1 with four line/record",
    )
    parser.add_argument(
        "-2",
        "--fastq2",
        required=True,
        help="Paired end Fastq file 2 with four line/record",
    )
    parser.add_argument(
        "-m",
        "--cutoff",
        type=int,
        default=4,
        help="minimum read count for each read cluster (default: 4)",
    )
    parser.add_argument(
        "-x",
        "--idxBase",
        type=int,
        default=13,
        help="how many base in 5' end as index? (default: 13)",
    )
    parser.add_argument(
        "-q",
        "--barcodeCutOff",
        type=int,
        default=30,
        help="Minimum base calling quality for barcode sequence (default=30)",
    )
    parser.add_argument(
        "-c",
        "--constant_region",
        default="CATCG",
        help="Constant sequence after tags (default: CATCG ,e.g. Douglas's index-R1R)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Threads to use (deflaut: 1)"
    )
    parser.add_argument(
        "-a",
        "--mismatch",
        type=int,
        default=1,
        help="Allow how many mismatch in constant region (deflaut: 1)",
    )
    parser.add_argument(
        "-r",
        "--read",
        required=True,
        choices=["read1", "read2"],
        help="barcode on first N bases of which read from pair end",
    )
    parser.add_argument(
        "-f",
        "--fraction",
        default=0.66,
        type=float,
        help="Fraction of base to call a concensus base",
    )


def clustering(
    outputprefix,
    inFastq1,
    inFastq2,
    idx_base,
    min_family_member_count,
    barcode_cut_off,
    constant,
    threads,
    allow_mismatch,
    which_side,
    fraction_threshold,
):
    json_file = outputprefix + ".json"
    barcode_dict = defaultdict(list)
    result = recordsToDict(
        outputprefix,
        inFastq1,
        inFastq2,
        idx_base,
        barcode_cut_off,
        constant,
        barcode_dict,
        allow_mismatch,
        which_side,
        programname,
    )
    barcode_dict, read_num, barcode_count, discarded_sequence_count = result
    logger.info("Extracted: %i barcode group" % (barcode_count + 1))
    logger.info("discarded: %i sequences" % (discarded_sequence_count))
    logger.info("Parsed:    %i seqeucnes" % (read_num))

    barcode_member_counts = list(
        map(lambda index: len(barcode_dict[index]), barcode_dict.keys())
    )
    p = plotBCdistribution(barcode_member_counts, outputprefix)
    dictToJson(barcode_dict, json_file)
    barcode_dict.clear()

    clustering = Clustering(
        outputprefix, min_family_member_count, json_file, threads, fraction_threshold
    )

    output_cluster_count, read1File, read2File = clustering.run()
    # all done!
    logger.info("Finished writing error free reads")
    logger.info("[Summary]                        ")
    logger.info("read1:                     %s" % (read1File))
    logger.info("read2:                     %s" % (read2File))
    logger.info("output clusters:           %i" % (output_cluster_count))
    logger.info(
        "Percentage retained:       %.3f"
        % (float(output_cluster_count) / read_num * 100)
    )
    return 0


def run(args):
    """
    main function:
        controlling work flow
        1. generate read clusters by reading from fq1 and fq2
        2. obtain concensus sequence from read clusters
        3. writing concensus sequence to files
    """
    start = time.time()
    outputprefix = args.outputprefix
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idx_base = args.idxBase
    min_family_member_count = args.cutoff
    barcode_cut_off = args.barcodeCutOff
    constant = args.constant_region
    threads = args.threads
    allow_mismatch = args.mismatch
    which_side = args.read
    fraction_threshold = args.fraction
    if fraction_threshold >= 1:
        raise SeqUtilsError("Fraction cannot > 1")

    # print out parameters
    logger.info("[Parameters] ")
    logger.info("indexed bases:                     %i" % (idx_base))
    logger.info("minimum coverage:                  %i" % (min_family_member_count))
    logger.info("min mean barcode quality:          %i" % (barcode_cut_off))
    logger.info("outputPrefix:                      %s" % (outputprefix))
    logger.info("threads:                           %i" % (threads))
    logger.info("using constant regions:            %s" % (constant))
    logger.info("allowed mismatches:                %i" % (allow_mismatch))
    logger.info("Fraction to call concnesus:        %.2f" % (fraction_threshold))

    # divide reads into subclusters
    clustering(
        outputprefix,
        inFastq1,
        inFastq2,
        idx_base,
        min_family_member_count,
        barcode_cut_off,
        constant,
        threads,
        allow_mismatch,
        which_side,
        fraction_threshold,
    )
    logger.info(
        "time lapsed:      %2.3f min" % (np.true_divide(time.time() - start, 60))
    )
