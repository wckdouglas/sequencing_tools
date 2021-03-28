#!/usr/bin/env python

import os
import filecmp
from collections import defaultdict
import logging
from sequencing_tools.fastq_tools import readfq
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(os.path.basename(__file__))

PROG_PREFIX = ""
if os.environ["_"].endswith("poetry"):
    PROG_PREFIX = "poetry run "
    logger.info('Using poetry')
else:
    logger.info('Using python')

test_data_path = os.path.dirname(os.path.realpath(__file__)) + "/data"


def test_bam():
    in_bam = test_data_path + "/test.bam"
    out_bed = test_data_path + "/out.bed"
    command = (
        "{PREFIX} seqtools bam2bed -i {in_bam} --primary  "
        "| sort -k1,1 -k2,2n -k3,3n "
        "| {PREFIX} seqtools dedup -i - "
        "> {out_bed}".format(in_bam=in_bam, out_bed=out_bed, PREFIX=PROG_PREFIX)
    )
    os.system(command)
    assert filecmp.cmp(out_bed, test_data_path + "/test.bed")
    os.remove(out_bed)


def test_multi():
    in_bam = test_data_path + "/multi.bam"
    out_bam = test_data_path + "/multi.out"
    command = "{PREFIX} seqtools filterMulti -i  {in_bam} -o - | samtools view > {out_bam}".format(
        in_bam=in_bam, out_bam=out_bam, PREFIX=PROG_PREFIX
    )
    os.system(command)
    assert filecmp.cmp(out_bam, test_data_path + "/multi.result")
    os.remove(out_bam)


def same_fq(fq1, fq2):
    id_dict1 = defaultdict(set)
    for seqid, seq, qual in readfq(fq1):
        id_dict1[seqid].add(seq + qual)

    id_dict2 = defaultdict(set)
    for seqid, seq, qual in readfq(fq2):
        id_dict2[seqid].add(seq + qual)
    return id_dict1 == id_dict2


def test_correct():
    in_bam = test_data_path + "/tag.bam"
    out_fq = test_data_path + "/tag.fq"
    command = "{PREFIX} seqtools demux -i  {in_bam} -o {out_fq} -c -t RX".format(
        in_bam=in_bam, out_fq=out_fq, PREFIX=PROG_PREFIX
    )
    os.system(command)
    assert same_fq(out_fq, test_data_path + "/corrected.conserve.fq")

    command = "{PREFIX} seqtools demux -i  {in_bam} -o {out_fq} -t RX".format(
        in_bam=in_bam, out_fq=out_fq, PREFIX=PROG_PREFIX
    )
    os.system(command)
    assert same_fq(out_fq, test_data_path + "/corrected.qual.fq")
    os.remove(out_fq)


def test_filter():
    in_bam = test_data_path + "/tag.bam"
    out_bam = test_data_path + "/filtered.out"
    command = "{PREFIX} seqtools filterSoftClip  --pe -s 0 -i {in_bam} -o - | samtools view > {out_bam}".format(
        in_bam=in_bam, out_bam=out_bam, PREFIX=PROG_PREFIX
    )
    os.system(command)
    assert filecmp.cmp(out_bam, test_data_path + "/clipped.result")
    os.remove(out_bam)


def test_stranded_base_count():
    golden_file = test_data_path + "/pileup.txt"
    command = (
        "{PREFIX} seqtools pileup -i {path}/MT_TF.bam "
        "-f {path}/MT_TF.fa -c 0 --min_coverage 0 -q 0  "
        "> {path}/test_pileup.txt".format(path=test_data_path, PREFIX=PROG_PREFIX)
    )

    os.system(command)
    assert filecmp.cmp(golden_file, test_data_path + "/test_pileup.txt")
