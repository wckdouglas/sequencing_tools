#!/usr/bin/env python

import os
import filecmp
from collections import defaultdict
from sequencing_tools.fastq_tools import readfq 

test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'

def test_bam():
    in_bam = test_data_path + '/test.bam' 
    out_bed = test_data_path + '/out.bed' 
    command = 'bam_to_bed.py -i {in_bam} | sort -k1,1 -k2,2n |  deduplicate_bed.py -i -  > {out_bed}'.format(in_bam = in_bam,
                                                                                                               out_bed = out_bed)
    os.system(command)
    assert(filecmp.cmp(out_bed, test_data_path + '/test.bed'))
    os.remove(out_bed)


#def test_multi():
#    in_bam = test_data_path + '/multi.bam'
#    out_bam = test_data_path + '/multi.out'
#    command = 'reduce_multi_reads.py -i  {in_bam} -o - > {out_bam}'.format(in_bam = in_bam,
#                                                                            out_bam = out_bam)
#    os.system(command)
#    assert(filecmp.cmp(out_bam, test_data_path + '/multi.result'))
#    os.remove(out_bam)


def same_fq(fq1, fq2):
    id_dict1 = defaultdict(set)
    for seqid, seq, qual in readfq(fq1):
        id_dict1[seqid].add(seq+qual)

    id_dict2 = defaultdict(set)
    for seqid, seq, qual in readfq(fq2):
        id_dict2[seqid].add(seq+qual)
    return id_dict1 == id_dict2
    

def test_correct():
    in_bam = test_data_path + '/tag.bam'
    out_fq = test_data_path + '/tag.fq'
    command = 'bam_read_cluster.py -i  {in_bam} -o {out_fq} -c -t RX'.format(in_bam = in_bam,
                                                                             out_fq = out_fq)
    os.system(command)
    assert(same_fq(out_fq, test_data_path + '/corrected.conserve.fq'))

    command = 'bam_read_cluster.py -i  {in_bam} -o {out_fq} -t RX'.format(in_bam = in_bam,
                                                                             out_fq = out_fq)
    os.system(command)
    assert(same_fq(out_fq, test_data_path + '/corrected.qual.fq'))
    os.remove(out_fq)

def test_filter():
    in_bam = test_data_path + '/tag.bam'
    out_bam = test_data_path + '/filtered.out'
    command = 'filter_soft_clip.py  --pe -s 0 -i {in_bam} -o - > {out_bam}'.format(in_bam = in_bam,
                                                                            out_bam = out_bam)
    os.system(command)
    assert(filecmp.cmp(out_bam, test_data_path + '/clipped.bam'))
    os.remove(out_bam)