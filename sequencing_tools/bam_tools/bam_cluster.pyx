from __future__ import print_function
from matplotlib import use as mpl_use
mpl_use('Agg')
import numpy as np
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool
import pysam
import array
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import defaultdict
import time
from sequencing_tools.bam_tools.read_cluster cimport readGroup
from sequencing_tools.fastq_tools.function_clip import hamming_distance


accpetable_flag = [99,147,
                163,83,
                77, 141,
                1123, 1171,
                1187, 1107]
cdef bool qualify_aln(AlignedSegment aln):
    '''
    Check if alignment is properly mapped
    '''

    return aln.flag in accpetable_flag

def add_data_member(member_count_list, member_count):
    '''
    make count table for member
    '''
    cdef:
        int m

    for m in member_count_list:
        member_count[m] += 1
    return member_count


def plot_bc_member(member_count, filename):
    '''
    Plotting Distribution of member count across all barcodes, as QC
    '''

    sns.set_style('white')
    figurename = filename.split('.')[0] + '.pdf'
    df = pd.DataFrame(list(member_count.items()),
                    columns=['members','counts']) \
        .assign(log_count = lambda d: np.log(d['counts']))
    with sns.plotting_context('paper', font_scale = 1.3):
        p = sns.FacetGrid(data = df)
    p.map(sns.barplot, 'members','log_count')
    p.set_xlabels('Member in barcode family')
    p.set_ylabels('Barcode count')
    p.savefig(figurename)
    print('Plotted: %s' %figurename, file = sys.stderr)


def cluster_bam(tag, bool conserved, AlignmentFile in_bam, out_fastq):
    '''
    Loop through alignments within in_bam file,
    extract barcode, and check if next alignment has the same barcode,
    if yes, put in as a list of alginments
    if no, use existing list to generate concensus sequence and output as fastq
    '''
    cdef:
        AlignedSegment aln
        int iter_count
        str name
        str barcode
        int group_count = 0
        int out_count = 0
        bool barcode_not_too_bad

    reference_list = in_bam.references
    member_count = defaultdict(int)
    for iter_count, aln in enumerate(in_bam):
        if qualify_aln(aln):
            if group_count == 0:
                read_group = readGroup(aln, tag)
                group_count += 1
            else:
                name = aln.query_name
                barcode = aln.get_tag(tag)
                #barcode_is_same = barcode == read_group.barcode
                barcode_not_too_bad =  hamming_distance(barcode, read_group.barcode) < 2
                if barcode_not_too_bad:
                    read_group.put_alignment(aln)
                else:
                    # conclude group and write
                    read_group.cluster(conserved)
                    read_group.generate_fastq_record()
                    out_fastq.write(read_group.fastq_record)
                    out_count += len(read_group.concensus_flag1)
                    member_count = add_data_member(read_group.member_count_list, member_count)
                    # reinitialize group
                    read_group = readGroup(aln, tag)
            if iter_count % 5000000 == 0 and iter_count != 0:
                print('Parsed %i alignments' %iter_count, file = sys.stderr)

    #spit out the rest in the stupid list
    read_group.cluster(conserved)
    read_group.generate_fastq_record()
    out_fastq.write(read_group.fastq_record)
    out_count += len(read_group.concensus_flag1)
    member_count = add_data_member(read_group.member_count_list, member_count)
    if out_fastq.name != '/dev/stdout':
        plot_bc_member(member_count, out_fastq.name)
    return out_count
