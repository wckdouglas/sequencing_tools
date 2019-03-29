from __future__ import print_function
import numpy as np
from matplotlib import use as mpl_use
mpl_use('Agg')  # Must be before importing matplotlib.pyplot or pylab
import matplotlib.pyplot as plt
import sys
import ujson
import re
from multiprocessing import Pool
from builtins import zip, map, range
from functools import partial
from cpython cimport bool
import six
from sequencing_tools.stats_tools import cy_mean, hamming_distance
from sequencing_tools.fastq_tools._fastq_tools cimport fastqRecord
from sequencing_tools.fastq_tools import readfq
from sequencing_tools.io_tools import xopen
from sequencing_tools.consensus_tools import ErrorCorrection


def plotBCdistribution(barcode_family_count, outputprefix):
    #plotting inspection of barcode distribution
    barcode_family_count = np.array(barcode_family_count, dtype=np.int64)
    hist, bins = np.histogram(barcode_family_count[barcode_family_count<50],bins=50)
    centers = (bins[:-1] + bins[1:]) / 2
    width = 0.7 * (bins[1] - bins[0])
    figurename = '%s.pdf' %(outputprefix)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(centers,hist,align='center',width=width)
    ax.set_xlabel("Number of occurence")
    ax.set_ylabel("Count of tags")
    ax.set_yscale('log',nonposy='clip')
    ax.set_title(outputprefix.split('/')[-1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(figurename)
    print('Plotted %s.' %figurename, file = sys.stderr)
    return 0



def dictToJson(barcode_dict, json_file):
    with open(json_file,'w') as f:
        [f.write(ujson.encode(items) + '\n') for items in six.iteritems(barcode_dict)]
    print('written %s' %(json_file) + '\n', file = sys.stderr)
    return 0


def writeSeqToFiles(read1, read2, output_cluster_count, result):
    if result!='No':
        print('@cluster%i_%s' %(output_cluster_count, result[0]), file = read1)
        print('@cluster%i_%s' %(output_cluster_count, result[1]), file = read2)
        return 1
    else:
        return 0


class Clustering():
    def __init__(self, outputprefix, min_family_member_count, json_file,
                        threads, fraction_threshold):
    # From index library, generate error free reads
    # using multicore to process read clusters

        self.read1File = outputprefix + '_R1_001.fastq.gz'
        self.read2File = outputprefix + '_R2_001.fastq.gz'
        self.min_family_member_count = min_family_member_count
        self.json_file = json_file
        self.threads = threads
        self.fraction_threshold = fraction_threshold
        self.correction_module = ErrorCorrection(mode = 'vote', 
                                threshold = fraction_threshold)
    
    def run(self):
        cdef:
            int counter = 0
            int output_cluster_count = 0

        pool = Pool(self.threads,maxtasksperchild=10000)
        with xopen(self.read1File,'w') as read1, \
                xopen(self.read2File,'w') as read2, 
                open(self.json_file,'r') as infile:
            write_func = partial(writeSeqToFiles,read1, read2)
            processes = pool.imap_unordered(self.__ErrorFreeReads__, infile, chunksize = 10000)
            #processes = map(error_func, infile) ## debug
            for result in processes:
                output_cluster_count += write_func(output_cluster_count, result)
                counter += 1
                if counter % 1000000 == 0:
                    print('Processed %i read clusters.' %(counter), file = sys.stderr)
        pool.close()
        pool.join()
        return output_cluster_count, read1File, read2File


    def __ErrorFreeReads__(self, json_record):
        """
        main function for getting concensus sequences from read clusters.
        return  a pair of concensus reads with a 4-line fastq format
        see functions: 1. filterRead,
                    2. concensusPairs,
                    3. calculateConcensusBase
        """
        # skip if not enough sequences to perform voting
        cdef:
            str index
            int member_count
            str sequence_left, quality_left
            str sequence_right, quality_right
            str left_record, right_record

        record = ujson.decode(json_record)
        index = str(record[0])
        table = np.array(record[1], dtype=str)
        member_count = table.shape[0]
        if member_count >= min_family_member_count:
            sequence_left, quality_left, sequence_right, quality_right = self.__ConcensusPairs__(table[:,0], 
                                                                                            table[:,2],
                                                                                            table[:,1]
                                                                                            table[:,3])
            left_record = '%s_%i_readCluster\n%s\n+\n%s' %(index, member_count, sequence_left, quality_left)
            right_record = '%s_%i_readCluster\n%s\n+\n%s' %(index, member_count, sequence_right, quality_right)
            return left_record, right_record
        else:
            return 'No'


    def __ConcensusPairs__(self, seq_left_list, qual_left_list, seq_right_list, qual_right_list):
        """ given a pair of reads as defined as the class: seqRecord
        return concensus sequence and mean quality of the pairs,
            as well as the number of reads that supports the concnesus pairs
        see function: concensusSeq, calculateConcensusBase
        """
        #extract table
        cdef:
            str sequence_left, quality_left
            str sequence_right, quality_right

        # get concensus left reads first
        sequence_left, quality_left = self.correction_module.Correct(seq_left_list, qual_left_list)
        # get concensus right reads
        sequence_right, quality_right = self.correction_module.Correct(seq_right_list, qual_right_list)
        return sequence_left, quality_left, sequence_right, quality_right

############### clustering #####################


cpdef int readClusteringR2(barcode_dict, idx_base, barcode_cut_off, constant,
                   constant_length, hamming_threshold, usable_seq, failed_file,
                   low_complexity_composition, read1, read2):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idx_bases/ of read 1
    """
    cdef:
        str id_left, seq_left, qual_left
        str id_right, seq_right, qual_right
        str barcode, constant_region
        int barcodeQualmean
        bool no_N_barcode, is_low_complexity_barcode, hiQ_barcode, accurate_constant


    id_left, seq_left, qual_left = read1.id, read1.seq, read1.qual
    id_right, seq_right, qual_right = read2.id, read2.seq, read2.qual
    assert id_left.split(' ')[0] == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode = seq_right[:idx_base]
    constant_region = seq_right[idx_base:usable_seq]
    barcodeQualmean = int(cy_mean(map(ord, qual_right[:idx_base])) - 33)

    no_N_barcode = 'N' not in barcode
    is_low_complexity_barcode = bool(low_complexity_composition.search(barcode))
    hiQ_barcode = barcodeQualmean > barcode_cut_off
    accurate_constant = hamming_distance(constant, constant_region) <= hamming_threshold
    min_qual_left = min(map(ord, qual_left))
    min_qual_right = min(map(ord, qual_right))
    qual_pass = np.min([min_qual_left,min_qual_right])  >= 53

    if no_N_barcode and hiQ_barcode and accurate_constant: #and not is_low_complexity_barcode):
        seq_right = seq_right[usable_seq:]
        qual_right = qual_right[usable_seq:]
        barcode_dict[barcode].append([seq_left,seq_right,qual_left, qual_right])
        return 0
    else:
        failed_file.write('\t'.join([id_left, seq_left, qual_left, seq_right, qual_right]) + '\n')
        return 1

cpdef int readClusteringR1(barcode_dict, idx_base, barcode_cut_off, constant,
                     constant_length, hamming_threshold, usable_seq, failed_file,
                     low_complexity_composition, read1, read2):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idx_bases/ of read 1
    """
    cdef:
        str id_left, seq_left, qual_left
        str id_right, seq_right, qual_right
        str barcode, constant_region
        int barcode_qual_min
        bool no_N_barcode, is_low_complexity_barcode, hiQ_barcode, accurate_constant


    id_left, seq_left, qual_left = read1.id, read1.seq, read1.qual
    id_right, seq_right, qual_right = read2.id, read2.seq, read2.qual
    assert id_left.split(' ')[0] == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode = seq_left[:idx_base]
    constant_region = seq_left[idx_base:usable_seq]
    barcode_qual_min = int(min(map(ord,qual_left[:idx_base])) - 33)

    no_N_barcode = 'N' not in barcode
    is_low_complexity_barcode = bool(low_complexity_composition.search(barcode))
    hiQ_barcode = barcode_qual_min >= barcode_cut_off
    accurate_constant = hamming_distance(constant, constant_region) <= hamming_threshold
    min_qual_left = min(map(ord, qual_left))
    min_qual_right = min(map(ord, qual_right))
    qual_pass = min([min_qual_left,min_qual_right]) >= 53


    if no_N_barcode and hiQ_barcode and accurate_constant and qual_pass: #and not low_complexity_barcode:
        seq_left = seq_left[usable_seq:]
        qual_left = qual_left[usable_seq:]
        barcode_dict[barcode].append([seq_left, seq_right, qual_left, qual_right])
        return 0
    return 1

def recordsToDict(str outputprefix, str inFastq1, str inFastq2, int idx_base, int barcode_cut_off,
                str constant, barcode_dict, int allow_mismatch, str which_side, str programname):

    cdef:
        int discarded_sequence_count = 0
        int constant_length = len(constant)
        int usable_seq = idx_base + constant_length
        int mul = 6
        str failed_reads
        int read_num
        fastqRecord read1, read2

    low_complexity_base = ['A' * mul,'C' * mul,'T' * mul,'G' * mul]
    low_complexity_composition = re.compile('|'.join(low_complexity_base))

    failed_reads = outputprefix + '-failed.tsv'
    with xopen(inFastq1,'r') as fq1, xopen(inFastq2,'r') as fq2, open(failed_reads,'w') as failed_file:


        if which_side == 'read2':
            cluster_reads = partial(readClusteringR2, barcode_dict, idx_base, barcode_cut_off,
                            constant, constant_length, allow_mismatch, usable_seq,
                            failed_file, low_complexity_composition)
        elif which_side == 'read1':
            cluster_reads = partial(readClusteringR1, barcode_dict, idx_base, barcode_cut_off,
                            constant, constant_length, allow_mismatch, usable_seq,
                            failed_file, low_complexity_composition)

        iterator = enumerate(zip(readfq(fq1), readfq(fq2)))
        for read_num, (read1,read2) in iterator:
            discarded_sequence_count += cluster_reads(read1, read2)
            if read_num % 10000000 == 0:
                print('[%s] Parsed: %i sequence' %(programname,read_num), file = sys.stderr)

    barcode_count = len(barcode_dict.keys())
    return barcode_dict, read_num, barcode_count, discarded_sequence_count


