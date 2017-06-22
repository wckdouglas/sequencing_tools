from itertools import izip, imap, product
from scipy.spatial.distance import hamming
from multiprocessing import Pool
from functools import partial
import numpy as np
import gzip
cimport numpy as np
from sys import stderr
from cpython cimport bool
import io
import os

# define fastq record type
cdef class fastqRecord:
    cdef:
        public str id
        public str seq
        public str qual

    def __init__(self, str id, str seq, str qual):
        self.id = id
        self.seq = seq
        self.qual = qual


def readfq(fp): # this is a generator function
    '''
        https://github.com/lh3/readfq/blob/master/readfq.py
    '''
    cdef:
        str l, name, seq

    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield fastqRecord(name, ''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield fastqRecord(name, seq, ''.join(seqs)) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


cpdef float hammingDistance(str expected_constant, str constant_region):
    '''
    Calculating hamming distance from two strings
    '''
    cdef float dist = hamming(list(expected_constant),list(constant_region))
    return dist


def clip_read1(barcode_cut_off, constant, constant_no_evaluation, prefix_split,
            idx_base, usable_seq, float hamming_threshold,
            fastqRecord read1, fastqRecord read2):
    """
    from each read1, clipped barcode and constant regions and make it as ID
    only when 1. barcode quality > cutoff and
              2. constant region matched
    """

    cdef:
        float barcode_mean_qual
        bool no_N_barcode, hiQ_barcode, accurate_constant
        str prefix, seq_record


    seq_name = read1.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read1.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read1.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(map(ord, barcode_qual)) - 33
    prefix = barcode[:prefix_split]

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual > barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hammingDistance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_left = read1.seq[usable_seq:]
        qual_left = read1.qual[usable_seq:]
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, seq_left, qual_left) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, read2.seq, read2.qual)
        return 1, prefix, seq_record
    else:
        return 0, None, None


def clip_read2(barcode_cut_off, constant, constant_no_evaluation, prefix_split,
            idx_base, usable_seq, float hamming_threshold,
            fastqRecord read1, fastqRecord read2):
    """
    from each read1, clipped barcode and constant regions and make it as ID
    only when 1. barcode quality > cutoff and
              2. constant region matched
    """

    cdef:
        float barcode_mean_qual
        bool no_N_barcode, hiQ_barcode, accurate_constant
        str prefix, seq_record


    seq_name = read2.id.split(' ')[0]
#    assert seq_name == id_right.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(id_right, id_left)
    barcode =  read2.seq[:idx_base]
    barcode_qual = read1.qual[:idx_base]
    constant_region = read2.seq[idx_base:usable_seq]
    barcode_mean_qual = np.mean(map(ord, barcode_qual)) - 33
    prefix = barcode[:prefix_split]

    no_N_barcode = 'N' not in barcode
    hiQ_barcode = barcode_mean_qual > barcode_cut_off
    accurate_constant = True if constant_no_evaluation else  hammingDistance(constant, constant_region) <= hamming_threshold

    if no_N_barcode and hiQ_barcode and accurate_constant:
        seq_right = read2.seq[usable_seq:]
        qual_right = read2.qual[usable_seq:]
        seq_record = '@%s_%s/1\n%s\n+\n%s\n' %(barcode, seq_name, read1.seq, read1.qual) +\
                    '@%s_%s/2\n%s\n+\n%s' %(barcode, seq_name, seq_right, qual_right)
        return 1, prefix, seq_record
    else:
        return 0, None, None


def open_files(output_prefix, prefix_split):
    '''
    making prefixed fastq files
    '''
    barcode_prefix = product(['A','C','T','G'],repeat=prefix_split)
    barcode_prefix = [''.join(prefix) for prefix in barcode_prefix]
    if prefix_split > 0:
        file_dict = {prefix: open('%s_%s.fq' %(output_prefix, prefix), 'w') for prefix in barcode_prefix}
    else:
        file_dict = {prefix: open('%s.fq' %(output_prefix), 'w') for prefix in barcode_prefix}
    return file_dict

def gzopen(filename, read_flag = 'rb'):
    if 'r' in read_flag:
        return os.popen('zcat '+ filename)
    elif 'w' in read_flag:
        return open(filename, read_flag)


def run_pairs(outputprefix, inFastq1, inFastq2, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname,
            prefix_split, read):
    '''
    Feed in paried end fastq file
    loop over pairs and generate interleaved fastq
    '''
    cdef:
        int usable_seq
        float hamming_threshold
        str out_R1, out_R2
        int out, count, out_count = 0
        fastqRecord read1, read2
        str prefix, seq_record

    constant_length = len(constant)
    constant_no_evaluation = constant_length < 1
    hamming_threshold = float(allow_mismatch)/constant_length if not constant_no_evaluation else 0
    constant_no_evaluation = constant_length < 2
    usable_seq = idx_base if constant_no_evaluation else idx_base + constant_length


    file_dict = open_files(outputprefix, prefix_split)
    with gzopen(inFastq1, read_flag = 'rb') as in1, gzopen(inFastq2, read_flag = 'rb') as in2:
        if read == 'read1':
            clipping = partial(clip_read1, barcode_cut_off, constant,
                            constant_no_evaluation, prefix_split, idx_base,
                            usable_seq, hamming_threshold)

        elif read == 'read2':
            clipping = partial(clip_read2, barcode_cut_off, constant,
                            constant_no_evaluation, prefix_split, idx_base,
                            usable_seq, hamming_threshold)


        for count, (read1, read2) in enumerate(izip(readfq(in1), readfq(in2))):
            out, prefix, seq_record = clipping(read1, read2)
            out_count += out
            if out == 1:
                file_dict[prefix].write(seq_record + '\n')
            if count % 10000000 == 0 and count != 0:
                stderr.write('[%s] Parsed %i records\n'%(programname, count))
    [f.close() for f in file_dict.values()]

    stderr.write('[%s] Parsed:           %i sequences\n' %(programname, count))
    stderr.write('[%s] Output:           %i sequences\n' %(programname, out_count))
    return 0


def run_pairs_stdout(inFastq1, inFastq2, idx_base,
            barcode_cut_off, constant, allow_mismatch, programname,
            prefix_split, read):
    '''
    Feed in paried end fastq file
    loop over pairs and generate interleaved fastq
    '''
    cdef:
        int usable_seq
        float hamming_threshold
        str out_R1, out_R2
        int out, count, out_count = 0
        fastqRecord read1, read2
        str prefix, seq_record

    constant_length = len(constant)
    constant_no_evaluation = constant_length < 1
    hamming_threshold = float(allow_mismatch)/constant_length if not constant_no_evaluation else 0
    usable_seq = idx_base if constant_no_evaluation else idx_base + constant_length


    with gzopen(inFastq1, read_flag = 'rb') as in1, gzopen(inFastq2, read_flag = 'rb') as in2:
        if read == 'read1':
            clipping = partial(clip_read1, barcode_cut_off, constant,
                            constant_no_evaluation, prefix_split, idx_base,
                            usable_seq, hamming_threshold)

        elif read == 'read2':
            clipping = partial(clip_read2, barcode_cut_off, constant,
                            constant_no_evaluation, prefix_split, idx_base,
                            usable_seq, hamming_threshold)

        for count, (read1, read2) in enumerate(izip(readfq(in1), readfq(in2))):
            out, prefix, seq_record = clipping(read1, read2)
            out_count += out
            if out == 1:
                print seq_record
            if count % 10000000 == 0 and count != 0:
                stderr.write('[%s] Parsed %i records\n'%(programname, count))

    stderr.write('[%s] Parsed:           %i sequences\n' %(programname, count))
    stderr.write('[%s] Output:           %i sequences\n' %(programname, out_count))
    return 0
