#!/usr/bin/env python

import sys
import os

def read_fastq(file_fq):
    """
    takes a fastq file as input
    yields idSeq, sequence and score
    for each fastq entry
    http://codereview.stackexchange.com/questions/32897/efficient-parsing-of-fastq
    """

    line_count = 0

    while True:

        line = file_fq.readline()
        line_count += 1
        #break if we hit the end of the file
        if not line:
            break

        if line_count == 1:
            idSeq = line.strip().lstrip('@')

        elif line_count == 2:
            sequence = line.strip()

        elif line_count == 4:
            score = line.strip()
            line_count = 0
            yield idSeq,sequence, score



def open_file(fastqfile):
    if fastqfile.endswith('gz'):
        return os.popen("zcat %s" %fastqfile)
    else:
        return open(fastqfile,'r')


def main():
    if len(sys.argv) != 4:
        sys.exit('usage: python %s <interleaved.fq> <f.fq> <r.fq>' %sys.argv[0])
    fastq_in = sys.argv[1]
    fastq_forward = sys.argv[2]
    fastq_reverse = sys.argv[3]
    r1_count, r2_count = 0, 0
    print 'Reading from %s ' %fastq_in
    print 'Writing to %s and %s' %(fastq_forward, fastq_reverse)
    with open(fastq_forward.rstrip('.gz'),'w') as r1, \
            open(fastq_reverse.rstrip('.gz'),'w') as r2:
        seq_id_1 = ''
        in_file = open_file(fastq_in) if fastq_in != '-' else sys.stdin
        for seq_id, seq, qual in read_fastq(in_file):
            _id = seq_id.split('/')[0]
            if seq_id.endswith('/1'):
                seq_id_1 = _id
                r1.write('@%s\n%s\n+\n%s\n' %(_id, seq, qual))
                r1_count += 1
            else:
                assert _id == seq_id_1, 'Wrong ID!!! %s || %s' %(_id, seq_id_1)
                r2.write('@%s\n%s\n+\n%s\n' %(_id, seq, qual))
                r2_count += 1
    assert r1_count == r2_count, 'Not equal reads!!!! %i != %i' %(r1_count, r2_count)
    if fastq_forward.endswith('.gz'):
        os.system('gzip -f %s' %(fastq_forward.rstrip('.gz')))
        os.system('gzip -f %s' %(fastq_reverse.rstrip('.gz')))
    print 'Written %i records' %(r1_count)


if __name__ == '__main__':
    main()
