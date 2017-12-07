#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from sequencing_tools.fastq_tools import read_interleaved


def main():
    if len(sys.argv) != 4:
        print('usage: python %s <interleaved.fq> <f.fq> <r.fq>\n' %sys.argv[0], file = sys.stderr)
        sys.exit(0)

    fastq_in = sys.argv[1]
    fastq_forward = sys.argv[2]
    fastq_reverse = sys.argv[3]
    r1_count, r2_count = 0, 0
    print('Reading from %s ' %fastq_in, file = sys.stdout)
    print('Writing to %s and %s' %(fastq_forward, fastq_reverse), file = sys.stdout)

    with open(fastq_forward.rstrip('.gz'),'w') as r1, \
            open(fastq_reverse.rstrip('.gz'),'w') as r2:
        seq_id_1 = ''
        in_file = open(fastq_in) if fastq_in != '-' else sys.stdin
        for R1, R2 in read_interleaved(in_file):

            print('@%s\n%s\n+\n%s' %(R1.id, R1.seq, R1.qual), file = r1)
            r1_count += 1
            print('@%s\n%s\n+\n%s' %(R2.id, R2.seq, R2.qual), file = r2)
            r2_count += 1

    assert r1_count == r2_count, 'Not equal reads!!!! %i != %i' %(r1_count, r2_count)
    if fastq_forward.endswith('.gz'):
        os.system('gzip -f %s' %(fastq_forward.rstrip('.gz')))
        os.system('gzip -f %s' %(fastq_reverse.rstrip('.gz')))
    print('Written %i records' %(r1_count), file = sys.stdout)


if __name__ == '__main__':
    main()
