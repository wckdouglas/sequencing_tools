import fileinput
from operator import itemgetter
import numpy as np


'''
usage: cat bed12 | python bed12_to_intron.py  > intron.bed
'''


class Bed12Record():
    def __init__(self, line):
        '''
        object for analyzing bed12 line
        '''
        fields = line.strip().split('\t')
        chrom, start, end, tid, exon_count, \
            strand, exon_starts, exon_ends = itemgetter(0,1,2,3,4,5,8,9)(fields)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.exon_count = int(exon_count)
        self.tid = tid
        self.exon_starts = np.array(list(map(int, exon_starts.strip(',').split(','))))
        self.exon_ends = np.array(list(map(int, exon_ends.strip(',').split(','))))
        self.exon_length = np.array([e-s for e,s in zip(self.exon_ends, self.exon_starts)])
        if self.strand == '-':
            self.exon_length = self.exon_length[::-1]
            self.reversed_exon_starts = self.exon_ends[::-1]
        self.cumulative_exon_starts = np.append([0], np.cumsum(self.exon_length)[:-1])
        self.transcript_length = self.end - self.start
        assert(len(self.exon_starts)==len(self.exon_ends))
        assert(len(self.exon_starts)== self.exon_count)

    def get_introns(self):
        '''
        generator for introns
        '''
        for i, (next_s, end) in enumerate(zip(self.exon_starts[1:], self.exon_ends)):
            intron_count = self.exon_count - i if self.strand == '-' else i + 1
            intron_line = '{chrom}\t{start}\t{end}\t{tid}\t{intron_count}\t{strand}'\
                    .format(chrom = self.chrom,
                            start = end,
                            end = next_s,
                            tid = self.tid,
                            intron_count = intron_count,
                            strand = self.strand)
            yield intron_line

    def genomic_position(self, tpos):
        '''
        translate transcriptome position to genomic position

        |-------------|-----------*---------|-----------|
        /    exon1   /\        exon2      / \    exon3  \
       /            /  \                 /   \           \
      |------------|xxxx|---------*-----|xxxxx|----------|  
                    intron              intron
        '''

        self.idx = self.cumulative_exon_starts < tpos # which exon the position is on
        shifted_pos = self.cumulative_exon_starts[self.idx][-1]
        offset_from_exon_start = tpos - shifted_pos  # how far is the position to the exon start
        if self.strand == "+":
            out_pos = offset_from_exon_start + self.exon_starts[self.idx][-1]
        else:
            out_pos = self.reversed_exon_starts[self.idx][-1] - offset_from_exon_start
        return out_pos
