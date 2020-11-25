import os
from operator import itemgetter
from collections import defaultdict
import numpy as np
import pysam


'''
usage: cat bed12 | python bed12_to_intron.py  > intron.bed
'''


class Bed12Record():
    """
    Parser for bed12 line
    """
    def __init__(self, line):
        '''
        Analyzing bed12 line

        Args:
            line (str): bed12 line 
        
        Example::

            with open('gene.bed12','r') as bed12:
                for line in bed12:
                    transcript = Bed12Record(line)
        '''
        fields = line.strip().split('\t')
        chrom, start, end, tid, exon_count, \
            strand, exon_starts, exon_sizes = itemgetter(0,1,2,3,9,5,11,10)(fields)
        self.chrom = chrom #: chromosome name
        self.start = int(start) #: genomic start site as in bed12 
        self.end = int(end) #: genomic end site as in bed12 
        self.strand = strand #: srand
        self.exon_count = int(exon_count) #: number of exons for the transcript
        self.tid = tid #: transcript ID
        self.exon_starts = np.array(list(map(int, exon_starts.strip(',').split(',')))) # list of exon starts
        self.exon_sizes = np.array(list(map(int, exon_sizes.strip(',').split(',')))) # list of exon sizes
        self.exon_ends = np.array([s+l for s,l in zip(self.exon_starts, self.exon_sizes)]) # list of exon ends
        if self.strand == '-':
            self.exon_sizes = self.exon_sizes[::-1]
            self.reversed_exon_starts = self.exon_ends[::-1]
        self.transcript_length = self.end - self.start 
        assert(len(self.exon_starts)==len(self.exon_ends))
        assert(len(self.exon_starts)== self.exon_count)

    def get_introns(self):
        '''
        generator for introns, for each transcript (bed12 line), return all its intron

        Returns:
            generator: intron_line, 6-columns bed for the introns
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

        Args:
            tpos: transcript position

        Returns:
            int: the corresponding genomic position

        Illustration::

                |-------------|-----------*---------|-----------|
                /    exon1   /\        exon2      / \    exon3  \
            /            /  \                 /   \           \
            |------------|xxxx|---------*-----|xxxxx|----------|  
                            intron              intron
        '''

        self.idx = self.exon_starts < tpos # which exon the position is on
        shifted_pos = self.exon_starts[self.idx][-1]
        offset_from_exon_start = tpos - shifted_pos  # how far is the position to the exon start
        if self.strand == "+":
            out_pos = offset_from_exon_start + self.exon_starts[self.idx][-1]
        else:
            out_pos = self.reversed_exon_starts[self.idx][-1] - offset_from_exon_start
        return out_pos



class GTFRecord():
    '''
    parsing a GTF record line
    '''
    def __init__(self, gtf_line):
        self.fields = gtf_line.split('\t')
        self.chrom = self.fields[0] #: chromosome name
        self.feature_type = self.fields[2] #: feature type of the transcript
        self.start = self.fields[3] #: genomic start
        self.end = self.fields[4] #: genomic end
        self.strand = self.fields[6] #: strand
        self.start, self.end = int(self.start), int(self.end)
        self.info = self.__parse_extra_fields__(self.fields[-1]) #: info field as a dictionary


    def __parse_extra_fields__(self, extra_field):
        info_fields = extra_field.split(';')
        info_dict = {}
        for info in info_fields[:-1]:
            row_fields = info.strip().split(' ')
            info_dict[row_fields[0]] = row_fields[1].strip('\'"')

        return info_dict


class Transcript():
    '''
    Transcript plotting
    '''
    def __init__(self, transcript_dict):
        '''
        Args:
        dictionary: transcript_dict[exons|trascript] = list
                    something like transcripts = defaultdict(list) [start and end]
        '''
        self.transcript = transcript_dict['transcript'][0]
        self.exons = transcript_dict['exon']
        self.strand = self.transcript.strand

    def plot(self, ax, plot_transcript=False, y = 0, 
             xs = None, xe = None, fontsize=15,
             gene_name = None):
        '''
        plot gene model
        '''
        if plot_transcript:
            ax.hlines(y = y + 1, 
                      xmin = ts[tid]['transcript'][0].start, 
                      xmax=ts[tid]['transcript'][0].end,
                     linewidth=10, color='darkblue')

        #plot intron/transcript
        start = max(xs, self.transcript.start) if xs else self.transcript.start
        end = min(xe, self.transcript.end) if xe else self.transcript.end

        ss = np.linspace(start, end, 20)
        arrow = "<-" if self.strand == '+' else '->'
        for s,e in zip(ss, ss[1:]):
            ax.annotate("", 
                        xy=(s,y), 
                        xytext = (e+5,y), 
                        arrowprops=dict(arrowstyle=arrow))


        #plot exons
        starts = list(map(lambda x: x.start, self.exons))
        ends = list(map(lambda x: x.end, self.exons))
        for start, end in zip(starts, ends):
            start = max(xs, start) if xs else start
            end = min(xe, end) if xe else end
            ax.hlines(y = y, 
                  xmin = start, 
                  xmax = end,
                  linewidth=15, 
                  color='darkblue')
            
            ss = np.linspace(start, end, 5)
            for s in ss:
                arrow = ">" if self.strand == '+' else '<'
                ax.text(s, y,arrow, fontsize=fontsize,va='center',
                        color='white')

        if gene_name:
            ax.text(end + (xe-xs)*0.02, y, self.transcript.info[gene_name], 
                    fontsize=fontsize, va='center', ha='left', color='darkblue')




class GeneModel():
    '''
    Exon and transcripts GTF
    '''
    def __init__(self, gtf_file):
        assert os.path.isfile(gtf_file + '.tbi'), 'Require indexed GTF'
        self.gtf = pysam.Tabixfile(gtf_file)
        self.transcripts = None

    def fetch_transcripts(self, chrom, start, end):
        self.transcripts = defaultdict(lambda: defaultdict(list))
        for line in self.gtf.fetch(chrom, start, end):
             g = GTFRecord(line)
             if g.feature_type != "gene":
                 self.transcripts[g.info['transcript_id']][g.feature_type].append(g)
        return self.transcripts



