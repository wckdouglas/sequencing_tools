import os
from operator import itemgetter
from collections import defaultdict
import numpy as np
import pysam
from .transcriptome import Exon, Transcript


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
        exon_starts = fields[11].strip(',').split(',')
        exon_sizes = fields[10].strip(',').split(',')
        exon_starts = int(fields[1]) + np.array(list(map(int, exon_starts))) # list of exon genomic starts
        exon_ends = exon_starts + np.array(list(map(int, exon_sizes))) # list of exon sizes
        keys = ['chrom','tid','strand','tx_start', 'tx_end', 
                'cds','cde', 'exon_count', 'exon_starts', 'exon_ends']
        values = [fields[0], fields[3], fields[5],fields[1], fields[2],
                fields[6], fields[7], fields[9], 
                ','.join(map(str, exon_starts)) + ',', 
                ','.join(map(str, exon_ends)) + ',']
        transcript = {k:v for k,v in zip(keys, values)}
        self.transcript = Transcript(transcript)


    def get_introns(self):
        '''
        generator for introns, for each transcript (bed12 line), return all its intron

        Returns:
            generator: intron_line, 6-columns bed for the introns
        '''
        exons = list(self.transcript.exons.values())
        for intron_count, (exon_3, exon_5) in enumerate(zip(exons[1:], exons)):
            if self.transcript.strand == '-':
                start, end = exon_3.end, exon_5.start
            else:
                start, end = exon_5.end, exon_3.start
            intron_line = '{chrom}\t{start}\t{end}\t{tid}\t{intron_count}\t{strand}'\
                    .format(chrom = self.transcript.chrom,
                            start = start + 1,
                            end = end,
                            tid = self.transcript.id,
                            intron_count = intron_count + 1,
                            strand = self.transcript.strand)
            yield intron_line

    def genomic_position(self, tpos):
        return self.transcript.genomic_position(tpos)



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


class TranscriptPlot():
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



