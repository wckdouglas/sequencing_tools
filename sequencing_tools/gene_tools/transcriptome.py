#!/usr/bin/env python

'''
contain transcript structure:

1. transcript
2. exon
'''

import sys
from collections import defaultdict
from pandas import read_csv
import logging
log = logging.getLogger(__name__)


class Exon():
    '''
    Exon info
    '''
    def __init__(self, start, end, exon_num, strand):
        self.start = int(start)
        self.end = int(end)
        self.exon_num = int(exon_num)
        self.strand = strand
        self.length = self.end - self.start
        self.contain_cds = 0
        self.contain_cde = 0
        self.after_cds = 0
        self.after_cde = 0
        self.coding_bases = 0
        self.cumulative_transcript_length = 0  #exon end coordincate - transcript start


    def __ContainWholeCD__(self):
        '''
        The exon contain the whole coding sequence
        '''
        return self.contain_cds and self.contain_cde


    def __IsCodingExon__(self):
        '''
        does the exon contains any coding bases
        '''
        return self.contain_cde or self.contain_cds or (not self.after_cde and self.after_cds)
    

class Transcript():
    '''
    transcript  info:
    '''
    def __init__(self, transcript):
        self.chrom = transcript['chrom']
        self.id = transcript['tid']
        self.strand = '+' if transcript['strand'] in ['+', 1] else '-'
        self.tx_start = transcript['tx_start']
        self.tx_end = transcript['tx_end']
        self.cds = transcript['cds']
        self.cde = transcript['cde']
        self.exon_count = transcript['exon_count']
        self.exon_starts = transcript['exon_starts'].split(',')[:-1]
        self.exon_ends = transcript['exon_ends'].split(',')[:-1]
        self.five_UTR_length = self.tx_start - self.cds
        if self.strand == '-': # reverse exon starts and ends for reverse strand
            self.exon_starts = self.exon_starts[::-1]
            self.exon_ends = self.exon_ends[::-1]
            self.five_UTR_length = self.tx_end - self.cde
            self.cds, self.cde = self.cde, self.cds

        self.cds_off_set =  -1
        self.exons = {} #(populate by (exon start, exon end, exon length))
        self.__MakeTranscripts__()




    def __MakeTranscripts__(self):
        """
        populate the exon list
        """
        cumulative_transcript_length = 0
        after_cds = 0 # walked passed cds??
        after_cde = 0
        for exon_num in range(self.exon_count):
            long_cds = 0
            exon = Exon(self.exon_starts[exon_num], 
                        self.exon_ends[exon_num],
                        exon_num + 1,
                        self.strand)
            
            if exon.start <= self.cds <= exon.end:
                after_cds += 1
                exon.contain_cds = 1
                distance_from_exon_start = self.cds - exon.start if exon.strand=="+" else exon.end - self.cds
                self.cds_off_set = distance_from_exon_start + cumulative_transcript_length

            if exon.start <= self.cde <= exon.end:
                after_cde += 1
                exon.contain_cde = 1
                distance_from_exon_end = self.cde - exon.start if exon.strand=="-" else exon.end - self.cde


            exon.after_cds = after_cds
            exon.after_cde = after_cde

            if exon.after_cds:
                if not exon.contain_cds and not exon.contain_cde:
                    '''
                    Exon:   5'---------------------3'
                    CDS: |->                         ->|
                    '''
                    exon.coding_bases = exon.length
            
                elif exon.contain_cds and not exon.contain_cde:
                    '''
                    Exon: 5'---------------------3'
                    CDS:        |->
                    '''
                    exon.coding_bases = exon.length - distance_from_exon_start

                elif not exon.contain_cds and exon.contain_cde:
                    '''
                    Exon: 5'---------------------3'
                    CDS:              ->|
                    '''
                    exon.coding_bases = exon.length - distance_from_exon_end

                elif exon.contain_cds and exon.contain_cde:
                    '''
                    Exon: 5'---------------------3'
                    CDS:      |->           ->|
                    '''
                    exon.coding_bases = exon.length - distance_from_exon_end - distance_from_exon_start

            assert(exon.length >= 0)
            cumulative_transcript_length += exon.length
            exon.cumulative_transcript_length = cumulative_transcript_length

            self.exons[exon_num + 1] = exon
    

    def FilterExon(self):
        '''
        get the largest exon after CDS within 500bp
        '''

        for exon in self.exons.values():
            if self.__CloseExon__(exon):
                # we don't allow CDS-CDE-containing the exon
                yield exon


    def __CloseExon__(self, exon):
        '''
        The whole exon is in 500bp downstream of CDS
        '''
#        return self.cds_off_set < exon.cumulative_transcript_length - exon.length < (500 + self.cds_off_set)
        return 50 <= exon.cumulative_transcript_length  <= 500 or \
             exon.cumulative_transcript_length - exon.length <= 500 <= exon.cumulative_transcript_length


    def __FirstCodingExon__(self):
        return list(filter(lambda ex: ex.contain_cds==1, self.exons.values()))[0]


class Transcriptome():
    """
    all annotated transcripts
    read in a refflat file
    """

    def __init__(self, refflat = None, coding_only=True):
        self.url = refflat or 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz'
        self.transcript_dict = self.MakeTranscriptome(coding_only=coding_only)

    
    def MakeTranscriptome(self, coding_only=False):
        '''
        build transcriptome database, index by gene name
        '''
        refflat = read_csv(self.url, sep='\t',
                            usecols = [0,1,2,3,4,5,6,7,8, 9,10],
                            names = ['Gene name', 'tid', 'chrom', 'strand', 
                                        'tx_start', 'tx_end',
                                        'cds','cde', 'exon_count', 
                                        'exon_starts','exon_ends'],
                            chunksize = 1000) 
 
        transcript_dict = defaultdict(lambda: defaultdict(Transcript))
        for tab in refflat:
            for i, transcript in tab.iterrows():
                if transcript['cds'] != transcript['cde'] or not coding_only:
                    transcript_dict[transcript['Gene name']][transcript['tid']] = Transcript(transcript)
        return transcript_dict



        

                
        
