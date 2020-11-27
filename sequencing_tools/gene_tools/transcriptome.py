#!/usr/bin/env python

"""
Modeul for manipulating transcriptome annotation data

"""

import sys
from collections import defaultdict
from pandas import read_csv, read_sql_query
import sqlite3
import logging
from ..utils import SeqUtilsError
logging.basicConfig(level = logging.INFO)
log = logging.getLogger('Transcriptome')


class Exon():
    '''
    Exon object storing the coordinate, exon rank and strand information

    Args:
        start (int): genomic start coordinate of the exon
        end (int): genomic end coordinate of the exon
        exon_num (int): exon rank
        strand (string): strandeness of the exon (either + or -)
    
    '''
    def __init__(self, start, end, exon_num, strand):
        self.start = int(start)  #: Exon start position on the chromosome
        self.end = int(end) #: Exon end position on the chromosome
        self.exon_num = int(exon_num) #: Exon rank
        self.strand = strand #: Exon strand
        self.length = self.end - self.start #: Exon size
        self.contain_cds = 0 #: 1 if this exon contains coding start site else 0
        self.contain_cde = 0 #: 1 if this exon contains coding end site else 0
        self.after_cds = 0 #: 1 if downstream of coding start site else 0
        self.after_cde = 0 #: 1 if downstream of coding end  site else 0
        self.coding_bases = 0 #: how many bases are within coding frame
        self.cumulative_transcript_length = None   #: the transcript position of last exon base
        self.transcript_start = None #: the transcript position of the first exon base
        self.transcript_end = None #: the transcript position of the last exon base

    def add_cumulative_length(self, cumulative_transcript_length):
        self.cumulative_transcript_length = cumulative_transcript_length
        self.transcript_start = self.cumulative_transcript_length - self.length
        self.transcript_end = self.cumulative_transcript_length

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
    
    def __contains__(self, transcript_position):
        '''
        dose this exon contain a certain transcript position?
        '''
        return self.transcript_start <= transcript_position <= self.cumulative_transcript_length
    
    def __len__(self):
        return self.length
    

class Transcript():
    def __init__(self, transcript):
        """
        Transcript object storing information of a transcript
            
            
        Args:
            transcript (dict): with the following keys:
                        1. chrom
                        2. tid
                        3. strand
                        4. tx_start
                        5. tx_end
                        6. cds
                        7. cde
                        8. exon_count
                        9. exon_starts: comma delimited
                        10. exon_ends: comma delimited

        
        """
        self.chrom = transcript['chrom'] #: chromosome name
        self.id = transcript['tid'] #: transcript ID
        self.strand = '+' if transcript['strand'] in ['+', 1] else '-' #: strand ('+' or '-')
        self.tx_start = int(transcript['tx_start']) #: transcription start site (genomic)
        self.tx_end = int(transcript['tx_end'])  #: transcription end site (genomnic)
        self.cds = int(transcript['cds']) #: coding start site (genomic)
        self.cde = int(transcript['cde']) #: coding end site (genomic)
        self.exon_count = int(transcript['exon_count']) #: number of exons in this transcript 
        self.exon_starts = transcript['exon_starts'].split(',')[:-1] #: str, comma list of exon start sites 
        self.exon_ends = transcript['exon_ends'].split(',')[:-1] #: str, comman list of exon end sites
        self.five_UTR_length = self.tx_start - self.cds #: how long is the 5' UTR? 5' as of genomic coordinate
        self.transcript_length = 0 #: transcript size
        if self.strand == '-': # reverse exon starts and ends for reverse strand
            self.exon_starts = self.exon_starts[::-1]
            self.exon_ends = self.exon_ends[::-1]
            self.five_UTR_length = self.tx_end - self.cde
            self.cds, self.cde = self.cde, self.cds

        self.cds_off_set =  -1
        self.exons = {} #: dictionary with key as exon rank, values are :class:`sequencing_tools.gene_tools.transcriptome.Exon`
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
            exon.add_cumulative_length(cumulative_transcript_length)

            self.exons[exon_num + 1] = exon
        self.transcript_length = cumulative_transcript_length
    
    def blocks(self, tstart_pos: int, tend_pos: int):
        '''
        given a start position and end position along the transcript,
        return the block starts and block sizes on the genome scale

        Example::

                                        exon 1                           exon 2
                Transcript: |===========|-----------------|=================|
                Amplicon:        |------>                 <--------|
                blocks:             block1                   block2
                Return:        [(a,     b)                (c,      d)   ]

        Args:
            tstart_pos (int): left position on the transcript
            tend_pos (int): right position on the transcript

        Returns:
            list: list of tuples containing the (start, end) of each block
        '''
        assert (tend_pos > tstart_pos)
        assert(tend_pos <= self.transcript_length)
        assert(tstart_pos >= 0)

        start_collecting = 0
        collected_all_exon = 0
        blocks = []

        for exon_number, exon in self.exons.items():
            '''
            Examples: 

                            exon 1                           exon 2
            Transcript: |===========|-----------------|=================|
            Amplicon:        |------>                 <--------|
            blocks:             block1                   block2

            
                            exon 1                           
            Transcript: |===========================|
            Amplicon:        |--------------|
            blocks:             block1                   


                            exon 1                exon 2                 exon3
            Transcript: |===========|---------|==========|--------|=================|
            Amplicon:        |------>         |----------|        <--------|
            blocks:            block1            block2            block3


            '''
            if tstart_pos in exon:
                start_collecting = 1
                if self.strand == "+":
                    block_start = exon.start + (tstart_pos - exon.transcript_start)
                else:
                    block_end = exon.end - (tstart_pos - exon.transcript_start)

                if  tend_pos in exon:
                    # example 2
                    if self.strand == '+':
                        block_end = exon.start + (tend_pos - exon.transcript_start)
                    else:
                        block_start = exon.end - (tend_pos - exon.transcript_start)
                    collected_all_exon = 1
                
                else:
                    # example 1 or 3
                    if self.strand == '+':
                        block_end = exon.end
                    else:
                        block_start = exon.start
                blocks.append((block_start, block_end))
            
            elif collected_all_exon == 0 and start_collecting == 1:
                if self.strand == '+':
                    block_start = exon.start
                else:
                    block_end = exon.end

                if tend_pos in exon:
                    # exon 2 from example 1, or exon 3 from example 3
                    if self.strand == '+':
                        block_end = exon.start + (tend_pos - exon.transcript_start)
                    else:
                        block_start = exon.end - (tend_pos - exon.transcript_start)
                    collected_all_exon = 1
                else:
                    # exon 2 from example 3
                    if self.strand == "+":
                        block_end = exon.end
                    else:
                        block_start = exon.start
                
                blocks.append((block_start, block_end))
        amplicon_size = tend_pos - tstart_pos
        assert( sum(map(lambda x: x[1] - x[0], blocks)) == amplicon_size) 
        return blocks

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
        assert( 0 <= tpos <= self.transcript_length )
        for exon in self.exons.values():
            if tpos in exon:
                offset = tpos - exon.transcript_start
                return exon.start - offset

    
    def __FirstCodingExon__(self):
        return list(filter(lambda ex: ex.contain_cds==1, self.exons.values()))[0]


class Transcriptome():
    """
    all annotated transcripts
    read in a refflat file

    Args: 
        refflat:  refflat file or `url <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz>`_
        sqldb: Can be found under `Annotationhub <https://annotationhub.bioconductor.org/package2/AHEnsDbs>`_ (if no refflat is supplied)
        coding_only: only index coding transcript?
    

    Usage::

        url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz'
        txome = Transcriptome(refflat=url, coding_only=False)
        gene = txome.get_gene('WASH7P')
        tx = gene['NR_024540']
        print(tx.exons[1].start, tx.exons[2].end) # 29320 24891
        print(tx.genomic_position(300)) #18171
        print(tx.blocks(0, 300))
        #[(29320, 29370), (24737, 24891), (18270, 18366)]

    Or::

        # wget http://s3.amazonaws.com/annotationhub/AHEnsDbs/v87/EnsDb.Hsapiens.v87.sqlite
        txome = Transcriptome(sqldb="EnsDb.Hsapiens.v87.sqlite", coding_only=False)

    """

    def __init__(self, refflat = None, sqldb=None, coding_only=True):
        self.coding_only = coding_only #: Only index coding genes (mRNA)?
        self.transcript_dict = defaultdict(lambda: defaultdict(Transcript)) #: Dict storing transcripts, with gene name as key and values are :class:`sequencing_tools.gene_tools.transcriptome.Transcript`
        self.transcript_count = 0 #: how many transcript was indexed?
        if sqldb:
            self.sqldb = sqldb
            self.sql_connection = sqlite3.connect(self.sqldb)
            self.__MakeTranscriptomeFromSqlite()
        else:
            self.url = refflat or 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz'
            self.__MakeTranscriptomeFromRefFlat()

    
    def __MakeTranscriptomeFromRefFlat(self):
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
 
        for tab in refflat:
            for i, transcript in tab.iterrows():
                self.transcript_count += 1
                if transcript['cds'] != transcript['cde'] or not self.coding_only:
                    self.transcript_dict[transcript['Gene name']][transcript['tid']] = Transcript(transcript)

    def __MakeTranscriptomeFromSqlite(self):
        '''
        indexing from sqlite db
        '''
        cursor = self.sql_connection.cursor()
        log.info('Reading database %s' %self.sqldb)

        #print info
        query = 'select * from metadata'
        for result in cursor.execute(query):
            if result[0] in {'Db type','Organism','genome_build', 'ensembl_version'}:
                log.info('\t{} = {}'.format(result[0], result[1]))
        ###############################################################################################
        
        query = '''
            select a.tx_id, a.tx_biotype, a.tx_seq_start, a.tx_seq_end, 
                    a.tx_cds_seq_start, a.tx_cds_seq_end,
                    c.exon_seq_start, c.exon_seq_end, b.exon_idx,
                    d.seq_strand, d.seq_name, d.gene_name
            from tx as a 
            inner join tx2exon as b ON a.tx_id = b.tx_id 
            inner join exon as c on b.exon_id = c.exon_id 
            inner join gene as d on a.gene_id = d.gene_id
        '''
        transcriptDB = read_sql_query(query, con = self.sql_connection) \
            .groupby(['tx_id','tx_biotype', 'tx_seq_start', 'tx_seq_end', 
                        'tx_cds_seq_start','tx_cds_seq_end', 'seq_strand',
                        'seq_name','gene_name'], 
                    as_index=False) \
            .agg({'exon_seq_start': lambda xs: ','.join(map(str, xs)) + ',',
                    'exon_seq_end': lambda xs: ','.join(map(str, xs)) + ',',
                    'exon_idx': 'count'})\
            .rename(columns = {'seq_name':'chrom',
                            'tx_id': 'tid',
                            'seq_strand': 'strand',
                            'tx_seq_start':'tx_start',
                            'tx_seq_end': 'tx_end',
                            'tx_cds_seq_start':'cds',
                            'tx_cds_seq_end':'cde',
                            'exon_seq_start':'exon_starts',
                            'exon_seq_end': 'exon_ends',
                            'exon_idx':'exon_count'}) 
        log.info('Loaded %i exons' %(transcriptDB.shape[0]))

        for i, transcript in transcriptDB.iterrows():
            if transcript['tx_biotype'] == 'protein_coding' or not self.coding_only:
                self.transcript_dict[transcript['gene_name']][transcript['tid']] = Transcript(transcript)

    def get_gene(self, gene_name):
        '''
        Get transcripts from gene

        Args:
            gene_name (str): 

        Returns:
            dict: dictionary with key as transcript id as key and values are :class:`sequencing_tools.gene_tools.transcriptome.Transcript`

        '''
        if gene_name not in self.transcript_dict.keys():
            raise SeqUtilsError('%s not in database' %gene_name)
        return self.transcript_dict[gene_name]

        

                
        
