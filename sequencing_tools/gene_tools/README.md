# sequencing_tools.gene_tools #

## Transcriptome ##

input:
- refflat file: (default: [hg19 refflat](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz) )
- coding_only: True if only protein coding genes are needed

Example usage:
```
from sequencing_tools.gene_tools import Transcriptome
transcriptome = Transcriptome(coding_only=True)
transcript_dict = transcriptome.MakeTranscriptome()
```

The output transcript_dict is:
- key: gene name
- value: dict
  - with key: transcript ID
  - value: Transcript object
    - chrom, start, end, strand,
    - exons:
      - start, end, exon_num, strand, exon length, contain_cds, contain_cde, after_cds, after_cde, coding_bases, cumulative_transcript_length


## Bed12 ##

input:
- bed12 line


Example usage:
```
from sequencing_tools.gene_tools import Bed12Record
with open('bedfile') as bed:
    for line in bed:
        bedline = Bed12Record(line)
        
        # generate intron coordinates
        for intron in bedline.get_introns()
            print(intron) # 6 columns bed for outputing intron coordinates
        
        '''
        translate transcriptome position to genomic position

        |-------------|-----------*---------|-----------|
        /    exon1   /\        exon2      / \    exon3  \
       /            /  \                 /   \           \
      |------------|xxxx|---------*-----|xxxxx|----------|  
                    intron              intron
        '''

        # lets say we want to find the genomic position for the 10th base on this transcript

        genomic_position = bedline.genomic_position(10)
```


## GTF ##
Parser for gtf file

input:
   - GTFline

output:
- object with attr:
    - chrom: str
    - start: int
    - end: int
    - strand: str
    - feature_type: str
    - info: dict
    - fields: list separated by '\t'

Example usage:
```
from sequencing_tools.gene_tools import GTFRecord
with open('GTF_file') as gtf:
    for line in gtf:
        if not line.startswith('#'):
            gtf_line = GTFRecord(line)

```
