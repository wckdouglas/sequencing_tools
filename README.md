
# Sequencing tools #

[![CI](https://github.com/wckdouglas/sequencing_tools/workflows/CI/badge.svg)](https://github.com/wckdouglas/sequencing_tools/actions)[![Documentation Status](https://readthedocs.org/projects/sequencing-tools/badge/?version=latest)](https://sequencing-tools.readthedocs.io/en/latest/?badge=latest)[![Docker Cloud Build Status](https://hub.docker.com/repository/docker/wckdouglas/sequencing_tools/general)](https://img.shields.io/docker/cloud/build/wckdouglas/sequencing_tools)


This repo stores different utils for TGIRT-seq/NGS pipelines. This consist of several [ready-to-use scripts](#scripts) and two [module](#modules) for operations on fastq files and BAM alignments.

---

## Installation ##


```
git clone https://github.com/wckdouglas/sequencing_tools.git
cd sequencing_tools
pip install -r requirements.txt
python setup.py install --user
```

Prerequisits:

* pysam>0.11.0
* numpy>=1.12.1
* cython>=0.25
* matplotlib>=2.0.0
* scipy>=0.19.0
* seaborn>=0.7.1
* networkx
* ujson
* six
* pytest

---

<h2 id='scripts'> Ready-to-use scripts </h2>

* [umi clipper](#clip)
* [Paired end fastq merger](#merge_pe)
* [Unique BAM splitter](#split_bam)
* [Multimap BAM reducer](#fix_multi)
* [BAM to BED converter](#b2b)
* [UMI tagger](#bam_tag)
* [BED demultiplexer](#dedup_bed)
* [FastQ deinterleave](#deinterleaved)
* [Poisson UMI](#poisson_umi)
* [BAM base counter](#base_count)

<h3 id='clip'>  UMI clipper </h3>

Clipping UMI from paired end fastq

```
usage: seqtools clipFQ [-h] [-o OUT_FILE] -1 FASTQ1 -2 FASTQ2 [-x IDXBASE]
                     [-q BARCODECUTOFF] [-a MISMATCH] [-r {read1,read2}]

Clip the barcode sequence and attached to the front of seq id

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_FILE, --out_file OUT_FILE
                        Interleaved Fastq files (default: -)
  -1 FASTQ1, --fastq1 FASTQ1
                        Paired end Fastq file 1 with four line/record
  -2 FASTQ2, --fastq2 FASTQ2
                        Paired end Fastq file 2 with four line/record
  -x IDXBASE, --idxBase IDXBASE
                        how many base in 5' end as index? (default:
                        XXXXXXXXXXXXX) X as umi bases, can also add constant
                        regions at the back, such as XXXXCATGC, if CATGC is
                        the constant region
  -q BARCODECUTOFF, --barcodeCutOff BARCODECUTOFF
                        Average base calling quality for barcode sequence
                        (default=20)
  -a MISMATCH, --mismatch MISMATCH
                        Allow how many mismatch in constant region (deflaut:
                        1)
  -r {read1,read2}, --read {read1,read2}
                        Which read is the UMI on?
```

<h3 id='merge_pe'> Paired end fastq merger </h3>

Merging paired end fastq to a fragment, only output reads that are successfully merged.

```
usage: seqtools mergepe [-h] [-i INTERLEAVED] [-o OUTFILE] [-m MIN_LEN]
                      [-e ERROR] [-a]

Merging interleaved, paired-end fastq file and output overlapped regions only
with error correction using cutadapt module to find overlapping regions

optional arguments:
  -h, --help            show this help message and exit
  -i INTERLEAVED, --interleaved INTERLEAVED
                        Interleaved Fastq files (default: -)
  -o OUTFILE, --outfile OUTFILE
                        Merged fastq file (default: -)
  -m MIN_LEN, --min_len MIN_LEN
                        Minimum length of sequence to output (default: 18)
  -e ERROR, --error ERROR
                        Maximum error rate of alignment (default: 0.1)
  -a, --all             Output all bases (default: only overlapping regions)
```

<h3 id='split_bam'> Unique bam splitter</h3>

This script split bowtie2/hisat2 bam file to uniquely mapped or multiple mapped alignments 


```
usage: seqtools splitBam [-h] [-a {bowtie2,hisat2}] [-i INBAM] -o OUTPREFIX

Splitting bam file to uniquely mapped and multiple mapped (only support
bowtie2 and hisat2)

optional arguments:
  -h, --help            show this help message and exit
  -a {bowtie2,hisat2}, --aligner {bowtie2,hisat2}
                        PRogram generating the bam file (default: bowtie2)
  -i INBAM, --inBam INBAM
                        bam file name, or stdin (default: stdin)
  -o OUTPREFIX, --outprefix OUTPREFIX
                        output prefix ($OUTPUTPREFIX.unique.bam,
                        $OUTPUTPREFIX.multi.bam)
```


<h3 id='fix_multi'> Multimap BAM reducer </h3>

This scripts process multiply mapped reads and output ribosomal reads or the shortest read pairs from the alignments
input bam file should be generated from split_uniq_bam.py

```
usage: seqtools filterMulti [-h] -i INFILE -o OUTFILE [-b] [-z]

Process multiply mapped bam and get either shortest insert size or ribosomal
reads

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input bam/sam file (use < - > for stdin)
  -o OUTFILE, --outfile OUTFILE
                        Output bam/sam file (use < - > for stdout)
  -b, --bam_in          Provide this flag if bam instead of sam is used as
                        input
  -z, --bam_out         Provide this flag if bam is needed for output
```


<h3 id='b2b'> BAM to bed converter</h3>

```
usage: seqtools bam2bed [-h] -i IN_BAM [-o OUT_BED] [-m MIN_SIZE] [-M MAX_SIZE]
                     [-t TAG] [-a]

Making paired-end bam into bed file for every fragment

optional arguments:
  -h, --help            show this help message and exit
  -i IN_BAM, --in_bam IN_BAM
                        BAM file name, or stdin (-) ** name sorted
  -o OUT_BED, --out_bed OUT_BED
                        BED file output (default: - )
  -m MIN_SIZE, --min_size MIN_SIZE
                        minimum fragment size to report
  -M MAX_SIZE, --max_size MAX_SIZE
                        minimum fragment size to report
  -t TAG, --tag TAG     tag to extract
  -a, --all             output supplementary as single fragment
```

<h3 id='bam_tag'> UMI tagger </h3>

Adding UMI sequence to RX tag in bam. Alignments must be from fastq files generated by [clip_fastq.py](#clip)

```
usage: seqtools bamTag [-h] -i IN_BAM [-o OUT_BAM] [-t TAG] [-d DELIM]
                      [-f FRAGMENT]

Putting UMI as RX tag in bam, enables picard Markduplicates with
BARCODE_TAG=RX

optional arguments:
  -h, --help            show this help message and exit
  -i IN_BAM, --in_bam IN_BAM
                        BAM file name, or stdin (-)
  -o OUT_BAM, --out_bam OUT_BAM
                        BAM file output (default: - )
  -t TAG, --tag TAG     Tag id (default: RX )
  -d DELIM, --delim DELIM
                        Deliminator separating read id and bc (default: _ )
  -f FRAGMENT, --fragment FRAGMENT
                        after splitting read name using {delim}, which fragment
                        is UMI? can use -1 as last piece (default: 0)
```


<h3 id='dedup_bed'> Deuplicate BED </h3>

Deduplicate fragments with UMI in a bed file.

```
usage: seqtools dedup [-h] [-i INFILE] [-o OUTFILE] [-t THRESHOLD]
                          [-d DELIM] [-f F]

Demultiplexing UMI bed file with the follwing columns: 1. chrom name 2. start
3. end 4. {$UMI}_{$READ_ID} 5. score 6.strand The program internally used
hamming distance matrix of the barcodes to generate connected network and
identified UMI clusters

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        input bedfile, can be "-" or "/dev/stdin" for stdin
                        (default: -).Format should be sorted BED: 1. chrom 2.
                        start 3. end 4. {barcode}_{id} 6. strand
  -o OUTFILE, --outfile OUTFILE
                        output bedfile, can be "-" or "/dev/stdout" for stdin
                        (default: -).
  -t THRESHOLD, --threshold THRESHOLD
                        How many error between barcodes can be tolerated?
                        (default = 1)
  -d DELIM, --delim DELIM
                        Deliminator separating read id and bc (default: _ )
  -f F                  after splitting read name using {delim}, which
                        fragmnet is UMI? can use -1 as last piece (default: 0)
```

<h3 id='deinterleaved'> Deinterleaved fastq </h3>

Deinterleaved an interleaved fastq file.

```
usage: seqtools deinterleave [-h] [-i INFQ] -1 READ1 [-2 READ2]

Deinterleaving an interleaved fastq file

optional arguments:
  -h, --help            show this help message and exit
  -i INFQ, --infq INFQ  interleaved fq (defualt: -)
  -1 READ1, --read1 READ1
                        read 1 fastq output
  -2 READ2, --read2 READ2
                        read 2 fastq output
```


<h3 id='poisson_umi'> Adjust fragment counts for UMI saturations  </h3>

UMI saturation can be a problem for highly-expressed geens, this tool adjust the fragment count by implementing a [poisson model](https://www.ncbi.nlm.nih.gov/pubmed/21562209). 

![](https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/umi_correct.png)


where *m* is the diversity of UMI (*4^nt*, *nt* is number of bases as UMI), *k* is the number of different UMI being detected for the fragment and *n* is the true number of fragments that we are interested at.

Input for this tool is a 6-columns BED file with each fragment named as {UMI}_{READ_ID}  


```
usage: seqtools calcUMI [-h] -i IN_BED [-o OUT_BED] [--umi UMI]

Adjusting fragment count for UMI saturation, see paper: Counting individual
DNA molecules by the stochastic attachment of diverse labels.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_BED, --in_bed IN_BED
                        BED file name, or stdin (-) ** name sorted
  -o OUT_BED, --out_bed OUT_BED
                        BED file output (default: - )
  --umi UMI             Number of nucleotide as umi (default: 6)
```

<h3 id='base_count'> Extracting base count from BAM along the genome/genomic regions  </h3>

```
usage: seqtools pileup [-h] -i BAM -f FASTA [-b BASES] [-q QUAL]
                              [-c CROP] [-r BED] [--no_indel]
                              [--min_coverage MIN_COVERAGE]
                              [--concensus_fasta CONCENSUS_FASTA]

Pileup whole genome, only output bases where coverage > 0

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     Input bam file (indexed)
  -f FASTA, --fasta FASTA
                        reference fasta file
  -b BASES, --bases BASES
                        number of bases to look at every iteration (default:
                        100000)
  -q QUAL, --qual QUAL  base quality to filter (defulat: 30)
  -c CROP, --crop CROP  Crop how many bases from ends (defulat: 0)
  -r BED, --bed BED     bed file for regions (default: whole genome)
  --no_indel            Not considering alignments with Indel
  --min_coverage MIN_COVERAGE
                        Minimum coverage to output
  --concensus_fasta CONCENSUS_FASTA
                        Generate concensus fasta (only work if bed file is
                        provided)
```


---

<h2 id='modules'> Modules </h2>

Detail usage can be found at [readthedocs](https://sequencing-tools.readthedocs.io/en/latest/)


### viz_tools.cor_plot ####

A specialized paired correlation plot

```
from sequencing_tools.viz_tools import cor_plot
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('white')

d = np.random.rand(100,10)
d = pd.DataFrame(d)
fig = plt.figure()
cor_plot(d, fig)
fig.savefig('cor.png',bbox_inches='tight')
```

![](https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/cor.png)
