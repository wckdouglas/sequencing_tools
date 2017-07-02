# TGIRT-seq tools #

This repo stores different utils for TGIRT-seq/NGS pipelines. This consist of several [ready-to-use scripts](#scripts) and two [module](#modules) for operations on fastq files and BAM alignments.

---

## Installation ##

```
pip install git+https://github.com/wckdouglas/tgirt_seq_tools.git
```

If this gives an error, try:

```
git clone https://github.com/wckdouglas/tgirt_seq_tools.git
cd tgirt_seq_tools
python setup.py install --user
```

Prerequisits:

* pysam>0.11.0
* numpy
* cython
* pyBigWig

```
pip install pysam==0.11.0
pip install numpy
pip install cython
pip install pyBigWig
```
---

<h2 id='scripts'> Ready-to-use scripts <h2>

<h3 id='clip'>  Clip UMI from fastq </h3>

```
usage: clip_fastq.py [-h] [-o OUTPUTPREFIX] -1 FASTQ1 -2 FASTQ2 [-x IDXBASE]
                     [-q BARCODECUTOFF] [-a MISMATCH] [-s {0,1,2,3,4}]
                     [-r {read1,read2}]

Clip the barcode sequence and attached to the front of seq id

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTPREFIX, --outputprefix OUTPUTPREFIX
                        Interleaved Fastq files (default: -)
  -1 FASTQ1, --fastq1 FASTQ1
                        Paired end Fastq file 1 with four line/record
  -2 FASTQ2, --fastq2 FASTQ2
                        Paired end Fastq file 2 with four line/record
  -x IDXBASE, --idxBase IDXBASE
                        how many base in 5' end as index? (default:
                        XXXXXXXXXXXXX) X as umi bases can add constant regions
                        add back, such as XXXXCATGC
  -q BARCODECUTOFF, --barcodeCutOff BARCODECUTOFF
                        Average base calling quality for barcode sequence
                        (default=20)
  -a MISMATCH, --mismatch MISMATCH
                        Allow how many mismatch in constant region (deflaut:
                        1)
  -s {0,1,2,3,4}, --prefix_split {0,1,2,3,4}
                        Using how many bases on the barcode to split the
                        fastq? A choice of 3 will generate 4^3 = 64 files
                        (deflaut: 4)
  -r {read1,read2}, --read {read1,read2}
                        Which read is the UMI on?

```

### Split unique bam ###

This script split bowtie2/hisat2 bam file to uniquely mapped or multiple mapped alignments


```
usage: split_uniq_bam.py [-h] [-a {bowtie2,hisat2}] [-i INBAM] -o OUTPREFIX

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


### Reduce multi reads ###

This scripts process multiply mapped reads and output ribosomal reads or the shortest read pairs from the alignments
input bam file should be generated from split_uniq_bam.py

```
usage: reduce_multi_reads.py [-h] -i INFILE -o OUTFILE [-b] [-z]

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


### BAM to bed ###
This script process name-sorted paired-end bam and output bed file storing fragment positions.

```
usage: bam_to_bed.py [-h] -i IN_BAM [-o OUT_BED] [-m MIN_SIZE] [-M MAX_SIZE]

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
```

### bam_umi_tag.py ###

Adding UMI sequence to RX tag in bam. Alignments must be from fastq files generated by [clip_fastq.py](#clip)

```
Putting UMI as RX tag in bam, enables picard Markduplicates with
BARCODE_TAG=RX

optional arguments:
  -h, --help            show this help message and exit
  -i IN_BAM, --in_bam IN_BAM
                        BAM file name, or stdin (-) ** name sorted
  -o OUT_BAM, --out_bam OUT_BAM
                        BAM file output (default: - )
```

---

<h2 id='modules'> Modules </h2>

### tgirt_Seq_tools.fastq_tools ###

<h4 id='fastq_record'> *class* tgirt_seq_tools.fastq_tools.fastqRecord(id, seq, qual) </h4>

Parameters:
- id - sequence name
- sequence -  actual sequence
- quality - base qualities of the sequence

Example:
```
from tgirt_seq_tools.fastq_tools import fastqRecord
seq_name = 'seq1'
sequence = 'AACCTTGG'
seq_qual = '!!!!!!!!'
record = fastqRecord(seq_name, sequence, seq_qual)
print record.id, record.seq, record.qual
```


#### *function* tgirt_seq_tools.fastq_tools.gzopen ####

**python** gzip library for reading is [slow](http://aripollak.com/pythongzipbenchmarks/).This function provide a faster way open gzip file, using **GNU** ```zcat``` for backend.

This can treat as normal open(filename, 'r') and return a file handle

```
fastq_tools.gzopen(filename, 'r')
```

return: file handle

#### *iterator* fastq_tools.readfq(file) ####

This is a fast [fastq iterator](https://github.com/lh3/readfq/blob/master/readfq.py) that returns a [fastqRecord](#fastq_record).

Parameter:

* fp -file handle of a fastq file

Return:

* fastqRecord object
    * name - sequence id
    * seq - sequence
    * qual - quality


### bam_tools ###
