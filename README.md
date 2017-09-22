[![Build Status](https://travis-ci.org/wckdouglas/sequencing_tools.svg?branch=master)](https://travis-ci.org/wckdouglas/sequencing_tools)

# Sequencing tools #

This repo stores different utils for TGIRT-seq/NGS pipelines. This consist of several [ready-to-use scripts](#scripts) and two [module](#modules) for operations on fastq files and BAM alignments.

---

## Installation ##

```
pip install git+https://github.com/wckdouglas/sequencing_tools.git
```

If this gives an error, try:

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
* python-cjson>=1.2.0
* scipy>=0.19.0
* seaborn>=0.7.1
* pyBigWig

---

<h2 id='scripts'> Ready-to-use scripts </h2>

* [umi clipper](#clip)
* [Paired end fastq merger](#merge_pe)
* [Unique BAM splitter](#split_bam)
* [Multimap BAM reducer](#fix_multi)
* [BAM to BED converter](#b2b)
* [UMI tagger](#bam_tag)

<h3 id='clip'>  UMI clipper </h3>

Clipping UMI from paired end fastq

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

<h3 id='merge_pe'> Paired end fastq merger </h3>

Merging paired end fastq to a fragment, only output reads that are successfully merged.

```
usage: pe_fq_merge.py [-h] [-i INTERLEAVED] [-o OUTFILE] [-m MIN_LEN]
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


<h3 id='fix_multi'> Multimap BAM reducer </h3>

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


<h3 id='b2b'> BAM to bed converter</h3>

This script process name-sorted paired-end bam and output bed file storing fragment positions. Outer ends from read1 and read2 are used as fragment ends. 

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

<h3 id='bam_tag'> UMI tagger </h3>

Adding UMI sequence to RX tag in bam. Alignments must be from fastq files generated by [clip_fastq.py](#clip)

```
usage: bam_umi_tag.py [-h] -i IN_BAM [-o OUT_BAM] [-t TAG]

Putting UMI as RX tag in bam, enables picard Markduplicates with
BARCODE_TAG=RX

optional arguments:
  -h, --help            show this help message and exit
  -i IN_BAM, --in_bam IN_BAM
                        BAM file name, or stdin (-) ** name sorted
  -o OUT_BAM, --out_bam OUT_BAM
                        BAM file output (default: - )
  -t TAG, --tag TAG     Tag id (default: RX )
```

---

<h2 id='modules'> Modules </h2>

* [fastq_tools](#fastq)
* [bam_tools](#bam_tools)

<h3 id='fastq'> tgirt_Seq_tools.fastq_tools </h3>

<h4 id='fastq_record'> <i>class</i> sequencing_tools.fastq_tools.fastqRecord(id, seq, qual) </h4>

Parameters:
- id - sequence name (string)
- sequence -  actual sequence (string)
- quality - base qualities of the sequence (string)

Example:
```
from sequencing_tools.fastq_tools import fastqRecord
seq_name = 'seq1'
sequence = 'AACCTTGG'
seq_qual = '!!!!!!!!'
record = fastqRecord(seq_name, sequence, seq_qual)
print record.id, record.seq, record.qual
```


#### *function* fastq_tools.gzopen ####

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

#### *iterator* fastq_tools.read_interleaved(file) ####

This is a interleaved fastq file parser, return a pairs of read

Parameter:

* fp: file handle of a interleaved fastq file

Return:

* R1: fastqRecord object
    * name: sequence id
    * seq: sequence
    * qual: quality

* R2: fastqRecord object
    * name: sequence id
    * seq: sequence
    * qual: quality

#### fastq_tools.complement ####

Find complement a sequence.

parameter:

* seq - sequence to be complemented (string)

return:

* out - complemented sequence (string)


#### fastq_tools.reverse_complement ####

Reverse complement a sequence.

parameter:

* seq - sequence to be reverse complemented (string)

return:

* out - reverse complemented sequence (string)

<h3 id='bam_tools'> sequencing_tools.bam_tools </h3>

#### bam_tools.cigar_to_str ####

cigar string to string, only extract cigar op == M, I or S

usage: cigar_to_str(cigar_string)

parameter:

* cigar_string: standard cigar string from BAM file

return:

* cigar_seq: same length string as sequence, with every character matched the aligned status of the base

example:
```
$ cigar_seq = cigar_to_str('3S5M1I3M')
$ print cigar_seq
'SSSMMMMMIMMM'
```

#### bam_tools.concordant_alignment ####

Check if alignment is properly mapped flag in [99, 147, 163, 83]

usage: concordant_alignment(alignment)

parameter:

* alignment - pysam alignment segment

return:

* is_concordant - boolean


#### bam_tools.concordant_pairs ####

Check if pair: read1.flag and read2.flag ==  (99, 147) or (83,163)

usage: concordant_pairs(read1, read2)

parameters:

* read1 - read1 of the pair (pysam alignment segment)
* read2 - read2 of the pair (pysam alignment segment)

return:

* Boolean - True if 83 and 163 or 99 and 147 for their flags otherwise False

#### bam_tools.fragment_ends ####

Get start and end position of a pair of reads


usage: fragment_ends(read1, read2)

parameters:

* read1 - read1 of the pair (pysam alignment segment)
* read2 - read2 of the pair (pysam alignment segment)

Return:

* start - leftmost positoin of the pair
* end - rightmost position of the pair

#### bam_tools.get_strand ####

Get strand of the paired fragment

usage: get_strand(alignment)

parameter:

* alignment - an aligned segment in bam (pysam alignment segment)

return:

* strand
 - "+" if it is reverse read2 or forward read1
 - "-" if it is reverse read1 or forward read2

#### bam_tools.make_cigar_seq ####

*Generator* convert number and operator into a sequence of aligned status of bases, see [split_cigar](#split_cigar)

usage: make_cigar_seq(cigar_numbers, cigar_operator)

parameter:

* cigar_numbers - list of numbers in string format
* cigar_operator - list of single character

return:
* *generator* cigar_base - sequence of cigar base (Ignored soft clip)


Example:

```
$ for c in make_cigar_seq(['3','5','1','3'],['S','M','I','M']):
$    print c

MMMMM
I
MMM
```

#### bam_tools.make_regions ####

*Generator* segment chromosome in to regions

usage: make_regions(chromosome_length, how_many_bases_every_time)

Parameter:

* chromosome_length - last base you want to look at
* how_many_bases_every_time - segment size

Return:

* tuple (start,end)
	1.  start - start of segment
    2.  end - end of segment


example:

```
$ for start, end in make_regions(100,10):
$    print start, end

0 10
10 20
20 30
30 40
40 50
50 60
60 70
70 80
80 90
90 100
```

#### bam_tools.read_ends ####

Get read end positions, output start and end position of a read

usage: read_ends(AlignedSegment)

Parameters:

* AlignedSegment -  a pysam alignment

Return:

* start:  leftmost positoin of the read
* end:    rightmost position of the read

#### bam_tools.remove_insert ####

*Generator* remove insertion base from aligned sequence

usage: remove_insert(sequence, quality_string, cigar_seq)

Parameters:

* sequence - DNA sequence from BAM
* quality_string - qual string from BAM
* cigar_seq - cigar seq from [cigar_to_str](#bam_tools)

Yield:

* base - base that are not annotated as insertion
* base_qual - BAQ associated with the base


<h4 id='split_cigar'> bam_tools.split_cigar </h4>

Split cigar string to numpy array

usage: split_cigar(cigar_string)

Parameters:

* cigar_string: cigar string, e.g. 63M

Return:

* *tuple* ([list of numbers],
			[list of cigar operators correspongs to the numbers])

Example:

```
$ split_cigar('63M')

[[63], [M]]
```

<h2 id='extra'> Extra </h2>

### sequencing_tools.douglas_colors ###

#### douglas_colors.douglas_palette ####

Automatic set color if seaborn is installed, otherwise return list of colors

usage: douglas_palette()

Example:

```
import matplotlib.pyplot as plt
import numpy as np
from sequencing_tools.douglas_colors import douglas_palette
douglas_palette()
plt.figure()
ax=plt.subplot();
for i in range(14):
    ax.plot(np.arange(10),np.arange(10) - i,label=i)
ax.legend(bbox_to_anchor=(1,1))
plt.savefig('palette.png',bbox_inches='tight')
```

![](https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/palette.png)

#### douglas_colors.cor_plot ####

A specialized paired correlation plot

```
from sequencing_tools.douglas_colors import cor_plot
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('white')

d = np.random.rand(100,10)
d = pd.DataFrame(d)
fig =cor_plot(d)
fig.savefig('cor.png',bbox_inches='tight')
```

![](https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/cor.png)
