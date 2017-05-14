# TGIRT-seq tools #

This repo stores different utils for tgirt-seq pipeline


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
pysam>0.10.0
numpy
cython

```
pip install pysam==0.9.0
pip install numpy
pip install cython
```
---

## Split unique bam ##

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

---

## Reduce multi reads ##

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

---

## BAM to bed ##
from pair-end BAM, get fragment coordinate as BED file

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
