# TGIRT-seq tools #

This repo stores different utils for tgirt-seq pipeline



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
