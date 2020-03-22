# sequencing_tools.fasta_tools #

## reading fasta file ##

Iterate over a fasta file

input:
    file handle of a fasta file

Example usage:
```
$ from sequencing_tools.fastq_tools import readfa
$ with open('test.fa') as fa:
$     for seqid, seq in readfa(fa):
$         print(seqid, seq)
```


## reading multiple alignment fasta file ##

Plotting and making consensus sequence from multiple alignments

inut:
    fa_file: fasta file name


Example usage:
```
    ma = multi_alignment("fasta file")

    ax = plt.subplot(111)
    ma.plot(ax = ax) # plotting multiple alignment
    
    consensus_seq, scores = ma.concensus() # making consensus sequence, score is shows the proportion of sequence having the consensus base

    matrix = ma.PairMatrix() # Pairwise hamming distance matrix computing for each pair of sequences
```

