# sequencing_tools.consensus_tools #

## calculate QUAL character  ##

input:
    prob: probability of the base call being an correct base call 

output:
    qual_char: a Phred score character

example usage:
```
$ from sequencing_tools.consensus_tools import prob_to_qual_string
$ prob_to_qual_string(0.9999999999)
'I'
$ ct.prob_to_qual_string(0.8)
"'"
```

## calculate base call probability from base quality  ##

input:
    base_qual: a list of illumina numeric Phred score  ([33 - 73])

output:
    iterator of base call error error probability

example usage:
```
$ from sequencing_tools.consensus_tools import qual_to_prob
$ phred = 'III++ABC'
$ phred = map(lambda x: ord(x) - 33, phred)
$ [q for q in qual_to_prob(pred)]
[0.0001,
 0.0001,
 0.0001,
 0.1,
 0.1,
 0.000630957344480193,
 0.0005011872336272725,
 0.00039810717055349735]
```

## calculate posterior probability for consensus base given a list of nucleotides ##

input:
    column_bases: a numpy list of nucleotide that, such as ['A','C','A','A','A']
    column_qualities: a numpy list of qual score resepective to the column bases, such as [43, 10, 25, 30, 39]
    guess_base: character [A, C, T or G] that matches one of the column)baess

output:
    log_posterior: log posterior probability of the 'guess_base' is the correct consensus base

example usage:
```
$ from sequencing_tools.consensus_tools import calculatePosterior
$ import numpy as np
$ bases = np.array(list('ACAAA'))
$ quals = np.array([43, 10, 25, 30, 39])
$ for b in ['A', 'C']:
$    print(calculatePosterior(bases, quals, b))
-3.405541190667555
-36.045225444348695
```


## Error correction in fastq records ##

input:
    mode:  prob or vote
    threshold: only consider for "vote" mode, as a cutoff for returning a "N" if not enough fraction of bases agree
    
output:
    seq: string of consensus nucleotides
    qual: string of quality score resepective to the consensus bases

example usage for using a posterior probability method that is suggested by (https://github.com/fulcrumgenomics/fgbio/wiki/Calling-Consensus-Reads):
```
$ from sequencing_tools.consensus_tools import ErrorCorrection
$ ec = ErrorCorrection(mode='prob')
$ seq, qual = ec.Correct(['AAACA','AAAAA','AAACA','AAACA'],
            ['IIIII','IIIAI','FFFFF', 'FFFFF'])
$ print( (seq, qual) )
('AAACA', 'IIIII')
```

Example usage for using a "voting" mechanism that is suggested by [SafeSeq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111315/):
```
$ from sequencing_tools.consensus_tools import ErrorCorrection
$ ec = ErrorCorrection(mode='vote')
$ seq, qual = ec.Correct(['AAACA','AAAAA','AAACA','AAACA'],
            ['IIIII','IIIAI','FFFFF', 'FFFFF'])
$ print( (seq, qual) )
('AAANA', 'III!I')
```


## Generate consensus sequencing from Bam file ##

input:
    bam: bam file path
    min_cov: minimum coverage to consider, any position lower than this will have "N" base

output:
    consensus sequence

usage:
```
from sequencing_tools.consensus_tools import ConsensusAlignments
ca = ConsensusAlignments(bam_file)
ca.concensus_seq('chrM',1000,1100)
```
