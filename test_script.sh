set -e

python --version
python -c "\
try:
    import sequencing_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    

seqtools bam2bed -h
seqtools demux -h
seqtools bamTag -h
seqtools bedpe2bed -h
seqtools clipFQ -h
seqtools dedup -h
seqtools deinterleave -h
seqtools cov2bw -h
seqtools correctFQ -h
seqtools filterSoftClip -h
seqtools filterUMI -h
seqtools mergepe -h
seqtools calcUMI -h
seqtools filterMulti -h
seqtools spliceBam -h
seqtools splitBam -h
seqtools pileup -h
