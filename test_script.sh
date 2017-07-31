set -e

python --version
python -c "\
try:
    import tgirt_seq_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    
which bam_to_bed.py
which bam_umi_tag.py
which bedpe_to_bed.py
which clip_fastq.py
which deinterleave_fastq.py
which depth_to_bigwig.py
which filterSoftClip.py
which reduce_multi_reads.py
which split_bam.py
which split_uniq_bam.py
which stranded_base_count.py
