set -e

python --version
python -c "\
try:
    import tgirt_seq_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    
bam_to_bed.py
bam_umi_tag.py
bedpe_to_bed.py
clip_fastq.py
deinterleave_fastq.py
depth_to_bigwig.py
filterSoftClip.py
reduce_multi_reads.py
split_bam.py
split_uniq_bam.py
stranded_base_count.py
