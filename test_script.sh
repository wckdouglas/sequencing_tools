set -e

python --version
python -c "\
try:
    import tgirt_seq_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    

bam_to_bed.py -h
bam_umi_tag.py -h
bedpe_to_bed.py -h
clip_fastq.py -h
depth_to_bigwig.py -h
deinterleave_fastq.py -h
filter_soft_clip.py -h
reduce_multi_reads.py -h
split_bam.py -h
split_uniq_bam.py -h
stranded_base_count.py -h
pe_fq_merge.py -h
filter_umi.py -h
