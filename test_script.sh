set -e

python --version
python -c "\
try:
    import tgirt_seq_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    

bin/bam_to_bed.py -h
bin/bam_umi_tag.py -h
bin/bedpe_to_bed.py -h
bin/clip_fastq.py -h
bin/deinterleave_fastq.py -h
bin/depth_to_bigwig.py -h
bin/filterSoftClip.py -h
bin/reduce_multi_reads.py -h
bin/split_bam.py -h
bin/split_uniq_bam.py -h
bin/stranded_base_count.py -h
