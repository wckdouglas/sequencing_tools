set -e

python --version
python -c "\
try:
    import tgirt_seq_tools
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
    
python bin/bam_to_bed.py -h
python bin/bam_umi_tag.py -h
python bin/bedpe_to_bed.py -h
python bin/clip_fastq.py -h
python bin/deinterleave_fastq.py -h
python bin/depth_to_bigwig.py -h
python bin/filterSoftClip.py -h
python bin/reduce_multi_reads.py -h
python bin/split_bam.py -h
python bin/split_uniq_bam.py -h
python bin/stranded_base_count.py -h
