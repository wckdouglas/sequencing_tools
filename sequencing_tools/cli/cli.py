import argparse

from sequencing_tools.cli import (bam_read_cluster, bam_to_bed, bam_umi_tag,
                                  bedpe_to_bed, clip_fastq, deduplicate_bed,
                                  deinterleave_fastq, depth_to_bigwig,
                                  fastq_pe_cluster, filter_soft_clip,
                                  filter_umi, pe_fq_merge,
                                  poisson_umi_adjustment, reduce_multi_reads,
                                  split_bam, split_uniq_bam,
                                  stranded_base_count)


def add_opt(subparsers):
    # add subparsers
    bam_to_bed.getopt(subparsers)  # bam2bed
    bam_read_cluster.getopt(subparsers)  # demux
    bam_umi_tag.getopt(subparsers)  # bamTag
    bedpe_to_bed.getopt(subparsers)  # bedpe2bed
    clip_fastq.getopt(subparsers)  # clipFQ
    deduplicate_bed.getopt(subparsers)  # dedup
    deinterleave_fastq.getopt(subparsers)  # deinterleave
    depth_to_bigwig.getopt(subparsers)  # cov2bw
    fastq_pe_cluster.getopt(subparsers)  # correctFQ
    filter_soft_clip.getopt(subparsers)  # filterSoftClip
    filter_umi.getopt(subparsers)  # filterUMI
    pe_fq_merge.getopt(subparsers)  # mergepe
    poisson_umi_adjustment.getopt(subparsers)  # calcUMI
    reduce_multi_reads.getopt(subparsers)  # filterMulti
    split_bam.getopt(subparsers)  # spliceBam
    split_uniq_bam.getopt(subparsers)  # splitBam
    stranded_base_count.getopt(subparsers)  # pileup


def main():
    parser = argparse.ArgumentParser(description="Tools for manipulating NGS data")
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True
    add_opt(subparsers)

    args = parser.parse_args()
    if args.subcommand == "bam2bed":
        bam_to_bed.run(args)

    elif args.subcommand == "demux":
        bam_read_cluster.run(args)

    elif args.subcommand == "bamTag":
        bam_umi_tag.run(args)

    elif args.subcommand == "bedpe2bed":
        bedpe_to_bed.run(args)

    elif args.subcommand == "clipFQ":
        clip_fastq.run(args)

    elif args.subcommand == "dedup":
        deduplicate_bed.run(args)

    elif args.subcommand == "deinterleave":
        deinterleave_fastq.run(args)

    elif args.subcommand == "cov2bw":
        depth_to_bigwig.run(args)

    elif args.subcommand == "correctFQ":
        fastq_pe_cluster.run(args)

    elif args.subcommand == "filterSoftClip":
        filter_soft_clip.run(args)

    elif args.subcommand == "filterUMI":
        filter_umi.run(args)

    elif args.subcommand == "mergepe":
        pe_fq_merge.run(args)

    elif args.subcommand == "calcUMI":
        poisson_umi_adjustment.run(args)

    elif args.subcommand == "filterMulti":
        reduce_multi_reads.run(args)

    elif args.subcommand == "spliceBam":
        split_bam.run(args)

    elif args.subcommand == "spitBam":
        split_uniq_bam.run(args)

    elif args.subcommand == "pileup":
        stranded_base_count.run(args)
