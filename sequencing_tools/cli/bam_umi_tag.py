from sequencing_tools.bam_tools.umi_tag import add_umi_tag


def getopt(subparsers):
    parser = subparsers.add_parser(
        name="bamTag",
        description="Putting UMI as RX tag in bam, enables picard Markduplicates with BARCODE_TAG=RX",
    )
    parser.add_argument(
        "-i", "--in_bam", required=True, help="BAM file name, or stdin (-)"
    )
    parser.add_argument(
        "-o", "--out_bam", default="-", help="BAM file output (default: - )"
    )
    parser.add_argument("-t", "--tag", default="RX", help="Tag id (default: RX )")
    parser.add_argument(
        "-d",
        "--delim",
        default="_",
        help="Deliminator separating read id and bc (default: _ )",
    )
    parser.add_argument(
        "-f",
        "--fragment",
        default="0",
        help="after splitting read name using {delim}, which"
        + "fragment is UMI? can use -1 as last piece (default: 0)",
    )


def run():
    args = getopt()
    add_umi_tag(args.in_bam, args.out_bam, args.tag, args.delim, int(args.fragment))


if __name__ == "__main__":
    main()
