import pysam
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import sys

def add_umi_tag(in_bam, out_bam):
    cdef:
        AlignmentFile inbam, outbam
        AlignedSegment aln
        int aln_count
        str id, umi

    print >> sys.stderr, 'Parsing from %s to %s' %(in_bam, out_bam)
    with pysam.Samfile(in_bam, 'rb') as inbam:
        with pysam.Samfile(out_bam,'wb', template=inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                id = aln.query_name
                splitted_id = id.split('_')
                umi = splitted_id[0]
                aln.query_name = splitted_id[1]
                aln.tags += [('RX',umi)]
                outbam.write(aln)
    print >> sys.stderr, 'Parsed %i alignments from %s to %s' %(aln_count, in_bam, out_bam)
    return 0
