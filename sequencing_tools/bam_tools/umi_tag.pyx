from __future__ import print_function

import re
import sys
from functools import partial

import pysam
from pysam.libcalignmentfilecimportAlignmentFile import AlignedSegment

cpdef add_umi_tag(str in_bam, str out_bam, str tag, str delim, int frag):
    cdef:
        AlignmentFile inbam, outbam
        AlignedSegment aln
        int aln_count
        str id, umi

    print('Parsing from %s to %s' %(in_bam, out_bam), file = sys.stderr)
    with pysam.Samfile(in_bam, 'rb') as inbam:
        with pysam.Samfile(out_bam,'wb', template=inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                id = aln.query_name
                splitted_id = id.split(delim)
                umi = splitted_id[frag]
                aln.query_name = id.replace(umi + delim, '').replace(delim+umi, '')
                aln.tags += [(tag,umi)]
                outbam.write(aln)
    print('Parsed %i alignments from %s to %s' %(aln_count, in_bam, out_bam), file = sys.stderr)
    return 0


cdef str get_umi_from_tag(str tag, AlignedSegment aln):
    return aln.get_tag(tag)

cdef str get_umi_from_name(AlignedSegment aln):
    cdef:
        str aln_id = aln.query_name
    return  aln_id.split('_')[0]

def filter_umi(inbam, outbam, threshold, tag, prefix):
    cdef:
        int out_count = 0
        int aln_count = 0
        str umi 
        AlignedSegment aln 

    threshold = int(threshold)
    bases = [(threshold + 1) * x for x in 'ACTG']
    nucleotide_runs = re.compile('|'.join(bases))

    get_umi = partial(get_umi_from_tag, tag) if tag else partial(get_umi_from_name)
    with pysam.Samfile(inbam) as inbam:
        with pysam.Samfile(outbam,'w', template = inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                umi = get_umi(aln)
                if not nucleotide_runs.search(umi):
                    outbam.write(aln)
                    out_count += 1
    print('Parsed %i alignments\noutput %i alignments' %(aln_count, out_count), file=sys.stderr)
    return 0
