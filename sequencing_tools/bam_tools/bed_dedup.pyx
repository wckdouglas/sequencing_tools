from __future__ import print_function
import sys
import fileinput 
from operator import itemgetter
from itertools import combinations, groupby
from functools import partial
from collections import Counter, defaultdict
from networkx import Graph, connected_components
from sequencing_tools.stats_tools import hamming_distance, levenshtein_distance
from sequencing_tools.bam_tools.umi_network import demultiplex_directional, demultiplex_adj
import sys
import six
from libc.stdint cimport uint32_t
long = six.integer_types[-1]

cdef class fragment_group:
    '''
    Data structure for storing fragments with same start, end positions and strands
    store barcode as a dictionary with values indicating numbers of same fragment
    '''

    cdef: 
        str chrom, start, end, strand
        long fragment_size
        list unique_barcodes
        str cigar
        set cigarlist
        dict barcodes_set


    def __init__(self, chrom, start, end, strand, bc, cigar):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.fragment_size = long(self.end) - long(self.start)

        self.barcodes_set = dict()
        self.unique_barcodes

        # add first record
        self.barcodes_set[cigar] = self.barcodes_set.setdefault(cigar, dict()) 
        self.barcodes_set[cigar][bc] = self.barcodes_set[cigar].setdefault(bc, 0) + 1

    def add_member(self, bc, cigar):
        '''
        add a member to the fragment group, add barcode to barcode dictionary
        '''
        self.barcodes_set[cigar] = self.barcodes_set.setdefault(cigar, dict()) 
        self.barcodes_set[cigar][bc] = self.barcodes_set[cigar].setdefault(bc, 0) + 1


    def demultiplexing_barcodes(self, int threshold):
        '''
        Demultiplex the barcode list
        '''

        cdef:
            str cigar
            dict barcodes_dict
            int _max_member_count
            long max_member_count = 0
            str _barcode_name
            list _barcodes
            list _temp_unique_barcodes 

        self.unique_barcodes = []
        for cigar, barcodes_dict in six.iteritems(self.barcodes_set):
            _temp_unique_barcodes = []
            if threshold > 0:
                if len(barcodes_dict.keys()) == 1: # for the singular fragment
                                                   # generate a phantom list with the single fragment record
                    _max_member_count = list(barcodes_dict.values())[0]
                    _barcode_name = '{barcode}_{count}_members'.format(barcode = list(barcodes_dict.keys())[0], 
                                                                    count = _max_member_count)
                    _temp_unique_barcodes.append(_barcode_name)

                else: # for more than 1 unique barcode
                    #_barcodes, _max_member_count = demultiplex_adj(barcodes_dict, threshold = threshold)
                    _barcodes, _max_member_count = demultiplex_directional(barcodes_dict, threshold = threshold)
                    _temp_unique_barcodes.extend(_barcodes)


            else: # if no error toleration, all unique UMI is their own cluster
                _barcodes = ['{barcode}_{count}_members'.format(barcode = barcode, count = count) \
                                        for barcode, count in six.iteritems(barcodes_dict)]
                _temp_unique_barcodes.extend(_barcodes)
                _max_member_count = max(six.itervalues(barcodes_dict))


            if cigar: # if no cigarstring, it will be empty string, return false in here
                _temp_unique_barcodes = list(map(lambda x: x + '_' + cigar, _temp_unique_barcodes))
            

            self.unique_barcodes.extend(_temp_unique_barcodes)
            max_member_count = max(max_member_count, _max_member_count)
        return max_member_count


    def output_bed_line(self, file=sys.stdout):
        '''
        output every unique fragments
        '''
        cdef:
            str barcode
            str template
            uint32_t i
            int member_count = 0

        for i, barcode in enumerate(self.unique_barcodes):
            member_count += int(barcode.split('_')[1])
            template = '{chrom}\t{start}\t{end}\t{barcode}\t{length}\t{strand}' \
                .format(chrom = self.chrom,
                    start = self.start,
                    end = self.end,
                    barcode = barcode,
                    length = self.fragment_size,
                    strand= self.strand)
            print(template, file = file)
        return i + 1, member_count
        

    def check_fragment(self, str chrom, str start, str end, str strand):
        '''
        check if new fragment belong to this fragment group
        '''
        cdef:
            bint chrom_same, start_end_same, strand_same

        chrom_same = (chrom == self.chrom)
        start_end_same = (start == self.start and end == self.end)
        strand_same = (strand == self.strand)
        return chrom_same and start_end_same and strand_same

    
    def get_unique_umi(self):
        return self.unique_barcodes


def dedup_bed(in_file_handle, out_file_handle, threshold, str delim, int f, int ct):
    cdef:
        str line, bc, read_name, chrom, start, end, strand
        str bc_line
        uint32_t in_count, _out_count, _member_count
        long out_count = 0, member_count
        fragment_group barcode_group
        str cigar = ''
        list fields
        long max_member_count = 0

    barcode_group = None
    assert((ct > 5 and isinstance(ct, int)) or ct == -1, 'Unacceptable cigar field, needs to be an integer > 5')
    for in_count, line in enumerate(in_file_handle):
        fields = line.strip().split('\t')
        read_name = fields[3]
        bc = read_name.split(delim)[f]
        chrom, start, end, strand = itemgetter(0,1,2,5)(fields)

        cigar = fields[ct] if ct > 0 else ''  # only make cigar field when ct is not -1

        if not barcode_group:
            '''
            if no initial barcode group, initilize one
            '''
            barcode_group = fragment_group(chrom, start, end, strand, bc, cigar)

        elif barcode_group.check_fragment(chrom, start, end, strand):
            '''
            add member if same chrom, start, end & strand
            '''
            barcode_group.add_member(bc, cigar)

        else:
            '''
            summarize previous group and
            initialize new fragment group if not same coordinate
            '''
            _max_member_count = barcode_group.demultiplexing_barcodes(threshold)
            _out_count, _member_count = barcode_group.output_bed_line(file = out_file_handle)
            member_count += _member_count
            out_count += _out_count
            max_member_count = max(max_member_count, _max_member_count)

            barcode_group = fragment_group(chrom, start, end, strand, bc, cigar)

    # output the last group
    _max_member_count = barcode_group.demultiplexing_barcodes(threshold)
    _out_count, _member_count = barcode_group.output_bed_line(file = out_file_handle)
    member_count += _member_count
    out_count += _out_count
    max_member_count = max(max_member_count, _max_member_count)
    assert(member_count == in_count + 1, 'Wrong output lines!')

    print('Iput %i lines, output %i lines with highest duplicate with %i members' %(in_count + 1, out_count, max_member_count), file=sys.stderr)
    return 0
