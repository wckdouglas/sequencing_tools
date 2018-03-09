from __future__ import print_function
import sys
import fileinput 
from operator import itemgetter
from itertools import combinations, groupby
from functools import partial
from collections import Counter, defaultdict
from networkx import Graph, connected_components
from sequencing_tools.stats_tools import hamming_distance, levenshtein_distance
import sys
import six


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
        self.unique_barcodes = []

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

        for cigar, barcodes_dict in six.iteritems(self.barcodes_set):
            if threshold > 0:
                if len(barcodes_dict.keys()) == 1: # for theingular fragment
                                                    # generate a phantom list with the single fragment record
                    self.unique_barcodes = ['{barcode}_{count}_members'.format(barcode = list(barcodes_dict.keys())[0], 
                                                                    count = list(barcodes_dict.values())[0])]
                else: # for more than 1 unique barcode
                    self.unique_barcodes = demultiplex(barcodes_dict, threshold = threshold)
                    if cigar: # if no cigarstring, it will be empty string, return false in here
                        self.unique_barcodes = map(lambda x: x + '_' + cigar, self.unique_barcodes)

            else: # if no error toleration, all unique UMI is their own cluster
                self.unique_barcodes = ['{barcode}_{count}_members'.format(barcode = barcode, count = count) \
                                        for barcode, count in six.iteritems(barcodes_dict)]

                if cigar: # if no cigarstring, it will be empty string, return false in here
                    self.unique_barcodes = map(lambda x: x + '_' + cigar, self.unique_barcodes)


    def output_bed_line(self):
        '''
        output every unique fragments
        '''
        cdef:
            str barcode
            str template

        for barcode in self.unique_barcodes:
            template = '{chrom}\t{start}\t{end}\t{barcode}\t{length}\t{strand}' \
                .format(chrom = self.chrom,
                    start = self.start,
                    end = self.end,
                    barcode = barcode,
                    length = self.fragment_size,
                    strand= self.strand)
            yield template

    def check_fragment(self, chrom, start, end, strand):
        '''
        check if new fragment belong to this fragment group
        '''
        chrom_same = (chrom == self.chrom)
        start_end_same = (start == self.start and end == self.end)
        strand_same = (strand == self.strand)
        return chrom_same and start_end_same and strand_same


cpdef int barcode_distance(barcode_pair):
    '''
    calculating hamming distance of a barcode pair
    inpurt is tuple of two barcodes
    '''
    cdef:
        str a, b
        str i, j

    a, b = barcode_pair
    return levenshtein_distance(a, b)

def make_graph(comparison, threshold):
    '''
    Using a graph to connect all umi with <= threshold mismatches
    '''
    G = Graph()
    for pair in comparison:
        if barcode_distance(pair) <= threshold:
            G.add_edge(pair[0],pair[1])
        else:
            G.add_node(pair[0])
            G.add_node(pair[1])
    return G

def unique_barcode_from_graph(graph, barcodes):
    '''
    Merging barcode families, using the pre-built network-of-barcode 
    '''
    unique_barcode = []
    for subgraph in connected_components(graph):
        subgraph = list(subgraph)
        if len(subgraph) == 1:
            barcode_id = subgraph[0]
            member_counts = barcodes[barcode_id]
            barcode_id = barcode_id + '_' + str(member_counts) + '_members'
        else:
            member_counts = sum(barcodes[bc] for bc in subgraph)
            barcode_id = subgraph[0] + '_' + str(member_counts) + '_members'
        unique_barcode.append(barcode_id)
    return unique_barcode

def demultiplex(barcodes, threshold=1):
    '''
    demultiplexing barcode families
    '''
    comparison = combinations(barcodes.keys(),r=2)
    graph = make_graph(comparison, threshold)
    unique_barcode =  unique_barcode_from_graph(graph, barcodes)
    return unique_barcode


def dedup_bed(in_file_handle, out_file_handle, threshold, str delim, int f, int ct):
    cdef:
        str line, bc, read_name, chrom, start, end, strand
        str bc_line
        int in_count
        int out_count = 0
        #fragment_group barcode_group
        str cigar = ''

    barcode_group = None
    for in_count, line in enumerate(in_file_handle):
        fields = line.strip().split('\t')
        read_name = fields[3]
        bc = read_name.split(delim)[f]
        chrom, start, end, strand = itemgetter(0,1,2,5)(fields)

        if ct > 5:
            cigar = fields[ct]

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
            barcode_group.demultiplexing_barcodes(threshold)
            for bc_line in barcode_group.output_bed_line():
                print(bc_line, file=out_file_handle)
                out_count += 1

            barcode_group = fragment_group(chrom, start, end, strand, bc, cigar)

    # output the last group
    barcode_group.demultiplexing_barcodes(threshold)
    for bc_line in barcode_group.output_bed_line():
        print(bc_line, file=out_file_handle)
        out_count += 1

    print('Iput %i lines, output %i lines' %(in_count + 1, out_count), file=sys.stderr)
    return 0
