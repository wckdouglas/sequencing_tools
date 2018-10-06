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

class fragment_group:
    '''
    Data structure for storing fragments with same start, end positions and strands
    store barcode as a dictionary with values indicating numbers of same fragment
    '''


    def __init__(self, coordinates, fragment_list, str umi_delim = '_', int umi_f = 0, int ct = -1):
        chrom, start, end, strand = coordinates
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.umi_f = umi_f
        self.umi_delim = umi_delim
        self.ct =ct
        self.fragment_list = fragment_list


        self.fragment_size = long(self.end) - long(self.start)
        self.fragment_count = 0
        self.barcodes_set = dict()
        self.unique_barcodes = []
        self.outline_template = '{chrom}\t{start}\t{end}\t{barcode}\t{length}\t{strand}' 


    def resolve_fragments(self):
        '''
        loop over iterator: fragment_list,
        extract cigar and umi

        '''
        cdef:
            str fragment
            list fields
            str cigar, umi
            int fragment_count

        for fragment_count, fragment in enumerate(self.fragment_list):
            fields = fragment.strip().split('\t')
            umi = fields[3].split(self.umi_delim)[self.umi_f]
            cigar = fields[self.ct]
            self.barcodes_set[cigar] = self.barcodes_set.setdefault(cigar, dict()) 
            self.barcodes_set[cigar][umi] = self.barcodes_set[cigar].setdefault(umi, 0) + 1
        self.fragment_count = fragment_count + 1


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
            str umi
            int member_count

        self.unique_barcodes = []

        for cigar, barcodes_dict in six.iteritems(self.barcodes_set):
            _temp_unique_barcodes = []

            if threshold > 0: 
                '''
                with toleration
                demultiplex or output single fragment
                '''
                if len(barcodes_dict.keys()) == 1: 
                    '''
                    for the singular fragment
                    generate a phantom list with the single fragment record
                    '''
                    _max_member_count = list(barcodes_dict.values())[0]
                    _barcode_name = '{barcode}_{count}_members'.format(barcode = list(barcodes_dict.keys())[0], 
                                                                    count = _max_member_count)
                    _temp_unique_barcodes.append(_barcode_name)

                else: # for more than 1 unique barcode
                    #_barcodes, _max_member_count = demultiplex_adj(barcodes_dict, threshold = threshold)
                    _barcodes, _max_member_count = demultiplex_directional(barcodes_dict, threshold = threshold)
                    _temp_unique_barcodes.extend(_barcodes)


            else: 
                ''' 
                if no error toleration, 
                all unique UMI is their own cluster
                '''
                _barcodes = ['{barcode}_{count}_members'.format(barcode = umi, count = member_count) \
                                        for umi, member_count in six.iteritems(barcodes_dict)]
                _temp_unique_barcodes.extend(_barcodes)
                _max_member_count = max(six.itervalues(barcodes_dict))


            if self.ct > 5:   
                '''
                if no cigarstring, 
                it will be empty string, 
                return false in here 
                '''
                _temp_unique_barcodes = list(map(lambda x: x + '_' + cigar, _temp_unique_barcodes))
            

            self.unique_barcodes.extend(_temp_unique_barcodes)
            max_member_count = max(max_member_count, _max_member_count)
        return max_member_count


    def print_demultiplexed(self, out_file=sys.stdout):
        '''
        output every unique fragments
        '''
        cdef:
            str barcode
            str outline
            uint32_t i
            int member_count = 0

        for i, barcode in enumerate(self.unique_barcodes):
            member_count += int(barcode.split('_')[1])
            outline = self.outline_template\
                .format(chrom = self.chrom,
                        start = self.start,
                        end = self.end,
                        barcode = barcode,
                        length = self.fragment_size,
                        strand= self.strand)
            print(outline, file = out_file)
        return i + 1, member_count

    
    def get_unique_umi(self):
        return self.unique_barcodes


cpdef fragment_coordinates(str bedline):
    fields = bedline.strip().split('\t')
    return  itemgetter(0,1,2,5)(fields)


def dedup_bed(in_file_handle, out_file_handle, threshold, str delim, int f, int ct):
    cdef:
        str bc, read_name, chrom, start, end, strand
        str bc_line
        uint32_t in_count = 0 , _out_count, _member_count
        long out_count = 0, total_member_count = 0
        str cigar = ''
        list fields
        long max_member_count = 0

    barcode_group = None
    assert((ct > 5 and isinstance(ct, int)) or ct == -1, 'Unacceptable cigar field, needs to be an integer > 5')

    for coordinates, fragment_list in groupby(in_file_handle, fragment_coordinates):
        barcode_group = fragment_group(coordinates, 
                                        fragment_list, 
                                        umi_delim = delim, 
                                        umi_f = f, 
                                        ct = ct)
        barcode_group.resolve_fragments()
        _max_member_count = barcode_group.demultiplexing_barcodes(threshold = threshold)
        _out_count, _member_count = barcode_group.print_demultiplexed(out_file = out_file_handle)
        total_member_count += _member_count
        out_count += _out_count
        max_member_count = max(max_member_count, _max_member_count)
        in_count += barcode_group.fragment_count

    assert(total_member_count == in_count + 1, 'Wrong output lines!')

    print('Iput %i lines, output %i lines with highest duplicate with %i members' %(in_count + 1, out_count, max_member_count), file=sys.stderr)
    return 0
