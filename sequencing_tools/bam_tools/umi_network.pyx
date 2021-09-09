from collections import defaultdict
from itertools import combinations

import six
from networkx import Graph, connected_components

from sequencing_tools.stats_tools import hamming_distance, levenshtein_distance


'''
Core code copied from umi_tools
https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/network.py
'''

def iter_nearest_neighbours(umis, substr_idx):
    '''
    Added by Matt 06/05/17
    use substring dict to get (approximately) all the nearest neighbours to
    each in a set of umis.
    '''
    cdef:
        int i 
        str u
        str nbr
        set neighbours 

    for i, u in enumerate(umis, 1):
        neighbours = set()
        for idx, substr_map in six.iteritems(substr_idx):
            u_sub = u[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.difference_update(umis[:i])
        for nbr in neighbours:
            yield u, nbr


def get_substr_slices(umi_length, idx_size):
    '''
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    '''
    cdef:
        int cs, r, s

    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    for s in sub_sizes:
        yield(offset, offset + s)
        offset += s


def build_substr_idx(umis, umi_length, min_edit):
    '''
    Build a dictionary of nearest neighbours using substrings, can be used
    to reduce the number of pairwise comparisons.
    '''

    cdef:
        str u

    substr_idx = defaultdict(
        lambda: defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx


def get_adj_list_directional(umis, counts, threshold=1):
    ''' identify all umis within the hamming distance threshold
    and where the counts of the first umi is > (2 * second umi counts)-1
    output:
        {'UMI node1': [sub node1, sub node2],
        'umi node2': []}
    '''

    cdef:
        str umi1, umi2

    adj_list = {umi: [] for umi in umis}
    if len(umis) > 25:
        umi_length = len(umis[0])
        substr_idx = build_substr_idx(umis, umi_length, threshold)
        iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
    else:
        iter_umi_pairs = combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if hamming_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2]*2)-1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1]*2)-1:
                adj_list[umi2].append(umi1)

    return adj_list


def breadth_first_search(node, adj_list):
    cdef:
        set searched = set()
        set queue = set()
        str next_node

    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def get_connected_components_adjacency(umis, graph, counts):
    ''' find the connected UMIs within an adjacency dictionary'''

    # TS: TO DO: Work out why recursive function doesn't lead to same
    # final output. Then uncomment below

    # if len(graph) < 10000:
    #    self.search = breadth_first_search_recursive
    # else:
    #    self.search = breadth_first_search

    cdef:
        set found = set()
        list components = list()
        str node

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            # component = self.search(node, graph)
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)
    return components



def group_directional(clusters, adj_list, counts):
    ''' return groups for directional method'''

    cdef:
        set observed = set()
        list groups = []
        str node
        list temp_cluster

    for cluster in clusters:
        if len(cluster) == 1:
            groups.append(list(cluster))
            observed.update(cluster)
        else:
            cluster = sorted(cluster, key=lambda x: counts[x],
                                reverse=True)
            # need to remove any node which has already been observed
            temp_cluster = []
            for node in cluster:
                if node not in observed:
                    temp_cluster.append(node)
                    observed.add(node)
            groups.append(temp_cluster)

    return groups



def demultiplex_directional(umi_counter, threshold=1):
    '''
    umi counter is a counter for each umi (dictionary),
    ''' 

    cdef:
        int max_member_count = 0
        list out_umis = []
        list sub_nodes
        str umi
        list umi_group

    umis = list(umi_counter.keys())
    adj_list = get_adj_list_directional(umis, umi_counter)
    clusters = get_connected_components_adjacency(umis, adj_list, umi_counter)
    umi_groups = group_directional(clusters, adj_list, umi_counter)

    out_umis = []
    for umi_group in umi_groups:
        member_count = sum(umi_counter[umi] for umi in umi_group)
        max_member_count = max(max_member_count, member_count)
        umi_id = umi_group[0] + '_' + str(member_count) + '_members'
        out_umis.append(umi_id)
    return out_umis, max_member_count



# this should be the adjacnecy method defined in umi_tools

def make_graph(comparison, int threshold):
    '''
    Using a graph to connect all umi with <= threshold mismatches
    '''
    cdef:
        str umi_1, umi_2

    G = Graph()
    for umi_1, umi_2 in comparison:
        if levenshtein_distance(umi_1, umi_2) <= threshold:
            G.add_edge(umi_1, umi_2)
        else:
            G.add_node(umi_1)
            G.add_node(umi_2)
    return G

def unique_barcode_from_graph(graph, barcodes):
    '''
    Merging barcode families, using the pre-built network-of-barcode 
    '''
    cdef:
        set subgraph
        list subgraph_list
        str bc, barcode_id
        int member_count
        list unique_barcode = []
        long max_member_count = 0

    for subgraph in connected_components(graph):
        subgraph_list = list(subgraph)
        if len(subgraph_list) == 1:
            barcode_id = subgraph_list[0]
            member_count = barcodes[barcode_id]
            barcode_id = barcode_id + '_' + str(member_count) + '_members'
        else:
            member_count = sum(barcodes[bc] for bc in subgraph)
            barcode_id = subgraph_list[0] + '_' + str(member_count) + '_members'
        unique_barcode.append(barcode_id)
        max_member_count = max(max_member_count, member_count)
    return unique_barcode, max_member_count

def demultiplex_adj(barcodes, threshold=1):
    '''
    demultiplexing barcode families
    '''

    comparison = combinations(barcodes.keys(),r=2)
    graph = make_graph(comparison, threshold)
    unique_barcode, max_member_count =  unique_barcode_from_graph(graph, barcodes)
    return unique_barcode, max_member_count
