import argparse
import logging
import sys
import numpy as np


def map_index_to_kmer(index, kmer_len):
    """Given some index, map it to a k-mer

    Args:
        index (int): [0, 4^k - 1] index that will map back to its string
        kmer_len (int): The kmer length used to calculate it
    """
    mapping = ['A', 'C', 'G', 'T']

    result = ''
    for _ in range(kmer_len):
        result = mapping[index % 4] + result
        # shift the base 4 string to the right
        index = index // 4

    return result

def build_graph(pdist, kmer_len):
    """Based on some probability distribution given, build the resultant
    debrujin graph.
    """

    # now we construct a graph with these kmer counts    
    graph = {}
    for i in range(len(pdist)):
        if pdist[i] == 0:
            continue
        kmer = map_index_to_kmer(i, kmer_len)
        kmer_from, kmer_to = kmer[0:kmer_len - 1], kmer[1: kmer_len]

        if kmer_from not in graph.keys():
            graph[kmer_from] = {}
        
        if kmer_to not in graph.keys():
            graph[kmer_to] = {}
        
        if kmer_to not in graph[kmer_from]:
            graph[kmer_from][kmer_to] = 0

        graph[kmer_from][kmer_to] += pdist[i]

    return graph


def count_edges(graph):
    """Given a graph, calculates the indegree and outdegree for each node."""
    edge_counts = {}    

    for node in graph.keys():
        edge_counts[node] = {
            'indeg': 0,
            'outdeg': sum(list(graph[node].values()))
        }

    for node in graph.keys():
        for key, count in graph[node].items():
            edge_counts[key]['indeg'] += count
    
    return edge_counts


# does a single pass of depth first search, stopping when we cannot go anymore
def dfs(graph, start):
    to_add = ''
    node_visit = start

    # do a dfs to find the path starting at start
    while True:
        next_node = None

        # find a node with a non-zero value
        for node, value in graph[node_visit].items():
            if value > 0:
                next_node = node
                break

        if next_node is None:
            break
        
        assert graph[node_visit][next_node] > 0
        graph[node_visit][next_node] -= 1
        to_add += next_node[-1]
        node_visit = next_node

    return to_add, graph


def find_eulerian_path(graph, start, k):
    """Given a graph, and the starting point, finds the Eulerian path.

    Args:
        graph (Dict[str, Dict[str, int]]): Graph with a set Eulerian path.
        start: The starting kmer.

    Returns:
        str: The final constructed genome.
    """
    # first grab the possible sequence
    to_add, graph = dfs(graph, start)
    seq = start + to_add

    # now we modify the start_seq by finding some node within the start sequence
    # that has an edge left
    while True:
        node_left_with_edge = None
        
        # these are the nodes already added to our sequence
        all_kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]

        possible_nodes = set()

        for kmer in all_kmers:
            possible_nodes.add(kmer[0:len(kmer) - 1])
            possible_nodes.add(kmer[1:len(kmer)])

        logging.debug(possible_nodes)

        for node in possible_nodes:
            assert node in graph

            for value in graph[node].values():
                if value > 0:
                    node_left_with_edge = node
                    break

            if node_left_with_edge is not None:
                break

        # stop if we have no more nodes left
        if node_left_with_edge is None:
            break

        # now we dfs on the node
        logging.debug(f'Num edges: {sum([sum(value for value in graph[node].values()) for node in graph.keys()])}')
        add_seq, graph = dfs(graph, node_left_with_edge)
        logging.debug(f'Current length to add: {len(add_seq)}, current length: {len(seq)}')
        seq = seq[:seq.find(node_left_with_edge) + k - 1] + add_seq + seq[seq.find(node_left_with_edge) + k - 1:]
        logging.debug(f'Final length: {len(seq)}')

    logging.debug(f'Num edges: {sum([sum(value for value in graph[node].values()) for node in graph.keys()])}')
    return seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, required=True)
    parser.add_argument('--L', type=int, required=True)
    parser.add_argument('--p', required=True)
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    # first we grab the number of each kmer that we want
    # then, we perform short assembly on it
    pdist_map = np.load(args.p)
    input_pdist = pdist_map[list(pdist_map.keys())[0]]

    # now, we find the exact number of kmers that we should expect by multiplying it
    # by the amount that would get us the total length which is l - k + 1
    f1_norm = args.L - args.k + 1
    # assume that rounding should be okay
    pdist = np.round(input_pdist * f1_norm)
    
    if sum(pdist) > f1_norm:
        # then remove some random kmer, as it won't matter as much -- should max be 1
        while sum(pdist) > f1_norm:
            sample = np.random.choice(len(pdist), p=input_pdist)
            if pdist[sample] > 0:
                pdist[sample] -= 1

    elif sum(pdist) < f1_norm:
        # add some random kmer, as it won't matter as much -- should max be 1
        samples = np.random.choice(len(pdist), f1_norm - int(sum(pdist)), p=input_pdist)
        for sample in samples:
            pdist[sample] += 1

    # handle the special case of k being 1, where you simply just print it all out
    # one by one
    if args.k == 1:
        output_str = ''
        for i in range(4):
            output_str += (map_index_to_kmer(i, 1) * int(pdist[i]))
        print(output_str)
        logging.debug(len(output_str))
        exit()

    graph = build_graph(pdist, args.k)

    edge_counts = count_edges(graph)
    
    # the number of out edges should be n - k + 1
    logging.debug(f'{sum([edge["outdeg"] for edge in edge_counts.values()])} number of edges')

    # now sort this in terms of the number of indeg - outdeg.
    # the smallest such will be the starting point, i.e. outdegs are more than the indegrees
    nodes = list(graph.keys())
    nodes = sorted(nodes, key=lambda x: edge_counts[x]['indeg'] - edge_counts[x]['outdeg'])

    path = find_eulerian_path(graph, nodes[0], args.k)
    # print the substring of it
    print(path)

    logging.debug(len(path))
