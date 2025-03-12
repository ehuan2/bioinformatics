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
    seq = start
    node_visit = start

    # do a dfs to find the path starting at start
    while any(value > 0 for value in graph[node_visit].values()):
        # sample from a probability distribution such that the values represent the probability mass
        total = sum(graph[node_visit].values())
        next_node = np.random.choice(
            list(graph[node_visit].keys()),
            p=[value / total for value in graph[node_visit].values()]
        )
        graph[node_visit][next_node] -= 1
        seq += next_node[-1]
        node_visit = next_node

    return seq, graph


def find_eulerian_path(graph, start, k):
    """Given a graph, and the starting point, finds the Eulerian path.

    Args:
        graph (Dict[str, Dict[str, int]]): Graph with a set Eulerian path.
        start: The starting kmer.

    Returns:
        str: The final constructed genome.
    """
    # first grab the possible sequence
    seq, graph = dfs(graph, start)

    # now we modify the start_seq by finding some node within the start sequence
    # that has an edge left
    node_left_with_edge = None

    while node_left_with_edge is not None:
        node_left_with_edge = None
        for node in [seq[i:i + k] for i in range(len(seq) - k + 1)]:
            for value in graph[node].values():
                if value > 0:
                    node_left_with_edge = node
                    break
            if node_left_with_edge is not None:
                break

        # now we dfs on the node
        add_seq, graph = dfs(graph, node_left_with_edge)
        seq = seq[:seq.find(add_seq)] + add_seq + seq[seq.find(add_seq) + k:]

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
    pdist = pdist_map[list(pdist_map.keys())[0]]

    # now, we find the exact number of kmers that we should expect by multiplying it
    # by the amount that would get us the total length which is l - k + 1
    f1_norm = args.L - args.k + 1
    pdist *= f1_norm

    graph = build_graph(pdist, args.k)

    edge_counts = count_edges(graph)

    # now sort this in terms of the number of indeg - outdeg.
    # the smallest such will be the starting point, i.e. outdegs are more than the indegrees
    nodes = list(graph.keys())
    nodes = sorted(nodes, key=lambda x: edge_counts[x]['indeg'] - edge_counts[x]['outdeg'])

    print(find_eulerian_path(graph, nodes[0], args.k))
