import argparse
from collections import deque
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
    dual_graph = {}
    for i in range(len(pdist)):
        if pdist[i] == 0:
            continue
        kmer = map_index_to_kmer(i, kmer_len)
        kmer_from, kmer_to = kmer[0:kmer_len - 1], kmer[1: kmer_len]

        if kmer_from not in graph.keys():
            graph[kmer_from] = {}
            dual_graph[kmer_from] = set()
        
        if kmer_to not in graph.keys():
            graph[kmer_to] = {}
            dual_graph[kmer_to] = set()
        
        if kmer_to not in graph[kmer_from]:
            graph[kmer_from][kmer_to] = 0

        graph[kmer_from][kmer_to] += pdist[i]
        
        # add edges going both ways
        dual_graph[kmer_from].add(kmer_to)
        dual_graph[kmer_to].add(kmer_from)

    return graph, dual_graph


def bfs(graph, start):
    """Perform bfs as a helper function to the number of components.

    Args:
        graph (Dict[str, List[str]]): The given graph
        start (str): The starting node.

    Returns:
        set: The set of visited graphs seen
    """
    if start not in graph:
        return []
    
    visited = set()
    queue = deque([start])
    visited.add(start)
    
    while len(queue) != 0:
        node = queue.popleft()
        
        for neighbor in graph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    return visited


def get_components(graph):
    """Given a graph, computes the components.

    Args:
        graph (Dict[str, List[str]]): The deBruijn graph.

    Returns:
        List[str]: The components as the list of keys.
    """
    comp_sizes = []
    nodes_visited = set()

    for key in graph.keys():
        if key in nodes_visited:
            # do the breadth first search with nodes_visited
            continue

        visited_from_key = bfs(graph, key)
        for visited_node in visited_from_key:
            nodes_visited.add(visited_node)
        nodes_visited.add(key)
        comp_sizes.append(visited_from_key)
    
    return comp_sizes


def build_components(graph, dual_graph):
    components = get_components(dual_graph)
    
    graph_subsets = []
    for component in components:
        graph_subsets.append({
            key: graph[key] for key in component
        })

    return graph_subsets


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
def dfs(graph, start, stop_at_self):
    to_add = ''
    node_visit = start

    # do a dfs to find the path starting at start
    while True:
        next_nodes = []
        # find all nodes with a non-zero value, and sample them
        # randomly
        for node, value in graph[node_visit].items():
            if value > 0:
                next_nodes.append((node, value))

        if len(next_nodes) == 0:
            break
        
        total = sum(value for _, value in next_nodes)
        next_node = np.random.choice(
            list([node for node, _ in next_nodes]),
            p=[value / total for _, value in next_nodes]
        )

        assert graph[node_visit][next_node] > 0
        graph[node_visit][next_node] -= 1
        to_add += next_node[-1]
        node_visit = next_node

        if stop_at_self and (next_node == start):
            break

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
    to_add, graph = dfs(graph, start, False)
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
            if node not in graph:
                continue

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
        add_seq, graph = dfs(graph, node_left_with_edge, True)
        logging.debug(f'Current length to add: {len(add_seq)}, current length: {len(seq)}')
        seq = seq[:seq.find(node_left_with_edge) + k - 1] + add_seq + seq[seq.find(node_left_with_edge) + k - 1:]
        logging.debug(f'Final length: {len(seq)}')

    logging.debug(f'Num edges: {sum([sum(value for value in graph[node].values()) for node in graph.keys()])}')
    return seq


def force_eulerian_graph(graph, edge_counts):
    # forces a eulerian graph by adding the insufficient outdegrees to the insufficient indegrees
    # first, we need to sort the edges into the insufficient outdegrees and the insufficient indegrees
    insufficient_indegs = list(filter(
        lambda node: edge_counts[node]['indeg'] < edge_counts[node]['outdeg'], list(graph.keys())
    ))
    insufficient_outdegs = list(filter(
        lambda node: edge_counts[node]['indeg'] > edge_counts[node]['outdeg'], list(graph.keys())
    ))

    if len(insufficient_indegs) == 0 or len(insufficient_outdegs) == 0:
        return graph

    def f1_norm(values):
        if all(value == 0 for value in values):
            return np.array(values)
        return values / np.sum(values)

    # now, for each of the insufficient outdegrees, we will create a probability distribution
    # over the indegrees and add those edge counts to the graph
    for node in insufficient_outdegs:
        # first generate the counts from the node to the indegrees
        out_to_in = []
        for indeg in insufficient_indegs:
            count = 0.0 if indeg not in graph[node] else graph[node][indeg]
            out_to_in.append(count)

        out_to_in = f1_norm(out_to_in)
        out_to_in *= np.abs(edge_counts[node]['indeg'] - edge_counts[node]['outdeg'])
        out_to_in = np.round(out_to_in)

        for indeg, add_edge_count in zip(insufficient_indegs, out_to_in):
            if indeg not in graph[node]:
                graph[node][indeg] = add_edge_count
            else:
                graph[node][indeg] += add_edge_count

    return graph


def score_components(components):
    components_edges_count = [
        sum(
            count
            for neighbour_dict in component.values()
            for count in neighbour_dict.values()
        )
        for component in components
    ]

    print(sorted(zip(components_edges_count, components), key=lambda x: x[0]))


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

    graph, dual_graph = build_graph(pdist, args.k)

    # TODO: Calculate the eulerian path in these different components, then merge them
    # We can merge them by simply adding them or by adding some ratio of them
    # depending on how much they matter according to the probability mass that
    # belongs to them
    # we can merge them by adding them straight up
    components = build_components(graph, dual_graph)

    # get a probability distribution of the components that matter most
    components_pdist = score_components(components)

    component_paths = []
    for component in components:
        edge_counts = count_edges(component)
        
        # the number of out edges should be n - k + 1
        logging.debug(f'{sum([edge["outdeg"] for edge in edge_counts.values()])} number of edges')

        # now sort this in terms of the number of indeg - outdeg.
        # the smallest such will be the starting point, i.e. outdegs are more than the indegrees
        nodes = filter(
            lambda key: (edge_counts[key]['indeg'] + edge_counts[key]['outdeg']) % 2 == 1,
            list(component.keys())
        )
        nodes = list(component.keys())
        nodes = sorted(nodes, key=lambda x: edge_counts[x]['indeg'] - edge_counts[x]['outdeg'])

        # now, we also modify the graph such that we make it a perfect Eulerian cycle
        # by adding extra edges between a lack of outdegrees to the lack of indegrees
        # according to the same probability distribution as the kmers we would expect
        # graph = force_eulerian_graph(graph, edge_counts)

        latest_path = find_eulerian_path(
            component,
            nodes[0] if len(nodes) > 0 else list(component.keys())[0],
            args.k
        )

        component_paths.append(latest_path)
        logging.debug(len(latest_path))

    # then based on this probability distribution, add on components path
    # until we hit our required L
