# assemble_genome.py.
# Author: Eric Huang
# Date: March 2nd 2025

"""
This is the solution to problem 4, the genome assembly problem.

As input it takes a list of sequences which needs to be assembled.
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import deque
from functools import reduce
from kmer_comp import get_kmer_frequency_array, map_kmer_to_index
from typing import List, Dict, Tuple
from tqdm import tqdm
import argparse
import logging
import numpy as np
import os
import sys

DEBUG = False


def get_sequences(input_file) -> List[SeqRecord]:
    """Parses the input file path to get the list of fasta sequences.

    Args:
        input_file (str): Path to input file.

    Returns:
        List[SeqRecord]: The list of fasta sequence records
    """
    sequences = []

    # parse the FASTA file
    with open(input_file, 'r') as handle:
        # Parse the FASTA file
        for record in SeqIO.parse(handle, 'fasta'):
            logging.debug(f'Record id: {record.id}')
            # return the first sequence we get
            sequences.append(record)

    return sequences


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


def error_correct_kmers(kmer_freqs: List[float | int], kmer_len: int, tol: int) -> List[float | int]:
    """Given a set of kmer sequences, return the error corrected list.

    Args:
        kmer_freqs (List[float | int]): Inputted list of kmer sequences.
        kmer_len (int): Length of a kmer
        tol (int): The tolerance of a kmer frequency, if it is below, it will add this somewhere else

    Returns:
        List[float | int]: Outputted error corrected list!
    """
    original_sum = sum(kmer_freqs)

    for i in tqdm(range(len(kmer_freqs))):
        if kmer_freqs[i] <= tol:
            # then we want to merge this with another kmer, the best matching one
            # which we will define as the frequency of a neighbour that has the highest one
            kmer = map_index_to_kmer(i, kmer_len)
            # now let's generate all 4^k pairs of kmers that are close to this
            possible_kmers = []
            bases = ['A', 'C', 'G', 'T']

            for j in range(len(kmer)):
                for base in bases:
                    # skip the same base
                    if kmer[j] == base:
                        continue
                    possible_kmers.append(
                        kmer[:j] + base + kmer[j+1:]
                    )

            candidates = [
                (candidate, map_kmer_to_index(candidate), kmer_freqs[map_kmer_to_index(candidate)])
                for candidate in possible_kmers
            ]

            candidates = sorted(candidates, key=lambda x: x[2], reverse=True)
            if kmer_freqs[candidates[0][1]] == 0:
                continue
            kmer_freqs[candidates[0][1]] += kmer_freqs[i]
            kmer_freqs[i] = 0

    assert original_sum == sum(kmer_freqs)
    return kmer_freqs


def bfs(graph: Dict[str, List[str]], start: str) -> set:
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


def get_components(graph: Dict[str, List[str]]) -> List[str]:
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


def check_eulerian_path(graph: Dict[str, Dict[str, int]]):
    """Checks for the Eulerian path condition and prints out all nodes that
        have a discrepancy between the indegree and outdegree.

    Args:
        graph (Dict[str, List[str]]): The given graph.
    """
    degree_counts = {
        key: {
            'indeg': 0,
            'outdeg': 0
        }
        for key in graph.keys()
    }

    for key in graph.keys():
        degree_counts[key]['outdeg'] += sum(graph[key].values())
        for neighbour, count in graph[key].items():
            degree_counts[neighbour]['indeg'] += count

    diff_counts = {
        key: value['indeg'] - value['outdeg']
        for key, value in degree_counts.items()
    }

    no_match = list(filter(lambda x: x != 0, diff_counts.values()))
    logging.debug(
        f'Min diff: {min(diff_counts.items(), key=lambda x: x[1])}, '
        f'Max diff: {max(diff_counts.items(), key=lambda x: x[1])}, '
        f'Average diff: {sum(map(lambda x: abs(x[1]), diff_counts.items())) / len(diff_counts)}'
    )
    logging.debug(f'Match count: {len(diff_counts) - len(no_match)}, Not match: {len(no_match)}')
    # logging.debug(diff_counts)
    return diff_counts


def get_debruijn_graph(kmer_freq: List[int | float], kmer_len: int) -> Dict[str, Dict[str, int]]:
    """Given a list of kmer frequencies, return the debrujin nodes who are (k - 1)-mers

    Args:
        kmer_freq (List[int | float]): kmer frequency
        kmer_len (int): Length of the kmer

    Returns:
        Dict[str, Dict[str, int]]: Mapping from each node to its neighbours, adjacency list
            which will contain the amount of edges it has per unique string
    """
    # given the kmer frequencies, build the graph that matches things up for us
    # first we need to recover the entire list of kmers
    graph = {}
    
    # dual graph will be computed both ways
    dual_graph = {}

    # for each kmer frequency we have, add it to the list
    for i in range(len(kmer_freq)):
        kmer = map_index_to_kmer(i, kmer_len)

        if kmer_freq[i] > 0:
            kmer_from, kmer_to = kmer[0:len(kmer) - 1], kmer[1: len(kmer)]

            if kmer_from not in graph.keys():
                graph[kmer_from] = {}
                dual_graph[kmer_from] = []
            if kmer_to not in graph.keys():
                graph[kmer_to] = {}
                dual_graph[kmer_to] = []

            if kmer_to not in graph[kmer_from]:
                graph[kmer_from][kmer_to] = kmer_freq[i]
            else:
                graph[kmer_from][kmer_to] += kmer_freq[i]
            
            # add edges going both ways
            dual_graph[kmer_from].append(kmer_to)
            dual_graph[kmer_to].append(kmer_from)

    # Now we want to analyze the number of nodes, edges and components
    components = get_components(dual_graph)
    if DEBUG:
        logging.debug(f'Number of nodes: {len(graph.keys())}')
        logging.debug(f'Number of edges: {sum([sum(edgelist.values()) for edgelist in graph.values()])}')
        logging.debug(f'Components and their sizes: {list(map(lambda x: len(x), components))}')
        # next we want to check whether the eulerian path condition holds -- i.e. the indegree = outdegree
    
    largest_component = sorted(components, key=lambda x: len(x), reverse=True)[0]
    logging.debug(f'Largest: {len(largest_component)}')

    # now we take the subgraph, which is only composed of the largest component
    subgraph = {
        key: graph[key]
        for key in largest_component
    }

    return subgraph

# form of edges of src to dest int number of times
Edge = Tuple[str, str, int]

def get_edges(diff_counts: Dict[str, int]) -> Tuple[List[Edge], List[Edge], List[str]]:
    """Given the difference counts of edges with different indegrees and outdegrees,
        calculates both the edges needed to add to the graph to balance.

        We can also get the other direction, where we remove the edges that affect the
        indegree and outdegree counts.

    Args:
        diff_counts (Dict[str, int]): The difference between the indegree and outdegree.
        
        Note that negative means there are more outdegrees.

    Returns:
        Tuple[List[Edge], List[Edge]]: The list of edges to add and to subtract.
    """
    edges_to_add = []
    for key in diff_counts.keys():
        if diff_counts[key] > 0:
            # now we find all other keys s.t. they are negative and possible to draw an edge to
            # i.e. key[1:len(key)] == other_key[0:len(key) - 1]
            for other_key in diff_counts.keys():
                if diff_counts[key] <= 0:
                    break

                if (key != other_key and
                        key[1:len(key)] == other_key[0:len(key) - 1] and
                        diff_counts[other_key] < 0):
                    edges_added = min(abs(diff_counts[other_key]), diff_counts[key])
                    edges_to_add.append((key, other_key, edges_added))
                    diff_counts[key] -= edges_added
                    diff_counts[other_key] += edges_added

    edges_to_remove = []
    for key in diff_counts.keys():
        if diff_counts[key] < 0:
            # now we find all other keys s.t. they are positive and possible to remove an edge to
            # i.e. key[1:len(key)] == other_key[0:len(key) - 1]
            for other_key in diff_counts.keys():
                if diff_counts[key] >= 0:
                    break

                if key != other_key and key[1:len(key)] == other_key[0:len(key) - 1] and diff_counts[other_key] > 0:
                    edges_added = min(diff_counts[other_key], abs(diff_counts[key]))
                    edges_to_remove.append((key, other_key, edges_added))
                    diff_counts[key] += edges_added
                    diff_counts[other_key] -= edges_added

    return edges_to_add, edges_to_remove

def find_eulerian_path(graph: Dict[str, List[str]]) -> str:
    """Given a graph, finds the Eulerian path.

    Args:
        graph (Dict[str, List[str]]): Graph with a set Eulerian path.

    Returns:
        str: The final constructed genome.
    """
    pass


def assemble_genome(seqs: List[SeqRecord], kmer_len: int, tol: int) -> str:
    """Given a list of sequences to assemble, return a string
    representing the final assembled genome.

    This will be tackled through short assembly.    

    Args:
        seqs (List[SeqRecord]): List of sequences to assemble
        k (int): The length of our kmers
        tol (int): The tolerance of the kmer frequency

    Returns:
        str: The assembled genome
    """
    # first grab all possible k-mers from the sequence of records
    kmer_freqs = [np.array(get_kmer_frequency_array(seq.seq, kmer_len)) for seq in seqs]
    # sum up across the rows
    total_kmer_freq = np.sum(kmer_freqs, axis=0)

    # Then we do some error correction on these k-mers
    total_kmer_freq = error_correct_kmers(total_kmer_freq, kmer_len, tol)

    # Then we want to build the graph based on the short read assembly
    graph = get_debruijn_graph(total_kmer_freq, kmer_len)

    # now, we try to make the Eulerian condition hold -- we'll try to connect all the nodes until
    # just one pair have a missing indegree and a missing outdegree
    diff_counts = check_eulerian_path(graph)

    # using these diff counts, we'll try to match outdegrees and indegrees, and create edges to add
    edges_to_add, edges_to_remove = get_edges(diff_counts)
    
    # next, we add these edges to the graph
    for src, dest, count in edges_to_add:
        if dest in graph[src]:
            graph[src][dest] += count
        else:
            graph[src][dest] = count

    # logging.debug(graph)
    _ = check_eulerian_path(graph)

    # and we remove them too
    for src, dest, count in edges_to_remove:
        if dest in graph[src]:
            graph[src][dest] -= count
            assert graph[src][dest] >= 0
    
    diff_counts = check_eulerian_path(graph)
    
    # now we take the node with the most out degrees as the starting path
    # continue going for as long as possible
    exit()
    # Finally, we want to calculate the Eulerian path!
    return find_eulerian_path(graph, start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='The input file with sequences to assemble')
    parser.add_argument('-o', '--output_dir', help='Output directory, which if set will write everything to there')
    parser.add_argument('-k', '--kmer_len', default=4, type=int, help='The kmer length')
    parser.add_argument('-t', '--tol', default=1, type=int, help='The kmer frequency tolerance')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
        DEBUG = True

    assembled_genome = assemble_genome(get_sequences(args.input), args.kmer_len, args.tol)

    # next handle the output directories
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    output_file_path = (
        os.path.join(args.output_dir, f'{os.path.basename(args.input)}_assembled_genome.fna')
        if args.output_dir else
        'assembled_genome.fna'
    )

    with open(output_file_path, 'w') as output_file:
        output_file.write('>seq0\n')
        output_file.write(str(assembled_genome))
