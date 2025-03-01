# build_deBruijn.py.
# Author: Eric Huang
# Date: February 25th 2025

"""
This is the solution to problem 2, the deBruijn graph building question.

As input it takes a FASTA file which represents the (k + 1)-mer which represents
the nodes, and then the ...
"""

from Bio import SeqIO
from Bio.Seq import Seq
from typing import List
import argparse
import logging
import os
import sys


def get_sequences(input_file) -> List[str]:
    """Parses the input file path to get the list of fasta sequences.

    Args:
        input_file (str): Path to input file.

    Returns:
        List[str]: The list of fasta sequences
    """
    sequences = []

    # parse the FASTA file
    with open(input_file, 'r') as handle:
        # Parse the FASTA file
        for record in SeqIO.parse(handle, 'fasta'):
            logging.debug(f'Record id: {record.id}, sequence: {record.seq}')
            # return the first sequence we get
            sequences.append(record.seq)

    return sequences

'''
def get_reverse_complements(seqs: List[str]) -> List[str]:
    """Given a list of sequences, get the list of reverse complements

    Args:
        seqs (List[str]): List of sequences

    Returns:
        List[str]: list of reverse complements
    """
    complement_mapping = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
    }
    
    reverse_complements = []

    for seq in seqs:
        # take the reverse, then find its complement by mapping the nucleotides
        # to their complement
        next_seq = ""
        for base in seq[::-1]:
            next_seq += complement_mapping[base]
        logging.debug(f'Reverse complement of {seq} is {next_seq}')
        reverse_complements.append(next_seq)
    return reverse_complements
'''
def get_reverse_complements(seqs: List[Seq]) -> List[Seq]:
    """Given a list of sequences, get the list of reverse complements

    Args:
        seqs (List[Seq]): List of sequences

    Returns:
        List[Seq]: list of reverse complements
    """
    return [seq.reverse_complement() for seq in seqs]


def get_debruijn_edges(seqs: List[Seq]) -> List[tuple]:
    """Given a list of (k+1)-mers that represent edges, return the adjacency list
    of the nodes which are the corresponding k-mers.

    Args:
        seqs (List[Seq]): List of strings.

    Returns:
        List[tuple]: List of output tuples representing the adjacency list, in sorted order.
    """
    edgeset = set()
    for seq in seqs:
        edgeset.add((str(seq[0:len(seq) - 1]), str(seq[1: len(seq)])))
    return sorted(list(edgeset))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='Input FASTA file path')
    parser.add_argument('-o', '--output_dir', help='Output directory, which if set will write to <input_filename>_k-mers.txt')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    input_file = args.input_file

    logging.debug(f'Input file: {input_file}')

    seqs = get_sequences(input_file)
    
    if args.debug:
        # the reverse complements of the reverse complements should be the exact same
        # as the original seqs given
        assert sorted(seqs) == sorted(get_reverse_complements(get_reverse_complements(seqs)))
        logging.debug(f'Reverse complements holds!')

    # now add its reverse complements too
    reverse_complements = get_reverse_complements(seqs)
    logging.debug(f'Reverse complements of {seqs} is {reverse_complements}')
    seqs.extend(reverse_complements)

    # then for this set of sequences, let's grab all of its edges
    debruijn_edges = get_debruijn_edges(seqs)

    # handle the output directory location
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    output_file_path = (
        os.path.join(args.output_dir, f'{os.path.basename(input_file)}_deBruijn.txt')
        if args.output_dir else
        'deBruijn.txt'
    )

    # write out to the appropriate file
    with open(output_file_path, 'w') as output_file:
        for edge in debruijn_edges:
            output_file.write(f'({edge[0]}, {edge[1]})\n')
            
            # good check when debug
            if args.debug:
                assert len(edge[0]) == len(edge[1])
                assert len(edge[0]) == len(seqs[0]) - 1
