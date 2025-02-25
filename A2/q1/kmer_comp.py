# kmer_comp.py.
# Author: Eric Huang
# Date: February 25th 2025

"""
This is the solution to problem 1, the k-mer composition question.

As input it takes a FASTA file alongside an integer k and outputs
a text file k-mers.txt with the k-mer frequency array sorted in lexicographic
order.

i.e. we have our frequency array F_s defined as a 4^k length array
where the ith element shows the number of times that the ith k-mer appears

e.g.:
ACGCGGCTCTGAAA
k = 2,
2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0

is the output since AA has 2, AC has 1, etc.
"""

from Bio import SeqIO
from typing import List
import argparse
import logging
import math
import os
import sys


def get_sequence(input_file) -> str:
    """Parses the input file path to get the fasta sequence.

    Args:
        input_file (str): Path to input file.

    Returns:
        str: The fasta sequence
    """
    # parse the FASTA file
    with open(input_file, 'r') as handle:
        # Parse the FASTA file
        for record in SeqIO.parse(handle, 'fasta'):
            logging.debug(f'Record id: {record.id}, sequence: {record.seq}')
            # return the first sequence we get
            return record.seq


def map_kmer_to_index(kmer):
    """Given some k-mer, maps it to the 4^k index.

    Args:
        kmer (str): Given some k-mer, maps it to its proper 4^k index.
    """
    
    # assuming we have AAAA = 0, etc., we can define this as:
    # A = 0, C = 1, G = 2, T = 3 and we have that each position represents
    # the power to multiply from, so if we have XYZ then this would result in
    # X * 4^2 + Y * 4^1 + Z * 4^0, etc.
    # e.g.: CGT = 1 * 4^2 + 2 * 4^1 + 3 * 4^0 = 27

    mapping = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

    # reverse the kmer which will help with finding the index more easily
    kmer = kmer[::-1]
    
    index = 0
    for i in range(len(kmer)):
        index += mapping[kmer[i]] * (4 ** i)

    return index


def gen_all_kmers(kmer) -> List[str]:
    """Given a kmer string, produce a list of kmers that remove its ambiguity.

    Args:
        kmer (str): The kmer string to map and produce all kmers

    Returns:
        List[str]: Reduces the ambiguous degenerative pairs to their proper pairing.
    """
   
    # this approach assumes that every single base is evenly likely given
    # this is taken from this wikipedia page on IUPAC notation:
    # https://en.wikipedia.org/wiki/Nucleic_acid_notation
    kmer_mapping = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'U': ['T'], # treat uracil exactly the same as thymine
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

    kmers = ['']

    for i in range(len(kmer)):
        kmers = [
            kmer_prefix + next_char
            for kmer_prefix in kmers
            for next_char in kmer_mapping[kmer[i]]
        ]

    return kmers


def get_kmer_frequency_array(seq, k) -> List[float | int]:
    """Given some DNA sequence, returns the kmer frequency array.

    Args:
        seq (str): The given DNA sequence seq.
        k (int): The kmer length.

    Returns:
        List[float] | List[int]: A list of frequencies that are either partial weights
            or if there are no ambiguities, all ints
    """
    # populate the frequencies of size 4^k
    kmer_frequencies = [0 for _ in range(4 ** k)]

    # we should iterate n - k + 1 times, from [0, n - k] as our start range
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]

        # then, transform each kmer to a list of kmers that could represent it
        # after which, we will evenly split the weighting of the kmers' frequency
        # across all possible kmers
        kmers_possible = gen_all_kmers(kmer)
        logging.debug(f'Found all possible kmers : {kmers_possible}')

        for kmer_possible in kmers_possible:
            kmer_index = map_kmer_to_index(kmer_possible)

            # weighting to be int if there is no ambiguity
            weight = 1 / len(kmers_possible) if len(kmers_possible) != 1 else 1

            kmer_frequencies[kmer_index] += weight
            logging.debug(
                f'{kmer_possible}: Its index at weight {weight}, {kmer_index}'
            )

    return kmer_frequencies


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='Input FASTA file path')
    parser.add_argument('-o', '--output_dir', help='Output directory, which if set will write to <input_filename>_k-mers.txt')
    parser.add_argument('-k', type=int, required=True, help='k-mer value')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    input_file = args.input_file
    k = args.k

    logging.debug(f'Input file: {input_file}')
    logging.debug(f'K-mer length: {k}')

    seq = get_sequence(input_file)
    kmer_frequency = get_kmer_frequency_array(seq, k)

    assert k <= len(seq), 'Cannot have a k-mer that is larger than the sequence'

    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    output_file_path = (
        os.path.join(args.output_dir, f'{os.path.basename(input_file)}_len_{k}_k-mers.txt')
        if args.output_dir else
        'k-mers.txt'
    )

    with open(output_file_path, 'w') as output_file:
        output_file.write(' '.join(map(str, kmer_frequency)))

    # if there is debug set, ensure that the total frequency is what you would expect it to be
    if args.debug:
        total_freq = 0
        for freq in kmer_frequency:
            total_freq += freq
        assert math.isclose(total_freq, len(seq) - k + 1) # the total frequency should be n - k + 1
        logging.debug('Total frequency matches properly!')
