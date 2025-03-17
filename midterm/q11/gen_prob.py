# gen_prob.py
# Author: Eric Huang
# Date: March 11th 2025

import argparse
import random
from Bio.Seq import Seq
import numpy as np


def map_kmer_to_index(kmer):
    """Given some k-mer, maps it to the 4^k index.

    Args:
        kmer (Seq): Given some k-mer, maps it to its proper 4^k index.
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


# This file is meant as a script to generate fasta file inputs of the specified size
if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    k = 4
    frequency = np.zeros(4**k, dtype=float)

    # let's generate some things
    mapping = {
        map_kmer_to_index('AATC'): 0.75,
        map_kmer_to_index('CTGA'): 0.25
    }

    for key, freq in mapping.items():
        frequency[key] = freq
        print(map_index_to_kmer(key, k))

    np.savez(f'{args.output}.npz', frequency)
