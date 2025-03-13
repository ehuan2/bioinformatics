# verify_q11.py
# Author: Eric Huang
# Date: March 11th 2025

"""
This script verifies that q11 output is right by:
1. taking the output and breaking it down to its kmers
2. verifies its length
3. verifies that breaking down to its kmers will match the input npz file
"""

import argparse
from Bio import SeqIO
import numpy as np

def get_sequence(input_file):
    """Parses the input file path to get the fasta sequence.

    Args:
        input_file (str): Path to input file.

    Returns:
        Seq: The fasta sequence
    """
    # parse the FASTA file
    with open(input_file, 'r') as handle:
        # Parse the FASTA file
        for record in SeqIO.parse(handle, 'fasta'):
            # return the first sequence we get
            return record.seq


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

if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--len', type=int, required=True)
    parser.add_argument('-k', type=int, required=True)
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-s', '--input_seq')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    
    # first parse the output fna sequence
    seq = get_sequence(args.output)
    print(len(seq))
    assert len(seq) == args.len

    # grab all k-mers of length k and update a frequency array
    frequency = np.zeros(4**args.k, dtype=float)
    for i in range(args.len - args.k + 1):
        kmer = seq[i:i + args.k]
        kmer_index = map_kmer_to_index(kmer)
        frequency[kmer_index] += 1
 
    if args.input_seq:
        input_seq = get_sequence(args.input_seq)
        assert len(seq) == len(input_seq)

        input_seq_freq = np.zeros(4**args.k, dtype=float)
        for i in range(args.len - args.k + 1):
            kmer = seq[i:i + args.k]
            kmer_index = map_kmer_to_index(kmer)
            
            input_seq_freq[
                map_kmer_to_index(input_seq[i:i + args.k])
            ] += 1

        if not np.all(frequency == input_seq_freq):
            print(frequency, input_seq_freq)
            raise RuntimeError(f'Did not pass test, frequencies do not match')
       
    for i in range(4**args.k):
        frequency[i] /= (args.len - args.k + 1)

    input_pdist = np.load(args.input_file)
    input_pdist = input_pdist[list(input_pdist.keys())[0]]

    # now let's compare this output frequency with the input frequency
    close_array = np.isclose(frequency, input_pdist)
    for i in range(len(close_array)):
        if not close_array[i]:
            print(frequency[i], input_pdist[i])

    assert all(close_array)
    print(f'Output {frequency} is good for the input {input_pdist}!')
