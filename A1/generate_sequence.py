# generate_sequence.py
# Author: Eric Huang
# Date: January 21st 2025

import argparse
import random

def gen_random_seq(seq_length):
    chars = ['A', 'T', 'C', 'G']
    return "".join(random.choice(chars) for _ in range(seq_length))

# This file is meant as a script to generate fasta file inputs of the specified size
if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_sequences', type=int, required=True)
    parser.add_argument('-l', '--sequence_lengths', nargs='+', type=int, required=True)
    args = parser.parse_args()

    num_sequences = args.num_sequences
    sequence_lengths = args.sequence_lengths

    # make sure it's well formatted
    assert num_sequences == len(sequence_lengths)

    for i in range(num_sequences):
        print(f'>seq{i}')
        print(gen_random_seq(sequence_lengths[i]))
