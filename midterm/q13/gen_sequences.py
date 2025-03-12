# gen_sequences.py
# Author: Eric Huang
# Date: March 11th 2025

import argparse
import random
from Bio.Seq import Seq
import numpy as np

def gen_random_seq(seq_length):
    chars = ['A', 'T', 'C', 'G']
    return "".join(random.choice(chars) for _ in range(seq_length))

def random_mutation(sequence, sequence_length):
    mutate_index = random.randint(0, sequence_length - 1)
    chars = set(['A', 'T', 'C', 'G'])
    chars.discard(sequence[mutate_index])
    new_nucleotide = random.choice(list(chars))
    return sequence[:mutate_index] + new_nucleotide + sequence[mutate_index + 1:]

# This file is meant as a script to generate fasta file inputs of the specified size
if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_sequences', type=int, required=True)
    parser.add_argument('-l', '--len', type=int, required=True)
    parser.add_argument('-m', '--mutations', type=int, default=1)
    parser.add_argument('-c', '--copies', type=int, default=2)
    args = parser.parse_args()

    num_sequences = args.num_sequences
    sequence_length = args.len

    counter = 0

    sequences = []
    copies_counter = []

    for i in range(num_sequences):
        seq = Seq(gen_random_seq(sequence_length))
        copies = random.randint(2, args.copies)
        counter += copies
        copies_counter.append(counter)
        for j in range(copies):
            if random.randint(1, 2) == 1:
                sequences.append(str(seq))
            else:
                sequences.append(str(seq.reverse_complement()))

    assert num_sequences >= args.mutations
    mutation_indices = np.random.choice(num_sequences, size=args.mutations, replace=False)
    # print(f'Mutation indices: {mutation_indices}, mutation indices real: {[copies_counter[index] - 1 for index in mutation_indices]}')

    for i in mutation_indices:
        # mutate one of these...
        sequences[copies_counter[i] - 1] = random_mutation(sequences[copies_counter[i] - 1], sequence_length)

    random.shuffle(sequences)

    for i in range(len(sequences)):
        print(f'>seq{i}')
        print(sequences[i])
