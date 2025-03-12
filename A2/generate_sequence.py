# generate_sequence.py
# Author: Eric Huang
# Date: January 21st 2025

import argparse
import datetime
import os
import random

def gen_random_seq(seq_length, ambiguous):
    chars = ['A', 'T', 'C', 'G']
    if ambiguous:
        chars.extend(['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'])
    return "".join(random.choice(chars) for _ in range(seq_length))


# This file is meant as a script to generate fasta file inputs of the specified size
if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_sequences', type=int, required=True)
    parser.add_argument('-l', '--sequence_lengths', nargs='+', type=int)
    parser.add_argument('-s', '--same_sequence_length', type=int)
    parser.add_argument('-o', '--output_dir', type=str)
    parser.add_argument('-a', '--ambiguous', action='store_true', help='Include ambiguous characters')
    args = parser.parse_args()

    num_sequences = args.num_sequences
    
    assert args.sequence_lengths or args.same_sequence_length
    sequence_lengths = (
        args.sequence_lengths
        if args.sequence_lengths else
        [args.same_sequence_length for _ in range(args.num_sequences)]
    )

    # make sure it's well formatted
    assert num_sequences == len(sequence_lengths)

    if args.output_dir:
        # Create the output directory if it doesn't exist
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        output_file_path = os.path.join(
            args.output_dir,
            f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.fna'
        )

        with open(output_file_path, 'w') as output_file:
            for i in range(num_sequences):
                # use the current timestamp as the example file name
                output_file.write(f'>seq{i}\n')
                output_file.write(gen_random_seq(sequence_lengths[i], args.ambiguous))
                output_file.write('\n')

    else:
        for i in range(num_sequences):
            print(f'>seq{i}')
            print(gen_random_seq(sequence_lengths[i], args.ambiguous))
