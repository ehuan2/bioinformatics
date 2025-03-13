import argparse
import numpy as np

# This file is meant as a script to generate fasta file inputs of the specified size
if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', type=int, required=True)
    parser.add_argument('--output')
    args = parser.parse_args()

    total_nums = 4**args.k
    x = np.random.rand(total_nums)
    # Normalize the array to ensure it sums to 1
    p = x / x.sum()

    print(p)

    np.savez(f'{args.output}.npz', p)
