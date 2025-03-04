# print_seq_length.py.
# Author: Eric Huang
# Date: March 2nd 2025

"""
This one here is simply a file to test out the length of each sequence read.
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from typing import List
import argparse
import logging
import numpy as np
import os
import sys


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
            print(input_file, len(record.seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir', default='./output_dir', help='Output directory, which if set will write everything to there')
    args = parser.parse_args()

    for file in os.listdir(args.output_dir):
        get_sequences(os.path.join(args.output_dir, file))
