# edit_distance.py
# Author: Eric Huang
# Date: March 11th 2025
import argparse
import os
import logging
import sys

def parse_fasta(input_file):
    """Parses the input FASTA file and returns an array of tuples of sequence names and sequences

    Args:
        input_file (str): Path to the FASTA file.
    
    Returns:
        seqs (array of dicts): array of dicts of sequence names and sequences
    """
    assert os.path.isfile(input_file)

    with open(input_file, 'r') as file:
        lines = [line.strip() for line in file.readlines()]

    seqs = []
    current_str = ''
    current_key = None

    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if current_key != None:
                seqs.append(current_str)
                current_str = ''
            current_key = line[1:]
        else:
            current_str += line

    seqs.append(current_str)
    return seqs


def reverse_complement(seq):
    mapping = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    reverse_seq = ''
    for char in seq:
        reverse_seq += mapping[char]
    return reverse_seq[::-1]


def gen_reverse_complement_pair(seq):
    # given a nucleotide sequence, generate the tuple representing the sequence and its reverse complement
    # this tuple will be unique because we will sort it lexicographically
    seq_list = [seq, reverse_complement(seq)]
    return tuple(sorted(seq_list))


def parse_reads(input_string):
    """
    Parse the input string into input and output reads.
    
    Args:
        input_string (str): The input string containing the reads.
    
    Returns:
        tuple: A tuple containing the input read and the output read.
    """
    start_input = input_string.find('[')
    end_input = input_string.find(']->[')
    start_output = end_input + 4
    end_output = input_string.find(']', start_output)
    
    input_read = input_string[start_input+1:end_input]
    output_read = input_string[start_output:end_output]
    
    return input_read, output_read


def parse_output(output):
    assert os.path.isfile(output)

    with open(output, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    
    mutations = []    
    
    for line in lines:
        seq_1, seq_2 = parse_reads(line)
        mutations.append((seq_1, seq_2))

    return mutations

def count_diff(seq1, seq2):
    assert len(seq1) == len(seq2)
    count = 0
    for i in range(len(seq1)):
        count += 1 if seq1[i] != seq2[i] else 0
    return count


if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    # first take the input file and get their counts
    # then modify those counts accordingly to the output file
    # then, verify that all the final counts are correct
    seqs = parse_fasta(args.input)

    seq_counts = {}

    for seq in seqs:
        if seq not in seq_counts.keys():
            seq_counts[seq] = 1
        else:
            seq_counts[seq] += 1

    for start, end in parse_output(args.output):
        # start has to 100% be in the input
        assert start in seq_counts.keys() 
        # but the end might be the reverse complement instead!
        assert end in seq_counts.keys() or reverse_complement(end) in seq_counts.keys()

        # decrement the count
        seq_counts[start] -= 1
        
        # increment the end 
        if end not in seq_counts.keys():
            seq_counts[end] = 1
        else:
            seq_counts[end] += 1

        # verify that start and end are only one nucleotide apart
        assert count_diff(start, end) == 1

    # now check that the total counts are correct
    for seq, count in seq_counts.items():
        reverse = reverse_complement(seq)
        total_count = (seq_counts[reverse] if reverse in seq_counts.keys() else 0) + count
        assert total_count == 0 or total_count > 1
 
    print(f'Solution was correct!')
 