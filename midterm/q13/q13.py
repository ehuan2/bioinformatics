# q13.py
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


# returns whether or not pairs are one nucleotide away or not
def is_single_nucleotide_apart(seq1, seq2):
    diffs = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            diffs += 1
            if diffs >= 2:
                return False

    return diffs == 1


# returns whether or not the given sequence is one nucleotide away from
# a different pair
def is_one_nucleotide_away(seq, pair_two):
    return (
        is_single_nucleotide_apart(seq, pair_two[0]) or
        is_single_nucleotide_apart(seq, pair_two[1])
    )


# finds which of the pair is one nucleotide away
def find_nucleotide_diff(seq1, pair2):
    if is_single_nucleotide_apart(seq1, pair2[0]):
        return pair2[0]
    elif is_single_nucleotide_apart(seq1, pair2[1]):
        return pair2[1]
    raise RuntimeError('Should not reach here...')


if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    # now parse the fasta file with two nucleotide strings
    seqs = parse_fasta(args.reads)
    logging.debug(seqs)

    count_mapping = {}

    # first generate a mapping from a key and its reverse complement to its count
    for seq in seqs:
        # if the count is 1, then include which original sequence it came from
        key = gen_reverse_complement_pair(seq)
        if key in count_mapping.keys():
            count_mapping[key]['count'] += 1
        else:
            count_mapping[key] = {
                'count': 1,
                'seq': set(),
                'neighbours': set()
            }
        
        # add the current sequence to the set
        count_mapping[key]['seq'].add(seq)

    logging.debug(count_mapping)

    # then, calculate the graph that can represent this
    keys = list(count_mapping.keys())
    count_one_keys = set()
    for i in range(len(keys)):
        if count_mapping[keys[i]]['count'] == 1:
            count_one_keys.add(keys[i])

        # then add the neighbours as necessary
        for j in range(len(keys)):
            if i == j:
                continue
            if any(
                is_one_nucleotide_away(seq, keys[j])
                for seq in count_mapping[keys[i]]['seq']
            ):
                count_mapping[keys[i]]['neighbours'].add(keys[j])

    logging.debug(count_mapping)
    logging.debug(count_one_keys)

    # now, we need to sort each node with count of 1 based on how we need to merge it
    # we'll sort it based on the number of neighbours
    # now for each of these one count, one neighbour's nodes we need to merge them first
    # and as we merge, we update the graph with its counts and neighbours
    while len(count_one_keys) > 0:
        count_one_list = sorted(
            list(count_one_keys),
            key=lambda key: len(count_mapping[key]['neighbours']),
            reverse=True # sort it such that the most neighbours appear first
        )
        logging.debug(f'Current keys and neighbours: {[(key, count_mapping[key]["neighbours"]) for key in count_one_list]}')

        # choose the first one and merge it to one of its neighbours
        # should always have a neighbour, or else error -- and this neighbour
        # should have only correct reads -- if we only merge to correct reads
        # then we should never force an errorneous read to be correct
        current_node = count_one_list[0]
        logging.debug(f'Current key: {current_node}, {count_mapping[current_node]["neighbours"]}')
        assert len(count_mapping[current_node]['neighbours']) > 0
        assert count_mapping[current_node]['count'] == 1

        merge_with_node = sorted(
            list(count_mapping[current_node]['neighbours']),
            key=lambda key: count_mapping[key]['count'],
            reverse=True
        )[0] # get the node with the most counts so far to merge with
        
        # now, we should merge these two nodes together, by merging
        # into the other node. We merge by updating the count_one_keys
        # as well as the neighbours
        logging.debug(f'Merging {current_node} into {merge_with_node}')
        current_seq = list(count_mapping[current_node]["seq"])[0]
        print(f'{current_seq}->{find_nucleotide_diff(current_seq, merge_with_node)}')

        count_mapping[merge_with_node]['count'] += 1
        count_mapping[current_node]['count'] = 0
        count_mapping[current_node]['neighbours'] = set()

        for neighbour in count_mapping[current_node]['neighbours']:
            count_mapping[neighbour]['neighbours'].remove(current_node)
        
        count_one_keys.remove(current_node)
        count_one_keys.discard(merge_with_node)


    # test cases: 
    # 1. the mutation is the very last one 
    # 2. the mutation is some earlier one and the very last one
    # 3. mutation adds up to some other one
    # 4. one mutation in one nucleotide leads to a change from one to another
    # and then another one changes it to a different one
    # e.g.: AG, CG, CG, AT
    # 5. Different types of graphs! We would expect any type of graph where:
    # a. Each node is either on its own with a count >= 2
    # b. Each node is connected to at least one other node

    # TODO: test with an empty case
