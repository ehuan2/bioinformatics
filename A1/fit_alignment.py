# fit_alignment.py
# Author: Eric Huang
# Date: January 20th 2025
import argparse
import os
from enum import Enum

# set this flag to true if you want to see the debug print statements
DEBUG = False


class SequenceStep(Enum):
    """Enum for declaring the possible alignment sequence steps we can take
    
    END: Signals the end of the sequence step, i.e. stop
    N1_TAKE: takes the nucleotide from the nucl_1 string
    N2_TAKE: takes the nucleotide from the nucl_2 string
    MISMATCH: matches or mismatches the string, i.e. takes from both
    """
    END = 0
    N1_TAKE = 1
    N2_TAKE = 2
    MISMATCH = 3


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
                seqs.append({
                    'key': current_key,
                    'seq': current_str
                })
                current_str = ''
            current_key = line[1:]
        else:
            current_str += line

    seqs.append({
        'key': current_key,
        'seq': current_str
    })
    return seqs


def print_array(nucl_1, nucl_2, arr):
    """Given two nucleotide strings and its edit distance 2D array, print all the rows and columns out formatted nicely

    Args:
        nucl_1 (str): The first nucleotide string
        nucl_2 (str): The second nucleotide string
        arr (2D int array): The array to print out and format nicely
    """
    print("  - ", end="")
    for char in nucl_2:
        print(f"{char} ", end="")
    print()

    for x in range(len(arr)):
        row = arr[x]
        print(f"{nucl_1[x - 1] if x != 0 else '-'} ", end="")
        for y in row:
            print(f"{y} ", end="")
        print()


def verify_alignment_score(score, n1, n2):
    """Verifies that the n1 and n2 alignments match up with the reported score.

    Args:
        score (int): The score calculated through the alignment result
        n1 (str): The first nucleotide string alignment found
        n2 (str): The second nucleotide string alignment found
    """
    assert len(n1) == len(n2)

    # calculate the score based on the alignments
    calc_score = 0
    for i in range(len(n1)):
        # make sure not both of them are gaps
        assert not (n1[i] == '-' and n2[i] == '-')
        if n1[i] == n2[i]:
            calc_score += 1
        else:
            calc_score -= 1
    assert calc_score == score
    print("Alignment Score is properly calculated!")


def verify_substring_alignment(nucl_i, n_i):
    """Verifies whether the found alignment string n_i is a valid substring
       of nucl_i and verifies it's a valid alignment.

    Args:
        nucl_i (str): The original nucleotide string
        n_i (str): The string to verify
    """
    # we can verify by getting rid of any gaps
    # and then searching for this substring
    assert n_i.replace('-', '') in nucl_i
    print("Substring is a proper alignment!")


def verify_alignment(nucl_i, n_i):
    """Verifies whether the alignment string n_i is a valid alignment of nucl_i

    Args:
        nucl_i (str): The original nucleotide string
        n_i (str): The string to verify
    """
    # we can verify by getting rid of any gaps
    # and then making sure it's equivalent to the original nucleotide string
    assert n_i.replace('-', '') == nucl_i
    print("String is a proper alignment!")


def fit_alignment(nucl_1, nucl_2):
    """Finds the best substring of nucl_1 to fit with nucl_2.

    Runs in O(vw) time because we're running a local alignment sequence and finding the best
    score and backtracking from that score.

    We can reword this problem to be the following:
        Given the alignment grid, what's the maximal alignment score we can get from alignments
        starting from the 0th column to the mth column.

    To accomplish this, we essentially run the exact same algorithm with the exception that 
    our starting point is anywhere alongside the first column, and our ending point is anywhere
    alongside the last column.
    
    To make this work, we just need to simply not penalize starting from anywhere in the first column,
    initializing everything to a 0, and marking it as an end.

    Args:
        nucl_1 (str): The string to substring. Must be of equal or greater length than nucl_2
        nucl_2 (str): The string to fit with.

    Returns:
        sequence_score: The final sequence score
        n1: The final sequence for nucl_1
        n2: The final sequence for nucl_2
    """
    assert len(nucl_1) >= len(nucl_2)

    # Tackle this through dynamic programming
    # We'll create an array matrix that's n x m size (or rather v x w size)
    # where n is the number of characters in nucl_1 and m is the number of characters in nucl_2
    # we differ in this case by not finding the best global alignment, rather we find the best alignment
    # that starts in the 0th column and ends in the 1st column

    n, m = len(nucl_1), len(nucl_2)
    dist_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # now start in the top left corner, at (0, 0) and fill out the columns and rows as needed
    for i in range(n + 1):
        dist_arr[i][0] = 0
        backtrack_arr[i][0] = (i - 1, 0, SequenceStep.END)

    for i in range(m + 1):
        dist_arr[0][i] = -i
        backtrack_arr[0][i] = (0, i - 1, SequenceStep.N2_TAKE)

    # set the base case
    backtrack_arr[0][0] = (-1, -1, SequenceStep.END)

    # now we need to iterate in the order of (1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (2, 2) where
    # the sum of the coordinates is what's coming up next.
    # we do this all the way until (n, m)
    for i in range(2, n + m + 1):
        # basically, we have i ranging from [2, n + m] inclusive
        # and we want that x + y = i, with 1 \leq x \leq n and 1 \leq y \leq m
        # so, we can set y = i - x, and we know that 1 \leq i - x \leq m
        # or, x is in the range [i - m, i - 1] and the range [1, n], or [max(i - m, 1), min(i - 1, n)]
        for x in range(max(i - m, 1), min(i, n + 1)):
            y = i - x
            
            # now we take the maximal edit distance by simply penalizing it if it's wrong
            # n?_take_cost means the cost of choosing a gap in the opposite string and choosing this string
            n1_take_cost = dist_arr[x - 1][y] - 1
            n2_take_cost = dist_arr[x][y - 1] - 1
            # mismatch cost means the cost of matching them up, rewarding +1 if it matches, -1 otherwise
            mismatch_cost = dist_arr[x - 1][y - 1] + (1 if nucl_1[x - 1] == nucl_2[y - 1] else -1)

            max_score = max(
                n1_take_cost,
                n2_take_cost,
                mismatch_cost
            )

            dist_arr[x][y] = max_score

            if mismatch_cost == max_score:
                # then this is the best path
                backtrack_arr[x][y] = (x - 1, y - 1, SequenceStep.MISMATCH)
            elif n1_take_cost == max_score:
                backtrack_arr[x][y] = (x - 1, y, SequenceStep.N1_TAKE)
            elif n2_take_cost == max_score:
                backtrack_arr[x][y] = (x, y - 1, SequenceStep.N2_TAKE)

    if DEBUG:
        print_array(nucl_1, nucl_2, dist_arr)
    
    # now let's create both sequences by following the backtrack matrix
    seq_1 = []
    seq_2 = []

    # first we need to find the highest value in the final column
    coord = None
    final_max_score = None
    for i in range(n + 1):
        if coord == None or final_max_score == None or final_max_score < dist_arr[i][m]:
            final_max_score = dist_arr[i][m]
            coord = (i, m)

    # the max number of steps should be n + m
    # keep track of the coordinates of where we are
    for i in range(n + m + 1):
        # stop once we get to the start
        if coord[1] == 0:
            break

        # for each coordinate we have, figure out its previous step and recreate the steps
        coord = backtrack_arr[coord[0]][coord[1]]
        if coord[2] == SequenceStep.N1_TAKE:
            seq_1.append(nucl_1[coord[0]])
            seq_2.append('-')
        elif coord[2] == SequenceStep.N2_TAKE:
            seq_1.append('-')
            seq_2.append(nucl_2[coord[1]])
        elif coord[2] == SequenceStep.MISMATCH:
            seq_1.append(nucl_1[coord[0]])
            seq_2.append(nucl_2[coord[1]])
        elif coord[2] == SequenceStep.END:
            # make sure we're in the first column now
            assert coord[1] == 0

    # create the final strings we need
    n1 = "".join(reversed(seq_1))
    n2 = "".join(reversed(seq_2))

    return final_max_score if final_max_score != None else 0, n1, n2


if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    DEBUG = args.debug

    # now parse the fasta file with two nucleotide strings
    seqs = parse_fasta(args.input)
    nucl_1, nucl_2 = seqs[0]['seq'], seqs[1]['seq']

    # and finally output the best fit alignment
    score, n1, n2 = fit_alignment(nucl_1, nucl_2)
    print(score)
    print(n1)
    print(n2)
    
    if DEBUG:
        verify_alignment_score(score, n1, n2)
        verify_substring_alignment(nucl_1, n1)
        verify_alignment(nucl_2, n2)

