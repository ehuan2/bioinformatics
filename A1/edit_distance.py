# edit_distance.py
# Author: Eric Huang
# Date: January 20th 2025
import argparse
import os
from enum import Enum

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
    """Parses the input FASTA file (expecting two lines) and returns the two strings

    Args:
        input_file (str): Path to the FASTA file.
    
    Returns:
        nucl_1 (str): First parsed nucleotide string.
        nucl_2 (str): Second parsed nucleotide string.
    """
    assert os.path.isfile(input_file)

    with open(input_file, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        return lines[1], lines[3]


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


def edit_distance(nucl_1, nucl_2):
    """Calculates and returns the minimal number of edit operations to transform nucl_1 to nucl_2

    Args:
        nucl_1 (str): First parsed nucleotide string.
        nucl_2 (str): Second parsed nucleotide string.
    
    Returns:
        edit_distance: The minimal number of edit operations.
    """
    # Tackle this through dynamic programming
    # We'll create an array matrix that's n x m size
    # where n is the number of characters in nucl_1 and m is the number of characters in nucl_2

    n, m = len(nucl_1), len(nucl_2)
    dist_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # now start in the top left corner, at (0, 0) and fill out the columns and rows as needed
    for i in range(n + 1):
        dist_arr[i][0] = i
        backtrack_arr[i][0] = (i - 1, 0, SequenceStep.N1_TAKE)

    for i in range(m + 1):
        dist_arr[0][i] = i
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
            
            # now we take the minimal edit distance by simply penalizing it if we need to edit
            # n?_take_cost means the cost of choosing a gap in the opposite string and choosing this string
            n1_take_cost = dist_arr[x - 1][y] + 1
            n2_take_cost = dist_arr[x][y - 1] + 1
            # mismatch cost means the cost of matching them up
            mismatch_cost = dist_arr[x - 1][y - 1] + (0 if nucl_1[x - 1] == nucl_2[y - 1] else 1)

            min_cost = min(
                n1_take_cost,
                n2_take_cost,
                mismatch_cost
            )

            dist_arr[x][y] = min(
                n1_take_cost,
                n2_take_cost,
                mismatch_cost
            )

            if mismatch_cost == min_cost:
                # then this is the best path
                backtrack_arr[x][y] = (x - 1, y - 1, SequenceStep.MISMATCH)
            elif n1_take_cost == min_cost:
                backtrack_arr[x][y] = (x - 1, y, SequenceStep.N1_TAKE)
            elif n2_take_cost == min_cost:
                backtrack_arr[x][y] = (x, y - 1, SequenceStep.N2_TAKE)

    print_array(nucl_1, nucl_2, dist_arr)
    
    # now let's print the two sequences
    seq_1 = []
    seq_2 = []

    # the max number of steps should be n + m
    # keep track of the coordinates of where we are
    coord = (n, m)
    for i in range(n + m):
        if coord[0] == 0 and coord[1] == 0:
            break

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
            assert coord[0] == coord[1] and coord[0] == 0

    print("".join(reversed(seq_1)))
    print("".join(reversed(seq_2)))

    return dist_arr[n][m]


if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    args = parser.parse_args()

    # now parse the fasta file with two nucleotide strings
    nucl_1, nucl_2 = parse_fasta(args.input)

    # and finally output the edit distnace between the two
    print(edit_distance(nucl_1, nucl_2))
