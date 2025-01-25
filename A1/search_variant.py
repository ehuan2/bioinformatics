# search_variant.py
# Author: Eric Huang
# Date: January 21st 2025
import argparse
import os
from enum import Enum
import time
from functools import wraps


DEBUG = False
VERBOSE = False


def timing_wrapper(func):
    """Helper function to build an annotation to print out timing

    Args:
        func (func): Function to time

    Returns: Wrapper that times it for us
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time

        if DEBUG or VERBOSE:
            print(f"Function '{func.__name__}' finished in {execution_time:.2f} seconds")
        return result
    return wrapper


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

    seqs = {}
    current_str = ''
    current_key = None

    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if current_key != None:
                seqs[current_key] = current_str.upper()
                current_str = ''
            current_key = line[1:]
        else:
            current_str += line

    seqs[current_key] = current_str.upper()
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


def generate_kmers(seq, k=10, gen_precise=True):
    """Given a nucleotide sequence, generate all k-mers

    Args:
        seq (str): String to generate k-mers from
        k (int): Size of k-mer string
        gen_precise (bool): Flag set false if we want to generate the variants to match against
    
    Returns:
        Set of all possible kmers and their relative positioning within the sequence
    """
    # the number of k-mers = number of nucleotides from first to last
    # = (l - k) - 0 + 1 = l - k + 1
    all_kmers = []
    
    for i in range(len(seq) - k + 1):
        # we'll perform a naive approach and add all possible kmers, instead of pattern matching later
        base_sequence = seq[i:i+k]
        if not gen_precise:
            for j in range(k):
                for char in ['A', 'T', 'C', 'G']:
                    all_kmers.append(
                        (base_sequence[:j] + char + base_sequence[j + 1:], i)
                    )
        else:
            all_kmers.append((base_sequence, i))

    return all_kmers


@timing_wrapper
def build_inverse_index(db, k=10, gen_precise=True):
    """Given a database of sequences, index its kmers

    Args:
        db (array of dicts): Array of dictionaries
        k (int): Size of k-mer
        gen_precise (bool): boolean passed to generate_kmers for whether we want to have more matches
    
    Returns:
        inverse_index (dictionary of string to list of tuples of strings and positions):
            Maps from kmers to a list of string ids
    """
    # first build the kmers
    inverse_index = {}

    for seq_id, seq in db.items():
        kmers = generate_kmers(seq, k, gen_precise)
        for kmer, pos in kmers:
            if kmer not in inverse_index.keys():
                inverse_index[kmer] = []
            inverse_index[kmer].append((seq_id, pos))
    
    return inverse_index


@timing_wrapper
def get_hits_from_index(inverse_index, query_kmers):
    """Return the hits of our query kmers by giving the keys from the inverse index

    Args:
        inverse_index (dict of string to tuple of string and int):
            Maps kmers to list of sequences IDs and positions
        query_kmers (array of tuple of strings and ids):
            The query kmers to search over of their string and position ids

    Return:
        seq_ids (list of 4-tuples):
            All the resulting hits, as 4-tuple, their database sequence id,
            their positioning in the database sequence,
            the positioning in the query itself,
            and finally the actual k-mer itself!
    """
    hits = []

    for query in query_kmers: 
        seq, query_pos = query

        if seq not in inverse_index.keys():
            continue
        
        for seq_id, seq_pos in inverse_index[seq]:
            hits.append((seq_id, seq_pos, query_pos, seq))

    return hits


@timing_wrapper
def get_top_sequences(query_hits, k=10):
    """From the query hits, get the top-k sequences based on frequency

    Args:
        query_hits (set of 4-tuples):
            List of tuples representing the sequence id and its position,
            the query position and the actual match
        k (int, optional): Top-k. Defaults to 10.

    Returns:
        top_query_hits: The top-k query hits based on frequency. It'll be a list of
        3-tuple, where we have the seq_id, the number of hits,
        and a list of 3-tuple of the location of hits, the query positioning and the matching sequence
    """
    # for each query hit, let's create an aggregate statistic
    seq_id_mapping = {}
    
    for seq_id, pos, query_pos, seq in query_hits:
        if seq_id not in seq_id_mapping.keys():
            seq_id_mapping[seq_id] = []
        seq_id_mapping[seq_id].append((pos, query_pos, seq))

    sequence_list = []
    for seq_id in seq_id_mapping.keys():
        sequence_list.append((seq_id, len(seq_id_mapping[seq_id]), seq_id_mapping[seq_id]))

    sequence_list.sort(key=lambda x : x[1], reverse=True)
    return sequence_list[:k]

@timing_wrapper
def get_local_alignment(nucl_1, nucl_2, blast_scoring=False):
    """Returns the best local alignment and its score between two nucleotide sequences

    Args:
        nucl_1 (str): First nucleotide sequence
        nucl_2 (str): Second nucleotide sequence
        blast_scoring (bool, default to False): Flag to use blast n scoring or not
    
    Return:
        score: The final alignment score
        n1: The alignment of nucl_1
        n2: The alignment of nucl_2

    Note: For local alignments, we return a subset of both of them
    """
    # Tackle this through dynamic programming
    # We'll create an array matrix that's n x m size (or rather v x w size)
    # where n is the number of characters in nucl_1 and m is the number of characters in nucl_2
    # we differ in this case by not finding the best global alignment, rather we find the best alignment
    # with the best score

    n, m = len(nucl_1), len(nucl_2)
    dist_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack_arr = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # now start in the top left corner, at (0, 0) and fill out the columns and rows as needed
    for i in range(n + 1):
        dist_arr[i][0] = 0
        backtrack_arr[i][0] = (-1, -1, SequenceStep.END)

    for i in range(m + 1):
        dist_arr[0][i] = 0
        backtrack_arr[0][i] = (-1, -1, SequenceStep.END)

    # keep track of what the final score will be
    coord = None
    final_max_score = None

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
            n1_take_cost = dist_arr[x - 1][y] + (
                # penalize the gap extension less than the gap opening
                (-2 if backtrack_arr[x - 1][y][2] == SequenceStep.N1_TAKE else -5)
                if blast_scoring else
                -1
            )
            n2_take_cost = dist_arr[x][y - 1] + (
                # penalize the gap extension less than the gap opening
                (-2 if backtrack_arr[x][y - 1][2] == SequenceStep.N1_TAKE else -5)
                if blast_scoring else
                -1
            )
            # mismatch cost means the cost of matching them up, rewarding +1 if it matches, -1 otherwise
            mismatch_cost = dist_arr[x - 1][y - 1] + (
                (2 if blast_scoring else 1)
                if nucl_1[x - 1] == nucl_2[y - 1] else
                (-3 if blast_scoring else -1)
            )

            max_score = max(
                n1_take_cost,
                n2_take_cost,
                mismatch_cost,
                0
            )

            if coord == None or final_max_score == None or final_max_score < max_score:
                final_max_score = max_score
                coord = (x, y)

            dist_arr[x][y] = max_score

            if mismatch_cost == max_score:
                # then this is the best path
                backtrack_arr[x][y] = (x - 1, y - 1, SequenceStep.MISMATCH)
            elif n1_take_cost == max_score:
                backtrack_arr[x][y] = (x - 1, y, SequenceStep.N1_TAKE)
            elif n2_take_cost == max_score:
                backtrack_arr[x][y] = (x, y - 1, SequenceStep.N2_TAKE)
            elif max_score == 0:
                backtrack_arr[x][y] = (-1, -1, SequenceStep.END)

    if VERBOSE:
        print_array(nucl_1, nucl_2, dist_arr)
    
    # now let's create both sequences by following the backtrack matrix
    seq_1 = []
    seq_2 = []

    # the max number of steps should be n + m
    # keep track of the coordinates of where we are
    for i in range(n + m + 1):
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
            # stop once we reach the end
            break

    # create the final strings we need
    n1 = "".join(reversed(seq_1))
    n2 = "".join(reversed(seq_2))

    return final_max_score if final_max_score != None else 0, n1, n2


@timing_wrapper
def get_alignments(database_seq, query_seq, query_hits, blast_scoring=False):
    """From a database sequence and list of ids that returns the list of alignments and their scores

    Args:
        database_seq (dictionary): Mapping of seq IDs to sequences
        query_seq (str): Query sequence
        query_hits (array of strings): Array of seq IDs
        blast_scoring (bool, defaults to False): Flag to use blast scoring or not

    Returns:
        alignments (list of tuples of score, str ID, string, string):
            1. Score
            2. String ID database match
            3. Alignment of the database
            4. Alignment of the query
    """
    final_scores = []

    for hit in query_hits:
        score, n1, n2 = get_local_alignment(database_seq[hit], query_seq, blast_scoring)
        final_scores.append((score, hit, n1, n2))

    final_scores.sort(key=lambda x: x[0], reverse=True)
    return final_scores


@timing_wrapper
def create_ranges(pos_tuple_list, k=10, k_multiplier=2):
    """Given a list of 3-tuples, merge k-mer matches such that we naturally extend
    to find the range of k-mers that will be the highest scoring pair

    Args:
        pos_tuple_list (list of 3-tuples):
            List of 3-tuples including the db position, the query position and the matching sequence
        k (int): The size of the k-mer
        k_multiplier (int): The range to look at for merging. Default is 2.

    Returns:
        List of 3-tuple, first is a 2-tuple of the range of the db position,
        second is the 2-tuple of the range of the query position,
        and the 
    """
    
    # first, we'll create ranges for every single 3-tuple
    # second, we'll merge these ranges if and only if they are right next to each other
        # i.e. we merge them one by one

    range_list = []
    for db_pos, query_pos, _ in pos_tuple_list:
        # we'll make the ranges [) for easier matching
        # note that both ranges will be with where the k-mer starts, not where it ends for the latter
        range_list.append(((db_pos, db_pos + 1), (query_pos, query_pos + 1)))

    # now for each range list we'll sort it based on the db_pos
    range_list.sort(key=lambda x: x[0][0])

    # now we'll go through constructing the new ranges
    # we'll do this by setting a lower pointer, and iterating on it while we can still do it
    # once we can no longer iterate, we add it to the final range list
    final_range_list = []

    lower_range = range_list[0]

    # we'll only merge if we can do a two-hit seeding
    def can_merge(r1, r2):
        r1_db_pos_tuple, r1_query_pos = r1
        r2_db_pos_tuple, r2_query_pos = r2
        return (
            r1_db_pos_tuple[1] + k_multiplier * k >= r2_db_pos_tuple[0] and
            r1_query_pos[1] + k_multiplier * k >= r2_query_pos[0]
        )

    def merge_range(r1, r2):
        r1_db_pos_tuple, r1_query_pos = r1
        r2_db_pos_tuple, r2_query_pos = r2
        return (
            (r1_db_pos_tuple[0], r2_db_pos_tuple[1]),
            (r1_query_pos[0], r2_query_pos[1])
        )

    for next_range in range_list[1:]:
        if can_merge(lower_range, next_range):
            lower_range = merge_range(lower_range, next_range)
        else:
            final_range_list.append(lower_range)
            lower_range = next_range

    return final_range_list


@timing_wrapper
def get_seq_score(n1, n2, blast_scoring=False):
    """Given two sequences, n1, n2, return the final score of the two

    Args:
        n1 (str): The first string
        n2 (str): The second string
        blast_scoring (bool, optional):
            Boolean for whether or not we should use blast scoring.
            Defaults to False.

    Returns:
        score:
            The final score of the alignment of the two.
            Assumes that they are the same size.
    """
    score = 0
    for i in range(min(len(n1), len(n2))):
        if n1[i] == '-':
            score += (
                (-2 if i - 1 >= 0 and n1[i - 1] == '-' else -5)
                if blast_scoring else
                -1
            )
        elif n2[i] == '-':
            score += (
                (-2 if i - 1 >= 0 and n2[i - 1] == '-' else -5)
                if blast_scoring else
                -1
            )
        elif n1[i] == n2[i]:
            score += 2 if blast_scoring else 1
        elif n1[i] != n2[i]:
            score += -3 if blast_scoring else -1

    # add gap penalities if need be at the very end
    if len(n1) != len(n2):
        score += -5 + (-2) * abs(len(n1) - len(n2) - 1)

    return score


@timing_wrapper
def get_hsp_alignments(database_seq, query_seq, query_hits, blast_scoring=False, k=10, k_multiplier=2, use_local_alignments=False):
    """From a database sequence and list of ids that returns the list of alignments and their scores from HSPs

    Args:
        database_seq (dictionary): Mapping of seq IDs to sequences
        query_seq (str): Query sequence
        query_hits (array of 3-tuple): See get_top_sequences for more information
            First is the sequence id
            Second is the number of matches
            Third is a list of 3-tuple with db position, query position and sequence itself
        blast_scoring (bool, defaults to False): Flag to use blast scoring or not
        k (int): The size of the k-mer
        k_multiplier (int): The range to look at for merging. Default is 2.
        use_local_alignments (bool):
            Decides whether or not to use gapped extensions and try to calculate the best local alignment
            or simply calculate the range's score. Defaults to False.

    Returns:
        alignments (list of tuples of score, str ID, string, string):
            1. Score
            2. String ID database match
            3. Alignment of the database
            4. Alignment of the query
    """
    final_scores = []

    for hit, _, pos_tuple_list in query_hits:
        # now, for each hit that we have, we need to process the tuple list
        # to find the extension ranges to score on
        ranges = create_ranges(pos_tuple_list)

        # now, we take our ranges and choose the longest one, and then we find the local alignment
        ranges.sort(key=lambda x: x[0][1] - x[0][0], reverse=True)
        if VERBOSE:
            print(ranges)

        # Run the local alignment on this section against the entire query sequence
        # To ensure best alignments, we should extend the database sequence we look at to be
        # at least the length of the query sequence. 
        query_seq_length = len(query_seq)

        # so we want to make sure that the total length is query_seq_length +- k * k_multiplier
        top_range = ranges[0][0]
        
        if use_local_alignments:
            length_to_add = query_seq_length - (top_range[1] + k - top_range[0]) + 2 * k * k_multiplier

            # and we'll add this length equally if we can, otherwise, skewed one way
            lower_bound = max(0, top_range[0] - int(length_to_add / 2))
            upper_bound = min(len(database_seq[hit]),
                            top_range[1] + k + int(length_to_add / 2))

            # but this should speed things up significantly as we're looking at a subset of candidates
            score, n1, n2 = get_local_alignment(
                database_seq[hit][lower_bound:upper_bound],
                query_seq,
                blast_scoring
            )
        else:
            query_seq_range = ranges[0][1]
            n1 = database_seq[hit][top_range[0]:top_range[1]]
            n2 = query_seq[query_seq_range[0]:query_seq_range[1]]
            score = get_seq_score(n1, n2)
            
        # then, we add this to the final score and return it
        final_scores.append((score, hit, n1, n2))

    final_scores.sort(key=lambda x: x[0], reverse=True)
    return final_scores


if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', required=True)
    parser.add_argument('--query', required=True)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--k', type=int, default=10)
    parser.add_argument('--top_k', type=int, default=10)
    parser.add_argument('--blast_scoring', action='store_true')
    parser.add_argument('--gen_close_match', action='store_true')
    parser.add_argument('--align_hsp', action='store_true')
    parser.add_argument('--k_multiplier', type=int, default=2)
    parser.add_argument('--use_local_alignments', action='store_true')
    args = parser.parse_args()

    DEBUG = args.debug
    VERBOSE = args.verbose
    k = args.k
    top_k = args.top_k
    blast_scoring = args.blast_scoring
    gen_precise = not args.gen_close_match
    align_hsp = args.align_hsp
    k_multiplier = args.k_multiplier
    use_local_alignments = args.use_local_alignments

    database_seq = parse_fasta(args.db)
    query_seq = [value for value in parse_fasta(args.query).values()][0]
    if DEBUG:
        print(f'Number of database keys: {len(database_seq.keys())}')

    # Step 1: Build an inverse index
    inverse_index = build_inverse_index(database_seq, k, gen_precise)
    # Step 2: Generate the query k-mers
    query_kmers = generate_kmers(query_seq, k)
    # Step 3: Get the hits from the inverse index of the query k-mers
    query_hits = get_hits_from_index(inverse_index, query_kmers)
    if DEBUG:
        print(f'Query hits: {len(query_hits)}')

    # Step 4: Aggregate the query hits and sort it based on number of hits, take the topk
    query_hits = get_top_sequences(query_hits, top_k)
    if DEBUG:
        print(f'Top-k hits based on matches is: {[seq_id for seq_id, _, _ in query_hits]}')
    
    # Step 5: Here we branch off into two alternative methods.
    # The first will be to use q2 as a subroutine, and simply run the alignments from the query hits
    # The second will be to use the HSPs. This one is a bit more complicated as we have multiple HSPs per
    # database sequence. We will be accomplishing this by 
    
    # A) From the query hits, generate and rank all the scores through local alignments
    if not align_hsp:
        # now let's adjust the query hits just for the sequence id since we're not doing hsp
        query_hits = [seq_id for seq_id, _, _ in query_hits]
        best_scores = get_alignments(database_seq, query_seq, query_hits, blast_scoring)
        print(best_scores[0][1])

    # B) From the query hits, generate and rank all the scores through HSPs, which should be a lot faster
    else:
        best_scores = get_hsp_alignments(
            database_seq,
            query_seq,
            query_hits,
            blast_scoring,
            k,
            k_multiplier,
            use_local_alignments
        )
        print(best_scores[0][1])

    if DEBUG:
        print(f'Order of scoring: {[(hit, score)  for score, hit, _, _ in best_scores]}')
