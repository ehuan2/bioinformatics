# build_newick_tree.py.
# Author: Eric Huang
# Date: March 1st 2025

"""
This is the solution to problem 3, the Newick Tree construction and analysis problem.

As input it takes a list of sequence IDs to download from NCBI, a FASTA file
which represents a candidate viral sequence to output a separate UPGMA tree
based on these methods.
"""

from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from kmer_comp import get_kmer_frequency_array
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
        List[str]: The list of fasta sequences
    """
    sequences = []

    # parse the FASTA file
    with open(input_file, 'r') as handle:
        # Parse the FASTA file
        for record in SeqIO.parse(handle, 'fasta'):
            logging.debug(f'Record id: {record.id}, sequence: {record.seq}')
            # return the first sequence we get
            sequences.append(record)

    return sequences


def build_upgma_tree(labels, distance_matrix):
    """Given a distance matrix and its corresponding labels, calculates the UPGMA tree and converts it to Newick format.

    Args:
        labels: The ordered labels of the distance matrix.
        distance_matrix: The symmetric matrix representing distances.
    """
    dm = DistanceMatrix(labels, distance_matrix)
    constructor = DistanceTreeConstructor()
    upgmatree = constructor.upgma(dm)
    return upgmatree


def build_tree(seqs: List[SeqRecord], dist_fn):
    """Given a list of sequence records, build the UPGMA trees.

    Args:
        seqs (List[SeqRecord]): Given sequences
        dist_fn: A function that gives a distance between two vectors.
    
    Returns: Returns the UPGMA tree based on the distance function and a list of sequences.
    """
    
    # first we consider processing the sequences and grabbing their k-mer frequency vectors
    labels = [seq.id for seq in seqs]
    frequency_vecs = [get_kmer_frequency_array(seq.seq) for seq in seqs]

    # next, we take these frequency vectors and generate distance matrices
    # according to: https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
    # we need to build the lower triangular matrix
    distance_matrix = [[dist_fn(frequency_vecs[i], frequency_vecs[j]) for j in range(i + 1)] for i in range(frequency_vecs)]

    # finally, we actually build the upgma tree
    return build_upgma_tree(labels, distance_matrix)


def euclidean_distance(vec1, vec2):
    # this defaults to the l2-norm
    return np.linalg.norm(np.array(vec1) - np.array(vec2))

def cosine_similarity(vec1, vec2):
    # the cosine similarity is simply just the dot product normalized by their l2-norms
    cosine_sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))

    # This however measures similarity and needs to be transformed into its own distance metric.

    # now because cosine similarity is the measure between two vectors and its angle,
    # and because no component is negative, its cosine similarity is between
    # 0 and 1, i.e. either exactly the same or at most orthogonal to each other.

    # So, it'll be within [0, 1] with 0 being more dissimilar and 1 being most similar.
    # We need to transform this to [0, \infty) with 0 being the least distance and most similar,
    # and \infty being the greatest distance and most dissimilar

    # To accomplish this, let's map 1 to 0 and 0 to \infty
    # this can be done easily by doing 1 / cosine_similarity - 1

    # In this case, we should with some small epsilon to prevent divide by 0
    return 1 / (cosine_sim + np.finfo(float).eps) - 1


def pearson_correlation(vec1, vec2):
    # it's simply the covariance normalized by its standard deviation
    # which is calculated with: https://numpy.org/devdocs/reference/generated/numpy.corrcoef.html
    # so we simply get the correlation coefficient between the two vectors
    pearson_corr = np.corrcoef(vec1, vec2)[0, 1]
    
    # Similarly to the cosine similarity, the pearson correlation also measures similarity
    # and needs to be transformed into a distance metric.

    # It should hold that the correlation between vec1 and vec2 represents how far
    # apart they are from each other, that with a negative correlation, vec1 and vec2
    # have frequencies that represent a negative correlation, and thus should be far away mutations
    # whereas if they are positive, then they are closer together.
    
    # So we need to map [-1, 1] to [0, \infty) with -1 being \infty and 1 being 0.
    # we can do this by using (1 - pearson) / (1 + pearson) by noting that
    # if pearson = -1, then we have \infty and if pearson = 1, then we have 0
    # if pearson = 0, then we have a distance of 1
    return (1 - pearson_corr) / (1 + pearson_corr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_viral_file', required=True, help='Input viral sequences FASTA file path')
    parser.add_argument('-c', '--input_candidate_file', required=True, help='Input candidate sequence FASTA file path')
    parser.add_argument('-o', '--output_dir', help='Output directory, which if set will write everything to there')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    input_viral_file = args.input_viral_file
    input_candidate_file = args.input_candidate_file

    seqs = get_sequences(input_viral_file)
    candidate_seq = get_sequences(input_candidate_file)[0]

    # handle the output directory location
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # first iterate by finding the upgma trees for viral seqs only
    euclidean_viral_tree = build_tree(seqs, euclidean_distance)
    cosine_viral_tree = build_tree(seqs, cosine_similarity)
    pearson_viral_tree = build_tree(seqs, pearson_correlation)

    # reiterate by finding the upgma trees for everything together
    seqs.append(candidate_seq)
    euclidean_all_seqs_tree = build_tree(seqs, euclidean_distance)
    cosine_all_seqs_tree = build_tree(seqs, cosine_similarity)
    pearson_all_seqs_tree = build_tree(seqs, pearson_correlation)

    # finally print them all out!
    Phylo.write(euclidean_viral_tree, "euclidean_viral.txt", "newick")
    Phylo.write(cosine_viral_tree, "cosine_viral.txt", "newick")
    Phylo.write(pearson_viral_tree, "pearson_viral.txt", "newick")

    Phylo.write(euclidean_all_seqs_tree, "Euclidean.txt", "newick")
    Phylo.write(cosine_all_seqs_tree, "Cosine.txt", "newick")
    Phylo.write(pearson_all_seqs_tree, "Pearson.txt", "newick")
