# build_newick_tree.py.
# Author: Eric Huang
# Date: March 1st 2025

"""
This is the solution to problem 3, the Newick Tree construction and analysis problem.

As input it takes a list of sequence IDs to download from NCBI, a FASTA file
which represents a candidate viral sequence to output a separate UPGMA tree
based on these methods.
"""

from Bio import SeqIO, Phylo, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from kmer_comp import get_kmer_frequency_array
from typing import List
import argparse
import logging
import numpy as np
import os
import sys


def download_viral_seqs(input_file_path: str) -> List[SeqRecord]:
    """Parses the input file path to get the list of fasta sequences.

    Args:
        input_file_path (str): Path to input file.

    Returns:
        List[SeqRecord]: The list of fasta sequences
    """
    
    viral_ids = []
    with open(input_file_path, 'r') as input_file:
        for line in input_file:
            # remove the first part before the colon, and split again, and remove
            # any whitespaces after that
            viral_ids.extend([
                viral_id.strip()
                for viral_id in line.split(':')[1].split(',')
            ])

    viral_dir = './viral_logs/'

    if not os.path.exists(viral_dir):
        os.makedirs(viral_dir)

    # now for each viral sequence, we get the https://bioperl.org/formats/sequence_formats/GenBank_sequence_format
    # genbank format first, then we parse this using BioPython again
    # to get the final sequence record
    for viral_id in viral_ids:
        viral_file_path = os.path.join(viral_dir, f'{viral_id}.gb')
        if os.path.exists(viral_file_path):
            logging.debug(f'{viral_file_path} is cached!')
            continue

        Entrez.email = 'e48huang@uwaterloo.ca'
        handle = Entrez.efetch(
            db="nucleotide",
            id=viral_id,
            rettype='gb',
            retmode='text'
        )
        
        with open(viral_file_path, 'w') as viral_file:
            viral_file.write(handle.read())
            logging.debug(f'Finish downloading {viral_id} to {viral_file_path}')

    # now we parse through all the files and get the sequences as we need
    sequences = []

    for viral_id in viral_ids:
        viral_file_path = os.path.join(viral_dir, f'{viral_id}.gb')
        # parse the FASTA file
        with open(viral_file_path, 'r') as handle:
            # Parse the FASTA file
            for record in SeqIO.parse(handle, 'gb'):
                logging.debug(f'Record id: {record.id}')
                sequences.append(record)

    return sequences


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


def build_tree(seqs: List[SeqRecord], dist_fn, kmer_len):
    """Given a list of sequence records, build the UPGMA trees.

    Args:
        seqs (List[SeqRecord]): Given sequences
        dist_fn: A function that gives a distance between two vectors.
    
    Returns: Returns the UPGMA tree based on the distance function and a list of sequences.
    """
    logging.debug(f'Building a tree for seqs: {[seq.id for seq in seqs]}, {dist_fn.__name__}')
    
    # first we consider processing the sequences and grabbing their k-mer frequency vectors from q1
    labels = [seq.id for seq in seqs]
    frequency_vecs = [
        get_kmer_frequency_array(seq.seq, kmer_len)
        for seq in seqs
    ]

    # next, we take these frequency vectors and generate distance matrices with our dist_fn
    # note: according to https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
    # we need to build the lower triangular matrix
    distance_matrix = [
        [dist_fn(frequency_vecs[i], frequency_vecs[j]) for j in range(i + 1)]
        for i in range(len(frequency_vecs))
    ]

    # finally, we actually build the upgma tree, given the labels
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
    return 1 / (
        cosine_sim + np.finfo(float).eps
        if cosine_sim == 0 else
        cosine_sim
    ) - 1


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
    return (1 - pearson_corr) / (
        1 + pearson_corr + np.finfo(float).eps
        if pearson_corr == -1 else
        1 + pearson_corr
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_viral_file', required=True, help='Input viral sequences text file path, writes to ./viral_seqs/')
    parser.add_argument('-c', '--input_candidate_file', required=True, help='Input candidate sequence FASTA file path')
    parser.add_argument('-o', '--output_dir', help='Output directory, which if set will write everything to there')
    parser.add_argument('-k', '--kmer_len', default=4, type=int, help='The kmer length')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    input_viral_file = args.input_viral_file
    input_candidate_file = args.input_candidate_file

    # First grab the viral sequences
    seqs = download_viral_seqs(input_viral_file)
    candidate_seq = get_sequences(input_candidate_file)[0]

    # next handle the output directories
    # handle the output directory location
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # first do the viral only analysis, i.e. the analysis without the candidate
    output_dir = args.output_dir if args.output_dir else './'
    viral_only_dir = os.path.join(output_dir, 'viral_only_trees/')
    if not os.path.exists(viral_only_dir):
        os.makedirs(viral_only_dir)

    # first iterate by finding the upgma trees for viral seqs only
    dist_fns = [euclidean_distance, cosine_similarity, pearson_correlation]

    map_dist_fn_to_filename = {
        euclidean_distance.__name__: 'Euclidean.nwk',
        cosine_similarity.__name__: 'Cosine.nwk',
        pearson_correlation.__name__: 'Pearson.nwk',
    }

    # finally, iterate over all distance functions and calculate the upgma tree
    for dist_fn in dist_fns:
        viral_tree = build_tree(seqs, dist_fn, args.kmer_len)
        Phylo.write(
            viral_tree,
            os.path.join(viral_only_dir, map_dist_fn_to_filename[dist_fn.__name__]),
            "newick"
        )

    # reiterate by finding the upgma trees for everything together
    seqs.append(candidate_seq)
    for dist_fn in dist_fns:
        viral_tree = build_tree(seqs, dist_fn, args.kmer_len)
        Phylo.write(
            viral_tree,
            os.path.join(output_dir, map_dist_fn_to_filename[dist_fn.__name__]),
            "newick"
        )
