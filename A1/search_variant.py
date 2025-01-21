# search_variant.py
# Author: Eric Huang
# Date: January 21st 2025
import argparse

def parse_fasta(db_fna, query_fna):
    """Parses the database and query FASTA files 

    Args:
        db_fna (_type_): _description_
        query_fna (_type_): _description_
    """

if __name__ == '__main__':
    # first parse for the input file
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', required=True)
    parser.add_argument('--query', required=True)
    args = parser.parse_args()

    database_seq, patient_seq = parse_fasta(args.db, args.query)
