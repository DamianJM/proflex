import pandas as pd
from collections import Counter
import argparse

def count_and_top_kmers(csv_file, column_index, k, top_n):
    df = pd.read_csv(csv_file, header=None)
    
    sequences = df.iloc[:, column_index].dropna()  # Drop any NaN values

    # Function to generate kmers of length k
    def generate_kmers(sequence, k):
        return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    
    # Counter to store kmer frequencies
    kmer_counter = Counter()
    
    # Iterate over each sequence and generate kmers
    for sequence in sequences:
        kmers = generate_kmers(sequence, k)
        kmer_counter.update(kmers)
    
    # Number of unique kmers
    num_unique_kmers = len(kmer_counter)
    
    # Get the top N kmers by abundance
    top_kmers = kmer_counter.most_common(top_n)
    
    return num_unique_kmers, top_kmers

def main():
    parser = argparse.ArgumentParser(description="Count unique kmers and get the top N kmers by abundance.")
    parser.add_argument("csv_file", type=str, help="Path to the CSV file")
    parser.add_argument("column_index", type=int, help="Index of the column to process (0-based index)")
    parser.add_argument("k", type=int, help="Length of the kmers")
    parser.add_argument("top_n", type=int, help="Number of top kmers to return")

    args = parser.parse_args()

    num_unique_kmers, top_kmers = count_and_top_kmers(args.csv_file, args.column_index, args.k, args.top_n)

    print(f'Number of unique {args.k}-mers: {num_unique_kmers}')
    print(f'Top {args.top_n} {args.k}-mers by abundance:')
    for kmer, count in top_kmers:
        print(f'{kmer}: {count}')

if __name__ == "__main__":
    main()

