import sys
from collections import defaultdict

def read_fasta(file_path):
    sequences = []
    with open(file_path, "r") as file:
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(''.join(sequence))
                    sequence = []
            else:
                sequence.append(line)
        if sequence:
            sequences.append(''.join(sequence))
    return sequences

def generate_kmers(sequence, k):
    """Generate k-mers from a single sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def count_unique_kmers(sequences, k):
    """Count unique k-mers in a list of sequences."""
    kmers = set()
    for seq in sequences:
        kmers.update(generate_kmers(seq, k))
    return len(kmers)

def main():
    if len(sys.argv) < 4:
        print("Usage: script.py <fasta_file> <k> <output_file>")
        return

    fasta_file_path = sys.argv[1]
    k = int(sys.argv[2])
    output_file_path = sys.argv[3]

    sequences = read_fasta(fasta_file_path)

    unique_kmer_count = count_unique_kmers(sequences, k)

    with open(output_file_path, "w") as output_file:
        output_file.write(f"Number of unique {k}-mers: {unique_kmer_count}\n")

    print(f"Number of unique {k}-mers: {unique_kmer_count}")

if __name__ == "__main__":
    main()
