import csv

def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary of sequences keyed by headers."""
    sequences = {}
    header = None
    sequence = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            sequences[header] = ''.join(sequence)
    
    return sequences

def write_csv(file_path, data):
    """Writes data to a CSV file."""
    with open(file_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Header', 'Sequence1', 'Sequence2'])
        for header, seqs in data.items():
            csvwriter.writerow([header] + seqs)

def main(fasta_file1, fasta_file2, output_csv):
    # Read the FASTA files
    sequences1 = read_fasta(fasta_file1)
    sequences2 = read_fasta(fasta_file2)
    
    # Match sequences based on headers
    matched_sequences = {}
    for header, seq1 in sequences1.items():
        if header in sequences2:
            matched_sequences[header] = [seq1, sequences2[header]]
    
    # Write matched sequences to CSV
    write_csv(output_csv, matched_sequences)

if __name__ == "__main__":
    fasta_file1 = 'pdb2fasta_output.fasta'
    fasta_file2 = 'TOTAL_global.fasta'
    output_csv = 'matched_sequences_global.csv'
    
    main(fasta_file1, fasta_file2, output_csv)
