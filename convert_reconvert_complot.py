import sys
import matplotlib.pyplot as plt

def read_sequences(file_path):
    with open(file_path, "r") as file:
        sequences = [line.strip() for line in file.readlines() if not line.startswith('>')]
    return sequences

def plot_sequences(original_sequence, converted_sequence, output_file_path):
    plt.figure(figsize=(10, 6))
    plt.plot(original_sequence, label='Original Sequence', color='blue')
    plt.plot(converted_sequence, label='Converted Sequence', color='orange')
    plt.xlabel('Position')
    plt.ylabel('Value')
    plt.title('Comparison of Original and Converted Sequences')
    plt.legend()
    plt.savefig(output_file_path)  # Save plot to the specified file
    plt.show()

def main():
    if len(sys.argv) < 5:
        print("Usage: script.py <original_file> <converted_file> <sequence_index> <output_file>")
        return

    original_file_path = sys.argv[1]
    converted_file_path = sys.argv[2]
    sequence_index = int(sys.argv[3])  # Index of the sequence of interest
    output_file_path = sys.argv[4]  # Output file for the plot

    # Read sequences from files
    try:
        original_sequences = read_sequences(original_file_path)
        converted_sequences = read_sequences(converted_file_path)
    except IOError as e:
        print(f"Error reading file: {e}")
        return

    if sequence_index < 0 or sequence_index >= len(original_sequences):
        print(f"Sequence index {sequence_index} is out of range.")
        return

    # Choose the sequence of interest
    original_sequence = original_sequences[sequence_index]
    converted_sequence = converted_sequences[sequence_index]

    # Convert sequences to numerical data
    try:
        original_sequence = [float(value) for value in original_sequence.split(',')]
        converted_sequence = [float(value) for value in converted_sequence.split(',')]
    except ValueError as e:
        print(f"Error converting sequences to float: {e}")
        return

    if len(original_sequence) != len(converted_sequence):
        print("Original and converted sequences are of different lengths.")
        return

    # Plot original sequence against converted sequence
    plot_sequences(original_sequence, converted_sequence, output_file_path)

if __name__ == "__main__":
    main()

