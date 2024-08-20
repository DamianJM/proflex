import sys
import numpy as np

def alpha(numbers, percentiles, alphabet='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'):
    n_bins = len(alphabet)
    bin_indices = np.digitize(numbers, percentiles) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    mapped_letters = [alphabet[idx] for idx in bin_indices]
    return mapped_letters

def reverse_alpha(letters, percentiles, alphabet='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'):
    n_bins = len(alphabet)
    bin_widths = np.diff(percentiles)
    bin_midpoints = percentiles[:-1] + bin_widths / 2
    letter_to_value = {alphabet[i]: bin_midpoints[i] for i in range(n_bins)}
    values = [letter_to_value[letter] for letter in letters]
    return values

def main():
    if len(sys.argv) < 5:
        print("Usage: script.py <input_file> <encoded_output_file> <decoded_output_file> <mode>")
        print("<mode> should be either 'global' or 'sequence'")
        return

    input_file_path = sys.argv[1]
    encoded_output_file_path = sys.argv[2]
    decoded_output_file_path = sys.argv[3]
    mode = sys.argv[4]

    alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    n_bins = len(alphabet)

    with open(input_file_path, "r") as input_file:
        lines = input_file.readlines()

    result_encoded = []
    result_decoded = []

    if mode == 'global':
        # Collect all values for global binning
        all_values = []
        for line in lines:
            if not line.startswith('>') and line.strip():  # Check for non-empty lines
                try:
                    values = list(map(float, line.strip().split(',')))
                    all_values.extend(values)
                except ValueError as e:
                    print(f"Skipping line due to error: {e}")
                    continue

        # Compute percentiles for global binning
        percentiles = np.percentile(all_values, np.linspace(0, 100, n_bins + 1))

        # Apply global binning and reverse it
        total_sequences = sum(1 for line in lines if not line.startswith('>') and line.strip())
        current_sequence = 0
        for line in lines:
            if not line.startswith('>') and line.strip():  # Check for non-empty lines
                try:
                    numbers = list(map(float, line.strip().split(',')))
                    encoded = alpha(numbers, percentiles)
                    decoded = reverse_alpha(encoded, percentiles)
                    result_encoded.append(encoded)
                    result_decoded.append(decoded)
                except ValueError as e:
                    print(f"Skipping line due to error: {e}")
                    result_encoded.append('ERROR')
                    result_decoded.append('ERROR')
                current_sequence += 1
                print(f"Processed sequence {current_sequence}/{total_sequences}")
            else:
                result_encoded.append(line)
                result_decoded.append(line)

    elif mode == 'sequence':
        # Sequence-specific binning
        total_sequences = sum(1 for line in lines if not line.startswith('>') and line.strip())
        current_sequence = 0
        for line in lines:
            if not line.startswith('>') and line.strip():  # Check for non-empty lines
                try:
                    numbers = list(map(float, line.strip().split(',')))
                    percentiles = np.percentile(numbers, np.linspace(0, 100, n_bins + 1))
                    encoded = alpha(numbers, percentiles)
                    decoded = reverse_alpha(encoded, percentiles)
                    result_encoded.append(encoded)
                    result_decoded.append(decoded)
                except ValueError as e:
                    print(f"Skipping line due to error: {e}")
                    result_encoded.append('ERROR')
                    result_decoded.append('ERROR')
                current_sequence += 1
                print(f"Processed sequence {current_sequence}/{total_sequences}")
            else:
                result_encoded.append(line)
                result_decoded.append(line)
    else:
        print(f"Invalid mode: {mode}. Please choose 'global' or 'sequence'.")
        return

    # Write encoded output
    with open(encoded_output_file_path, "w") as output_file:
        for line in result_encoded:
            if isinstance(line, list):
                output_file.write(','.join(map(str, line)) + '\n')
            else:
                output_file.write(line)

    # Write decoded output
    with open(decoded_output_file_path, "w") as output_file:
        for line in result_decoded:
            if isinstance(line, list):
                output_file.write(','.join(map(str, line)) + '\n')
            else:
                output_file.write(line)

if __name__ == "__main__":
    main()

