import sys
import numpy as np
from scipy.stats import wilcoxon

def read_values(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    
    values = []
    for line in lines:
        if line.startswith('>'):
            continue  # Skip header lines
        values.append([float(num) for num in line.strip().split(',')])
    
    return values

def calculate_p_values(original_values, reconverted_values):
    p_values = []
    for orig_seq, recon_seq in zip(original_values, reconverted_values):
        if all(o == r for o, r in zip(orig_seq, recon_seq)):
            p_values.append(1.0)  # No difference between original and reconverted
        else:
            _, p_value = wilcoxon(orig_seq, recon_seq, zero_method='zsplit')
            p_values.append(p_value)
    
    return p_values

def main():
    original_values_path = sys.argv[1]
    reconverted_values_path = sys.argv[2]
    output_p_values_path = sys.argv[3]

    # Read original and reconverted values
    original_values = read_values(original_values_path)
    reconverted_values = read_values(reconverted_values_path)

    # Calculate p-values
    p_values = calculate_p_values(original_values, reconverted_values)

    # Write p-values to output file
    with open(output_p_values_path, "w") as output_file:
        for p_value in p_values:
            output_file.write(f"{p_value}\n")

if __name__ == "__main__":
    main()

