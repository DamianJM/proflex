import sys
import numpy as np
import matplotlib.pyplot as plt

def read_p_values(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    
    p_values = [float(line.strip()) for line in lines]
    
    return p_values

def plot_p_value_bins(p_values, output_plot_path):
    # Define bins for p-values
    bins = np.logspace(np.log10(1e-10), 1, 200)  # Use logarithmic scale for bins
    
    # Create histogram-like plot with bins colored by number of occupants
    plt.figure(figsize=(10, 6))
    plt.hist(p_values, bins=bins, color='skyblue', edgecolor='black')
    plt.title('Histogram of p-values with Colored Bins (Log Scale)')
    plt.xlabel('p-value')
    plt.ylabel('Number of Occupants')
    plt.xscale('log')  # Use logarithmic scale for x-axis
    plt.savefig(output_plot_path)
    plt.show()

def main():
    p_values_path = sys.argv[1]
    output_plot_path = sys.argv[2]

    # Read p-values from input file
    p_values = read_p_values(p_values_path)

    # Plot p-values into specific bins colored by number of occupants
    plot_p_value_bins(p_values, output_plot_path)

if __name__ == "__main__":
    main()

