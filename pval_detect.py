import sys

def read_p_values(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    
    # Convert each line to a float and store in a list
    p_values = [float(line.strip()) for line in lines]
    
    return p_values

def count_p_values_above_threshold(p_values, threshold=0.05):
    # Count the number of p-values greater than the threshold
    count_above_threshold = sum(1 for p in p_values if p > threshold)
    return count_above_threshold

def main():
    if len(sys.argv) < 2:
        print("Usage: script.py <p_values_file>")
        return
    
    p_values_path = sys.argv[1]

    # Read p-values from input file
    p_values = read_p_values(p_values_path)

    # Count the number of p-values above 0.05
    count_above_005 = count_p_values_above_threshold(p_values, 0.05)

    print(f"Number of p-values above 0.05: {count_above_005}")

if __name__ == "__main__":
    main()
