import rpy2.robjects as robjects
import sys, os, shutil

robjects.r.source("/media/damian/NMA_run.R")

def scale_nma(nma_file):
    with open(nma_file, 'r') as f:
        numbers = [float(line.strip()) for line in f.readlines()]

    min_num = min(numbers)
    max_num = max(numbers)
    scaled_numbers = [(num - min_num) / (max_num - min_num) for num in numbers]
    
    # Map scaled numbers to alphabet
    alphabet = 'abcdefghijklmnopqrst'.upper()
    mapped_letters = [alphabet[int(num * 19)] for num in scaled_numbers]
    
    return mapped_letters

def output_nma(nma_file):
    # read numbers from file
    with open(nma_file, "r") as f:
        numbers = [str(line.strip()) for line in f.readlines()]
    
    return numbers


def main():
    pdb_dir = sys.argv[1]
    output_file = sys.argv[2]

    # Get the total number of PDB files in the directory
    pdb_files = [file for file in os.listdir(pdb_dir) if file.endswith('.pdb')]
    total_files = len(pdb_files)

    # Counter for tracking progress
    progress_counter = 0

    with open(output_file, 'a') as df_out:
        # Loop over PDB files in the directory
        for pdb_file in pdb_files:
            pdb_path = os.path.join(pdb_dir, pdb_file)
            header = f"{pdb_file},"
            generate_nma = robjects.r['generate_nma']
            result = generate_nma(pdb_path, "temp.txt")
            output = output_nma("temp.txt")
            df_out.write(header + ",".join(output) + '\n')
            shutil.move(f'PDB_DB/{pdb_file}', "PDB_DB/FINISHED/")
            # Increment progress counter
            progress_counter += 1

            # Print progress
            print(f"Processed {progress_counter} out of {total_files} files.")


if __name__ == "__main__":
    main()








