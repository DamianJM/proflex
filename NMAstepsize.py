# Process data set for NMA
# imports

import numpy as np
import sys

def load_data(filename:str):
    output_array = []
    with open(filename, 'r') as f:
        lines = [line.rstrip() for line in f]
        for line in lines:  
              temp = [float(i) for i in line.split(',')]
              output_array.append(temp)

    return output_array


def step_size(data: list):
    """Output average step size from list"""
    temp = np.array([])
    for i, item in enumerate(data):
        if i < len(data)-1:
            temp = np.append(temp, (abs(float(item) - float(data[i+1]))))
    
    return np.mean(temp)

def writeout_step(arr: list):
    return [step_size(i) for i in arr]

def main():
    data = sys.argv[1]
    arr = load_data(data)
    # step size calculation
    output = writeout_step(arr)
    with open("stepsize.txt", "w") as f:
        for i, item in enumerate(output):
            f.write(f'{i+1},{item}\n')


if __name__ == '__main__':
    main()
