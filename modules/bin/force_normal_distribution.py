#!/usr/bin/env python 

import sys
from scipy.special import ndtri
import numpy as np
import pandas as pd


def force_normalise(x, axis):
    '''Function that forces a normal distribution on a matrix'''
    return ndtri((pd.DataFrame(x).rank(axis=axis, ascending=True) - 0.5) / x.shape[axis])


if __name__ == '__main__':

    # Check if the expected number of input parameters is defined
    if len(sys.argv) < 3:
        print("Usage: infile axis")
        sys.exit()

    # Assign input parameters to variables
    infile=sys.argv[1]
    norm_axis=int(sys.argv[2])

    # Read in dataframe
    matrix = pd.read_csv(infile, sep='\t')

    # Set the index to the first '-' column
    matrix.set_index(matrix.columns[0], inplace=True)

    # Normalize the data
    normalized = force_normalise(matrix.values, norm_axis)

    # Convert matrix back to dataframe and add columns and an index
    final_matrix = pd.DataFrame(normalized)
    final_matrix.columns = matrix.columns
    final_matrix.index = matrix.index

    # Create out file name
    out_file_name = f"{infile.split('/')[-1].replace('.txt.gz', '')}.geneForceNormal.txt.gz"

    # Write dataframe out to file
    final_matrix.to_csv(out_file_name, compression='gzip', sep='\t')