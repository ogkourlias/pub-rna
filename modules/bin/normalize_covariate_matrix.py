#!/usr/bin/env python 

import pandas as pd
from scipy.special import ndtri
import numpy as np
import argparse

def force_normalise(x, axis):
    return ndtri((pd.DataFrame(x).rank(axis=axis, ascending=True) - 0.5) / x.shape[axis])


if __name__ == '__main__':

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_path", help="Path to the covariate matrix (REQUIRED)")
    parser.add_argument("--cohort_name", help="Name of the cohort (REQUIRED)")
    args = parser.parse_args()

    # Get input directory
    matrix_path = args.matrix_path
    if matrix_path == None:
        raise ValueError('Please provide a path to the covariate matrix')

    # Get cohort name
    cohort_name = args.cohort_name
    if cohort_name == None:
        raise ValueError('Please provide a cohort name')

    # Read dataframe
    matrix = pd.read_csv(matrix_path, sep='\t')

    # Rename first empty column
    matrix.rename(columns={'Unnamed: 0': 'metric'}, inplace=True)

    # Set metric name as row index
    matrix.set_index('metric', inplace=True)

    # Remove rows without variance
    matrix = matrix[matrix.apply(lambda row: row.nunique() > 1, axis=1)]

    # Normalize the data
    normalized = force_normalise(matrix.values, 1)

    # Convert the normalized data back to a pandas df and replace columns and row indices
    final_matrix = pd.DataFrame(normalized)
    final_matrix.columns = matrix.columns
    final_matrix['metric'] = matrix.index
    final_matrix.set_index('metric', inplace=True)

    # Write the dataframe to csv
    final_matrix.to_csv(f'{cohort_name}_covariates_wgs_pcs_normalized.txt', sep='\t')