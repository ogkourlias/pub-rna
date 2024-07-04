#!/usr/bin/env python 

import pandas as pd
import argparse

if __name__ == '__main__':

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--covariate_pcs", help="Path to the covariate PCs file (REQUIRED)")
    parser.add_argument("--covariates_explained_variance", help="Path to the covariates explained variance file (REQUIRED)")
    parser.add_argument("--explained_variance_threshold", help="The threshold of explained variance of the number of PCs to select (REQUIRED)")
    args = parser.parse_args()

    # Get covariates matrix path
    covariate_pcs = args.covariate_pcs
    if covariate_pcs == None:
        raise ValueError('Please provide a valid path to the covariates matrix')

    # Get covariates matrix path
    covariates_explained_variance = args.covariates_explained_variance
    if covariates_explained_variance == None:
        raise ValueError('Please provide a valid path to the covariates explained variance file')

    # Get input directory
    explained_variance_threshold = float(args.explained_variance_threshold)
    if explained_variance_threshold == None:
        raise ValueError('Please provide a valid explained variance threshold')

    # Load covariates explained variance dataframe    
    explained_variance_df = pd.read_csv(covariates_explained_variance, sep='\t')

    print(explained_variance_df.head())

    # Get the number of PCs that explain the threshold of variance
    print(explained_variance_df[explained_variance_df['CumulativeExplainedVariance'] > explained_variance_threshold])
    num_pcs = int(explained_variance_df[explained_variance_df['CumulativeExplainedVariance'] > explained_variance_threshold].iloc[0]['PC'][2:]) + 1

    # Load covariate PCs
    pc_df = pd.read_csv(covariate_pcs, sep='\t')

    # Filter the number of PCs 
    pc_df_filtered = pc_df[pc_df.columns[:num_pcs]] 

    # Set 'Sample' column as index of the dataframe
    pc_df_filtered.set_index('Sample',inplace=True)

    # Transpose the dataframe
    pc_df_transposed = pc_df_filtered.transpose()

    # Write out filtered and transposed dataframe to file
    pc_df_transposed.to_csv('covariates-pca_PCs-filtered-transpose.txt',sep='\t')