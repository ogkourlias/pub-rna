#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib.pyplot as plt 

if __name__ == '__main__':

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--pc_file", help="Path to the file containing the first two PCs (REQUIRED)")
    parser.add_argument("--z_score_threshold", help="Z-score threshold (REQUIRED)")
    args = parser.parse_args()

    # Get input directory
    pc_file = args.pc_file
    if pc_file == None:
        raise ValueError('Please provide a valid path to file with the first two PCs')

    # Get cohort name
    z_score_threshold = int(args.z_score_threshold)
    if z_score_threshold == None:
        raise ValueError('Please provide a valid z-score')

    # Load PCs
    df = pd.read_csv('pc1_2.txt', sep='\t')

    # Calculate z-score for first column
    df['z1'] = abs((df['Comp1'] - df['Comp1'].mean()) / df['Comp1'].std())
    df['z2'] = abs((df['Comp2'] - df['Comp2'].mean()) / df['Comp2'].std())

    # Filter out samples with z-score above threshold
    df_no_outliers = df[(df['z1'] < z_score_threshold) & (df['z2'] < z_score_threshold)]

    # Create outlier plot
    plt.scatter(df['Comp1'], df['Comp2'], c='red')
    plt.scatter(df_no_outliers['Comp1'], df_no_outliers['Comp2'], c='black')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.axvline(df['Comp1'].mean() - (3 * df['Comp1'].std()), color='red', ls='--')
    plt.axvline(df['Comp1'].mean() + (3 * df['Comp1'].std()), color='red', ls='--')
    plt.axhline(df['Comp2'].mean() - (3 * df['Comp2'].std()), color='red', ls='--')
    plt.axhline(df['Comp2'].mean() + (3 * df['Comp2'].std()), color='red', ls='--')
    plt.savefig('outliers.png')
    plt.close()

    # Write filtered df to file
    df_no_outliers.to_csv('pc1_2_no_outliers.txt', sep='\t', index=False)

    # Create list of remaining samples
    remaining_samples = df_no_outliers['Sample'].tolist()
    remaining_samples = [f'{s}\n' for s in remaining_samples]

    # Write list of remaning samples out to file
    with open('passed_samples.txt', 'w') as remaining_samples_file:
        remaining_samples_file.writelines(remaining_samples)

