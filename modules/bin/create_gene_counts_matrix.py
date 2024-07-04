#!/usr/bin/env python 

import pandas as pd
import argparse, os

if __name__ == '__main__':

    res_df = None

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", help="Path to the input dir (REQUIRED)")
    parser.add_argument("--cohort_name", help="Name of the cohort (REQUIRED)")
    args = parser.parse_args()

    # Get input directory
    input_dir = args.input_dir
    if input_dir == None:
        raise ValueError('Please provide a valid input directory')

    # Get cohort name
    cohort_name = args.cohort_name
    if cohort_name == None:
        raise ValueError('Please provide a cohort name')

    # # For every sample, read in gene counts and add to result dataframe
    # for sample in os.listdir(f'{input_dir}/star'):
    #     df = pd.read_csv(f'{input_dir}/star/{sample}/{sample}_ReadsPerGene.out.tab.gz', sep='\t', header=None)

    #     df = df.iloc[4:]
    #     df = df[[0, 3]]
    #     df.rename(columns={0: 'Gene', 3: sample}, inplace=True)

    #     try:
    #         res_df = pd.merge(res_df, df, on='Gene')
    #     except:
    #         res_df = df

    for sample in os.listdir(f'{input_dir}/'):
        if "qc" not in sample and "genotypes" not in sample and "exp" not in sample and "QC" not in sample:
            df = pd.read_csv(f'{input_dir}/{sample}/star/{sample}ReadsPerGene.out.tab.gz', sep='\t', header=None)

            df = df.iloc[4:]
            df = df[[0, 3]]
            df.rename(columns={0: 'Gene', 3: sample}, inplace=True)

            try:
                res_df = pd.merge(res_df, df, on='Gene')
            except:
                res_df = df


    # Write out resulting dataframe to out file
    res_df.to_csv(f'{cohort_name}_gene_counts.txt', sep='\t', index=False)