#!/usr/bin/env python 

import pandas as pd
import sys
import numpy as np

# Check if the expected number of input parameters is defined
if len(sys.argv) < 5:
    print("Usage: covariate_matrix wgs_matrix, gte_file, cohort_name")
    sys.exit()

# Store command arguments in variables
covariate_path = sys.argv[1]
wgs_path = sys.argv[2]
gte_path = sys.argv[3]
cohort_name = sys.argv[4]

# Load GTE file and create GTE dict
with open(gte_path, 'r') as gte_file:
    gte_dict = {l.split('\t')[0]: l.split('\t')[1].strip() for l in gte_file.readlines()}

# Preprocess covariate dataframe
covariate_df = pd.read_csv(covariate_path, sep='\t')
covariate_df.rename(columns={'Unnamed: 0': 'sample'}, inplace=True)
covariate_df.set_index('sample', inplace=True)
covariate_df = covariate_df.transpose()
covariate_df.reset_index(inplace=True)

# Define function to link genotype samples to RNA-seq samples
def link_genotype_to_rna(gt_sample):
    if gt_sample in gte_dict:
        return gte_dict[gt_sample]
    else:
        return np.nan

# Preprocess WGS dataframe
wgs_df = pd.read_csv(wgs_path, sep='\t')
wgs_df['index'] = wgs_df['sample'].apply(link_genotype_to_rna)
print(wgs_df.head())

wgs_df = wgs_df.dropna(subset=['index'])
wgs_df.drop(columns=['sample'], inplace=True)

# Merge covariate and WGS dataframe and wrangle to final dataframe
merged_df = pd.merge(covariate_df, wgs_df, on='index', how='inner')
merged_df.rename(columns={'index': 'metric'}, inplace=True)
merged_df.set_index('metric', inplace=True)
merged_df = merged_df.transpose()

# Save dataframe to storage
merged_df.to_csv(f'{cohort_name}_covariates_wgs_pcs.txt', sep='\t')