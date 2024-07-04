#!/usr/bin/env python 

import pandas as pd
import sys

# Check if the expected number of input parameters is defined
if len(sys.argv) < 4:
    print("Usage: cov_wgs_matrix residuals_pc_matrix")
    sys.exit()

cov_wgs_path=sys.argv[2]
res_pc_path=sys.argv[1]
cohort_name=sys.argv[3]

cov_wgs_df = pd.read_csv(cov_wgs_path, sep='\t')
res_pc_df = pd.read_csv(res_pc_path, sep='\t')


cov_wgs_df.rename(columns={'Unnamed: 0': 'metric'}, inplace=True)
cov_wgs_df.set_index('metric', inplace=True)
cov_wgs_df = cov_wgs_df.transpose()
cov_wgs_df.reset_index(inplace=True)
cov_wgs_df.rename(columns={'index': 'sample'}, inplace=True)

res_pc_df.rename(columns={'Sample': 'sample'}, inplace=True)
res_pc_df.reset_index(inplace=True)
res_pc_df.drop(columns=['index'], inplace=True)
res_pc_df = res_pc_df[res_pc_df.columns[:11]]
# res_pc_df.set_index('sample', inplace=True)
# res_pc_df = res_pc_df.transpose()

merged_df = pd.merge(cov_wgs_df, res_pc_df, on='sample', how='inner')
merged_df.set_index('sample', inplace=True)
merged_df = merged_df.transpose()
merged_df.rename(columns={'sample': 'covariate'})

merged_df.to_csv(f'{cohort_name}_covariates_wgs_pcs_res_pcs.txt', sep='\t')