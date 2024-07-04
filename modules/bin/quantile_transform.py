#!/usr/bin/env python

from sklearn.preprocessing import quantile_transform
import sys
import pandas as pd

# Check if the expected number of input parameters is defined
if len(sys.argv) < 2:
    print("Usage: infile")
    sys.exit()

# Assign input parameters to variables
input_file=sys.argv[1]

# Load gene count matrix
df = pd.read_csv(input_file, sep='\t')

# Set the Gene columns as index
df.set_index('-', inplace=True)

# Do the quantile transformation
res = quantile_transform(df.values)

# Reattach columns and index to resulting dataframe
res_df = pd.DataFrame(res, index=df.index, columns=df.columns)

# Create sample out file name
out_file_name = f"{input_file.split('/')[-1].replace('.txt.gz', '')}.quantile_transform.txt.gz"

# Save res dataframe
res_df.to_csv(out_file_name, sep='\t')


