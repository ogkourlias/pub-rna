import pandas as pd
import sys

# Check if the expected number of input parameters is defined
if len(sys.argv) < 3:
    print("Usage: infile sample_file")
    sys.exit()

# Assign input parameters to variables
infile=sys.argv[1]
sample_file=sys.argv[2]

# Load list of samples to include
with open(sample_file, 'r') as sample_file:
    samples_to_include = [s.strip() for s in sample_file.readlines()]

# Load matrix
df = pd.read_csv(infile, sep='\t')

# Create list of columns to include
included_cols = [df.columns[0]] + samples_to_include
included_cols = [col for col in included_cols if col in df.columns]

# Filter the dataframe
df = df[included_cols]

# Create output file name
out_file_name = f"{infile.split('/')[-1].replace('.txt.gz', '')}.samplesFiltered.txt.gz"

# Save dataframe to output file
df.to_csv(out_file_name, sep='\t', index=False)
