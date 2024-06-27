#!/usr/bin/env python3

"""
File:         heterozygosityCheck.py
Created:      2024/04/15
Last Changed: 2024/04/18
Author:       Peter Riesebos
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file, output_folder):
    # Load the .het file into a DataFrame
    het = pd.read_table(input_file, sep="\s+")

    # Calculate Het_Rate column
    het["Het_Rate"] = (het["OBS_CT"] - het["O(HOM)"]) / het["OBS_CT"]

    # Calculate mean and standard deviation of het_rate
    mean_het_rate = het["Het_Rate"].mean()
    std_dev_het_rate = het["Het_Rate"].std()

    # Return samples that deviate too much from the mean. To be thrown away.
    het_fail_samples = het[(het["Het_Rate"] < mean_het_rate - 3 * std_dev_het_rate) | (het["Het_Rate"] > mean_het_rate + 3 * std_dev_het_rate)]
    het_fail_samples_only = het_fail_samples["#IID"]
    het_fail_samples.to_csv(output_folder + 'HeterozygosityFailed.txt', sep='\t', index=False)
    het_fail_samples_only.to_csv(output_folder + 'HeterozygosityFailedSamplesOnly.txt', sep='\t', index=False, header=None)

    # Plot histogram
    plt.figure(figsize=(10, 6), facecolor="white")
    plt.hist(het["Het_Rate"], color="#000000", alpha=0.5, bins=20, edgecolor='black')

    # Add vertical lines for mean and ±3 standard deviations
    plt.axvline(mean_het_rate, color='blue', linestyle='--', label='Mean Het Rate')
    plt.axvline(mean_het_rate + 3 * std_dev_het_rate, color='red', linestyle='--', label='Mean + 3*SD')
    plt.axvline(mean_het_rate - 3 * std_dev_het_rate, color='green', linestyle='--', label='Mean - 3*SD')

    # Set labels and title
    plt.xlabel("Heterozygosity rate")
    plt.ylabel("Frequency")
    plt.title("Histogram of Heterozygosity Rate")
    plt.legend()

    # Save the plot
    plt.grid(True)
    plt.savefig(output_folder + 'HetCheck.png', bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process heterozygosity data and plot histogram.')
    parser.add_argument('input_file', type=str, help='Path to the input .het file')
    parser.add_argument('output_folder', type=str, help='Path to the output folder')
    args = parser.parse_args()

    main(args.input_file, args.output_folder)
