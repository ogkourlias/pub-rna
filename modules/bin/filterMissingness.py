#!/usr/bin/env python3

"""
File:         filterMissingness.py
Created:      2024/04/15
Last Changed: 2024/04/18
Author:       Peter Riesebos
"""

import argparse

def get_samples_with_low_missing_rate(file_path, threshold):
    samples = []
    with open(file_path, 'r') as file:
        # Skip header
        next(file)
        for line in file:
            print(line)
            sample_data = line.strip().split('\t')
            sample_id = sample_data[1]
            missing_rate = float(sample_data[4])
            if missing_rate < threshold:
                samples.append(sample_id)
    return samples

def write_samples_to_file(samples, output_file):
    with open(output_file, 'w') as file:
        for sample in samples:
            file.write(sample + '\n')

def main():
    parser = argparse.ArgumentParser(description='Filter samples based on missing rate in a .smiss file')
    parser.add_argument('file_path', type=str, help='Path to the .smiss file')
    parser.add_argument('output_file', type=str, help='Path to the output text file')
    parser.add_argument('--threshold', type=float, default=0.5, help='Missing rate threshold (default: 0.5)')
    args = parser.parse_args()

    samples_below_threshold = get_samples_with_low_missing_rate(args.file_path, args.threshold)

    write_samples_to_file(samples_below_threshold, args.output_file)
    print(f"Samples with missing rate < {args.threshold * 100:.2f}% written to {args.output_file}")

if __name__ == "__main__":
    main()
