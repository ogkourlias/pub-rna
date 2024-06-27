#!/usr/bin/env python3

"""
    usage:
        ./dotrevive.py -i [INPUT VCF FILE] -o [OUTPUT VCF FILE]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import gzip
from pathlib import Path
import os

# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.
    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="dotrevive")
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-o", "--output", type=str)
    return parser.parse_args()

def pop(poptxt, out):

    """
    Function to parse arguments.
    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """
    with open(poptxt, "r") as pop_f, open(out, "w") as out_f:
        pop_list = []  
        sample_dict = {}
        for i,line in enumerate(pop_f):
            line = line.strip("\n")
            if i == 0:
                for pop in line.split("\t")[1:]:
                    pop_list.append(pop)
                    sample_dict[pop] = []
                out_f.write(line + "\tassigned\n")
            else:
                shortest = 0
                assigned = "nan"
                split = line.split("\t")
                sample = split[0]
                pops = split[1:]
                for i,pop_float in enumerate(pops):
                    if float(shortest) == 0 or float(pop_float) < float(shortest):
                        shortest = pop_float
                        assigned = pop_list[i]
                sample_dict[assigned].append(sample)
                out_f.write(line + "\t" + assigned + "\n")
    for key in sample_dict.keys():
        Path(f"pops/{key}").mkdir(parents=True, exist_ok=True)
        with open(f"pops/{key}/{key}.txt", "w") as key_f:
            line_list = [f"0\t{val}" for val in sample_dict[key]]
            key_f.write("\n".join(line_list))
# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    # Perform comparison and write output.
    pop(args.input, args.output)
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
