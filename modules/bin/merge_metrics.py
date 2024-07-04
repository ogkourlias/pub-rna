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

def read_stats(file):
    with open(file, "r") as 


# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    # Perform comparison and write output.
    dot_revive = DotRevive(args.input, args.output)
    dot_revive.walk()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
