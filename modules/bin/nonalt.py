#!/usr/bin/env python3

"""
    usage:
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse



# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.

    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file idd as input chromosome.
    :return: Argparse idspace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="create-studies")
    parser.add_argument("-i", "--scount", type=str)
    parser.add_argument("-o", "--keepfile", type=str)
    return parser.parse_args()

def create_keep(scount, keepfile):
    with open(scount, 'r') as scount, open (keepfile, 'w') as keepfile:
        next(scount)
        for line in scount:
            line = line.strip().split("\t")
            obs_ct = sum(int(x) for x in line[1:])
            nonalt_ct = int(line[1]) + int(line[4])
            if nonalt_ct/obs_ct <= 0.9:
                keepfile.write(f"{line[0]}\n")  

            # MAIN
def main(args):
    """Main function"""
    args = arg_parse()
    create_keep(args.scount, args.keepfile)
    # FINISH
    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))