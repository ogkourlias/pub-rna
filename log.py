#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -f [INPUT VCF FILE] -r [VCF COMPARISON FILE]
        -ih [INPUT HEADERS TEXT FILE] -ch [COMPARISON HEADER TEXT FILE]
        -chr [CHROMOSOME TEXT FILE] -n [CHUNKSIZE (Variant amount per run)]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import os
import subprocess

def log():
    sessions = subprocess.run(["tmux", "list-sessions"], capture_output=True).stdout.decode("utf-8").split("\n")
    for session in sessions:
        session = session.split(":")[0]
        os.makedirs("logs", exist_ok=True)
        with open (f"logs/{session}.log", "w") as log:
            log.write(subprocess.run(["tmux", "capture-pane", "-pt", session], capture_output=True).stdout.decode("utf-8"))
    return 0

# MAIN
def main(args):
    """Main function"""
    # Get args.
    log()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

