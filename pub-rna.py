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
import argparse
import gzip
import cmd
import os
import csv
import subprocess

# File > Tissue > Study > Sample

# CLASSES
class CLI(cmd.Cmd):
    prompt = 'Pub-RNA>> '  # Change the prompt text
    intro = """

    Welcome to Pub-RNA CLI. Type help or ? to list commands.
    csv [csv-path] : Load CSV file. REQUIRED
    output [output-path] : Set base output path.
    exit : Exit the CLI

    Data Queries:
    tissues : List all tissues
    studies : List all studies
    samples : List all samples
    finished : List all finished samples


    """

    def __init__(self):
        """Loading CSV."""
        super(CLI, self).__init__()
        with open ("pub-rna.config", "r") as config:
            for line in config:
                if "csv_file" in line:
                    self.file = File(line.split("\"")[1])
                    self.file.open_csv()
                if "out_dir" in line:
                    self.output = line.split("\"")[1]
        print(f"Config loaded. \n Output path: {self.output}")
    
    def do_tissues(self, line):
        """Print all tissues."""
        print(self.file.tissues.keys())
    
    def do_studies(self, line):
        """Print all studies."""
        for study in line.split(" "):
            if not os.path.exists(f"{self.output}/{study}"):
                print(f"Study {study} has not been processed.")
            else:
                print(f"Study {study} has been processed.")

    def do_stop(self, line):
        """Stop study(s)."""
        for study in line.split(" "):
            os.system(f"tmux kill-session -t {study}")
            os.system(f"rm -r {study}")
    
    def do_start(self, line):
        """Start all studies."""
        abs_wd = os.path.abspath("../pub-rna/")
        for study in line.split(" "):
            if not os.path.exists(f"{self.output}/{study}"):
                print(f"Creating output directory for {study}.")
                os.mkdir(f"{self.output}/{study}")

            print(f"Creating {study}.config")
            with open (f"{abs_wd}/configs/{study}.config", "w") as study_config, open("pub-rna.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/{study}\"\n")
                    elif "input_txt" in line:
                        study_config.write(f"\tinput_txt = \"{abs_wd}/configs/{study}.txt\"\n")
                    else:
                        study_config.write(line)

                print(f"Creating sample list for {study}.")
                with open (f"{abs_wd}/configs/{study}.txt", "w") as sample_list:
                    for tissue in self.file.tissues:
                        if study in self.file.tissues[tissue].studies:
                            for sample in self.file.tissues[tissue].studies[study].samples:
                                sample_list.write(f"{self.file.tissues[tissue].studies[study].samples[sample].id}\n")
            
            if not os.path.exists(f"{abs_wd}/{study}"):
                os.mkdir(f"{abs_wd}/{study}")

            os.system(f"tmux new -d -s {study} -c {abs_wd}/{study}; tmux send-keys -t {study} \'ml Java && ../nextflow run ../main.nf -c {abs_wd}/configs/{study}.config\' Enter")


    
    def do_view(self, line):
        """Print samples."""
        ...
    
    def do_output(self, line):
        """Set output path."""
        self.output = line
    
    def do_finished_study(self, line):
        """Print all finished samples."""
        ...
        # for tissue in self.file.tissues:
        #     for study in self.file.tissues[tissue].studies:
        #         for sample in self.file.tissues[tissue].studies[study].samples:
        #             print(f"{tissue} : {study} : {sample}")


    def do_exit(self, line):
        """Exit the CLI."""
        return True

class File:
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, path):
        self.path = path

    def open_csv(self):
        self.tissues = {}
        with open(self.path, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            next(spamreader)
            for row in spamreader:
                if row[7] not in self.tissues:
                    self.tissues[row[7]] = Tissue(row[7])
                    self.tissues[row[7]].add_study(row[5], row[0])
                else:
                    self.tissues[row[7]].add_study(row[5], row[0])
        

class Tissue(File):
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, id): 
        self.id = id
        self.studies = {}
    
    def add_study(self, study_id, sample_id):
        if study_id not in self.studies:
            self.studies[study_id] = Study(study_id)
            self.studies[study_id].add_sample(sample_id)
        else:
            self.studies[study_id].add_sample(sample_id)
            


class Study(Tissue):
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, id):
        self.id = id
        self.samples = {}
    
    def add_sample(self, sample_id):
        if sample_id not in self.samples:
            self.samples[sample_id] = Sample(sample_id)
        else:
            print(f"Sample {sample_id} already exists")

class Sample(Study):
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, id):
        self.id = id

# FUNCTIONS
# def arg_parse():

#     """
#     Function to parse arguments.

#     :file: Input VCF file.
#     :reference: VCF file to be compared against.
#     :input_header: Text file containing tab seperated headers.
#     :comp_header: Text file containing tab seperated headers for compare file.
#     :chr: Test file idd as input chromosome.
#     :return: Argparse idspace object.
#     """

#     # Parser init
#     parser = argparse.ArgumentParser(prog="create-studies")
#     parser.add_argument("-s", "--studies", nargs='+', default=[])
#     parser.add_argument("-i", "--indir", type=str)
#     parser.add_argument("-o", "--outdir", type=str)
#     return parser.parse_args()


# MAIN
def main(args):
    """Main function"""
    CLI().cmdloop()
    # FINISH
    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))