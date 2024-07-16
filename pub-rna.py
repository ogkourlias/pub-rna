#!/usr/bin/env python3

"""
    usage:
    python3 pub-rna.py
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import glob
import argparse
import gzip
import cmd
import os
import csv
import subprocess
import pandas as pd
from functools import reduce

# File > Tissue > Study > Sample

# CLASSES
class CLI(cmd.Cmd):
    prompt = 'Pub-RNA>> '  # Change the prompt text
    intro = """

    Welcome to Pub-RNA CLI. Type help or ? to list commands.
    Before contniuing, make sure configs in modules and pub-rna.config are configured correctly.

    exit : Exit the CLI
    Queries:
    gt [study] [study] ... : Start genotyping pipeline for study(s)
    qc [study] [study] ... : Start genotype QC pipeline for study(s)
    exp [study] [study] ... : Start expression pipeline for study(s)
    impute [study] [study] ... : Start imputation pipeline for study(s)
    check : Check all studies and create todo lists in work directory
    stats : Create stats file for specified tissue(s)
    join_exp [tissue] [tissue] ... : Join expression files for all studies specified tissue(s)
    meta [tissue] [tissue] ... : Prepare meta-analysis files for specified tissue(s)

    """

    def __init__(self):
        super(CLI, self).__init__()
        """Loading CSV."""
        self.pipeline_dir = os.path.abspath("../pub-rna/")
        with open ("pub-rna.config", "r") as config:
            for line in config:
                if "csv_file" in line:
                    self.file = File(line.split("\"")[1])
                    self.file.open_csv()
                if "out_dir" in line:
                    self.output = line.split("\"")[1]
                if "work_dir" in line:
                    self.work_dir = line.split("\"")[1]
        print(f"Config loaded. \n Output path: {self.output}")
        print(f"Work path: {self.work_dir}")

        # Todo check
        self.todo = {}
        for tissue in self.file.tissues:
            self.todo[tissue] = []
            if os.path.exists(f"{self.output}/{tissue}"):
                print(f"{tissue} exists")
                for study in self.file.tissues[tissue].studies:
                    if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                        if os.path.exists(f"{self.output}/{tissue}/{study}"):
                            if os.path.exists(f"{self.output}/{tissue}/{study}/genotypes"):
                                if len(os.listdir(f"{self.output}/{tissue}/{study}/genotypes")) >= 22:
                                    continue
                                else:
                                    self.todo[tissue].append(study)
                            else:
                                self.todo[tissue].append(study)
                            #print(os.listdir(f"{self.output}/{tissue}/{study}"))
                        else:
                            self.todo[tissue].append(study)

    def do_help(self, line):
        print(self.intro)

    def do_gt(self, line):
        """Start genotyping pipeline for study(s)."""
        print(f"Creating gt_configs directory at {self.work_dir}/gt_configs")
        os.makedirs(f"{self.work_dir}/gt_configs", exist_ok=True)
        for study in line.split(" "):
            for tissue in self.file.tissues:
                if study in self.file.tissues[tissue].studies:
                    study_tissue = tissue
            
            # Create output directory
            print(f"Creating output directory for {study} at {self.output}/{study_tissue}/{study}")
            os.makedirs(f"{self.output}/{study_tissue}/{study}", exist_ok=True)
            # Create config file
            print(f"Creating {study}.config at {self.work_dir}/gt_configs/{study}.config")
            with open (f"{self.work_dir}/gt_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/align_gt.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/{study_tissue}/{study}\"\n")
                        continue
                    if "input_txt" in line:
                        study_config.write(f"\tinput_txt = \"{self.work_dir}/gt_configs/{study}.txt\"\n")
                        continue
                    study_config.write(line)

                print(f"Creating sample list for {study} at {self.work_dir}/gt_configs/{study}.txt")
                with open (f"{self.work_dir}/gt_configs/{study}.txt", "w") as sample_list:
                    for tissue in self.file.tissues:
                        if study in self.file.tissues[tissue].studies:
                            for sample in self.file.tissues[tissue].studies[study].samples:
                                sample_list.write(f"{self.file.tissues[tissue].studies[study].samples[sample].id}\n")

            os.makedirs(f"{self.work_dir}/gt/{study}", exist_ok=True)
            os.system(f"tmux new -d -s gt_{study} -c {self.work_dir}/gt/{study}; tmux send-keys -t gt_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/align_gt.nf -c {self.work_dir}/gt_configs/{study}.config\' Enter")

    def do_qc(self, line):
        """Start genotype QC pipeline for study(s)."""
        print(f"Creating qc_configs directory at {self.work_dir}/qc_configs")
        os.makedirs(f"{self.work_dir}/qc_configs", exist_ok=True)
        for study in line.split(" "):
            print(f"Creating output directory for {study}/qc. at {self.output}/whole-blood/{study}/qc")
            os.makedirs(f"{self.output}/whole-blood/{study}/qc", exist_ok=True)

            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/qc_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/qc.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/whole-blood/{study}/qc\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/whole-blood/{study}/genotypes\"\n")
                        continue
                    if "cohortName" in line:
                        study_config.write(f"\tcohortName = \"{study}\"\n")
                        continue
                    study_config.write(line)

            os.makedirs(f"{self.work_dir}/qc/{study}", exist_ok=True)
            os.system(f"tmux new -d -s qc_{study} -c {self.work_dir}/qc/{study}; tmux send-keys -t qc_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/qc.nf -c {self.work_dir}/qc_configs/{study}.config\' Enter")

    def do_exp(self, line):
        """Start all studies."""
        print(f"Creating exp_configs directory at {self.work_dir}/exp_configs")
        os.makedirs(f"{self.work_dir}/exp_configs", exist_ok=True)

        for study in line.split(" "):
            os.makedirs(f"{self.output}/whole-blood/{study}/exp", exist_ok=True)
            print(f"Creating {study}.config at {self.work_dir}/exp_configs/{study}.config")
            with open (f"{self.work_dir}/exp_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/exp.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/whole-blood/{study}/exp\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/whole-blood/{study}\"\n")
                        continue
                    if "cohortName" in line:
                        study_config.write(f"\tcohortName = \"{study}\"\n")
                        continue
                    study_config.write(line)
                    
            os.makedirs(f"{self.work_dir}/exp/{study}", exist_ok=True)
            os.system(f"tmux new -d -s exp_{study} -c {self.work_dir}/exp/{study}; tmux send-keys -t exp_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/exp.nf -c {self.work_dir}/exp_configs/{study}.config\' Enter")

    def do_impute(self, line):
        print(f"Creating impute_configs directory at {self.work_dir}/impute_configs")
        os.makedirs(f"{self.work_dir}/impute_configs", exist_ok=True)
        for study in line.split(" "):
            print(f"Creating impute_configs directory at {self.work_dir}/impute_configs")
            os.makedirs(f"{self.work_dir}/impute_configs", exist_ok=True)
            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/impute_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/impute.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/whole-blood/{study}/impute\"\n")
                        continue
                    if "in_dir" in line:
                        study_config.write(f"\tin_dir = \"{self.output}/whole-blood/{study}/qc/final\"\n")
                        continue
                    if "cohort_name" in line:
                        study_config.write(f"\tcohort_name = \"{study}\"\n")
                        continue
                    study_config.write(line)
            
            print(f"Creating work directory for {study} at {self.work_dir}/impute/{study}")
            os.makedirs(f"{self.work_dir}/impute/{study}", exist_ok=True)
            os.system(f"tmux new -d -s impute_{study} -c {self.work_dir}/impute/{study}; tmux send-keys -t impute_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/impute.nf -c {self.work_dir}/impute_configs/{study}.config\' Enter")
    
    def do_check(self, line):
        with open(f"{self.work_dir}/gt_todo.txt", "w") as gt_todo, open(f"{self.work_dir}/qc_todo.txt", "w") as qc_todo, open(f"{self.work_dir}/exp_todo.txt", "w") as exp_todo, open(f"{self.work_dir}/impute_todo.txt", "w") as impute_todo:
            for tissue in self.file.tissues:
                for study in self.file.tissues[tissue].studies:
                    if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                        if not os.path.exists(f"{self.output}{tissue}/{study}"):
                            gt_todo.write(f"{tissue} {study}\n")
                        else:
                            if not os.path.exists(f"{self.output}{tissue}/{study}/genotypes"):
                                gt_todo.write(f"{tissue} {study}\n")
                            else:
                                if len(os.listdir(f"{self.output}{tissue}/{study}/genotypes")) < 22:
                                    gt_todo.write(f"{tissue} {study}\n")
                                    continue
                                if not os.path.exists(f"{self.output}/{tissue}/{study}/qc"):
                                    qc_todo.write(f"{tissue} {study}\n")
                                else:
                                    if not os.path.exists(f"{self.output}{tissue}/{study}/qc/final"):
                                        qc_todo.write(f"{tissue} {study}\n")
                                        continue
                                    else:
                                        if not os.path.exists(f"{self.output}{tissue}/{study}/exp"):
                                            exp_todo.write(f"{tissue} {study}\n")
                                        else:
                                            if not os.path.exists(f"{self.output}{tissue}/{study}/impute"):
                                                impute_todo.write(f"{tissue} {study}\n")
    def do_stats(self, line):
        """Create stats file for specified tissue(s)."""
        for tissue in self.file.tissues:
            if os.path.exists(f"{self.output}/{tissue}"):
                with open(f"{self.output}/{tissue}/stats.txt", "w") as stats, open(f"{self.output}/{tissue}/stats_lower_30.txt", "w") as stats_lower:
                    stats.write(f"study\tsample_pre_filter\tsample_post_filter\tvariant_pre_filter\tvariant_post_filter\tvariant_post_impute\n")
                    for study in self.file.tissues[tissue].studies:
                        variant_pre_impute = "N/A"
                        sample_post_filter = 'N/A'
                        variant_post_filter = 'N/A'
                        variant_post_impute = 'N/A'
                        variant_pre_filter = 'N/A'
                        if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/final"):
                                with open(f"{self.output}/{tissue}/{study}/qc/final/{study}.beagle.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples" in line:
                                            sample_post_filter = line.split(" ")[0]
                                        if "variants" in line:
                                            variant_post_filter = line.split(" ")[0]

                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/variant_filter"):
                                variant_pre_filter = 0
                                for chr in range(1, 23):
                                    with open(f"{self.output}/{tissue}/{study}/qc/variant_filter/chr{chr}.prefilter.stats", "r") as qc_stats:
                                        for line in qc_stats:
                                            if "SN" in line:
                                                if "number of records" in line:
                                                    variant_pre_filter += int(line.split("\t")[-1])
                        else:
                            stats_lower.write(f"{study}\t{len(self.file.tissues[tissue].studies[study].samples)}\n")

                                                    # if os.path.exists(f"{self.output}/{tissue}/{study}/impute"):
                        #     with open(f"{self.output}/{tissue}/{study}/impute/impute.log", "r") as imp_stats:
                        #         for line in imp_stats:
                        #             if "variants" in line:
                        #                 variant_post_impute = line.split(" ")[0]
                        
                        stats.write(f"{study}\t{len(self.file.tissues[tissue].studies[study].samples)}\t{sample_post_filter}\t{variant_pre_filter}\t{variant_post_filter}\t{variant_post_impute}\n")
    
    def do_join_exp(self, line):
        for tissue in line.split(" "):
            df_list = []
            for study in self.file.tissues[tissue].studies:
                if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                    if os.path.exists(f"{self.output}/whole-blood/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz"):
                        with gzip.open(f"{self.output}/whole-blood/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", "r") as exp:
                            df = pd.read_csv(exp, sep="\t")
                            df_list.append(df)

                    else:
                        print(f"{study} not finished")
            if len(df_list) > 0:
                df = reduce(lambda x, y: pd.merge(x, y, on="-", how='outer'), df_list)
                df.fillna(0, inplace=True)
                df.to_csv(f"{self.output}/{tissue}/{tissue}.txt.gz", compression='gzip', sep="\t", index=False)

    def do_meta(self, line):
        for tissue in line.split(" "):
            with open(f"{self.output}/{tissue}-meta/gte.txt", "w") as gte:
                df_list = []
                for study in self.file.tissues[tissue].studies:
                    if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                        # Check if study has expression file
                        if os.path.exists(f"{self.output}/whole-blood/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz"):
                            with gzip.open(f"{self.output}/whole-blood/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", "r") as exp:
                                df = pd.read_csv(exp, sep="\t")
                                df['-'] = df['-'].str.split(".").str[0]
                                df_list.append(df)

                        # Create gte file
                            for sample in self.file.tissues[tissue].studies[study].samples:
                                gte.write(f"{sample}\t{sample}\t{study}\n")
                
                if len(df_list) > 0:
                    df = reduce(lambda x, y: pd.merge(x, y, on="-", how='outer'), df_list)
                    df.fillna(0, inplace=True)
                    df.to_csv(f"{self.output}/{tissue}-meta/{tissue}-exp.txt.gz", compression='gzip', sep="\t", index=False)
                    df.to_csv(f"{self.output}/{tissue}-meta/{tissue}-genes.txt.gz", columns = "-", header = False, index=False, compression='gzip')


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
                if row[7].replace(" ", "-").lower() not in self.tissues:
                    self.tissues[row[7].replace(" ", "-").replace("/", "-").lower()] = Tissue(row[7])
                    self.tissues[row[7].replace(" ", "-").replace("/", "-").lower()].add_study(row[5], row[0])
                else:
                    self.tissues[row[7].replace(" ", "-").replace("/", "-").lower()].add_study(row[5], row[0])
        

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
