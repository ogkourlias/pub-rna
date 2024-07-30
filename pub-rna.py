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
import seaborn as sns
import matplotlib.pyplot as plt
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
    gt_merge [tissue] [tissue] ... : Merge genotyping files for specified tissue(s)
    meta_prep [tissue] [tissue] ... : Prepare meta analysis files for specified tissue(s)
    meta [tissue] [tissue] ... : Start meta analysis pipeline for specified tissue(s)

    Suggested workflow:
    1. gt [study] [study] ...
    2. qc [study] [study] ...
    3. exp [study] [study] ...
    4. impute [study] [study] ...
    5. gt_merge [tissue] [tissue] ...
    6. meta_prep [tissue] [tissue] ... 
    7. meta [tissue] [tissue] ...
    
    """

    def __init__(self):
        super(CLI, self).__init__()
        """Loading CSV."""
        self.pipeline_dir = os.path.abspath("../pub-rna/")
        self.dir_blacklist = ["exp", "genotypes", "impute", "qc"]
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

    def do_gt_selective(self, line):
        """Start genotyping pipeline for study(s)."""
        print(f"Creating gt_configs directory at {self.work_dir}/gt_configs")
        os.makedirs(f"{self.work_dir}/gt_configs", exist_ok=True)
        line = line.split(" ")
        study_tissue = line[0].lower()
        studies = line[1:]
        for study in studies:
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
                        if study in self.file.tissues[tissue].studies and tissue == study_tissue:
                            for sample in self.file.tissues[tissue].studies[study].samples:
                                sample_list.write(f"{self.file.tissues[tissue].studies[study].samples[sample].id}\n")

            os.makedirs(f"{self.work_dir}/gt/{study}", exist_ok=True)
            os.system(f"tmux new -d -s gt_{study} -c {self.work_dir}/gt/{study}; tmux send-keys -t gt_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/align_gt.nf -c {self.work_dir}/gt_configs/{study}.config\' Enter")


    def do_qc(self, line):
        """Start genotype QC pipeline for study(s)."""
        print(f"Creating qc_configs directory at {self.work_dir}/qc_configs")
        os.makedirs(f"{self.work_dir}/qc_configs", exist_ok=True)
        for study in line.split(" "):
            for tissue in self.file.tissues:
                if study in self.file.tissues[tissue].studies:
                    study_tissue = tissue
            
            print(f"Creating output directory for {study}/qc. at {self.output}/{study_tissue}/{study}/qc")
            os.makedirs(f"{self.output}/{study_tissue}/{study}/qc", exist_ok=True)

            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/qc_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/qc.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/{study_tissue}/{study}/qc\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/{study_tissue}/{study}/genotypes\"\n")
                        continue
                    if "cohortName" in line:
                        study_config.write(f"\tcohortName = \"{study}\"\n")
                        continue
                    study_config.write(line)

            os.makedirs(f"{self.work_dir}/qc/{study}", exist_ok=True)
            os.system(f"tmux new -d -s qc_{study} -c {self.work_dir}/qc/{study}; tmux send-keys -t qc_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/qc.nf -c {self.work_dir}/qc_configs/{study}.config\' Enter")
        
    def do_qc_selective(self, line):
        """Start genotype QC pipeline for study(s)."""
        print(f"Creating qc_configs directory at {self.work_dir}/qc_configs")
        os.makedirs(f"{self.work_dir}/qc_configs", exist_ok=True)
        line = line.split(" ")
        study_tissue = line[0].lower()
        studies = line[1:]
        for study in studies:
            print(f"Creating output directory for {study}/qc. at {self.output}/{study_tissue}/{study}/qc")
            os.makedirs(f"{self.output}/{study_tissue}/{study}/qc", exist_ok=True)

            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/qc_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/qc.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/{study_tissue}/{study}/qc\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/{study_tissue}/{study}/genotypes\"\n")
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
            for tissue in self.file.tissues:
                if study in self.file.tissues[tissue].studies:
                    study_tissue = tissue

            os.makedirs(f"{self.output}/{study_tissue}/{study}/exp", exist_ok=True)
            print(f"Creating {study}.config at {self.work_dir}/exp_configs/{study}.config")
            with open (f"{self.work_dir}/exp_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/exp.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/{study_tissue}/{study}/exp\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/{study_tissue}/{study}\"\n")
                        continue
                    if "cohortName" in line:
                        study_config.write(f"\tcohortName = \"{study}\"\n")
                        continue
                    study_config.write(line)
                    
            os.makedirs(f"{self.work_dir}/exp/{study}", exist_ok=True)
            os.system(f"tmux new -d -s exp_{study} -c {self.work_dir}/exp/{study}; tmux send-keys -t exp_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/exp.nf -c {self.work_dir}/exp_configs/{study}.config\' Enter")

    def do_exp_selective(self, line):
        print(f"Creating exp_configs directory at {self.work_dir}/exp_configs")
        os.makedirs(f"{self.work_dir}/exp_configs", exist_ok=True)
        line = line.split(" ")
        study_tissue = line[0].lower()
        studies = line[1:]
        for study in studies:
            os.makedirs(f"{self.output}/{study_tissue}/{study}/exp", exist_ok=True)
            print(f"Creating {study}.config at {self.work_dir}/exp_configs/{study}.config")
            with open (f"{self.work_dir}/exp_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/exp.config", "r") as config:
                for line in config:
                    if "outDir" in line:
                        study_config.write(f"\toutDir = \"{self.output}/{study_tissue}/{study}/exp\"\n")
                        continue
                    if "inputDir" in line:
                        study_config.write(f"\tinputDir = \"{self.output}/{study_tissue}/{study}\"\n")
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
            for tissue in self.file.tissues:
                if study in self.file.tissues[tissue].studies:
                    study_tissue = tissue

            print(f"Creating impute_configs directory at {self.work_dir}/impute_configs")
            os.makedirs(f"{self.work_dir}/impute_configs", exist_ok=True)
            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/impute_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/impute.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/{study_tissue}/{study}/impute\"\n")
                        continue
                    if "in_dir" in line:
                        study_config.write(f"\tin_dir = \"{self.output}/{study_tissue}/{study}/qc/final\"\n")
                        continue
                    if "cohort_name" in line:
                        study_config.write(f"\tcohort_name = \"{study}\"\n")
                        continue
                    study_config.write(line)
            
            print(f"Creating work directory for {study} at {self.work_dir}/impute/{study}")
            os.makedirs(f"{self.work_dir}/impute/{study}", exist_ok=True)
            os.system(f"tmux new -d -s impute_{study} -c {self.work_dir}/impute/{study}; tmux send-keys -t impute_{study} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/impute.nf -c {self.work_dir}/impute_configs/{study}.config\' Enter")
        
    def do_impute_selective(self, line):
        print(f"Creating impute_configs directory at {self.work_dir}/impute_configs")
        os.makedirs(f"{self.work_dir}/impute_configs", exist_ok=True)
        line = line.split(" ")
        study_tissue = line[0].lower()
        studies = line[1:]
        for study in studies:
            print(f"Creating impute_configs directory at {self.work_dir}/impute_configs")
            os.makedirs(f"{self.work_dir}/impute_configs", exist_ok=True)
            print(f"Creating {study}.config")
            with open (f"{self.work_dir}/impute_configs/{study}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/impute.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/{study_tissue}/{study}/impute\"\n")
                        continue
                    if "in_dir" in line:
                        study_config.write(f"\tin_dir = \"{self.output}/{study_tissue}/{study}/qc/final\"\n")
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
        rows = []
        for tissue in line.split(" "):
            if os.path.exists(f"{self.output}/{tissue}"):
                with open(f"{self.output}/{tissue}/stats.txt", "w") as stats, open(f"{self.output}/{tissue}/stats_lower_30.txt", "w") as stats_lower:
                    col_list = [
                        "study",
                        "variant_pre_filter",
                        "variant_post_filter",
                        "variant_pre_impute",
                        "variant_post_impute",
                        "variant_post_missing",
                        "variant_post_het",
                        "variant_post_alt",
                        "variant_post_related",
                        "samples_input",
                        "sample_post_missing",
                        "sample_post_het",
                        "sample_post_alt",
                        "sample_post_related"
                    ]
                    stats.write('\t'.join(col_list) + "\n")
                    studies = [study for study in os.listdir(f"{self.output}/{tissue}") if ".txt" not in study]
                    for study in studies:
                        variant_pre_filter = 'nan'
                        variant_post_filter = 'nan'
                        variant_pre_impute = "nan"
                        variant_post_impute = 'nan'
                        variant_post_missing = 'nan'
                        variant_post_het = 'nan'
                        variant_post_alt = 'nan'
                        variant_post_related = 'nan'
                        samples_input = 'nan'
                        sample_post_missing = 'nan'
                        sample_post_het = 'nan'
                        sample_post_alt = 'nan'
                        sample_post_related = 'nan'
                        if os.path.exists(f"{self.output}/{tissue}/{study}"):
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/final"):
                                with open(f"{self.output}/{tissue}/{study}/qc/final/{study}.beagle.log", "r") as qc_stats:
                                    for line in qc_stats:
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
                            
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/missingness"):
                                sample_post_missing = 0
                                variant_post_missing = 0
                                with open(f"{self.output}/{tissue}/{study}/qc/missingness/{study}.keep.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples remaining" in line:
                                            sample_post_missing = line.split(" ")[1]
                                with open(f"{self.output}/{tissue}/{study}/qc/missingness/{study}.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "founders" in line:
                                            samples_input = line.split(" ")[0]
                            
                            elif os.path.exists(f"{self.output}/{tissue}/{study}/qc/missingness/null.log"):
                                sample_post_missing = 0
                                variant_post_missing = 0
                                with open(f"{self.output}/{tissue}/{study}/qc/missingness/null.keep.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples remaining" in line:
                                            sample_post_missing = line.split(" ")[1]
                                with open(f"{self.output}/{tissue}/{study}/qc/missingness/null.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "founders" in line:
                                            samples_input = line.split(" ")[0]
                            
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/het"):
                                sample_post_het = 0
                                variant_post_het = 0
                                with open(f"{self.output}/{tissue}/{study}/qc/het/{study}.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples remaining" in line:
                                            sample_post_het = line.split(" ")[1]
                            
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/alt_freq"):
                                sample_post_alt = 0
                                variant_post_alt = 0
                                with open(f"{self.output}/{tissue}/{study}/qc/alt_freq/{study}.over10.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples remaining" in line:
                                            sample_post_alt = line.split(" ")[1]
                            
                            if os.path.exists(f"{self.output}/{tissue}/{study}/qc/relatedness"):
                                sample_post_related = 0
                                variant_post_related = 0
                                with open(f"{self.output}/{tissue}/{study}/qc/relatedness/{study}.RelatednessCheck.log", "r") as qc_stats:
                                    for line in qc_stats:
                                        if "samples remaining" in line:
                                            sample_post_related = line.split(" ")[1]
                                            
                                                    # if os.path.exists(f"{self.output}/{tissue}/{study}/impute"):
                        #     with open(f"{self.output}/{tissue}/{study}/impute/impute.log", "r") as imp_stats:
                        #         for line in imp_stats:    
                        #             if "variants" in line:
                        #                 variant_post_impute = line.split(" ")[0]
                        
                        #stats.write(f"{study}\t{len(self.file.tissues[tissue].studies[study].samples)}\t{sample_post_filter}\t{variant_pre_filter}\t{variant_post_filter}\t{variant_post_impute}\n")
                        rows.append([study,
                            variant_pre_filter,
                            variant_post_filter,
                            variant_pre_impute,
                            variant_post_impute,
                            variant_post_missing,
                            variant_post_het,
                            variant_post_alt,
                            variant_post_related,
                            samples_input,
                            sample_post_missing,
                            sample_post_het,
                            sample_post_alt,
                            sample_post_related])

        df = pd.DataFrame(rows, columns=col_list)
        df.set_index('study', inplace=True)
        df = df.apply(pd.to_numeric, args=('coerce',))
        df.to_csv(f"{self.output}/stats.txt", sep="\t")
        sample_df = df[['samples_input','sample_post_missing', 'sample_post_het', 'sample_post_alt', 'sample_post_related']]
        #sample_df = sample_df.T
        print(sample_df)
        labs = ['Input','Missingness','Heterzygosity', 'Alt-Call Rate', 'Relatedness']
        colors = ['C0','C1', 'C2', 'C3', 'C4']
        # plt.figure()
        # for i in range(0, len(sample_df.columns)):
        #     plt.subplot(int(len(sample_df.columns) / 3), 3,i+1)
        #     sample_df.iloc[:,i].plot(color=colors, kind='bar', legend=False, xticks=[])
        #     plt.title(sample_df.columns[i])
        # fig = plt.gcf()
        # fig.legend(labs, loc='lower right', bbox_to_anchor=(1,-0.1), ncol=len(labs), bbox_transform=fig.transFigure, edgecolor=colors)
        n_cols = 3
        n_rows = int(len(sample_df.index) / n_cols)

        palette = sns.color_palette("Set2", n_colors=sample_df.shape[1])

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows), sharey=False)

        for i, (index, row) in enumerate(sample_df.iterrows()):
            row_idx = i // n_cols
            col_idx = i % n_cols
            sns.barplot(x=row.index, y=row.values, ax=axes[row_idx, col_idx], palette=palette)
            axes[row_idx, col_idx].set_title(f'{index}')
            axes[row_idx, col_idx].set_xlabel('')
            axes[row_idx, col_idx].set_xticklabels([])
            if col_idx == 0:
                axes[row_idx, col_idx].set_ylabel('Samples')

        handles = [plt.Line2D([0], [0], color=palette[j], lw=4) for j in range(sample_df.shape[1])]
        labels = sample_df.columns.tolist()
        fig.legend(handles, labs, loc='upper right', ncol=1, bbox_to_anchor=(0.5, 1.05),)

        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        plt.show()
        plt.savefig(f"{self.output}/sample_stats.png", dpi=300, bbox_inches='tight')
        plt.clf()
        plt.cla()
        plt.close()
        # make total
        total = sample_df.sum()
        total.plot(kind='bar', color=palette, legend=False, xticks=[],
         title='Samples remaining after GT filter steps.', ylabel='Samples', xlabel='Filter Steps', figsize=(10,5) )
        for i in range(len(total.values)):
            plt.text(i,total.values[i],total.values[i],  ha='center', va='bottom')
        plt.legend(handles, labs, loc=4, ncol=1, bbox_to_anchor=(0.5, 1.05),)
        plt.savefig(f"{self.output}/sample_stats_total.png", dpi=300, bbox_inches='tight')

    def do_join_exp(self, line):
        for tissue in line.split(" "):
            df_list = []
            for study in self.file.tissues[tissue].studies:
                if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                    if os.path.exists(f"{self.output}/{study_tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz"):
                        with gzip.open(f"{self.output}/{tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", "r") as exp:
                            df = pd.read_csv(exp, sep="\t")
                            df_list.append(df)

                    else:
                        print(f"{study} not finished")
            if len(df_list) > 0:
                df = reduce(lambda x, y: pd.merge(x, y, on="-", how='outer'), df_list)
                df.fillna(0, inplace=True)
                df.to_csv(f"{self.output}/{tissue}/{tissue}.txt.gz", compression='gzip', sep="\t", index=False)

    def do_meta_prep(self, line):
        for tissue in line.split(" "):
            os.makedirs(f"{self.output}/{tissue}/meta", exist_ok=True)
            with open(f"{self.output}/{tissue}/meta/gte.txt", "w") as gte:
                df_list = []
                for study in self.file.tissues[tissue].studies:
                    if len(self.file.tissues[tissue].studies[study].samples) >= 30:
                        # Check if study has expression file
                        if os.path.exists(f"{self.output}/{tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz"):
                            with gzip.open(f"{self.output}/{tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", "r") as exp:
                                df = pd.read_csv(exp, sep="\t")
                                df['-'] = df['-'].str.split(".").str[0]
                                df_list.append(df)

                        # Create gte file
                            for sample in self.file.tissues[tissue].studies[study].samples:
                                gte.write(f"{sample}\t{sample}\t{study}\n")
                
                if len(df_list) > 0:
                    df = reduce(lambda x, y: pd.merge(x, y, on="-", how='outer'), df_list)
                    df.fillna(0, inplace=True)
                    df.to_csv(f"{self.output}/{tissue}/meta/{tissue}-exp.txt.gz", compression='gzip', sep="\t", index=False)
                    df.to_csv(f"{self.output}/{tissue}/meta/{tissue}-genes.txt.gz", columns = "-", header = False, index=False, compression='gzip')

    def do_meta_prep_selective(self, line):
        meta_tissue = line.split(" ")[0]
        tissues = line.split(" ")[1:]
        os.makedirs(f"{self.output}/{meta_tissue}/meta", exist_ok=True)
        with open(f"{self.output}/{meta_tissue}/meta/gte.txt", "w") as gte:
            df_list = []
            for tissue in tissues:
                for study in os.listdir(f"{self.output}/{tissue}/"):
                    if os.path.exists(f"{self.output}/{tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz"):
                            with gzip.open(f"{self.output}/{tissue}/{study}/exp/10_final_covariate_correction/{study}_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz", "r") as exp:
                                df = pd.read_csv(exp, sep="\t")
                                df['-'] = df['-'].str.split(".").str[0]
                                df_list.append(df)
                    for sample in os.listdir(f"{self.output}/{tissue}/{study}/"):
                        if sample not in self.dir_blacklist:
                            gte.write(f"{sample}\t{sample}\t{study}\n")
                            
            if len(df_list) > 0:
                df = reduce(lambda x, y: pd.merge(x, y, on="-", how='outer'), df_list)
                df.fillna(0, inplace=True)
                df.to_csv(f"{self.output}/{meta_tissue}/meta/{meta_tissue}-exp.txt.gz", compression='gzip', sep="\t", index=False)
                df.to_csv(f"{self.output}/{meta_tissue}/meta/{meta_tissue}-genes.txt.gz", columns = "-", header = False, index=False, compression='gzip')


    def do_meta(self, line):
        for tissue in line.split(" "):
            print(f"Creating meta_configs directory at {self.work_dir}/meta_configs")
            os.makedirs(f"{self.work_dir}/meta_configs", exist_ok=True)
            print(f"Creating {tissue}.config")
            with open (f"{self.work_dir}/meta_configs/{tissue}.config", "w") as study_config, open(f"{self.pipeline_dir}/configs/meta.config", "r") as config:
                for line in config:
                    if "out_dir" in line:
                        study_config.write(f"\tout_dir = \"{self.output}/{tissue}/meta\"\n")
                        continue
                    if "in_dir" in line:
                        study_config.write(f"\tin_dir = \"{self.output}/{tissue}/genotypes\"\n")
                        continue
                    if "name" in line:
                        study_config.write(f"\nname = \"{tissue}\"\n")
                        continue
                    if "exp" in line:
                        study_config.write(f"\nexp = \"{self.output}/{tissue}/meta/{tissue}-exp.txt.gz\"\n")
                        continue
                    if "gte" in line:
                        study_config.write(f"\ngte = \"{self.output}/{tissue}/meta/gte.txt\"\n")
                        continue
                    if "genelist" in line:
                        study_config.write(f"\ngenelist = \"{self.output}/{tissue}/meta/{tissue}-genes.txt.gz\"\n")
                        continue
                    study_config.write(line)
            
                print(f"Creating work directory for {tissue} at {self.work_dir}/meta/{tissue}")
                os.makedirs(f"{self.work_dir}/meta/{tissue}", exist_ok=True)
                os.system(f"tmux new -d -s meta_{tissue} -c {self.work_dir}/meta/{tissue}; tmux send-keys -t meta_{tissue} \'ml Java && {self.pipeline_dir}/nextflow run {self.pipeline_dir}/modules/meta.nf -c {self.work_dir}/meta_configs/{tissue}.config\' Enter")
        
    def do_gt_merge(self, line):
        if len(line.split(" ")) >= 1:
            for tissue in line.split(" "):
                # for chr in range(1, 23):
                #     os.system(f"realpath {self.output}{tissue}/*/genotypes/chr{chr}.filtered.dose.vcf.gz >> {self.output}{tissue}/genotypes/chr{chr}-list.txt")
                #     os.system(f"sbatch --export=CHR={chr},IN_PATH={self.output}{tissue} ./modules/bin/combine.sh")

                os.system(f"for x in {{1..22}}; do realpath {self.output}{tissue}/*/impute/postimpute/chr$x.filtered.dose.vcf.gz >> {self.output}{tissue}/genotypes/chr$x-list.txt; done")
                os.system(f"for x in {{1..22}}; do sbatch --export=CHR=$x,IN_PATH={self.output}{tissue} ./modules/bin/combine.sh; done")
    def do_gt_merge_selective(self, line):
        meta_tissue = line.split(" ")[0]
        os.makedirs(f"{self.output}/{meta_tissue}/genotypes", exist_ok=True)
        if len(line.split(" ")) >= 1:
            for tissue in line.split(" ")[1:]:
                    os.system(f"for x in {{1..22}}; do realpath {self.output}{tissue}/*/impute/postimpute/chr$x.filtered.dose.vcf.gz >> {self.output}{meta_tissue}/genotypes/chr$x-list.txt; done")
            os.system(f"for x in {{1..22}}; do sbatch --export=CHR=$x,IN_PATH={self.output}{meta_tissue} ./modules/bin/combine.sh; done")

    def do_index(self, line):
        for tissue in line.split(" "):
            for study in os.listdir(f"{self.output}/{tissue}/"):
                for file in os.listdir(f"{self.output}/{tissue}/{study}/impute/postimpute/"):
                    if file.endswith(".vcf.gz"):
                        os.system(f"sbatch --export=FILE={self.output}/{tissue}/{study}/impute/postimpute/{file} ./modules/bin/index.sh")
    def do_gunzip(self, line):
        for tissue in line.split(" "):
            for study in os.listdir(f"{self.output}/{tissue}/"):
                for file in os.listdir(f"{self.output}/{tissue}/{study}/impute/postimpute/"):
                    if file.endswith(".vcf.gz"):
                        os.system(f"sbatch --export=FILE={self.output}/{tissue}/{study}/impute/postimpute/{file} ./modules/bin/gunzip.sh")

    def do_bgzip(self, line):
        for tissue in line.split(" "):
            for study in os.listdir(f"{self.output}/{tissue}/"):
                for file in os.listdir(f"{self.output}/{tissue}/{study}/impute/postimpute/"):
                    if file.endswith(".vcf"):
                        os.system(f"sbatch --export=FILE={self.output}/{tissue}/{study}/impute/postimpute/{file} ./modules/bin/bgzip.sh")
    
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
