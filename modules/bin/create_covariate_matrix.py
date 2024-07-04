#!/usr/bin/env python 

import pandas as pd
import gzip
import os
import argparse


def get_alignment_metrics(file_path, res_columns, res_values):

    included_metrics = [
            "ALIGNMENT_METRICS_BAD_CYCLES",
            "ALIGNMENT_METRICS_MEAN_READ_LENGTH",
            "ALIGNMENT_METRICS_PCT_ADAPTER",
            "ALIGNMENT_METRICS_PCT_CHIMERAS",
            "ALIGNMENT_METRICS_PCT_PF_READS",
            "ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED",
            "ALIGNMENT_METRICS_PCT_PF_READS_IMPROPER_PAIRS",
            "ALIGNMENT_METRICS_PCT_READS_ALIGNED_IN_PAIRS",
            "ALIGNMENT_METRICS_PF_ALIGNED_BASES",
            "ALIGNMENT_METRICS_PF_HQ_ALIGNED_BASES",
            "ALIGNMENT_METRICS_PF_HQ_ALIGNED_Q20_BASES",
            "ALIGNMENT_METRICS_PF_HQ_ALIGNED_READS",
            "ALIGNMENT_METRICS_PF_HQ_ERROR_RATE",
            "ALIGNMENT_METRICS_PF_HQ_MEDIAN_MISMATCHES",
            "ALIGNMENT_METRICS_PF_INDEL_RATE",
            "ALIGNMENT_METRICS_PF_MISMATCH_RATE",
            "ALIGNMENT_METRICS_PF_NOISE_READS",
            "ALIGNMENT_METRICS_PF_READS",
            "ALIGNMENT_METRICS_PF_READS_ALIGNED",
            "ALIGNMENT_METRICS_PF_READS_IMPROPER_PAIRS",
            "ALIGNMENT_METRICS_READS_ALIGNED_IN_PAIRS",
            "ALIGNMENT_METRICS_STRAND_BALANCE",
            "ALIGNMENT_METRICS_TOTAL_READS",
        ]

    with gzip.open(file_path, 'r') as alignment_metrics:
        content = alignment_metrics.read().decode()
        data = content.split('\n\n')[1].split('\n')
        columns = data[1].split('\t')
        values = data[2].split('\t')


    for i, c in enumerate(columns):
        
        if f'ALIGNMENT_METRICS_{c}' in included_metrics:
            res_columns.append(f'ALIGNMENT_METRICS_{c}')
            res_values.append(values[i])

    return res_columns, res_values

def get_insert_metrics(file_path, res_columns, res_values):

    included_metrics = [
        "INSERT_METRICS_MAX_INSERT_SIZE",
        "INSERT_METRICS_MEAN_INSERT_SIZE",
        "INSERT_METRICS_MEDIAN_ABSOLUTE_DEVIATION",
        "INSERT_METRICS_MEDIAN_INSERT_SIZE",
        "INSERT_METRICS_MIN_INSERT_SIZE",
        "INSERT_METRICS_MODE_INSERT_SIZE",
        "INSERT_METRICS_PAIR_ORIENTATION",
        "INSERT_METRICS_READ_PAIRS",
        "INSERT_METRICS_STANDARD_DEVIATION",
        "INSERT_METRICS_WIDTH_OF_10_PERCENT",
        "INSERT_METRICS_WIDTH_OF_20_PERCENT",
        "INSERT_METRICS_WIDTH_OF_30_PERCENT",
        "INSERT_METRICS_WIDTH_OF_40_PERCENT",
        "INSERT_METRICS_WIDTH_OF_50_PERCENT",
        "INSERT_METRICS_WIDTH_OF_60_PERCENT",
        "INSERT_METRICS_WIDTH_OF_70_PERCENT",
        "INSERT_METRICS_WIDTH_OF_80_PERCENT",
        "INSERT_METRICS_WIDTH_OF_90_PERCENT",
        "INSERT_METRICS_WIDTH_OF_95_PERCENT",
        "INSERT_METRICS_WIDTH_OF_99_PERCENT",
    ]

    if os.path.exists(f'{file_path}'):
        with gzip.open(file_path, 'r') as insert_metrics:
            content = insert_metrics.read().decode()
            data = content.split('\n\n')[1].split('\n')
            columns = data[1].split('\t')
            values = data[2].split('\t')
        for i, c in enumerate(columns):
            
            if f'INSERT_METRICS_{c}' in included_metrics:
                res_columns.append(f'INSERT_METRICS_{c}')
                res_values.append(values[i])       
    else:
        for i, c in enumerate(included_metrics):
                res_columns.append(f'{c}')
                res_values.append(0)    

    return res_columns, res_values 

def get_rna_seq_metrics(file_path, res_columns, res_values):
    included_metrics = [
        "RNASEQ_METRICS_CODING_BASES",
        "RNASEQ_METRICS_CORRECT_STRAND_READS",
        "RNASEQ_METRICS_IGNORED_READS",
        "RNASEQ_METRICS_INCORRECT_STRAND_READS",
        "RNASEQ_METRICS_INTERGENIC_BASES",
        "RNASEQ_METRICS_INTRONIC_BASES",
        "RNASEQ_METRICS_MEDIAN_3PRIME_BIAS",
        "RNASEQ_METRICS_MEDIAN_5PRIME_BIAS",
        "RNASEQ_METRICS_MEDIAN_5PRIME_TO_3PRIME_BIAS",
        "RNASEQ_METRICS_MEDIAN_CV_COVERAGE",
        "RNASEQ_METRICS_NUM_R1_TRANSCRIPT_STRAND_READS",
        "RNASEQ_METRICS_NUM_R2_TRANSCRIPT_STRAND_READS",
        "RNASEQ_METRICS_NUM_UNEXPLAINED_READS",
        "RNASEQ_METRICS_PCT_CODING_BASES",
        "RNASEQ_METRICS_PCT_CORRECT_STRAND_READS",
        "RNASEQ_METRICS_PCT_INTERGENIC_BASES",
        "RNASEQ_METRICS_PCT_INTRONIC_BASES",
        "RNASEQ_METRICS_PCT_MRNA_BASES",
        "RNASEQ_METRICS_PCT_R1_TRANSCRIPT_STRAND_READS",
        "RNASEQ_METRICS_PCT_R2_TRANSCRIPT_STRAND_READS",
        "RNASEQ_METRICS_PCT_RIBOSOMAL_BASES",
        "RNASEQ_METRICS_PCT_USABLE_BASES",
        "RNASEQ_METRICS_PCT_UTR_BASES",
        "RNASEQ_METRICS_PF_ALIGNED_BASES",
        "RNASEQ_METRICS_PF_BASES",
        "RNASEQ_METRICS_RIBOSOMAL_BASES",
        "RNASEQ_METRICS_UTR_BASES",
    ]


    with gzip.open(file_path, 'r') as rnaseq_metrics:
        content = rnaseq_metrics.read().decode()
        data = content.split('\n\n')[1].split('\n')
        columns = data[1].split('\t')
        values = data[2].split('\t')
        
    for i, c in enumerate(columns):
        
        if f'RNASEQ_METRICS_{c}' in included_metrics:
            res_columns.append(f'RNASEQ_METRICS_{c}')
            res_values.append(values[i])
    
    return res_columns, res_values

def get_star_log_metrics(file_path, res_columns, res_values):
    values = []

    included_metrics = [
        "STAR_LOG_Average_input_read_length",
        "STAR_LOG_Average_mapped_length",
        "STAR_LOG_Deletion_average_length",
        "STAR_LOG_Deletion_rate_per_base",
        "STAR_LOG_Insertion_average_length",
        "STAR_LOG_Insertion_rate_per_base",
        "STAR_LOG_Mismatch_rate_per_base_PCT",
        "STAR_LOG_Number_of_chimeric_reads",
        "STAR_LOG_Number_of_input_reads",
        "STAR_LOG_Number_of_reads_mapped_to_multiple_loci",
        "STAR_LOG_Number_of_reads_mapped_to_too_many_loci",
        "STAR_LOG_Number_of_splices_AT/AC",
        "STAR_LOG_Number_of_splices_Annotated_(sjdb)",
        "STAR_LOG_Number_of_splices_GC/AG",
        "STAR_LOG_Number_of_splices_GT/AG",
        "STAR_LOG_Number_of_splices_Non-canonical",
        "STAR_LOG_Number_of_splices_Total",
        "STAR_LOG_PCT_of_chimeric_reads",
        "STAR_LOG_PCT_of_reads_mapped_to_multiple_loci",
        "STAR_LOG_PCT_of_reads_mapped_to_too_many_loci",
        "STAR_LOG_PCT_of_reads_unmapped_other",
        "STAR_LOG_PCT_of_reads_unmapped_too_many_mismatches",
        "STAR_LOG_PCT_of_reads_unmapped_too_short",
        "STAR_LOG_Uniquely_mapped_reads_PCT",
        "STAR_LOG_Uniquely_mapped_reads_number"
    ]


    with gzip.open(file_path, 'r') as star_log_metrics:
        content = star_log_metrics.read().decode()
        data = content.split('\n\n')[1].split('\n')


    for line in data:
        metric = 'STAR_LOG_' + line.strip().split('|')[0].replace('%', 'PCT').replace(',', '').replace(':', '').strip().replace(' ', '_')

        try:
            value = line.strip().split('|')[1].strip().replace('%', '')
        except:
            continue
        
        if metric in included_metrics:
            res_columns.append(metric)
            res_values.append(value)
        
    return res_columns, res_values


def get_star_count_metrics(file_path, res_columns, res_values):
    fh = gzip.open(file_path, 'r')

    sampledata = {}
    mapped = [0,0,0]

    for line in fh:
            line = line.decode()
            elems = line.strip().split("\t")
            phen = elems[0]
            if not phen.startswith("ENSG"):
                for i in range(1,4):
                    ct = elems[i]
                    phenoname = elems[0]
                    if i == 1:
                        phenoname = "STAR_TAB_"+phenoname + "_sum"
                    elif i == 2:
                        phenoname = "STAR_TAB_"+phenoname + "_strandA"
                    else:
                        phenoname = "STAR_TAB_"+phenoname + "_strandB"
                    sampledata[phenoname] = ct
                    res_columns.append(phenoname)
                    res_values.append(ct)
            else:
                # sum up the counts for each column
                for i in range(1,4):
                    ct = int(elems[i])
                    mapped[i-1] += ct

    for i in range(0,len(mapped)): 
        if i == 0:
            phenoname = "STAR_TAB_N_mapped_sum"
        elif i == 1:
            phenoname = "STAR_TAB_N_mapped_strandA"
        else:
            phenoname = "STAR_TAB_N_mapped_strandB"
        sampledata[phenoname] = str(mapped[i])
        res_columns.append(phenoname)
        res_values.append(str(mapped[i]))

    fh.close()

    return res_columns, res_values


def get_all_metrics(sample, input_dir):
    columns = []
    values = []
    columns, values = get_alignment_metrics(f'{input_dir}/multiple_metrics/{sample}.alignment_summary_metrics.gz', columns, values)
    columns, values = get_insert_metrics(f'{input_dir}/multiple_metrics/{sample}.insert_size_metrics.gz', columns, values)
    columns, values = get_rna_seq_metrics(f'{input_dir}/rna_seq_metrics/{sample}_rnaseqmetrics.gz', columns, values)
    columns, values = get_star_log_metrics(f'{input_dir}/star/{sample}Log.final.out.gz', columns, values)
    columns, values = get_star_count_metrics(f'{input_dir}/star/{sample}ReadsPerGene.out.tab.gz', columns, values)
    return columns, values

if __name__ == '__main__':

    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", help="The output directory resulting from the alignment pipeline (REQUIRED)")
    parser.add_argument("--cohort_name", help="Name of the cohort (REQUIRED)")
    args = parser.parse_args()

    # Get input directory
    input_directory = args.input_dir
    if input_directory == None:
        raise ValueError('Please provide a valid input directory')

    # Get cohort name
    cohort_name = args.cohort_name
    if cohort_name == None:
        raise ValueError('Please provide a cohort name')

    columns = []
    values = []
    samples = []

    # Open error file
    error_file = open('failed_samples.txt', 'w')

    for sample in os.listdir(f'{input_directory}/'):
        if "QC" not in sample and "genotypes" not in sample and "exp" not in sample and "qc" not in sample:
                col, val = get_all_metrics(sample, input_directory + '/' + sample)
                columns = col
                values.append(val)
                samples.append(sample)
            # get_all_metrics(sample, input_directory)
            # print(sample)
            # try:
                # col, val = get_all_metrics(sample, input_directory + '/' + sample)
                # columns = col
                # values.append(val)
                # samples.append(sample)
            # except:
            #     error_file.write(f'Error processing sample {sample}\n')
            #     print(f'Error processing sample {sample}')

    # Close error file
    error_file.close()

    # Write out matrix to file
    df = pd.DataFrame(values, columns=columns)
    df = df.transpose()
    df.columns = samples
    df.to_csv(f'{cohort_name}_covariates.txt', sep='\t')