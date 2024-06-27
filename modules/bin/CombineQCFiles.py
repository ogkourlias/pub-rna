#!/usr/bin/env python3

"""
Author: Harm-Jan Westra, Modified by: Peter Riesebos
"""

import glob
import sys
import gzip
import math
import os


if len(sys.argv) < 3:
    print("Usage: outputDir outfile.txt.gz")
    sys.exit(-1)

outputDir = sys.argv[1]
outfile = sys.argv[2]

def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file,'rt')
    else:
        return open(file,'r')

def getfho(file):
    if file.endswith(".gz"):
        return gzip.open(file,'wt')
    else:
        return open(file,'w')

def getSample(file,removeStr):
    sample = file.split("/")[-1]
    sample = sample.replace(removeStr,"")
    sample = sample.replace(".gz","")
    return sample    

def getSamples(files,removeStr):
    print("{} files".format(len(files)))
    samples = set()
    for file in files:
        sample = getSample(file,removeStr)
        samples.add(sample)
    nrsamples = len(samples)
    print(f"{nrsamples} samples with {removeStr}")
    return samples

def toSortedArr(s):
    d = []
    for i in s:
        d.append(i)
    d.sort()
    return d

def getSamplesFromDir(files):
    samples = set()
    for filename in files:
        sample = getSampleFromDir(filename)
        samples.add(sample)
    return samples

def getSampleFromDir(filename):
    sample = filename.split("/")[-2]
    return sample

# inventorize samples
files = glob.glob(outputDir+"/multiple_metrics/{SRR*,ERR*}/multiple_metrics.alignment_summary_metrics.gz")
samples = getSamplesFromDir(files)
print("{} samples loaded sofar ".format(len(samples)))

files = glob.glob(outputDir+"/multiple_metrics/{SRR*,ERR*}/multiple_metrics.insert_size_metrics.gz")
samples = getSamplesFromDir(files)
print("{} samples loaded sofar ".format(len(samples)))

files = glob.glob(outputDir+"/rna_seq_metrics/*_rnaseqmetrics.gz")
samples.update(getSamples(files,"_rnaseqmetrics.gz"))
print("{} samples loaded sofar ".format(len(samples)))

files = glob.glob(outputDir+"/star/SRR*/*_ReadsPerGene.out.tab.gz")
samples.update(getSamples(files,"_ReadsPerGene.out.tab.gz"))
print("{} samples loaded sofar ".format(len(samples)))

files = glob.glob(outputDir+"/star/SRR*/*_Log.final.out.gz")
samples.update(getSamples(files,"_Log.final.out.gz"))
print("{} samples loaded total ".format(len(samples)))



# collection bins
data = {}
metrics = set()

## collectMultipleMetrics_QC
# alignment_summary_metrics - skip lines w/ #, skip empty lines
# header line starts with CATEGORY
# line with PAIR has average numbers
files = []
for pattern in ["SRR*", "ERR*"]:
    files.extend(glob.glob(os.path.join(outputDir, "multiple_metrics", pattern + "/multiple_metrics.alignment_summary_metrics.gz")))
for file in files:
    print("Parsing: "+file)
    sample = getSampleFromDir(file)
    # print(file+"\t"+sample)
    # sys.exit(0)
    fh = getfh(file)
    header = []
    sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("CATEGORY"):
            elems = line.strip().split("\t")
            for i in range(1,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    header.append("ALIGNMENT_METRICS_"+elems[i])
                    metrics.add("ALIGNMENT_METRICS_"+elems[i])
        elif line.startswith("PAIR"):
            elems = line.strip().split("\t")
            for i in range(1,len(elems)):
                sampledata[header[i-1]] = elems[i]
            break # no need to read the rest
        elif line.startswith("UNPAIRED"):
            elems = line.strip().split("\t")
            for i in range(1,len(elems)):
                sampledata[header[i-1]] = elems[i]
            break # no need to read the rest
    fh.close()    
    data[sample] = sampledata        

# insert_size_metrics - skip lines w/ #, skip empty lines; first line begins with; second line has stats 
# MEDIAN_INSERT_SIZE
files = []
for pattern in ["SRR*", "ERR*"]:
    files.extend(glob.glob(os.path.join(outputDir, "multiple_metrics", pattern + "/multiple_metrics.insert_size_metrics.gz")))
for file in files:
    print("Parsing: "+file)
    sample = getSampleFromDir(file)
    fh = getfh(file)
    header = []
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("MEDIAN_INSERT_SIZE"):
            elems = line.strip().split("\t")
            for i in range(0,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    dta = sampledata.get(col)
                    col = "INSERT_METRICS_"+col
                    header.append(col)
                    metrics.add(col)

            dataline = fh.readline()
            elems = dataline.strip().split("\t")
            for i in range(0,len(elems)):
                sampledata[header[i]] = elems[i]
            break
    fh.close()
    data[sample] = sampledata

## collectRnaSeqMetrics_QC
# *.rna_metrics.log
# header starts with ## METRICS CLASS; next line has column names, starting with PF_BASES
# line afterwards values
files = []
for pattern in ["SRR*", "ERR*"]:
    files.extend(glob.glob(os.path.join(outputDir, "rna_seq_metrics", pattern + "_rnaseqmetrics.gz")))
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,"_rnaseqmetrics.gz")
    fh = getfh(file)
    header = []
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("PF_BASES"):
            elems = line.strip().split("\t")
            for i in range(0,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    dta = sampledata.get(col)
                    col = "RNASEQ_METRICS_"+col
                    header.append(col)
                    metrics.add(col)

            dataline = fh.readline()
            elems = dataline.strip().split("\t")

            for i in range(0,len(elems)):
                sampledata[header[i]] = elems[i]
            break
    fh.close()
    data[sample] = sampledata

# STAR count metrics
files = []
for pattern in ["SRR*", "ERR*"]:
    files.extend(glob.glob(os.path.join(outputDir, "star", pattern + "/*_ReadsPerGene.out.tab.gz")))
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,"_ReadsPerGene.out.tab.gz")
    samples.add(sample)
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}
    fh = getfh(file)
    mapped = [0,0,0]
    
    for line in fh:
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
                metrics.add(phenoname)
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
            metrics.add(phenoname)

    fh.close()
    data[sample] = sampledata


# STAR final logs
allowedp = [
    "Number_of_input_reads",
    "Average_input_read_length",
    "Uniquely_mapped_reads_number",
    "Uniquely_mapped_reads_PCT",
    "Average_mapped_length",
    "Number_of_splices_Total",
    "Number_of_splices_Annotated_(sjdb)",
    "Number_of_splices_GT/AG",
    "Number_of_splices_GC/AG",
    "Number_of_splices_AT/AC",
    "Number_of_splices_Non-canonical",
    "Mismatch_rate_per_base_PCT",
    "Deletion_rate_per_base",
    "Deletion_average_length",
    "Insertion_rate_per_base",
    "Insertion_average_length",
    "Number_of_reads_mapped_to_multiple_loci",
    "PCT_of_reads_mapped_to_multiple_loci",
    "Number_of_reads_mapped_to_too_many_loci",
    "PCT_of_reads_mapped_to_too_many_loci",
    "PCT_of_reads_unmapped_too_many_mismatches",
    "PCT_of_reads_unmapped_too_short",
    "PCT_of_reads_unmapped_other",
    "Number_of_chimeric_reads",
    "PCT_of_chimeric_reads"
]
allowedpset = set()
for p in allowedp:
    allowedpset.add(p)

files = []
for pattern in ["SRR*", "ERR*"]:
    files.extend(glob.glob(os.path.join(outputDir, "star", pattern + "/*_Log.final.out.gz")))
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,"_Log.final.out.gz")
    samples.add(sample)
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    fh = getfh(file)
    for line in fh:
        line = line.strip()
        if "|" in line:
            line = line.replace("|","").strip()
            
            elems = line.split("\t")
            if len(elems) == 2:
                p = elems[0].strip().replace(" ","_")
                p = p.replace(":","")
                p = p.replace(",","")
                p = p.replace("%","PCT")
                # print(p)
                if p in allowedpset:
                    v = elems[1].replace("%","")
                    p = "STAR_LOG_"+p
                    sampledata[p] = v
                    metrics.add(p)
    fh.close()
    # sys.exit()
    data[sample] = sampledata


variant_count = glob.glob(outputDir+"/qc_logs/variant_count.txt")[0]
smiss_file_pre_filter = glob.glob(outputDir+"/PRE_FILTER.smiss")[0]
smiss_file_post_filter = glob.glob(outputDir+"/POST_FILTER.smiss")[0]

# get variant count for pre and post filter
try:
    with open(variant_count, 'r') as f:
        lines = f.readlines()
        first_row = lines[0].strip().split('\t')
        variant_count_pre = first_row[1]

        second_row = lines[1].strip().split('\t')
        variant_count_post = second_row[1]

except FileNotFoundError:
    print(f"Warning: variant_count file not found - {variant_count}")
    variant_count_pre, variant_count_post = 0

# Add total_lines as a metric for each sample
for sample in samples:
    sampledata = data.get(sample, {})
    sampledata['TOTAL_VARIANTS'] = str(variant_count_pre)
    sampledata['TOTAL_FILTERED_VARIANTS'] = str(variant_count_post)
    metrics.add('TOTAL_VARIANTS')
    metrics.add('TOTAL_FILTERED_VARIANTS')
    data[sample] = sampledata


def calculate_missingness(imiss_file, variant_count, metric_name, data_key, column_index):
    try:
        with getfh(imiss_file) as imiss_fh:
            next(imiss_fh)  # Skip the first row (headers)
            for line in imiss_fh:
                elems = line.strip().split()
                sample = elems[column_index].replace('-splitreads', '')
                missing_value = elems[-1].strip()
                sampledata = data.get(sample, {})
                sampledata[metric_name] = str(missing_value)
                metrics.add(metric_name)
                data[sample] = sampledata
                additional_row_value = (1 - float(missing_value)) * float(variant_count)
                rounded_down_value = math.floor(additional_row_value)
                sampledata[f'VARIANTS_WITH_GENOTYPE_CALL_{data_key}'] = str(rounded_down_value)
                metrics.add(f'VARIANTS_WITH_GENOTYPE_CALL_{data_key}')
                data[sample] = sampledata
    except FileNotFoundError:
        print(f"Warning: .smiss file not found - {imiss_file}")

# Calculate missingness and the additional row for each sample
calculate_missingness(smiss_file_pre_filter, variant_count_pre, 'MISSINGNESS_PRE_FILTER', 'PRE_FILTER', 0)
# Calculate missingness and the additional row for the filtered VCF file
calculate_missingness(smiss_file_post_filter, variant_count_post, 'MISSINGNESS_AFTER_FILTER', 'AFTER_FILTER', 1)

# currently hardcoded, needs to be changed
samples_txt_file = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/samplesWithPrediction_16_09_22_noOutliers.txt"

try:
    with open(samples_txt_file, 'r') as samples_txt_fh:
        next(samples_txt_fh)
        for line in samples_txt_fh:
            elems = line.strip().split('\t')
            if len(elems) > 2:  # Assuming at least three columns: sample, study, sra.library_layout
                sample = elems[0]
                study = elems[5]
                layout = elems[6]

                if sample in data:
                    sampledata = data[sample]
                    sampledata['STUDY'] = study
                    sampledata['LAYOUT'] = layout
                    metrics.add('STUDY')
                    metrics.add('LAYOUT')
                    data[sample] = sampledata
except FileNotFoundError:
    print(f"Warning: samples.txt file not found - {samples_txt_file}")


print("{} metrics loaded ".format(len(metrics)))
print("{} samples loaded ".format(len(samples)))

metrics = toSortedArr(metrics)
samples = toSortedArr(samples)

# for sample in samples:
#     print("{}".format(sample))

print("Writing: "+outfile)
fho = getfho(outfile)
header = "-\t"+"\t".join(samples)+"\n"
fho.write(header)
for metric in metrics:
    outln = metric
    for sample in samples:
        sampledata = data.get(sample)
        if sample is None:
            outln+="\tnan"
        else:
            v = sampledata.get(metric)
            if v is None:
                print("metric {} is missing for sample: {}".format(metric, sample))
                # sys.exit(-1)
                outln +="\tnan"
            else:
                outln +="\t"+v
    
    fho.write(outln+"\n")
fho.close()