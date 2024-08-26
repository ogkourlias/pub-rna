# PUB-RNA
This repository contains all binaries for the RNA-Seq based genotyping and meta-analysis pipeline. The pipeline is designed to take in RNA-Seq experiment RS IDs from a CSV and perform genotyping and meta-analysis on the samples. The pipeline is designed to be run on a high-performance computing cluster and is written with Nextflow, Python and R.
## Requirements and installation

To run the pipeline, you will need to have access to a high-performance computing cluster with the following software installed:

- Nextflow
- Apptainer/Singularity
- Python 3.10

The pipeline is made with the usage of a docker image in mind. To download the image and run the pipeline, you will need to pull a SIF file from docker Hub. This may be done with the following command:

``` 
apptainer pull docker://ogkourlias/pub-rna:latest
```

Once the .sif file is downloaded, you will need to fill in the 'pub-rna.config' file with the appropriate paths to the csv file containing the RS IDs (as formatted in te test_data dir), the output directory, and the path to the work dir, like this:

``` 
csv_file = "/path/to/file.csv"
out_dir = "/path/to/outdir/"
work_dir = "/path/to/workdir/"
```

Next, there are additional configuration files that need to be filled in. They can be located in configs directory. These are all nextflow configuration files and are used to set the parameters for the pipeline. Fill in everything as indicated in the template files.

## Usage
The pipeline is controlled by a python wrapper CLI script called 'pub-rna.py'. This script is used to run the pipeline and can be run with the following command:

```
python3 pub-rna.py
```

Once the CLI is active, you may type help to see the available commands. The pipeline steps can be run in the following order:

```
gt {RS_ID1 or RS_ID2, RS_ID3, ...}
exp {RS_ID1 or RS_ID2, RS_ID3, ...}
qc {RS_ID1 or RS_ID2, RS_ID3, ...}
impute {Tissue1 or Tissue2, Tissue3, ...}
meta {RS_ID1 or RS_ID2, RS_ID3, ...}
```

Additional options are available for each command and can be seen by help. The pipeline is designed to be run in the order of the commands above, but can be run in any order if desired. The pipeline will automatically check for the presence of the necessary files and will not rerun steps that have already been completed.

## License
This project is licensed under the GPL3.0 License - see the LICENSE.md file for details