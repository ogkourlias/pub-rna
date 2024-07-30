#!/bin/bash
#for x in {1..22}; do sbatch --export=CHR=$x combine.sh; done
#SBATCH --job-name=tabix
#SBATCH --error=logs/tabix.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --output=./logs/%j.out

ml HTSlib
bgzip ${FILE}
