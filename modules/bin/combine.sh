#!/bin/bash
#for x in {1..22}; do sbatch --export=CHR=$x combine.sh; done
#SBATCH --job-name=gt_merge
#SBATCH --error=logs/gt_merge.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --output=./logs/%j.out

ml BCFtools
ml HTSlib
bcftools merge -m none --force-samples -l ${IN_PATH}/genotypes/chr${CHR}-list.txt -Oz -o ${IN_PATH}/genotypes/chr${CHR}.vcf.gz
tabix ${IN_PATH}/genotypes/chr${CHR}.vcf.gz