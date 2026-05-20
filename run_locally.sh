#!/bin/bash -l

#SBATCH -A naiss2025-22-471
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 0-05:00:00
#SBATCH -J nf-varcall
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

## load nextflow module
#ml nextflow

# only used for testing!


nextflow run main.nf --input ../testfiles/input.csv -profile standard -resume \
    --reference ../testfiles/reference/GCF_003339765.1_Mmul_10_4chroms.fa \
    --reads_dir ../testfiles/reads/supersmall \
    --scaffold_list ../testfiles/test_4chroms.txt \
    --variant_caller freebayes \
    --outdir test/freebayes_new
