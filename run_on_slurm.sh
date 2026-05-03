#!/bin/bash -l

#SBATCH -A naiss2025-22-471
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 0-05:00:00
#SBATCH -J nf-varcall
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load nextflow module
ml nextflow

# nextflow run main.nf --input testfiles/input.csv -profile dardel -resume
nextflow run main.nf --input ../mapping_testfiles/input_2.csv -profile dardel_test -resume \
    --reference ../mapping_testfiles/reference/GCF_003339765.1_Mmul_10_4chroms.fa \
    --reads_dir ../mapping_testfiles/reads \
    --scaffold_list ../mapping_testfiles/scaffolds.txt \
    --name ../map_and_call_test/test_2 \
    --slurm_account naiss2025-22-471