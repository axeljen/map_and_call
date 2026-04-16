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
nextflow run main.nf --input ../mapping_testfiles/input.csv -profile dardel -resume \
    --popfile ../mapping_testfiles/popfile.txt \
    --reference ../mapping_testfiles/reference/GCF_003339765.1_Mmul_10_NC_041770.1.fna \
    --reads_dir ../mapping_testfiles/reads \
    --name ../map_and_call_test/freebayes_test \
    --scaffold_list ../mapping_testfiles/scaffolds.txt \
    --slurm_account naiss2025-22-471