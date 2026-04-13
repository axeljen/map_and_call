#!/bin/bash -l

#SBATCH -A naiss2026-4-127
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH -J nf-varcall
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# write your code here

# ml nextflow

# nextflow run main.nf --input testfiles/input.csv -profile dardel -resume
nextflow run main.nf --input testfiles/input.csv -profile standard -resume