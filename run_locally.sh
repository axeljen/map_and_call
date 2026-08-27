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

nextflow run main.nf \
    --input ../mapping_testfiles/sample_sheet_tutorial.csv \
    -profile standard \
    --reference ../mapping_testfiles/data/reference/reference_genome.fa.gz \
    --reads_dir ../mapping_testfiles/data/reads/ \
    --scaffold_list ../mapping_testfiles/data/reference/scaffold_list.txt \
    --variant_caller freebayes \
    --outdir test/freebayes_new_sing \
    -work-dir with_singularity \
    --with_singularity 
