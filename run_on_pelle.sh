#!/bin/bash -l

#SBATCH -A <NAISS_COMPUTE_PROJECT>
#SBATCH -p shared
#SBATCH -c 4
#SBATCH -t 5-00:00:00
#SBATCH -J map-and-call
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# load nextflow and modules
ml PDC/24.11 nextflow miniconda3

# path to input file with sample and read information, see readme for details
INPUT_CSV=/path/to/input.csv
# path to reference genome
REFERENCE=/path/to/reference_genome.fa
# path to reads directory. this argument can be omitted if full paths to reads are provided in the input csv (if so, remove the --reads_dir argument from the command below)
READS_DIR=/path/to/reads
# Optionally provide a list of scaffols upon which to perform variant calling. If not provided, also remove the arguemnt from the command, and the pipeline will be run for the entire genome
SCAFFOLD_LIST=/path/to/scaffolds.txt
# variant caller to use, either freebayes or bcftools
VARIANT_CALLER=freebayes
# chunk size in megabases for parallelizing variant calling. If many samples are included, this could be increased a bit to avoid too many jobs getting stuck in the queue perhaps.
CHUNK_SIZE=20
# output directory 
OUTDIR=path/to/output_directory
# account to use for slurm job submission
SLURM_ACCOUNT=<NAISS_COMPUTE_PROJECT>


# run the full pipeline
nextflow run main.nf --input $INPUT_CSV \
    -profile pelle -resume \
    --reference $REFERENCE \
    --reads_dir $READS_DIR \
    --scaffold_list $SCAFFOLD_LIST \
    --outdir $OUTDIR \
    --slurm_account $SLURM_ACCOUNT \
    --variant_caller $VARIANT_CALLER \
    --chunk_size $CHUNK_SIZE

##### Some additional, potentially useful flags:


## --skip_premapping_dedup
# By default, a a pcr deduplication step is run first on the raw reads, and then again on the mapped reads in the bam file.
# in modern data with reasonably low levels of duplication rates (say, < ~ 30 %), this is most likely redundant and could potentially even
# lead to over-deduplication and loss of data. So for modern, high quality data, this flag is recommended.

## --conda_cachedir
# by default, the conda environments created by the piepline will be stored in the .envs directory in the current working directory.
# If intending to run the pipeline several times, you can save some time and space by setting this to a dedicated location that
# you can then reuse for subseqeunt runs in different directories. Then, the environments doesn't need to be recreated.

## --use_mamba
# if you have a local installation of mamba, this argument will save some time when creating the conda environments. 
# mamba is not available as a module on dardel, so the default is to use conda.

