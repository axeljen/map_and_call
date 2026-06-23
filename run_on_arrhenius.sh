#!/bin/bash -l

#SBATCH -A <NAISS_COMPUTE_PROJECT>
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH -t 3-00:00:00
#SBATCH -J map-and-call
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# Load modules for nextflow and conda/mamba. Module names as per 2026-06-17.
module load buildtool-easybuild/5.2.1-hpca3ef7d197
module load Miniforge3/26.3.2-2
module load Nextflow/25.10.6-eb
module load Java/25.0.3-bdist

# Path to input file with sample and read information, see readme for details.
INPUT_CSV='/path/to/input.csv'

# Path to reference genome (fasta format).
REFERENCE='/path/to/reference_genome.fna'

# Path to reads directory. this argument can be omitted if full paths to reads
# are provided in the input csv (if so, remove the --reads_dir argument from
# the command below).
READS_DIR='/path/to/reads_directory'

# Optionally provide a list of scaffols upon which to perform variant calling.
# If not provided, also remove the arguemnt from the command, and the pipeline
# will be run for the entire genome.
SCAFFOLD_LIST='/path/to/scaffolds.txt'

# Variant caller to use, either 'freebayes' or 'bcftools'.
VARIANT_CALLER='freebayes'

# Chunk size in megabases for parallelizing variant calling. If many samples
# are included, this could be increased a bit to avoid too many jobs getting
# stuck in the queue perhaps.
CHUNK_SIZE=20

# Output directory (created if not present).
OUTDIR='/path/to/output_directory'

# Conda tmp/cache files. Enable to save file-disk quota. Suitable place  for
# these folders is in your project directory.
#export NXF_CONDA_CACHEDIR="/path/to/NXF_CONDA_CACHEDIR"
#export CONDA_PKGS_DIRS="/path/to/CONDA_PKGS_DIRS"
#export CONDA_ENVS_PATH="/path/to/CONDA_ENVS_PATH"
#mkdir -p "$NXF_CONDA_CACHEDIR" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_PATH"

# Run the full pipeline.
# Change to full path to main.nf if run from outside the map_and_call folder.
nextflow run main.nf \
    -profile 'arrhenius' \
    -resume \
    --use_mamba \
    --input "$INPUT_CSV" \
    --reference "$REFERENCE" \
    --reads_dir "$READS_DIR" \
    --scaffold_list "$SCAFFOLD_LIST" \
    --outdir "$OUTDIR" \
    --slurm_account "$SLURM_JOB_ACCOUNT" \
    --variant_caller "$VARIANT_CALLER" \
    --chunk_size "$CHUNK_SIZE"

##### Some additional, potentially useful flags:

## --skip_premapping_dedup
# By default, a PCR deduplication step is run first on the raw reads, and then
# again on the mapped reads in the bam file. In modern data with reasonably low
# levels of duplication rates (say, < ~ 30 %), this is most likely redundant and
# could potentially even lead to over-deduplication and loss of data. So for
# modern, high quality data, this flag is recommended.

## --conda_cachedir
# By default, the conda environments created by the piepline will be stored in
# the .envs directory in the current working directory.  If intending to run
# the pipeline several times, you can save some time and space by setting this
# to a dedicated location that you can then reuse for subseqeunt runs in
# different directories. Then, the environments doesn't need to be recreated.
# See also the NXF_CONDA_CACHEDIR variable above.

## -work-dir
# To save disk-file quota, one may try to write the work directory to a
# temporary folder.  For example: `-work-dir "$SNIC_TMP"/map_and_call`
