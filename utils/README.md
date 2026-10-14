This folder contains some helper scripts/utils that can be useful for converting/working with the results from map_and_call. Each script is described below.

## crams2bams.sh

If the pipeline was run with the --store_crams option, but you later end up wanting to work with bam files, this script is a quick wrapper around samtools view to convert the cram files back to bam format. Assuming that you're running the script from inside the utils folder, and have your results stored in a child directory of map_and_call, you can run the script as follows on a slurm cluster:

```bash
# set a compute account to use for slurm (this needs to be changed)
SLURM_ACCOUNT=naiss2026-4-1517  # change this to your own slurm account

# and point to the directory where the output from your map_and_call run is located
OUTPUT_DIR=/cfs/klemming/projects/supr/naiss2025-23-567/dev/mapcall_storage_usage/output_tutorial

PARTITION=shared # this works on dardel, but need to be changed for other clusters

# path to the reference genome
REFERENCE=$(find $OUTPUT_DIR/00_input_data/00_reference_genome -name "*.f*a.gz")

# path to all cramfiles
CRAMFILES=$(find $OUTPUT_DIR/02_bamfiles -name "*.cram")

# submit a slurm job for each file
for CRAMFILE in $CRAMFILES; do
    echo $CRAMFILE
    sbatch -A $SLURM_ACCOUNT -p $PARTITION -j $(basename $CRAMFILE .cram) crams2bams.sh $CRAMFILE $REFERENCE
done
```

