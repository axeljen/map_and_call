#!/bin/bash

#SBATCH --job-name=crams2bams
#SBATCH --output=crams2bams_logs/%j.out
#SBATCH --error=crams2bams_logs/%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

# on dardel load samtools
ml samtools

# check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found, please install."
    exit 1
fi

CRAM=$1
REFERENCE=$2

DELETE_CRAM=true

BAM="${CRAM%.cram}.bam"

samtools view -b -T "$REFERENCE" -o "$BAM" "$CRAM"
if [ "$DELETE_CRAM" = true ]; then
    samtools quickcheck "$BAM" && rm -f "$CRAM"
fi