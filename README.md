# map\_and\_call

Mapping and variant calling pipeline developed to handle everything from raw
fastq read input files, to filtered SNPs and indels ready for analysis. The
pipeline is built using [Nextflow](https://www.nextflow.io/) and handles all
dependencies internally using [Conda](https://conda.org/) environments. It was
designed specifically for smooth running on
[Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system/about-the-dardel-system-1.1053338)
(the compute cluster currently most frequently used at the [Swedish Museum of
Natural History](https://www.nrm.se)), but should be easy enough to adapt to
other environments.

Main steps include

- filtering fastq files using [fastp](https://github.com/opengene/fastp)
- short-read alignment using [bwa mem](https://github.com/lh3/bwa)
- variant calling using [bcftools](https://samtools.github.io/bcftools/) or
  [freebayes](https://github.com/freebayes/freebayes)
- summary statistics using [multiqc](https://seqera.io/multiqc/) and
  [qualimap](http://qualimap.conesalab.org/)
- *and more...*

## Quick start on dardel.pdc.kth.se

### 1. Clone the repository to a suitable place in your dardel project, and navigate to the directory

```
$ git clone https://github.com/axeljen/map_and_call.git
$ cd map_and_call
```

### 2. Prepare an input sample sheet with one row per sequence pair, and five columns with headers

```
sample_id;library;data_type;read_1;read_2
sample_1;lib1;1;sample_1_R1.fq.gz;sample_1_R2.fq.gz
sample_2;lib1;2;sample_2_R1.fq.gz;sample_2_R2.fq.gz
```

Where:

**sample_id** is a unique identifier for each sample.

**library** is used to differentiate between different libraries sequenced from
the same sample. These well be merged prior to deduplication. If the same
library was sequenced across different lanes, simply add one row per read pair
with the same library name, and the pipeline will handle merging per library
after mapping.

**data_type** is either `1` for modern sequencing data, or `2` for historical dna
(expecting shorter reads and more damage).

**read_1/read_2** points to the paths for the fastq files for this sequencing
run. Either specify the full path to the reads, or -- to keep the input file a
bit cleaner -- put all reads (or links to them) in a common directory, and point
to this directory with the `--reads_dir` argument when running the pipeline.
For example:

Start by creating symbolic links from all reads to a common directory:

```
$ mkdir reads
$ for read in $(find /dir/with/raw_data -name "*.fq.gz"); do
    ln -s "$read" reads/
  done
```

then use the basedir of the reads when running the pipeline, with `--reads_dir
reads`.

### 3. Edit the relevant variables in the [`run_on_dardel.sh`](run_on_dardel.sh) slurm script

As a minimum, you need to provide the CPU-account number (replace `<NAISS_COMPUTE_PROJECT>`),
and the paths to your sample sheet (`INPUT_CSV`), reference genome (`REFERENCE`).

### 4. Submit the pipeline to slurm using `sbatch`

```
$ sbatch --test-only run_on_dardel.sh
$ sbatch run_on_dardel.sh
```

## Output

If all goes well, the output directory should look something like:

    .
    ├── 00_input_data
    │   └── 00_reference_genome
    ├── 01_reports
    │   ├── 00_fastqc
    │   ├── 01_qualimap
    │   ├── 02_variantstats
    │   └── 03_damage_profiles
    ├── 02_bamfiles
    │   └── dedup_metrics
    ├── 03_genotypes
    │   ├── 00_raw_variants
    │   ├── 01_filtered_variants
    │   └── 02_maskfiles
    └── pipeline_info

### 00\_input\_data

Contains the index reference genome

### 01\_reports

Contains a number different QC reports for reads, mapped bam files and variants.

### 02\_bamfiles

Contains the final, mapped and processed bam/cram files for each sample, as
well as a deduplication metrics file for each sample.

### 03\_genotypes

#### 00\_raw\_variants

Contains the raw variants in vcf format.

#### 01\_filtered\_variants

Filtered SNPs and indels in vcf format, ready for downstream analyses.

#### 02\_maskfiles

Contains three bedfiles per sample:

- `<sample_id>_mappability_mask.bed` --- Callability mask across the genome,
  that is, this file contains genomic regions where we're confident in our
  ability to call genotypes
- `<sample_id>_homref_invariants.bed` --- This file contains all intervals in
  the reference genomes with sufficient read coverage for variant calling, but
  where no variants were called. That is, one can assume that these sites are
  homozygous reference for the particular sample.
- `<sample_id>_snp_mask.bed` --- Callability mask for SNPs: that is, this file
  contains genomic regions where we're confident in our ability to call SNPs if
  present. Any sites with indels will be excluded in this file.
