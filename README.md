# map\_and\_call

1. [Description](#description)
2. [Quick start on Dardel (TL;DR)](#quick-start-on-dardel.pdc.kth.se)
3. [Input](#input)
3. [Output](#output)
4. [Working example](#working-example)
5. [Workflow parameters](#workflow-parameters)

## Description

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

### 2. Prepare an input sample sheet

The sample sheet contains one row per sequencing read pair to include, and must
have five columns with headers: `sample_id`, `library`, `data_type`, `read_1`,
and `read_2`. Recommended format is tab-separated, but any whitespace,
semicolon, or comma should also work will also be detected automatically.

| sample_id | library | data_type | read_1              | read_2              |
|-----------|---------|-----------|----------------------|----------------------|
| sample_1  | lib1    | 1         | sample_1_R1.fq.gz    | sample_1_R2.fq.gz    |
| sample_2  | lib1    | 2         | sample_2_R1.fq.gz    | sample_2_R2.fq.gz    |

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

## Input

Input files:

1. **Reference genome** (nt sequences) in FASTA format. Can be compressed (gzip) or not.
2. **Paired-end sequencing reads** in FASTQ format. Can be compressed (gzip) or
not. Files need to have the suffix `.fq`, `.fastq`, `.fq.gz`, or `.fastq.gz`.

Furthermore, information about **sample names**, **library names**, and
**sample type** ("modern" or "historical") needs to be provided. See [Prepare
an input sample sheet](#2-prepare-an-input-sample-sheet) for details.

In addition, workflow parameters can be changed by editing the config file
([nextflow.config](nextflow.config)) or by using command-line options (see
[Workflow parameters](#workflow-parameters)).

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

## Working example

I have prepared a fully working example including a toy dataset. This can be
useful to test the pipeline in your environment, to ensure that everything
works as expected.

See also <https://github.com/axeljen/mapcall_tutorial> for more details.

To follow the working example, download the dataset:

```bash
  # download and unarchive the example data
  curl -L https://osf.io/download/d8gqw/ -o example_data.tar
  tar -xf example_data.tar
  rm example_data.tar
```

Now `example_data` contains everything you need for a complete run from reads
to variants, which should take ~10-20 minutes. See the included readme file for
instructions.

## Workflow parameters

Workflow parameters can be changed by editing the configuration file
([nextflow.config](nextflow.config)) or by using command-line options.

### Required inputs

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--input` | path | `null` | Path to the input CSV. Expected columns: `sample_id;lane;datatype;library;read_1;read_2`. |
| `--reference` | path | `null` | Path to the reference FASTA file. |

### Run naming and output

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--name` | string | `output` | Run name used for output directory and file naming. |
| `--outdir` | path | `output` | Output directory for workflow results. |

### Input-related optional parameters

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--popfile` | path | `null` | Optional population file (`sample_id;population`) for joint genotyping. |
| `--reads_dir` | path | `null` | Optional directory prepended to read paths in the input CSV. |
| `--bamfiles` | path | `null` | Text file listing BAM paths, one per line, for `call_variants.nf`. Sample IDs are extracted from read groups. |
| `--cleanreads_dir` | path | `null` | Directory containing clean reads for `read_qc.nf`. Supports glob patterns like `${sample_id}*{_R1,_R2,_1,_2,_paired,_merged}*.{fq,fastq}.gz`. |

### Workflow mode

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--mapper` | string | `bwa_mem` | Read mapper to use. Currently only `bwa_mem` is supported. |
| `--variant_caller` | string | `bcftools` | Variant caller to use. Supported values: `freebayes`, `bcftools`. |
| `--skip_variant_calling` | boolean | `false` | Stop after BAM/CRAM processing and skip variant calling and filtering. |
| `--map_historical_pairs` | boolean | `true` | Map paired-end reads from historical samples in addition to collapsed/singleton reads. |

### Preprocessing

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--min_readlength` | integer | `30` | Minimum read length after trimming. |
| `--trim_front` | integer | `0` | Bases to hard-trim from the start of reads with `fastp`. |
| `--trim_tail` | integer | `0` | Bases to hard-trim from the end of reads with `fastp`. |
| `--store_cleanreads` | boolean | `false` | Store trimmed/filtered reads in addition to final BAMs. Useful for debugging and downstream analyses of unmapped reads. |
| `--skip_premapping_dedup` | boolean | `false` | Skip pre-mapping deduplication. |
| `--premapping_dedup` | boolean | derived | Runs BBMap `clumpify` before mapping unless `--skip_premapping_dedup` is set. |
| `--skip_postmapping_dedup` | boolean | `false` | Skip post-mapping deduplication. |
| `--postmapping_dedup` | boolean | derived | Runs `samtools markdup` after mapping unless `--skip_postmapping_dedup` is set. |
| `--store_crams` | boolean | `true` | Store mapped reads in CRAM format instead of BAM. |
| `--min_mapqual` | integer | `30` | Minimum mapping quality. |
| `--min_basequal` | integer | `20` | Minimum base quality. |
| `--store_sample_fastqc` | boolean | `true` | Store FastQC reports for individual samples in addition to the MultiQC report. |

### Sex chromosome configuration

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--x_scaffolds` | string or list | `null` | Scaffold name(s) for the X chromosome. |
| `--y_scaffolds` | string or list | `null` | Scaffold name(s) for the Y chromosome. |
| `--z_scaffolds` | string or list | `null` | Scaffold name(s) for the Z chromosome. |
| `--w_scaffolds` | string or list | `null` | Scaffold name(s) for the W chromosome. |
| `--sex_assignment_lower_threshold` | float | `0.25` | Coverage ratio below which a sample is called hemizygous. |
| `--sex_assignment_upper_threshold` | float | `0.75` | Coverage ratio above which a sample is called homozygous. |

### Scaffold selection and masking

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--scaffold_list` | path | `null` | File with one scaffold per line. If unset, all scaffolds are used. |
| `--reference_mask` | BED path | `null` | BED file with regions excluded from callable sites, e.g. repeats or low-complexity regions. |

### Mapping configuration

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--modern_mapper` | string | `bwa_mem` | Mapper used for modern samples. Currently `bwa`-based. |
| `--short_reads_threshold` | integer | `70` | Reads of this length or shorter are mapped with `bwa aln`; longer reads with `bwa mem`. Ignored unless the historical mapper is split. |
| `--historical_mapper` | string | `bwa_mem` | Historical read mapper mode. Supported values include `split`, `bwa_mem`, and `bwa_aln`. |
| `--bwa_aln_flags` | string | `''` | Optional flags passed to `bwa aln` for historical/ancient samples. Ignored if `bwa_mem` is used. |

### Variant calling

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--chunk_size` | integer | `20` | Genomic interval size for parallel calling, in Mb. |
| `--ploidy` | integer | `2` | Expected ploidy. |
| `--min_base_quality` | integer | `20` | Minimum base quality used during variant calling. |
| `--store_raw_vcf` | boolean | `true` | Keep the unfiltered VCF in the output. |

### Variant filtering

Depth-related thresholds can be expressed as either:

- a float, interpreted as a fraction of the sample’s mean autosomal depth
- an integer, interpreted as an absolute read depth

If a fractional value falls below the minimum sensible depth, it is clamped automatically.

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--min_depth` | float or integer | `0.3` | Minimum depth to retain a genotype call. |
| `--max_depth` | float or integer | `2.0` | Maximum depth to retain a genotype call. |
| `--maf_threshold` | float | `0.0` | Minor allele frequency filter. `0.0` disables it. |
| `--fmiss_threshold` | float | `1.0` | Maximum fraction of missing genotypes. `1.0` disables it. |
| `--min_allele_balance` | float | `0.15` | Minimum allele balance for heterozygous calls. |
| `--min_global_dp` | integer | `null` | Minimum total read depth across all genotypes at a site to keep it. |
| `--max_global_dp` | integer | `null` | Maximum total read depth across all genotypes at a site to keep it. |
| `--snp_filter_expression` | string | `null` | Custom hard-filter expression for SNPs. If unset, caller-specific defaults are used. |
| `--indel_filter_expression` | string | `null` | Custom hard-filter expression for indels. If unset, caller-specific defaults are used. |

### Downsampling

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--downsample_bams` | boolean | `false` | Downsample BAMs before variant calling. |
| `--downsample_bams_coverage` | integer | `-1` | Target coverage for downsampling. `-1` means downsample to the shallowest sample in the dataset. |

### Historical DNA

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--damageprofiler_rescale` | boolean | `false` | Rescale base qualities using mapDamage for historical samples. |

### Advanced caller-specific defaults

These are internal defaults used when custom filter expressions are not provided.

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--snp_filter_expression_bcftools` | string | `QUAL<30 \|\| MQ<30` | Default SNP filter for `bcftools`. |
| `--indel_filter_expression_bcftools` | string | `QUAL<30 \|\| MQ<30` | Default indel filter for `bcftools`. |
| `--snp_filter_expression_freebayes` | string | `QUAL<30 \|\| MQM<30 \|\| MQMR<30` | Default SNP filter for `freebayes`. |
| `--indel_filter_expression_freebayes` | string | `QUAL<30 \|\| MQM<30 \|\| MQMR<30` | Default indel filter for `freebayes`. |
| `--snp_filter_expression_gatk` | string | `QD < 2.0 \|\| FS > 60.0 \|\| MQ < 40.0 \|\| MQRankSum < -12.5 \|\| ReadPosRankSum < -8.0` | Default SNP filter for GATK-style filtering. |
| `--indel_filter_expression_gatk` | string | `QD < 2.0 \|\| FS > 200.0 \|\| ReadPosRankSum < -20.0` | Default indel filter for GATK-style filtering. |

### HPC and execution tuning

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--slurm_account` | string | `null` | SLURM account string for cluster execution profiles. |
| `--max_memory` | string | `512.GB` | Hard memory ceiling across all processes. |
| `--max_cpus` | integer | `64` | Hard CPU ceiling across all processes. |
| `--max_time` | string | `96.h` | Hard time ceiling across all processes. |
| `--bwa_threads` | integer | `30` | CPU threads for BWA MEM alignment. |
| `--clumpify_threads` | integer | `8` | CPU threads for BBMap `clumpify` deduplication. |

### Conda and environment management

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--use_mamba` | boolean | `false` | Use Mamba instead of Conda when supported by the selected profile. |
| `--conda_cachedir` | path | `"$projectDir/.envs"` | Conda cache directory. Also used for the package cache. |

### Output and cleanup behavior

| Option | Type | Default | Description |
|---|---:|---:|---|
| `--progressive_cleanup` | boolean | `false` | Progressively clean work directories while the workflow runs. Experimental and breaks resumability. |
| `--use_scratch` | boolean | `false` | Use scratch space on compute nodes if available. |
| `--boost_cleanup` | boolean | `false` | Enable nf-boost cleanup of the work directory during execution. |
| `--boost_cleanup_interval` | duration | `180s` | Cleanup interval used by nf-boost. |

### Notes

- The workflow uses `manifest.nextflowVersion = '>=23.04.0'`.
- The default output directory is controlled by `--outdir`.
- Some parameters are derived internally, such as:
  - `--premapping_dedup`
  - `--postmapping_dedup`

  These are toggled indirectly by the corresponding `--skip_*` flags.

