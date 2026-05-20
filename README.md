# map_and_varcall

Mapping and variant calling pipeline developed to handle everything from raw fastq read input files, to filtered SNPs and indels ready for analysis. The pipeline is built using Nextflow and handles all dependencies internally using conda environments. It was designed specifically for smooth running on Dardel (the compute cluster currently most frequenctly used at the Swedish Museum of Natural History), but should be easy enough to adapt to other environments.

# Quick start on dardel

1. Clone the repository to a suitable place in your dardel project, and navigate to the directory:

    ```
    git clone XZY
    # navigate to the directory
    cd map_and_call
    ```
    
2. Prepare an input samplesheet with one row per sequence pair, and up to five columns:
    sample_id;library;datatype;fastq1.fq.gz;fastq2.fq.gz,
    where datatype is either 1 for modern sequencing data, or 2 for historical dna (expecting shorter reads and more damage). If one sample was sequenced more than once (several lanes or libraries), simply add multiple rows for the same sample_id, with the different fastq files, and the pipeline will handle merging per sample after mapping. An example of a samplesheet row:
    Sample1;lib1;1;reads/Sample1_lane1_R1.fastq.gz;reads/Sample1_lane1_R2.fastq.gz

    (library and datatype columns are optional, will default to a unique library per sample row and datatype 1, respectively, if not omitted).

    Paths to input reads may either be absolute path, or – preferable – just the basenames of the reads together with the --reads_dir argument when running the pipeline (keeps the samplesheet cleaner). If using --reads_dir, make sure that all reads are accessible in the same directory, for example by symlinking all reads to a common directory:

        mkdir reads
        for read in $(find /dir/with/raw_data -name "*.fq.gz");
        do
            ln -s $read reads/
        done


3. Edit the relevant parameters in the dardel_wrapper.sh slurm script:
    - Project ID to use for submitting jobs to slurm
    - Path to the reference genome (plain or gzipped) to be used for mapping and variant calling
    - Path to where to store the results
    - (A bunch of more "advanced" parameters are customizable in the nextflow.config file)

4. Submit the pipeline to slurm using the dardel_wrapper.sh script:

    sbatch dardel_wrapper.sh


# Output

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

## 02_maskfiles

Contains three bedfiles per sample:
- <sample_id>_mappability_mask.bed
Callability mask across the genome, that is, this file contains genomic regions where we're confident in our ability to call genotypes
- <sample_id>_homref_invariants.bed
This file contains all intervals in the reference genomes with sufficient read coverage for variant calling, but where no variants were called. That is, one can assume that these sites are homozygous reference for the particular sample.
- <sample_id>_snp_mask.bed
Callability mask for SNPs: that is, this file contains genomic regions where we're confident in our ability to call SNPs if present. Any sites with indels will be excluded in this file.


