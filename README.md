# map_and_varcall

Mapping and variant calling pipeline developed to handle everything from raw fastq read input files, to filtered SNPs and indels ready for analysis. The pipeline is built using Nextflow and handles all dependencies internally using conda environments. It was designed specifically for smooth running on Dardel (the compute cluster currently most frequenctly used at the Swedish Museum of Natural History), but should be easy enough to adapt to other environments.

# Quick start on dardel

1. Clone the repository to a suitable place in your dardel project, and navigate to the directory:

    ```
    git clone https://github.com/axeljen/map_and_call.git
    # navigate to the directory
    cd map_and_call
    ```
    
2. Prepare an input samplesheet with one row per sequence pair, and five columns with headers, like so:
    
        sample_id;library;data_type;read_1;read_2
        sample_1;lib1;1;sample_1_R1.fq.gz;sample_1_R2.fq.gz
        sample_2;lib1;2;sample_2_R1.fq.gz;sample_2_R2.fq.gz


    #### Where:  
    **sample_id** is a unique identifier for each sample.   

    **library** is used to differentiate between different libraries sequenced from the same sample. These well be merged prior to deduplication. If the same library was sequenced across different lanes, simply add one row per read pair with the same library name, and the pipeline will handle merging per library after mapping.  
    
    **data_type** is either 1 for modern sequencing data, or 2 for historical dna (expecting shorter reads and more damage).  
    
    **read_1/read_2** points to the paths for the fastq files for this sequencing run. Either specify the full path to the reads, or - to keep the input file a bit cleaner - put all reads (or links to them) in a common directory, and point to this directory with the --reads_dir argument when running the pipeline. For example:

    Symlink all reads to a common directory:


        mkdir reads
        for read in $(find /dir/with/raw_data -name "*.fq.gz");
        do
            ln -s $read reads/
        done


    And use the basedir of the reads when running the pipeline, with --reads_dir reads

3. Edit the relevant variables in the run_on_dardel.sh slurm script.

4. Submit the pipeline to slurm using the run_on_dardel.sh script.

        sbatch run_on_dardel.sh


# Output

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

## 00_input_data
Contains the index reference genome

## 01_reports
Contains a number different QC reports for reads, mapped bamfiles and variants

## 02_bamfiles
Contains the final, mapped and processed bam/cram files for each sample, as well as a deduplication metrics file for each sample.


## 03_genotypes

### 00_raw_variants

Contains the raw variants in vcf format.

### 01_filtered_variants

Filtered SNPs and indels in vcf format, ready for downstream analyses. 

### 02_maskfiles

Contains three bedfiles per sample:
- <sample_id>_mappability_mask.bed
Callability mask across the genome, that is, this file contains genomic regions where we're confident in our ability to call genotypes
- <sample_id>_homref_invariants.bed
This file contains all intervals in the reference genomes with sufficient read coverage for variant calling, but where no variants were called. That is, one can assume that these sites are homozygous reference for the particular sample.
- <sample_id>_snp_mask.bed
Callability mask for SNPs: that is, this file contains genomic regions where we're confident in our ability to call SNPs if present. Any sites with indels will be excluded in this file.


