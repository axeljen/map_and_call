#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//         CALL VARIANTS WORKFLOW - Variant calling from existing BAMs
//
// Purpose: Start from existing BAM files and run variant calling and filtering
// Use case: Run variant calling on existing alignments, or re-run with different parameters
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Subworkflow imports
// ─────────────────────────────────────────────────────────────────────────────
include { INDEX_REFERENCE } from '../subworkflows/index_reference'
include { VARIANT_CALLING } from '../subworkflows/variant_calling'
include { VARIANT_FILTERS } from '../subworkflows/variant_filters'

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { extract_sample_from_bam } from '../modules/samtools/extract_sample_from_bam'
include { samtools_dp } from '../modules/samtools/samtools_dp_process'
include { parse_region_depths } from '../modules/samtools/parse_region_depths'
include { samtools_downsample } from '../modules/samtools/samtools_downsample'
include { callable_regions } from '../modules/bedtools/callable_regions'
include { parse_summary_stats } from '../modules/summary_stats/parse_summary_stats'
include { combine_summary_tables } from '../modules/summary_stats/combine_summary_files'

// ─────────────────────────────────────────────────────────────────────────────
// Utility imports
// ─────────────────────────────────────────────────────────────────────────────
// import WorkflowUtils

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/*
 * Parse BAM file list and extract sample IDs from read groups
 * Input: path to text file with one BAM path per line
 * Returns: channel with tuples [bam, bai]
 */
def parse_bamfile_list(bamfile_list) {
    return channel
        .fromPath(bamfile_list)
        .splitText()
        .map { line -> line.trim() }
        .filter { line -> line && !line.startsWith('#') }
        .map { bam_path ->
            def bam = file(bam_path, checkIfExists: true)
            // Try common BAI naming patterns
            def index_candidates = [
            file("${bam_path}.bai"),
            file("${bam_path}.crai"),
            file(bam_path.replaceAll(/\.bam$/, '.bai')),
            file(bam_path.replaceAll(/\.cram$/, '.crai'))
                ].unique()
            def bai = index_candidates.find { bai -> bai.exists() }
            if (!bai) {
                error "BAI index file not found for ${bam_path}. Tried: ${index_candidates.collect{it.toString()}.join(', ')}"
            }
            return tuple(bam, bai)
        }
}

/*
 * Setup sex chromosome system and identify sex-linked contigs
 * Returns map with: sex_chrom_system, sex_linked_list, sex_limited_list, non_sex_limited_list
 */
def setup_sex_chromosome_system() {
    def result = [
        sex_chrom_system: 'unknown',
        sex_linked_list: [],
        sex_limited_list: [],
        non_sex_limited_list: []
    ]
    
    if (params.x_scaffolds || params.y_scaffolds) {
        if (params.z_scaffolds || params.w_scaffolds) {
            error "Please only specify either X and/or Y, OR Z and/or W scaffolds, not both."
        }
        result.sex_chrom_system = 'XY'
        result.sex_linked_list = [params.x_scaffolds, params.y_scaffolds].flatten()
        result.sex_limited_list = [params.y_scaffolds].flatten()
        result.non_sex_limited_list = [params.x_scaffolds].flatten()
    }
    else if (params.z_scaffolds || params.w_scaffolds) {
        result.sex_chrom_system = 'ZW'
        result.sex_linked_list = [params.z_scaffolds, params.w_scaffolds].flatten()
        result.sex_limited_list = [params.w_scaffolds].flatten()
        result.non_sex_limited_list = [params.z_scaffolds].flatten()
    }
    
    return result
}


// ─────────────────────────────────────────────────────────────────────────────
// Helper function: Calculate sample depth statistics and sex assignments
// ─────────────────────────────────────────────────────────────────────────────
def calculate_depth_and_sex(depth_avg_ch, sex_limited_list, non_sex_limited_list, 
                            sex_assignment_lower_threshold, sex_assignment_upper_threshold) {
    return depth_avg_ch
        .flatMap {
            sample, file ->
                file.splitCsv(header: false, sep: '\t').collect { row ->
                    def chrom = row[0]
                    def avg_depth = row[2].toDouble()
                    return tuple(sample, chrom, avg_depth)
                }
        }
        // Categorize chromosomes
        .map {
            sample, chrom, avg_depth ->
                def category = 'autosomes'
                if (sex_limited_list && sex_limited_list.contains(chrom)) {
                    category = 'sex_limited'
                }
                if (non_sex_limited_list && non_sex_limited_list.contains(chrom)) {
                    category = 'non_sex_limited'
                }
                return tuple(sample, chrom, avg_depth, category)
        }
        // Group by sample and category, calculate average depth per category
        .groupTuple(by: [0, 3])
        .map { sample, chroms, values, category ->
            def avg_depth = values.sum() / values.size()
            return tuple(sample, [(category): avg_depth])
        }
        .groupTuple(by: 0)
        .map { sample, category_maps ->
            // Merge all category maps into a single map per sample
            def depth_map = [:]
            category_maps.each { map ->
                depth_map.putAll(map)
            }
            tuple(sample, depth_map)
        }
        // Calculate ratio and assign sex
        .map { sample, depth_map ->
            def autosome_depth = depth_map.get('autosomes', 0.0)
            def non_sex_limited_depth = depth_map.get('non_sex_limited', 0.0)
            def sex_limited_depth = depth_map.get('sex_limited', 0.0)
            
            def ratio = (autosome_depth > 0 && non_sex_limited_depth > 0) ? 
                        non_sex_limited_depth / autosome_depth : 0.0
            
            def sex = 'unknown'
            if (ratio > 0) {
                if (ratio < sex_assignment_lower_threshold) {
                    sex = 'unknown'
                } else if (ratio >= 1.5) {
                    sex = 'unknown'
                } else if (sex_assignment_lower_threshold <= ratio && ratio <= sex_assignment_upper_threshold) {
                    sex = 'hemizygous'
                } else {
                    sex = 'homozygous'
                }
            }
            
            tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
        }
}


// ═══════════════════════════════════════════════════════════════════════════════
//                           MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════════

workflow {
    main:
    
    println "Running CALL_VARIANTS workflow: variant calling from existing BAM files"
    
    // Validate required parameters
    if (!params.bamfiles) {
        error "ERROR: --call_variants workflow requires --bamfiles parameter (path to text file listing BAM paths)"
    }
    
    // Validate variant caller parameter
    if (!['gatk_joint','gatk_haplotypecaller','freebayes','bcftools'].contains(params.variant_caller)) {
        error "Unsupported variant caller specified: ${params.variant_caller}. Must be one of gatk_joint, gatk_haplotypecaller, freebayes or bcftools."
    }
    
    // Setup sex chromosome system
    sex_config = setup_sex_chromosome_system()
    def sex_chrom_system = sex_config.sex_chrom_system
    def sex_linked_list = sex_config.sex_linked_list
    def sex_limited_list = sex_config.sex_limited_list
    def non_sex_limited_list = sex_config.non_sex_limited_list
    
    println "Sex-limited contigs: ${sex_limited_list}"
    println "Non-sex-limited contigs: ${non_sex_limited_list}"
    println "Inferred sex chromosome system: ${sex_chrom_system}"
    
    // Create channels for sex-linked contigs
    sex_linked_contigs = channel.value(sex_linked_list)
    sex_limited_contigs = channel.value(sex_limited_list)
    non_sex_limited_contigs = channel.value(non_sex_limited_list)
    
    // Reference genome
    ch_reference = channel.fromPath(params.reference, checkIfExists: true)
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     INDEX REFERENCE
    // ═══════════════════════════════════════════════════════════════════════════════
    
    INDEX_REFERENCE(
        ch_reference,
        params.scaffold_list,
        params.chunk_size,
        params.x_scaffolds ?: [],
        params.y_scaffolds ?: [],
        params.z_scaffolds ?: [],
        params.w_scaffolds ?: []
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PARSE BAM FILES AND EXTRACT SAMPLE IDS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    // Parse BAM file list
    bam_input_ch = parse_bamfile_list(params.bamfiles)

    // Extract sample IDs from BAM read groups
    extract_sample_from_bam(bam_input_ch)
    
    // Combine extracted sample IDs with BAM/BAI files
    bams_with_samples = extract_sample_from_bam.out.sample_bams
        .map {bam, bai, sample_id_file -> 
            def sample_id = sample_id_file.text.trim()
            tuple(sample_id, bam, bai)
        }
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     DEPTH CALCULATION AND OPTIONAL DOWNSAMPLING
    // ═══════════════════════════════════════════════════════════════════════════════
    
    // Calculate depth per region
    region_depths = samtools_dp(bams_with_samples
        .combine(INDEX_REFERENCE.out.refintervals))
        .groupTuple(by: 0)
    
    // Parse region depths and calculate sample depths
    sample_depths = calculate_depth_and_sex(
        parse_region_depths(region_depths, INDEX_REFERENCE.out.reference_fai).sample_depth_avg,
        sex_limited_list,
        non_sex_limited_list,
        params.sex_assignment_lower_threshold,
        params.sex_assignment_upper_threshold
    )
    
    // If downsampling is enabled, run downsampling
    if (params.downsample_bams) {
        // Calculate downsampling fractions
        if (params.downsample_bams_coverage != -1) {
            downsample_ch = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex ->
                    def fraction = (autosome_depth > 0) ? params.downsample_bams_coverage / autosome_depth : 0.0
                    fraction = fraction > 1 ? 1.0 : fraction
                    tuple(sample, fraction)
                }
        } else {
            // Downsample to minimum coverage
            min_autosome_depth = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex -> autosome_depth }
                .min()
            
            downsample_ch = sample_depths
                .combine(min_autosome_depth)
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex, target_min_depth ->
                    def fraction = (autosome_depth > 0) ? target_min_depth / autosome_depth : 0.0
                    fraction = fraction > 1 ? 1.0 : fraction
                    tuple(sample, fraction)
                }
        }
        
        // Downsample BAMs
        bams_for_calling = samtools_downsample(bams_with_samples
            .combine(downsample_ch, by: 0))
    } else {
        bams_for_calling = bams_with_samples
    }
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PREPARE SAMPLE STATS AND DEPTH CUTOFFS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    // Prepare sample stats channel
    sample_stats = sample_depths
        .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment ->
            if (sex_chrom_system == 'unknown') {
                return tuple(sample_id, autosomal_dp, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')
            } else {
                return tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, 'NA', 'NA', 'NA')
            }
        }
    
    // Calculate depth cutoffs for callable regions
    depth_cutoffs = sample_depths
        .map { sample_id, autosomal_dp, _non_sex_limited_dp, _sex_limited_dp, _ratio, sex_assignment ->
            def min_dp
            def max_dp
            
            if (params.min_depth instanceof Integer){
                min_dp = params.min_depth
            }
            else {
                min_dp = (autosomal_dp * params.min_depth)
            }
            if (params.max_depth instanceof Integer){
                max_dp = params.max_depth
            }
            else {
                max_dp = (autosomal_dp * params.max_depth)
            }
            if (min_dp < 1) {
                log.warn "Sample ${sample_id}: dynamic min_depth evaluated to ${min_dp} (autosomal depth: ${autosomal_dp}). Setting min_depth to 1."
                min_dp = 1
            }
            if (max_dp < 2) {
                log.warn "Sample ${sample_id}: dynamic max_depth evaluated to ${max_dp} (autosomal depth: ${autosomal_dp}). Setting max_depth to 2."
                max_dp = 2
            }
            tuple(sample_id, min_dp, max_dp, sex_assignment)
        }
        // Add the sample bedfile to this
        .combine(parse_region_depths.out.sample_depth_beds, by: 0)
    
    // Calculate callable regions
    callable_regions_input = depth_cutoffs
        .combine(sex_limited_contigs.toList())
        .combine(non_sex_limited_contigs.toList())
        .combine(INDEX_REFERENCE.out.reference_fai)

    callable_regions_out = callable_regions(callable_regions_input)
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     VARIANT CALLING
    // ═══════════════════════════════════════════════════════════════════════════════
    
    VARIANT_CALLING(
        bams_for_calling,
        ch_reference,
        INDEX_REFERENCE.out.reference_fai,
        INDEX_REFERENCE.out.reference_gzi,
        INDEX_REFERENCE.out.bwa_index,
        INDEX_REFERENCE.out.refintervals,
        params.variant_caller,
        params.popfile,
        sample_stats,
        params.store_raw_vcf
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     VARIANT FILTERING
    // ═══════════════════════════════════════════════════════════════════════════════
    
    VARIANT_FILTERS(
        VARIANT_CALLING.out.raw_vcfs,
        ch_reference,
        INDEX_REFERENCE.out.reference_fai,
        INDEX_REFERENCE.out.reference_gzi,
        INDEX_REFERENCE.out.refintervals,
        callable_regions_out,
        depth_cutoffs,
        sex_limited_contigs,
        params.variant_caller,
        params.snp_filter_expression,
        params.indel_filter_expression,
        params.snp_filter_expression_bcftools,
        params.snp_filter_expression_freebayes,
        params.snp_filter_expression_gatk,
        params.indel_filter_expression_bcftools,
        params.indel_filter_expression_freebayes,
        params.indel_filter_expression_gatk
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SUMMARY STATISTICS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    parse_summary_stats(
        VARIANT_CALLING.out.pops
            .combine(VARIANT_FILTERS.out.invariant_calls_bed, by: 0)
            .combine(VARIANT_FILTERS.out.callable_regions_bed, by: 0)
            .combine(VARIANT_FILTERS.out.snpable_regions_bed, by: 0)
            .combine(sample_stats, by: 0)
            .combine(VARIANT_FILTERS.out.filtered_snps_stats.toList())
            .combine(VARIANT_FILTERS.out.filtered_indel_stats.toList()),
        sex_chrom_system
    )

    combine_summary_tables(parse_summary_stats.out.summary_statistics.collect())
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PUBLISH OUTPUTS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    publish:
    reference_genome = INDEX_REFERENCE.out.bwa_index
    raw_vcf = params.store_raw_vcf ? VARIANT_CALLING.out.raw_vcf : channel.empty()
    filtered_snps = VARIANT_FILTERS.out.filtered_snps
    filtered_indels = VARIANT_FILTERS.out.filtered_indels
    callable_regions = VARIANT_FILTERS.out.callable_regions_bed
    snpable_regions = VARIANT_FILTERS.out.snpable_regions_bed
    invariant_calls = VARIANT_FILTERS.out.invariant_calls_bed
    raw_variant_stats = VARIANT_CALLING.out.raw_variant_stats
    raw_vcf_stats_plot = VARIANT_CALLING.out.raw_vcf_stats_plot
    filtered_snps_stats = VARIANT_FILTERS.out.filtered_snps_stats
    filtered_snps_stats_plot = VARIANT_FILTERS.out.filtered_snps_stats_plot
    filtered_indel_stats = VARIANT_FILTERS.out.filtered_indel_stats
    filtered_indel_stats_plot = VARIANT_FILTERS.out.filtered_indel_stats_plot
    summary_statistics = combine_summary_tables.out.table
}

// ═══════════════════════════════════════════════════════════════════════════════
//                             OUTPUT DEFINITIONS
// ═══════════════════════════════════════════════════════════════════════════════

output {
    reference_genome {
        path "00_input_data/00_reference_genome"
    }
    raw_vcf {
        enabled params.store_raw_vcf
        path "03_genotypes/00_raw_variants"
    }
    raw_variant_stats {
        path "01_reports/02_variantstats/00_raw_variants"
    }
    raw_vcf_stats_plot {
        path "01_reports/02_variantstats/00_raw_variants"
    }
    filtered_snps_stats {
        path "01_reports/02_variantstats/01_filtered_snps"
    }
    filtered_snps_stats_plot {
       path "01_reports/02_variantstats/01_filtered_snps"
    }
    filtered_indel_stats {
        path "01_reports/02_variantstats/02_filtered_indels"
    }
    filtered_indel_stats_plot {
        path "01_reports/02_variantstats/02_filtered_indels"
    }
    callable_regions {
        path "03_genotypes/02_maskfiles"
    }
    snpable_regions {
        path "03_genotypes/02_maskfiles"
    }
    invariant_calls {
        path "03_genotypes/02_maskfiles"
    }
    filtered_snps {
        path "03_genotypes/01_filtered_variants"
    }
    filtered_indels {
        path "03_genotypes/01_filtered_variants"
    }
    summary_statistics {
        path "01_reports"
    }
}
