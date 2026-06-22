#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                          SUBWORKFLOW IMPORTS
//
// Pipeline is modularized into focused subworkflows for maintainability.
// ═══════════════════════════════════════════════════════════════════════════════

include { INDEX_REFERENCE } from './subworkflows/index_reference'
include { PREPROCESS_MODERN } from './subworkflows/preprocess_modern'
include { PREPROCESS_HISTORICAL } from './subworkflows/preprocess_historical'
include { MAP_MODERN } from './subworkflows/map_modern'
include { MAP_HISTORICAL } from './subworkflows/map_historical'
include { PROCESS_BAMS } from './subworkflows/process_bams'
include { VARIANT_CALLING } from './subworkflows/variant_calling'
include { VARIANT_FILTERS } from './subworkflows/variant_filters'

// ─────────────────────────────────────────────────────────────────────────────
// Summary Statistics Modules (used in main workflow)
// ─────────────────────────────────────────────────────────────────────────────
include { parse_summary_stats } from './modules/summary_stats/parse_summary_stats'
include { combine_summary_tables } from './modules/summary_stats/combine_summary_files'

// ═══════════════════════════════════════════════════════════════════════════════
//                      UTILITY FUNCTIONS & HELPERS
//
// These helper functions provide reusable logic for common workflow operations.
// ═══════════════════════════════════════════════════════════════════════════════

/*
 * Parse input metadata CSV
 * Expected format: sample_id;taxon;read_1;read_2
 */
def parse_input(metadata_file) {
    return channel
        .fromPath(metadata_file)
        .splitCsv(header: true, sep: ';')
        .filter { row -> !row.sample_id.startsWith('#') }
        .map { row ->
            def data_type = row.data_type
            def sample_id = row.sample_id
            def library = row.library

            // Construct read paths based on params.reads_dir
            def read_1_path = params.reads_dir ? "${params.reads_dir}/${row.read_1}" : row.read_1
            def read_2_path = params.reads_dir ? "${params.reads_dir}/${row.read_2}" : row.read_2

            def read_1 = file(read_1_path, checkIfExists: true)
            def read_2 = file(read_2_path, checkIfExists: true)
            return [sample_id, data_type, library, read_1, read_2]
        }
        // add lane information if multiple entries per sample_id
        .groupTuple(by: 0)
        .flatMap { sample_id, data_type, library, reads_1, reads_2 ->
            if (data_type.size() == 1) {
                // Single lane: return one tuple with 'single_lane' as lane identifier
                return [[sample_id, 'single_lane', data_type[0], library[0], reads_1[0], reads_2[0]]]
            } else {
                // Multiple lanes: return list of tuples with lane numbers
                return (0..<data_type.size()).collect { idx ->
                    [sample_id, "${idx + 1}", data_type[idx], library[idx], reads_1[idx], reads_2[idx]]
                }
            }
        }
}

/*
 * Build a channel with one tuple per (sample cram/bam, reference bundle, interval).
 * Inputs:
 *   - bam_ch: tuples or values containing a cram/bam file (sample metadata allowed)
 *   - ref_ch: tuple/value holding [reference_genome, index_files]
 *   - intervals_ch: values like 'chr1:1-100000'
 * Output tuple format:
 *   [sample_id, taxon, cram, [reference_genome, index_files], interval]
 */
def expand_by_intervals(bam_ch, ref_ch, intervals_ch) {
    bam_ch
        .combine(ref_ch)
        .combine(intervals_ch)
        .map { sample_id, cram, reference, index_files, interval ->
            return tuple(sample_id, cram, [reference, index_files], interval)
        }
}

/*
 * Calculate sample depth statistics and sex assignments
 * Inputs:
 *   - depth_avg_ch: channel from parse_region_depths.sample_depth_avg
 *   - sex_limited_list: list of sex-limited scaffolds (Y/W)
 *   - non_sex_limited_list: list of non-sex-limited sex chromosomes (X/Z)
 * Output:
 *   channel with tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
 */
def calculate_depth_and_sex(depth_avg_ch, sex_limited_list, non_sex_limited_list) {
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
                if (sex_limited_list.contains(chrom)) {
                    category = 'sex_limited'
                }
                if (non_sex_limited_list.contains(chrom)) {
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
                if (ratio < params.sex_assignment_lower_threshold) {
                    sex = 'unknown'
                } else if (ratio >= 1.5) {
                    sex = 'unknown'
                } else if (params.sex_assignment_lower_threshold <= ratio && ratio <= params.sex_assignment_upper_threshold) {
                    sex = 'hemizygous'
                } else {
                    sex = 'homozygous'
                }
            }

            tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
        }
}

/*
 * Get filter expressions for variant calling based on variant caller type
 * Returns map with snp_filter_expr and indel_filter_expr
 */
def get_filter_expressions(caller_type, custom_snp=null, custom_indel=null) {
    // If custom expressions are provided, use those
    if (custom_snp && custom_indel) {
        return [snp_filter_expr: custom_snp, indel_filter_expr: custom_indel]
    }

    // Otherwise, return defaults based on variant caller
    if (caller_type in ['gatk_joint', 'gatk_haplotypecaller']) {
        return [snp_filter_expr: params.snp_filter_expression_gatk, indel_filter_expr: params.indel_filter_expression_gatk]
    }
    else if (caller_type == 'freebayes') {
        return [snp_filter_expr: params.snp_filter_expression_freebayes, indel_filter_expr: params.indel_filter_expression_freebayes]
    }
    else if (caller_type == 'bcftools') {
        return [snp_filter_expr: params.snp_filter_expression_bcftools, indel_filter_expr: params.indel_filter_expression_bcftools]
    }
    else {
        error "Unknown variant caller for filter expression: ${caller_type}"
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

/*
 * Main workflow - orchestrates subworkflows
 */
workflow {
    main:

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SETUP: Configuration and Input Parsing
    // ═══════════════════════════════════════════════════════════════════════════════

    // Setup sex chromosome system based on user parameters
    sex_config = setup_sex_chromosome_system()
    def sex_chrom_system = sex_config.sex_chrom_system
    def sex_linked_list = sex_config.sex_linked_list
    def sex_limited_list = sex_config.sex_limited_list
    def non_sex_limited_list = sex_config.non_sex_limited_list
    println "Sex-limited contigs: ${sex_limited_list}"
    println "Non-sex-limited contigs: ${non_sex_limited_list}"

    // Validate variant caller parameter
    if (!['gatk_joint','gatk_haplotypecaller','freebayes','bcftools'].contains(params.variant_caller)) {
        error "Unsupported variant caller specified: ${params.variant_caller}. Must be one of gatk_joint, gatk_haplotypecaller, freebayes or bcftools."
    }

    println "Inferred sex chromosome system based on user specifications: ${sex_chrom_system}"

    // Create channels for sex-linked contigs
    sex_linked_contigs = channel.value(sex_linked_list)
    sex_limited_contigs = channel.value(sex_limited_list)
    non_sex_limited_contigs = channel.value(non_sex_limited_list)

    // Input channel from metadata file
    ch_input = parse_input(params.input)
        .branch {
            sample_id, lane, data_type, library, r1, r2 ->
            modern: data_type == '1'
            historical: data_type == '2'
        }

    // Reference genome
    ch_reference = channel.fromPath(params.reference, checkIfExists: true)

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SUBWORKFLOW 1: INDEX REFERENCE
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
    //                     SUBWORKFLOW 2-3: PREPROCESS READS
    // ═══════════════════════════════════════════════════════════════════════════════

    PREPROCESS_MODERN(
        ch_input.modern,
        params.premapping_dedup
    )

    PREPROCESS_HISTORICAL(
        ch_input.historical,
        params.premapping_dedup
    )

    // Combine multiqc reports
    multiqc_rawreads_report = PREPROCESS_MODERN.out.multiqc_raw_report
        .mix(PREPROCESS_HISTORICAL.out.multiqc_raw_report)
        .collect()

    multiqc_cleanreads_report = PREPROCESS_MODERN.out.multiqc_clean_report
        .mix(PREPROCESS_HISTORICAL.out.multiqc_clean_report)
        .collect()

    clean_reads = PREPROCESS_MODERN.out.clean_reads
        .mix(PREPROCESS_HISTORICAL.out.clean_paired)
        .mix(PREPROCESS_HISTORICAL.out.clean_merged)

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SUBWORKFLOW 4-5: MAP READS
    // ═══════════════════════════════════════════════════════════════════════════════

    MAP_MODERN(
        PREPROCESS_MODERN.out.clean_reads,
        INDEX_REFERENCE.out.bwa_index,
        params.mapper
    )

    MAP_HISTORICAL(
        PREPROCESS_HISTORICAL.out.clean_paired,
        PREPROCESS_HISTORICAL.out.clean_merged,
        INDEX_REFERENCE.out.bwa_index,
        params.map_historical_pairs,
        params.historical_mapper
    )

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SUBWORKFLOW 6: ESTIMATE DEPTH & POST-PROCESS BAMS
    // ═══════════════════════════════════════════════════════════════════════════════

    PROCESS_BAMS(
        MAP_MODERN.out.bam,
        MAP_HISTORICAL.out.bam,
        ch_reference,
        INDEX_REFERENCE.out.reference_fai,
        INDEX_REFERENCE.out.reference_gzi,
        INDEX_REFERENCE.out.refintervals,
        sex_limited_list,
        non_sex_limited_list,
        sex_limited_contigs,
        non_sex_limited_contigs,
        params.postmapping_dedup,
        params.damageprofiler_rescale,
        params.downsample_bams,
        params.downsample_bams_coverage,
        params.min_depth,
        params.max_depth,
        params.sex_assignment_lower_threshold,
        params.sex_assignment_upper_threshold
    )

    // Prepare sample stats channel with consistent tuple size
    sample_stats = PROCESS_BAMS.out.sample_depths
        .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment ->
            if (sex_chrom_system == 'unknown') {
                return tuple(sample_id, autosomal_dp, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')
            } else {
                return tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, 'NA', 'NA', 'NA')
            }
        }

    // Add downsampled stats if downsampling was performed
    if (params.downsample_bams) {
        sample_stats = PROCESS_BAMS.out.sample_depths
            .combine(PROCESS_BAMS.out.sample_depths_downsampled, by: 0)
            .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment,
                   ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp, ds_ratio, ds_sex_assignment ->
                if (sex_chrom_system == 'unknown') {
                    return tuple(sample_id, autosomal_dp, 'NA', 'NA', 'NA', 'NA', ds_autosomal_dp, 'NA', 'NA')
                } else {
                    return tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment,
                                 ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp)
                }
            }
    }

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     SUBWORKFLOW 7-8: VARIANT CALLING & FILTERING
    // ═══════════════════════════════════════════════════════════════════════════════

    if (!params.skip_variant_calling) {

        VARIANT_CALLING(
            PROCESS_BAMS.out.final_bam,
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

        VARIANT_FILTERS(
            VARIANT_CALLING.out.raw_vcfs,
            ch_reference,
            INDEX_REFERENCE.out.reference_fai,
            INDEX_REFERENCE.out.reference_gzi,
            INDEX_REFERENCE.out.refintervals,
            PROCESS_BAMS.out.callable_regions,
            PROCESS_BAMS.out.depth_cutoffs,
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

        // ═══════════════════════════════════════════════════════════════════════════
        //                     SUMMARY STATISTICS GENERATION
        // ═══════════════════════════════════════════════════════════════════════════

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

    } // End of variant calling section

    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PUBLISH OUTPUTS
    // ═══════════════════════════════════════════════════════════════════════════════

    publish:
    // FastQC / MultiQC reports on raw and trimmed reads
    fastqc_raw = PREPROCESS_MODERN.out.fastqc_raw
        .mix(PREPROCESS_HISTORICAL.out.fastqc_rawreads)
    fastqc_clean = PREPROCESS_MODERN.out.fastqc_clean
        .mix(PREPROCESS_HISTORICAL.out.fastqc_cleanreads)

    multiqc_rawreads_report   = multiqc_rawreads_report
    multiqc_cleanreads_report = multiqc_cleanreads_report

    // clean reads
    clean_reads = clean_reads

    // Per-sample CRAM QC reports
    qualimap_reports = PROCESS_BAMS.out.qualimap_reports
    qualimap_downsampled_reports = params.downsample_bams ? PROCESS_BAMS.out.qualimap_downsampled_reports : channel.empty()

    // Raw CRAM files and deduplicated metrics
    bamfiles = PROCESS_BAMS.out.raw_bams
    bam_metrics = PROCESS_BAMS.out.cram_metrics
    mapping_depths = PROCESS_BAMS.out.mapping_depths
    cramfiles = params.store_crams ? PROCESS_BAMS.out.raw_crams : channel.empty()
    downsampled_cramfiles = (params.downsample_bams && params.store_crams) ? PROCESS_BAMS.out.downsampled_crams : channel.empty()

    // Reference genome
    reference_genome = INDEX_REFERENCE.out.bwa_index

    // Downsampled BAM files
    downsampled_bamfiles = params.downsample_bams ? PROCESS_BAMS.out.final_bam : PROCESS_BAMS.out.final_bam

    // Variant calling outputs (conditional on skip_variant_calling)
    raw_vcf = !params.skip_variant_calling ? VARIANT_CALLING.out.raw_vcf : channel.empty()
    filtered_snps = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_snps : channel.empty()
    filtered_indels = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_indels : channel.empty()

    // Callable regions
    callable_regions = !params.skip_variant_calling ? VARIANT_FILTERS.out.callable_regions_bed : channel.empty()
    snpable_regions = !params.skip_variant_calling ? VARIANT_FILTERS.out.snpable_regions_bed : channel.empty()
    invariant_calls = !params.skip_variant_calling ? VARIANT_FILTERS.out.invariant_calls_bed : channel.empty()

    // Variant statistics
    raw_variant_stats = !params.skip_variant_calling ? VARIANT_CALLING.out.raw_variant_stats : channel.empty()
    raw_vcf_stats_plot = !params.skip_variant_calling ? VARIANT_CALLING.out.raw_vcf_stats_plot : channel.empty()
    filtered_snps_stats = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_snps_stats : channel.empty()
    filtered_snps_stats_plot = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_snps_stats_plot : channel.empty()
    filtered_indel_stats = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_indel_stats : channel.empty()
    filtered_indel_stats_plot = !params.skip_variant_calling ? VARIANT_FILTERS.out.filtered_indel_stats_plot : channel.empty()

    // Summary statistics
    summary_statistics = !params.skip_variant_calling ? combine_summary_tables.out.table : channel.empty()

    // Damage profiles
    damage_profiles = PROCESS_BAMS.out.damage_reports
}

output {
    reference_genome {
        path "00_input_data/00_reference_genome"
    }
    clean_reads {
        enabled params.store_cleanreads
        path "00_input_data/01_clean_reads"
    }
    fastqc_raw {
        enabled params.store_sample_fastqc
        path "01_reports/00_fastqc/00_raw_reads/sample_fastqc_reports"
    }
    fastqc_clean {
        enabled params.store_sample_fastqc
        path "01_reports/00_fastqc/01_clean_reads/sample_fastqc_reports"
    }
    multiqc_rawreads_report {
        path "01_reports/00_fastqc/00_raw_reads"
    }
    multiqc_cleanreads_report {
        path "01_reports/00_fastqc/01_clean_reads"
    }
    qualimap_reports {
        path "01_reports/01_qualimap/01_qualimap_post_markdups"
    }
    qualimap_downsampled_reports {
        enabled params.downsample_bams
        path "01_reports/01_qualimap/02_qualimap_downsampled"
    }
    bamfiles {
        enabled !params.store_crams
        path "02_bamfiles"
    }
    downsampled_bamfiles {
        enabled params.downsample_bams && !params.store_crams
        path "02_bamfiles/downsampled"
    }
    mapping_depths {
        path "02_bamfiles/mapping_depths"
    }
    bam_metrics {
        enabled params.postmapping_dedup
        path "02_bamfiles/dedup_metrics"
    }
    cramfiles {
        enabled params.store_crams
        path "02_bamfiles"
    }
    downsampled_cramfiles {
        enabled params.downsample_bams && params.store_crams
        path "02_bamfiles/downsampled"
    }
    raw_vcf {
        path "03_genotypes/00_raw_variants"
        enabled params.store_raw_vcf && !params.skip_variant_calling
    }
    raw_variant_stats {
        path "01_reports/02_variantstats/00_raw_variants"
        enabled !params.skip_variant_calling
    }
    raw_vcf_stats_plot {
        path "01_reports/02_variantstats/00_raw_variants"
        enabled !params.skip_variant_calling
    }
    filtered_snps_stats {
        path "01_reports/02_variantstats/01_filtered_snps"
        enabled !params.skip_variant_calling
    }
    filtered_snps_stats_plot {
       path "01_reports/02_variantstats/01_filtered_snps"
        enabled !params.skip_variant_calling
    }
    filtered_indel_stats {
        path "01_reports/02_variantstats/02_filtered_indels"
        enabled !params.skip_variant_calling
    }
    filtered_indel_stats_plot {
        path "01_reports/02_variantstats/02_filtered_indels"
        enabled !params.skip_variant_calling
    }
    callable_regions {
        path "03_genotypes/02_maskfiles"
        enabled !params.skip_variant_calling
    }
    snpable_regions {
        path "03_genotypes/02_maskfiles"
        enabled !params.skip_variant_calling
    }
    invariant_calls {
        path "03_genotypes/02_maskfiles"
        enabled !params.skip_variant_calling
    }
    filtered_snps {
        path "03_genotypes/01_filtered_variants"
        enabled !params.skip_variant_calling
    }
    filtered_indels {
        path "03_genotypes/01_filtered_variants"
        enabled !params.skip_variant_calling
    }
    summary_statistics {
        path "01_reports"
        enabled !params.skip_variant_calling
    }
    damage_profiles {
        path "01_reports/03_damage_profiles"
    }
}
