#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: ESTIMATE_DEPTH
//
// Purpose: BAM post-processing, depth profiling, and sex assignment
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { samtools_merge } from '../modules/samtools/mergebams'
include { samtools_markdups } from '../modules/samtools/samtools_markdups'
include { samtools_bam2cram } from '../modules/samtools/samtools_bam2cram'
include { damage_profiler } from '../modules/mapdamage/damageprofiler'
include { damage_profiler_rescale } from '../modules/mapdamage/damageprofiler_rescale'
include { qualimap } from '../modules/qualimap/qualimap'
include { qualimap as qualimap_downsampled } from '../modules/qualimap/qualimap'
include { samtools_dp } from '../modules/samtools/samtools_dp_process'
include { samtools_dp as samtools_dp_downsampled } from '../modules/samtools/samtools_dp_process'
include { parse_region_depths } from '../modules/samtools/parse_region_depths'
include { parse_region_depths as parse_region_depths_downsampled } from '../modules/samtools/parse_region_depths'
include { samtools_downsample } from '../modules/samtools/samtools_downsample'
include { callable_regions } from '../modules/bedtools/callable_regions'

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

workflow PROCESS_BAMS {
    take:
    rawbam_modern          // tuple: [sample_id, library, datatype, bam, bai]
    rawbam_historical      // tuple: [sample_id, library, datatype, bam, bai]
    ch_reference           // path: reference FASTA
    reference_fai          // path: reference .fai index
    reference_gzi          // path: reference .gzi index
    refintervals_ch        // tuple: [region_id, interval]
    sex_limited_list       // list: sex-limited chromosome names (Y/W)
    non_sex_limited_list   // list: non-sex-limited sex chromosome names (X/Z)
    sex_limited_contigs    // channel value: sex-limited chromosomes
    non_sex_limited_contigs // channel value: non-sex-limited chromosomes
    postmapping_dedup      // val: boolean flag for post-mapping deduplication
    damageprofiler_rescale // val: boolean flag for damage rescaling
    downsample_bams        // val: boolean flag for downsampling
    downsample_coverage    // val: target coverage for downsampling (-1 for min coverage)
    min_depth              // val: minimum depth threshold (integer or fraction)
    max_depth              // val: maximum depth threshold (integer or fraction)
    sex_assignment_lower   // val: lower threshold for sex assignment
    sex_assignment_upper   // val: upper threshold for sex assignment
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Merge BAM files per sample (if multiple libraries exist)
    // ─────────────────────────────────────────────────────────────────────────────
    per_sample_bams = rawbam_modern
        .mix(rawbam_historical)
        .groupTuple(by: 0)
        .map { sample_id, _libraries, datatypes, bam_paths, bam_indices ->
            tuple(sample_id, datatypes[0], bam_paths.flatten(), bam_indices.flatten())
        }
        .branch { _sample_id, _datatype, bam_paths, _bam_indices ->
            multi: bam_paths.size() > 1
            single: bam_paths.size() == 1
        }

    merged_bams = samtools_merge(per_sample_bams.multi
        .map { sample_id, datatype, bam_paths, _bam_indices -> tuple(sample_id, datatype, bam_paths) })

    all_sample_bams = merged_bams
        .mix(per_sample_bams.single.map { sample_id, datatype, bam_paths, bam_indices ->
            tuple(sample_id, datatype, bam_paths[0], bam_indices[0])
        })
        .branch { _sample_id, datatype, _bam, _bai ->
            modern:     datatype == '1'
            historical: datatype == '2'
        }

    // ─────────────────────────────────────────────────────────────────────────────
    // Post-mapping deduplication (optional)
    // ─────────────────────────────────────────────────────────────────────────────
    if (postmapping_dedup) {
        dedup_bams = samtools_markdups(all_sample_bams.historical.mix(all_sample_bams.modern))
        all_sample_bams = dedup_bams.bam.branch { _sample_id, datatype, _bam, _bai ->
            modern:     datatype == '1'
            historical: datatype == '2'
        }
        cram_metrics = dedup_bams.metrics
    } else {
        cram_metrics = channel.empty()
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Historical DNA: Damage profiling and rescaling
    // ─────────────────────────────────────────────────────────────────────────────
    if (damageprofiler_rescale) {
        println "Running damage profiling and rescaling on historical samples..."
        rescaled_bams = damage_profiler_rescale(
            all_sample_bams.historical,
            ch_reference
        )
        final_bam_ch = all_sample_bams.modern
            .map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) }
            .mix(rescaled_bams.rescaled_bam)
        damage_reports = rescaled_bams.damage_reports
    } else {
        println "Damage profiling and rescaling is disabled. Skipping this step and using original BAMs for downstream analyses."
        rescaled_bams = damage_profiler(
            all_sample_bams.historical.map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) },
            ch_reference
        )
        final_bam_ch = all_sample_bams.modern
            .map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) }
            .mix(all_sample_bams.historical.map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) })
        damage_reports = rescaled_bams.damage_reports
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Mapping QC: BAM quality assessment
    // ─────────────────────────────────────────────────────────────────────────────
    qualimap(final_bam_ch)

    // ─────────────────────────────────────────────────────────────────────────────
    // Depth profiling and sex assignment
    // ─────────────────────────────────────────────────────────────────────────────
    final_bam_ch
        .combine(ch_reference)
        .combine(reference_fai)
        .combine(reference_gzi)
        .set { bams2cram_input }
    if (params.store_crams) {
        crams = samtools_bam2cram(bams2cram_input)
    }
    bams_for_calling_ch = final_bam_ch

    dp_input_ch = final_bam_ch
        .combine(refintervals_ch)

    region_depths = samtools_dp(dp_input_ch)
        .groupTuple(by: 0)

    sample_depths = calculate_depth_and_sex(
        parse_region_depths(region_depths, reference_fai).sample_depth_avg,
        sex_limited_list,
        non_sex_limited_list,
        sex_assignment_lower,
        sex_assignment_upper
    )

    // ─────────────────────────────────────────────────────────────────────────────
    // Downsample BAMs (optional)
    // ─────────────────────────────────────────────────────────────────────────────
    if (downsample_bams) {
        if (downsample_coverage != -1){
            // If a coverage is specified, calculate the fraction for each sample
            downsample_ch = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex ->
                    def fraction = (autosome_depth > 0) ? downsample_coverage / autosome_depth : 0.0
                    tuple(sample, fraction)
                }
        }
        else {
            // Find the sample with the lowest average autosomal depth
            min_autosome_depth = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex -> autosome_depth }
                .min()
            downsample_ch = sample_depths
                .combine(min_autosome_depth)
                .map { sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex, target_min_depth ->
                    def fraction = (autosome_depth > 0) ? target_min_depth / autosome_depth : 0.0
                    tuple(sample, fraction)
                }
        }
        
        // If the fraction is greater than 1, set it to 1 (no downsampling)
        downsample_ch = downsample_ch
            .map { sample, fraction ->
                if (fraction > 1) {
                    println "Warning: calculated downsampling fraction for sample ${sample} is greater than 1 (${fraction}). Setting fraction to 1 (no downsampling)."
                    fraction = 1.0
                }
                tuple(sample, fraction)
            }
            .set { downsample_fractions }

        // Downsample BAMs
        bams_for_downsampling = final_bam_ch
            .combine(downsample_fractions, by: 0)
        bams_for_calling_ch = samtools_downsample(bams_for_downsampling)
        
        // Run qualimap on the downsampled bams
        qualimap_downsampled(bams_for_calling_ch)

        // convert downsampled bams to cram if requested
        if (params.store_crams) {
            downsampled_crams = samtools_bam2cram(bams_for_calling_ch
                .combine(ch_reference)
                .combine(reference_fai)
                .combine(reference_gzi))
        }

        // Fetch the per-base coverage for the downsampled bams
        region_depths_downsampled = samtools_dp_downsampled(bams_for_calling_ch
            .combine(refintervals_ch))
            .groupTuple(by: 0)
        
        // Calculate depth and sex assignments for downsampled data
        sample_depths_downsampled = calculate_depth_and_sex(
            parse_region_depths_downsampled(region_depths_downsampled, reference_fai).sample_depth_avg,
            sex_limited_list,
            non_sex_limited_list,
            sex_assignment_lower,
            sex_assignment_upper
        )
        
        qualimap_downsampled_reports = qualimap_downsampled.out.qualimap_report
        sample_depths_for_filters = sample_depths_downsampled
    } else {
        qualimap_downsampled_reports = channel.empty()
        sample_depths_for_filters = sample_depths
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Calculate depth cutoffs for variant filters
    // ─────────────────────────────────────────────────────────────────────────────
    depth_cutoffs = sample_depths_for_filters
        .map { sample_id, autosomal_dp, _non_sex_limited_dp, _sex_limited_dp, _ratio, sex_assignment ->
            def min_dp
            def max_dp
            
            if (min_depth instanceof Integer){
                min_dp = min_depth
            }
            else {
                min_dp = (autosomal_dp * min_depth)
            }
            if (max_depth instanceof Integer){
                max_dp = max_depth
            }
            else {
                max_dp = (autosomal_dp * max_depth)
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

    // Prepare input for callable_regions by adding scaffold lists and reference
    callable_regions_input = depth_cutoffs
        .combine(sex_limited_contigs.toList())
        .combine(non_sex_limited_contigs.toList())
        .combine(reference_fai)

    mappability_mask = callable_regions(callable_regions_input)

    emit:
    final_bam = bams_for_calling_ch
    raw_bams = final_bam_ch
    raw_crams = params.store_crams ? crams.cram : channel.empty()
    downsampled_crams = (params.downsample_bams && params.store_crams) ? downsampled_crams.cram : channel.empty()
    sample_depths = sample_depths
    sample_depths_downsampled = downsample_bams ? sample_depths_downsampled : channel.empty()
    mapping_depths = parse_region_depths.out.sample_depth_avg
    depth_cutoffs = depth_cutoffs
    callable_regions = mappability_mask.callable
    qualimap_reports = qualimap.out.qualimap_report
    qualimap_downsampled_reports = qualimap_downsampled_reports
    cram_metrics = cram_metrics
    damage_reports = damage_reports
}
