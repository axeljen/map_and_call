#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                          MODULE IMPORTS
//
// Modules are organized by functional category for clarity and maintainability.
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// QC & Preprocessing Modules
// ─────────────────────────────────────────────────────────────────────────────
include { fastqc as fastqc_rawreads } from './modules/fastqc/fastqc_process'
include { fastqc as fastqc_cleanreads } from './modules/fastqc/fastqc_process'
include { multiqc_fastqc as multiqc_rawreads } from './modules/multiqc/multiqc_fastqc'
include { multiqc_fastqc as multiqc_cleanreads} from './modules/multiqc/multiqc_fastqc'
include { fastp } from './modules/fastp/trimming'
include { adapterremoval } from './modules/adapterremoval/adapterremoval'

// ─────────────────────────────────────────────────────────────────────────────
// Read Mapping & Alignment Modules
// ─────────────────────────────────────────────────────────────────────────────
include { bwa_index } from './modules/bwa/bwa_index'
include { bwa_mem } from './modules/bwa/bwa_mem'
include { bwa_mem as map_historical } from './modules/bwa/bwa_mem'
include { bwa_mem_singlereads as map_merged } from './modules/bwa/bwa_mem'
include { bwa_mem_singlereads as map_singletons } from './modules/bwa/bwa_mem'

// ─────────────────────────────────────────────────────────────────────────────
// Mapping QC & Analysis Modules
// ─────────────────────────────────────────────────────────────────────────────
include { damage_profiler } from './modules/mapdamage/damageprofiler'
include { qualimap } from './modules/qualimap/qualimap'
include { qualimap as qualimap_pre_dedup } from './modules/qualimap/qualimap'
include { qualimap as qualimap_downsampled } from './modules/qualimap/qualimap'
include { cramqc } from './modules/samtools/cramqc'
include { multiqc_cram as cram_multiqc } from './modules/multiqc/multiqc_cram'

// ─────────────────────────────────────────────────────────────────────────────
// BAM Processing & Analysis Modules
// ─────────────────────────────────────────────────────────────────────────────
include { samtools_index } from './modules/samtools/index_reference'
include { samtools_merge } from './modules/samtools/mergebams'
include { samtools_merge as merge_historical_bams } from './modules/samtools/mergebams'
include { samtools_markdups } from './modules/samtools/samtools_markdups'
include { samtools_stats } from './modules/samtools/samtools_stats'
include { samtools_dp } from './modules/samtools/samtools_dp_process'
include { samtools_dp as samtools_dp_downsampled } from './modules/samtools/samtools_dp_process'
include { samtools_dp as xlinked_dp } from './modules/samtools/samtools_dp_process'
include { samtools_dp as ylinked_dp } from './modules/samtools/samtools_dp_process'
include { parse_region_depths } from './modules/samtools/parse_region_depths'
include { parse_region_depths as parse_region_depths_downsampled } from './modules/samtools/parse_region_depths'
include { samtools_downsample } from './modules/samtools/samtools_downsample'
include { gatk_mark_duplicates } from './modules/gatk/gatk_mark_duplicates'

// ─────────────────────────────────────────────────────────────────────────────
// Reference Preparation Modules
// ─────────────────────────────────────────────────────────────────────────────
include { dochunks } from './modules/reference_intervals/dochunks'

// ─────────────────────────────────────────────────────────────────────────────
// Variant Calling Modules
// ─────────────────────────────────────────────────────────────────────────────
include { mpileup } from './modules/bcftools/bcftools_mpileup'
include { freebayes } from './modules/freebayes/freebayes'
// NOTE: GATK HaplotypeCaller modules are available but currently unused in workflow:
// include { create_sequence_dict } from './modules/gatk/create_sequence_dict'
// include { haplotype_caller } from './modules/gatk/haplotype_caller'
// include { haplotype_caller_erc } from './modules/gatk/haplotype_caller_erc'
// include { combine_gvcfs } from './modules/gatk/combine_gvcfs'
// include { genotype_gvcfs_combined as genotype_gvcfs } from './modules/gatk/genotype_gvcfs_combined'

// ─────────────────────────────────────────────────────────────────────────────
// Variant Filtering & Transformation Modules
// ─────────────────────────────────────────────────────────────────────────────
include { bcftools_merge } from './modules/bcftools/bcftools_merge'
include { bcftools_merge as bcftools_merge_snps } from './modules/bcftools/bcftools_merge'
include { bcftools_merge as bcftools_merge_indels } from './modules/bcftools/bcftools_merge'
include { bcftools_concat as bcftools_concat_raw } from './modules/bcftools/bcftools_concat'
include { bcftools_concat as bcftools_concat_snps } from './modules/bcftools/bcftools_concat'
include { bcftools_concat as bcftools_concat_indels } from './modules/bcftools/bcftools_concat'
include { bcftools_concat as bcftools_concat_final } from './modules/bcftools/bcftools_concat'
include { bcftools_norm } from './modules/bcftools/bcftools_norm'
include { select_snps } from './modules/bcftools/select_snps'
include { select_indels } from './modules/bcftools/select_indels'
include { select_invariant } from './modules/bcftools/select_invariant'
include { bcftools_filter_gatkindels } from './modules/bcftools/bcftools_filter_gatkindels'
include { bcftools_filter as bcftools_filter_snps } from './modules/bcftools/bcftools_filter'
include { bcftools_filter as bcftools_filter_indels } from './modules/bcftools/bcftools_filter'
include { bcftools_filter_fmiss_maf as bcftools_filter_fmiss_maf_snps } from './modules/bcftools/bcftools_filter_fmiss_maf'
include { bcftools_filter_fmiss_maf as bcftools_filter_fmiss_maf_indels } from './modules/bcftools/bcftools_filter_fmiss_maf'
include { ab_dp_filter as ab_dp_filter_snps } from './modules/custom_variant_filters/ab_dp_filter'
include { ab_dp_filter as ab_dp_filter_indels } from './modules/custom_variant_filters/ab_dp_filter'

// ─────────────────────────────────────────────────────────────────────────────
// Region Masking & Callable Site Modules
// ─────────────────────────────────────────────────────────────────────────────
include { callable_regions } from './modules/bedtools/callable_regions'
include { combine_bedfiles } from './modules/bedtools/combine_bedfiles'
include { combine_bedfiles as combine_homref } from './modules/bedtools/combine_bedfiles'
include { mappability_filter } from './modules/bedtools/mappability_filter'
include { mappability_filter as indel_mappability_filter } from './modules/bedtools/mappability_filter'
include { finalize_masks } from './modules/bedtools/finalize_masks'
include { combine_bedfiles as combine_homref_invariants } from './modules/bedtools/combine_bedfiles'
include { combine_bedfiles as combine_mappability_masks } from './modules/bedtools/combine_bedfiles'
include { combine_bedfiles as combine_mappability_masks_snps } from './modules/bedtools/combine_bedfiles'

// ─────────────────────────────────────────────────────────────────────────────
// Variant & Summary Statistics Modules
// ─────────────────────────────────────────────────────────────────────────────
include { vcf_stats as raw_vcf_stats } from './modules/variant_stats/vcf_stats_process'
include { combine_stats as combine_raw_vcf_stats } from './modules/variant_stats/combine_stats'
include { vcf_stats as filtered_snp_stats } from './modules/variant_stats/vcf_stats_process'
include { vcf_stats as filtered_indel_stats } from './modules/variant_stats/vcf_stats_process'
include { combine_stats as combine_filtered_snps_stats } from './modules/variant_stats/combine_stats'
include { combine_stats as combine_filtered_indels_stats } from './modules/variant_stats/combine_stats'
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
        .filter { row -> !row.data_type.startsWith('#') }
        .map { row ->
            def data_type = row.data_type
            def sample_id = row.sample_id
            def read_1 = file(row.read_1, checkIfExists: true)
            def read_2 = file(row.read_2, checkIfExists: true)
            return [sample_id, data_type, read_1, read_2]
        }
        // add lane information if multiple entries per sample_id
        .groupTuple(by: 0)
        .flatMap { sample_id, data_type, reads_1, reads_2 ->
            if (data_type.size() == 1) {
                // Single lane: return one tuple with 'single_lane' as lane identifier
                return [[sample_id, 'single_lane', data_type[0], reads_1[0], reads_2[0]]]
            } else {
                // Multiple lanes: return list of tuples with lane numbers
                return (0..<data_type.size()).collect { idx ->
                    [sample_id, "${idx + 1}", data_type[idx], reads_1[idx], reads_2[idx]]
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
 * Main workflow
 */
workflow {
    main:
    def sex_chrom_system = 'unknown'
    def sex_linked_contigs = []
    def sex_limited_contigs = []
    def non_sex_limited_contigs = []
    def sex_linked_list = []
    def sex_limited_list = []
    def non_sex_limited_list = []
    // Some quick input checks
    if (!['gatk_joint','gatk_haplotypecaller','freebayes','bcftools'].contains(params.variant_caller)) {
        error "Unsupported variant caller specified: ${params.variant_caller}. Must be one of gatk_joint, gatk_haplotypecaller, freebayes or bcftools."
    }
    if (params.x_scaffolds || params.y_scaffolds ){
        if (params.z_scaffolds || params.w_scaffolds){
            error "Please only specify either X and/or Y, OR Z and/or W scaffolds, not both. "
        }
        else {
            sex_chrom_system = 'XY'
            sex_linked_list = [params.x_scaffolds, params.y_scaffolds].flatten()
            sex_limited_list = [params.y_scaffolds].flatten()
            non_sex_limited_list = [params.x_scaffolds].flatten()
            sex_linked_contigs = channel.value(sex_linked_list)
            sex_limited_contigs = channel.value(sex_limited_list)
            non_sex_limited_contigs = channel.value(non_sex_limited_list)
        }
        }
    else if (params.z_scaffolds || params.w_scaffolds){
        sex_chrom_system = 'ZW'
        sex_linked_list = [params.z_scaffolds, params.w_scaffolds].flatten()
        sex_limited_list = [params.w_scaffolds].flatten()
        non_sex_limited_list = [params.z_scaffolds].flatten()
        sex_linked_contigs = channel.value(sex_linked_list)
        sex_limited_contigs = channel.value(sex_limited_list)
        non_sex_limited_contigs = channel.value(non_sex_limited_list)
    }
    else {
        sex_chrom_system = 'unknown'
        sex_linked_contigs = channel.value([])
        sex_limited_contigs = channel.value([])
        non_sex_limited_contigs = channel.value([])
    }
    println "Inferred sex chromosome system based on user specifications: ${sex_chrom_system}"

    // Input channel from metadata file
    ch_input = parse_input(params.input)
        .branch {
            sample_id, lane, data_type, r1, r2 ->
            modern: data_type == '1'
            historical: data_type == '2'
        }

    // Reference genome
    ch_reference = channel.fromPath(params.reference, checkIfExists: true)

    // =========================================================================
    // RAW READS FASTQC
    // =========================================================================

    // fastqc on raw reads
    ch_raw_reads_for_qc = ch_input.modern.map { sample_id, lane, _datatype, r1, r2 ->
            [sample_id, lane, [r1, r2], 'raw_reads']
        }
        // push historical reads through the same fastqc process
        .mix(ch_input.historical.map { sample_id, lane, _datatype, r1, r2 ->
            [sample_id, lane, [r1, r2], 'raw_reads']
        })
        

    fastqc_rawreads(ch_raw_reads_for_qc)
    // Collect all FastQC outputs and run MultiQC
    ch_multiqc_input = fastqc_rawreads.out.zip
        .map { _sample_id, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'raw_reads') }

    multiqc_rawreads(ch_multiqc_input)

    // =========================================================================
    // PREPROCESSING MODERN SAMPLES
    // =========================================================================

    fastp(ch_input.modern)

    // =========================================================================
    // PREPROCESSING HISTORICAL SAMPLES
    // =========================================================================

    // Push through adapterremoval
    adapterremoval(ch_input.historical)


    // =========================================================================
    // POST-PROCESSING FASTQC
    // =========================================================================

    // qc the clean reads too
    fastqc_cleanreads(
        fastp.out.reads
            .map { sample_id, lane, _datatype, r1, r2 ->
                [sample_id, lane, [r1, r2], 'clean_reads']
            }
        .mix(
            adapterremoval.out.trimmed_reads
                .map { sample_id, lane, r1, r2, collapsed, singletons ->
                    [sample_id, lane, [r1, r2, collapsed, singletons], 'clean_reads']
                }
            )
        )

    multiqc_cleanreads(fastqc_cleanreads.out.zip
        .map { _sample_id, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'clean_reads') }
        )

    // =========================================================================
    // INDEX REFERENCE GENOME AND PROCESS SCAFFOLD LIST AND CREATE INTERVALS
    // ========================================================================

    // First index the reference genome
    bwa_index_ch = bwa_index(ch_reference)
    faidx_and_chunks_ch = samtools_index(ch_reference)
    reference_fai = faidx_and_chunks_ch.reference_fai.first()
    reference_gzi = faidx_and_chunks_ch.reference_gzi.first()
   
    
    // if user provides a list of scaffolds upon which to call variants, fetch this
    if (params.scaffold_list) {
        scaffolds_ch = channel.fromPath(params.scaffold_list, checkIfExists: true)
            .splitCsv(header: false)
            .map { row -> row[0] }
    }
    else {
        // if no scaffold list is provided, we just use the reference fai index to get the scaffold names
        scaffolds_ch = reference_fai
            .map { fai ->
                def scaffolds = []
                fai.eachLine { line ->
                    def scaffold = line.split('\t')[0]
                    scaffolds << scaffold
                }
                return scaffolds
            }
            .collect()
    }

    // now generate reference intervals - these will be used for parallelization in the variant calling 
    // and filtering step, so not used in the mapping/preprocessing
    reference_intervals = dochunks(reference_fai, params.chunk_size, scaffolds_ch)

    // define a list of autosomes and potentially xz/yw chromosomes if given
    scaffolds = scaffolds_ch
            .flatten()
            .branch { scaffold ->
                autosomes: [params.y_scaffolds, params.w_scaffolds, params.x_scaffolds, params.z_scaffolds].flatten().contains(scaffold) == false
                sex_limited: [params.y_scaffolds, params.w_scaffolds].flatten().contains(scaffold)
                non_sex_limited: [params.x_scaffolds, params.z_scaffolds].flatten().contains(scaffold)
            }


    // =========================================================================
    // MAP MODERN SAMPLES
    // =========================================================================

    
    // combine reads with indexed reference
    mapping_ch = fastp.out.reads
        .map { sample_id, lane, _datatype, r1, r2 -> tuple(sample_id, lane, r1, r2) }
        .combine(bwa_index_ch.reference)
    // Choose mapping approach based on user parameter
    if (params.mapper == 'bwa_mem') {
        rawbam_ch = bwa_mem(mapping_ch)
    }
    else {
        error "Unsupported mapper specified: ${params.mapper}. Currently only 'bwa_mem' is supported."
    }

    // =========================================================================
    // Handle reference genome
    // =========================================================================

    refintervals_ch = reference_intervals
        .map { row -> row?.trim() }
        .flatMap { row -> row ? row.split(/\r?\n/) as List : [] }
        .filter { row -> row }
        .collect()
        .flatMap { intervals ->
            intervals.withIndex().collect { interval, idx -> tuple(idx + 1, interval) }
        }

    // =========================================================================
    // MAP HISTORICAL SAMPLES
    // =========================================================================

    // same for historical samples
    map_historical_pairs_ch = adapterremoval.out.trimmed_reads
        .map { sample_id, lane, r1, r2, _collapsed, _singletons -> tuple(sample_id, lane, r1, r2) }
        .combine(bwa_index_ch.reference)
    historical_pairbams_ch = map_historical(map_historical_pairs_ch )
    
    map_historical_collapsed_ch = adapterremoval.out.trimmed_reads
        .map { sample_id, lane, _r1, _r2, collapsed, _singletons -> tuple(sample_id, lane, collapsed) }

    historical_collapsed_bams_ch = map_merged(map_historical_collapsed_ch, 'collapsed', bwa_index_ch.reference)

    map_historical_singletions_ch = adapterremoval.out.trimmed_reads
        .map { sample_id, lane, _r1, _r2, _collapsed, singletons -> tuple(sample_id, lane, singletons) }
        .combine(bwa_index_ch.reference)
    
    historical_singleton_bams_ch = map_singletons(map_historical_singletions_ch, 'singletons', bwa_index_ch.reference)

    // collect historical bamfiles and merge
    historical_pairbams_ch
        .mix(historical_collapsed_bams_ch)
        .mix(historical_singleton_bams_ch)
        .groupTuple(by: [0,1])
        .map { sample_id, _lane, bam_paths, _bam_indices ->
            return tuple(sample_id, bam_paths)
        }
        .set {historical_bam_ch}

    historical_bams_merged_ch = merge_historical_bams(historical_bam_ch)

    // ========================================================================
    //  if any sample was split across lanes, merge bams here
    // ========================================================================

    // clean up bamchannel a bit and combine with historical bam channel
    datatypes_ch = ch_input.modern
        .mix(ch_input.historical)
        .map { sample_id, _lane, data_type, _r1, _r2 -> return tuple(sample_id, data_type) }
        .groupTuple(by: [0,1])

    rawbam_ch
        .groupTuple(by: [0])
        .map { sample_id, _lanes, bam_paths, _bam_indices ->
            return tuple(sample_id, bam_paths)
        }
        .mix(historical_bams_merged_ch
          .groupTuple(by: [0])
        )
        // combine with input channel to add back the information on datatype
        .combine(datatypes_ch, by: 0)
        .set {rawbam_ch} 

    // pull out those with multiple bams
    rawbam_ch
        .filter({ _sample_id, bam_paths, datatype -> bam_paths.size() > 1 })
        .set( {bams2merge_ch} )
    
    // merge those bams
    mergebams_out = samtools_merge(bams2merge_ch.map { sample_id, bam_paths, _datatype -> tuple(sample_id, bam_paths) })
    
    // and concatenate with those that did not need merging
    rawbam_ch
        .filter({ _sample_id, bam_paths, _datatype -> bam_paths.size() == 1 })
        .map { sample_id, bam_paths, datatype ->
            return tuple(sample_id, bam_paths[0], datatype)
        }
        .set( {singles_ch} )

    mergebams_out
        // add datatype back in 
        .combine(datatypes_ch, by: 0)
        .mix(singles_ch)
        .set( {final_bam_ch} )
    
    // =========================================================================
    // Deduplicate bam files
    // =========================================================================

    // add reference path to final_bam_ch
    final_bam_ch = final_bam_ch
        .combine(ch_reference)
        .map { sample_id, bam, datatype, reference ->
            return tuple(sample_id, bam, datatype, reference)
        }
    qualimap_pre_dedup(final_bam_ch.map { sample_id, bam, _datatype, reference -> tuple(sample_id, bam, reference) })
    
    dedup_bam_ch = samtools_markdups(final_bam_ch.map { sample_id, bam, _datatype, reference -> tuple(sample_id, bam, reference) })

    // =========================================================================
    // If historical samples, run these through damage profiler with rescaling
    // =========================================================================

    // add datatype back to dedupped maps
    historical_bams = dedup_bam_ch.bam
        .combine(datatypes_ch, by: 0)
        .filter {
            _sample_id, _cram, _crai, datatype -> datatype == "2"
        }

    rescaled_bams = damage_profiler(historical_bams.map { sample_id, cram, crai, _datatype -> tuple(sample_id, cram, crai) }, ch_reference)

    // =========================================================================
    // Combine the rescaled bam with the deduplicated, modern ones, which will be our final crams
    // =========================================================================

    dedup_bam_ch.bam
        .combine(datatypes_ch, by: 0)
        .filter{ _sample_id, _cram, _crai, datatype -> datatype == "1" }
        .map { sample_id, bam, crai, _datatype -> tuple(sample_id, bam, crai) }
        .mix(rescaled_bams.rescaled_bam)
        .set { final_bam_ch }

    // =========================================================================
    // Bam/cram qc
    // =========================================================================
    
    qualimap(final_bam_ch)

    // =========================================================================
    // Depth calculation for downstream filtering and possibly sex assignment
    // =========================================================================
    // Calculate average mapping depth per sample 
    // =========================================================================

    bams_for_calling_ch = final_bam_ch

    dp_input_ch = final_bam_ch
        .combine(refintervals_ch)

    region_depths = samtools_dp(dp_input_ch)
        .groupTuple(by: 0)

    sample_depths = calculate_depth_and_sex(
        parse_region_depths(region_depths, reference_fai).sample_depth_avg,
        sex_limited_list,
        non_sex_limited_list
    )

    sample_stats = sample_depths
        .map { sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex ->
        if (sex_chrom_system == 'unknown') {
            return tuple(sample, autosome_depth, 'NA', 'NA', 'NA', 'NA')
        }
        else {
            return tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex)
        }
        }
    

    // If downsample is set, estimate the fraction to downsample for each sample
    if (params.downsample_bams) {
        if (params.downsample_bams_coverage != -1){
            // if a coverage is specified, calculate the fraction for each sample based on the average autosomal depth
            downsample_ch = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex ->
                    def fraction = (autosome_depth > 0) ? params.downsample_bams_coverage / autosome_depth : 0.0
                    tuple(sample, fraction)
                }
        }
        else {
            // if no coverage is specified, find the sample with the lowest average autosomal depth and calculate fractions based on that
            min_autosome_depth = sample_depths
                .map { sample, autosome_depth, _non_sex_limited_depth, _sex_limited_depth, _ratio, _sex -> autosome_depth }
                .min()
            downsample_ch = sample_depths
                .combine(min_autosome_depth)
                .map { sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex, min_depth ->
                    def fraction = (autosome_depth > 0) ? min_depth / autosome_depth : 0.0
                    tuple(sample, fraction)
                }
        }
        // If the fraction is greater than 1, set it to 1 (no downsampling), and print a warning, because this obviously shouldn't be the case
        downsample_ch = downsample_ch
            .map { sample, fraction ->
                if (fraction > 1) {
                    println "Warning: calculated downsampling fraction for sample ${sample} is greater than 1 (${fraction}). Setting fraction to 1 (no downsampling). Please check the average depth for this sample and the specified downsampling coverage."
                    fraction = 1.0
                }
                tuple(sample, fraction)
            }
            .set { downsample_fractions }

        // send it off to downsampling
        bams_for_downsampling = final_bam_ch
            .combine(downsample_fractions, by: 0)
        bams_for_calling_ch = samtools_downsample(bams_for_downsampling)
        // run qualimap on the downsampled bams too
        qualimap_downsampled(bams_for_calling_ch)

        // fetch the per-base coverage also for the downsampled bams
        region_depths_downsampled = samtools_dp_downsampled(bams_for_calling_ch
            .combine(refintervals_ch))
            .groupTuple(by: 0)
        
        // Calculate depth and sex assignments for downsampled data
        sample_depths_downsampled = calculate_depth_and_sex(
            parse_region_depths_downsampled(region_depths_downsampled, reference_fai).sample_depth_avg,
            sex_limited_list,
            non_sex_limited_list
        )
        // add stats to stats channel
        sample_stats = sample_stats
            .combine(sample_depths_downsampled, by: 0)
            .map { sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex_assignment, autosome_depth_downsampled, non_sex_limited_depth_downsampled, sex_limited_depth_downsampled, ratio_downsample, sex_assignment_downsampled ->
                if (sex_chrom_system == 'unknown') {
                    return tuple(sample, autosome_depth, 'NA', 'NA', 'NA', 'NA', autosome_depth_downsampled, 'NA', 'NA', 'NA', 'NA')
                }
                else {
                    return tuple(sample, autosome_depth, non_sex_limited_depth, sex_limited_depth, ratio, sex_assignment, autosome_depth_downsampled, non_sex_limited_depth_downsampled, sex_limited_depth_downsampled, ratio_downsample, sex_assignment_downsampled)
                }
            }
        // And hereafter, use the downasmpled coverage for filtering
        sample_depths = sample_depths_downsampled
    }

    bams_output = bams_for_calling_ch

    // depending on the user parameter settings, set the coverage thresholds for custom DP filters
    depth_cutoffs = sample_depths
        .map { sample_id, autosomal_dp, _non_sex_limited_dp, _sex_limited_dp, _ratio, sex_assignment ->
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
            tuple(sample_id, min_dp, max_dp, sex_assignment)
        }
        // add the sample bedfile to this one
        .combine(parse_region_depths.out.sample_depth_beds, by: 0)

    mappability_mask = callable_regions(depth_cutoffs, sex_limited_contigs, non_sex_limited_contigs, reference_fai)



    // =========================================================================
    // Parse populations for joint calling if specified
    // =========================================================================
    // if a popfile was provided, get this info
    if (params.popfile) {
        pops = channel.fromPath(params.popfile, checkIfExists: true)
            .splitCsv(header: false, sep: ";")
            .map { sample_id, population ->
                tuple(sample_id, population)
            }
    }
    else {
        pops = sample_stats
            .map { sample_id, _autosomal_dp, _non_sex_limited_dp, _sex_limited_dp, _ratio, _sex_assignment, _ds_autosomal_dp, _ds_non_sex_limited_dp, _ds_sex_limited_dp ->
                tuple(sample_id, sample_id)
            }
    }

    // if bcftools is our caller, we want to do separate callings per population
    if (params.popfile && params.variant_caller == 'bcftools') {
        bams_for_calling_ch = channel.fromPath(params.popfile, checkIfExists: true)
            .splitCsv(header: false, sep: ';')
            .map { row -> tuple(row[0], row[1]) }
            .combine(bams_for_calling_ch, by: 0)
            .map { _sample_id, population, cram, crai ->
                tuple(population, cram, crai)
            }
            .groupTuple(by: 0)
    }
    else if (!params.popfile && params.variant_caller == 'bcftools') {
        // otherwise just assign a unique population to each sample, and call them independently
        bams_for_calling_ch = bams_for_calling_ch
            .map { sample_id, cram, crai ->
                tuple(sample_id, cram, crai)
            }
    } else if (params.variant_caller == 'freebayes') {
        // if freebayes is our caller, we want to do a joint calling across all samples, so we don't need to modify the bam channel, but we do need to add the population info to the bam channel so that we can feed it into freebayes later on
        bams_for_calling_ch = bams_for_calling_ch
            .map { _sample_id, cram, crai ->
                tuple('joint', cram, crai)
            }
            .groupTuple(by: 0)
    }
    // else if freebays is our caller, combine all bams into a single channel and send the pops along to freebayes later on
    
    
    // // ========================================================================= \\
    // // ######################################################################### \\
    // //                              Variant calling
    // // ######################################################################### \\
    // // ========================================================================= \\

    // // Variant calling can be made with four different routes, see below.        \\
    // // The output of all of these routes will be the same: One vcf file per      \\
    // // that has been combined across all samples. This is what we will later     \\
    // // feed into the filter steps.                                               \\

    if ( !['gatk_joint','gatk_haplotypecaller','freebayes','bcftools'].contains(params.variant_caller) ) {
        error "Unsupported variant caller specified: ${params.variant_caller}. Must be one of gatk_joint, gatk_haplotypecaller, freebayes or bcftools."
    }

    // ref bundle with all files potentially needed for variant calling
    ref_index_ch = bwa_index_ch.reference
        .map { _reference, index_files -> index_files }
    ref_bundle_ch = ch_reference
        .combine(samtools_index.out.reference_fai)
        .combine(samtools_index.out.reference_gzi)
        .combine(ref_index_ch)
        .map { row ->
            def reference = row[0]
            def fai = row[1]
            def gzi = row[2]
            def bwa_indices = row[3..-1]
            tuple(reference, [fai, gzi] + bwa_indices)
        }

    varcall_ch = 
        bams_for_calling_ch.combine(ref_bundle_ch)
        .combine(refintervals_ch)

    // // ========================================================================= \\
    // // Bcf-tools mpileup
    // // ========================================================================= \\

    if (params.variant_caller == 'bcftools') {
        mpileup(varcall_ch)
        pop_vcfs_ch = mpileup.out.vcf
    }

    else if (params.variant_caller == 'freebayes') {
        // collect pop assignments
        freebayes_pops = pops
            .map {sample, pop -> "${sample}=${pop}"}
            .collect().toList()
        varcall_ch = varcall_ch
            .combine(freebayes_pops)
        freebayes(varcall_ch)
        pop_vcfs_ch = freebayes.out.vcf
    }

    // Combine vcfs across populations/samples if there are more than one population
    // def npops = 
    npops = bams_for_calling_ch.groupTuple(by: 0)
        .map { population, _crams, _crais -> population }
        .collect()
        .map {
            populations -> 
            populations.size()
        }
        .branch { n ->
             multiple_pops: n > 1
             single_pop: n <= 1
        }
    // npops.view()
    pop_vcfs_ch.groupTuple(by: 0)
        // make sure that populations are always sorted in a consistent order
        .map { region_id, populations, vcfs, idxs ->
            def zipped = [populations, vcfs, idxs].transpose()
                .sort { a, b -> a[0] <=> b[0] } // sort by population name
            def (populations_sorted, vcfs_sorted, idxs_sorted) = zipped.transpose()
            tuple(region_id, populations_sorted, vcfs_sorted, idxs_sorted)
        }
        .map {region_id, pops, vcfs, idxs ->
            tuple(region_id, pops, vcfs, idxs)
        }
    npops.multiple_pops.combine(pop_vcfs_ch.groupTuple(by: 0)
        // make sure that populations are always sorted in a consistent order
        .map { region_id, populations, vcfs, idxs ->
            def zipped = [populations, vcfs, idxs].transpose()
                .sort { a, b -> a[0] <=> b[0] } // sort by population name
            def (populations_sorted, vcfs_sorted, idxs_sorted) = zipped.transpose()
            tuple(region_id, populations_sorted, vcfs_sorted, idxs_sorted)
        }
        .map {region_id, pops, vcfs, idxs ->
            tuple(region_id, pops, vcfs, idxs)
        }
        )
        .map { _npops, region_id, pops, vcfs, idxs ->
            tuple(region_id, pops, vcfs, idxs)
        }
    /////// FIX HERE: MULTIPOP SHOULD ONLY EMIT IF MULTUPLE POPULATIONS, 
    multipop_merge_in = npops.multiple_pops.combine(pop_vcfs_ch.groupTuple(by: 0)
            // make sure that populations are always sorted in a consistent order
            .map { region_id, populations, vcfs, idxs ->
                def zipped = [populations, vcfs, idxs].transpose()
                    .sort { a, b -> a[0] <=> b[0] } // sort by population name
                def (populations_sorted, vcfs_sorted, idxs_sorted) = zipped.transpose()
                tuple(region_id, populations_sorted, vcfs_sorted, idxs_sorted)
            }
            .map { region_id, pops, vcfs, idxs ->
                tuple(region_id, pops, vcfs, idxs)
            }
        )
        .map { _npops, region_id, pops, vcfs, idxs ->
            tuple(region_id, pops, vcfs, idxs)
        }
        // Broadcast singleton reference files across all regions
        .combine(samtools_index.out.reference_fasta)
        .combine(samtools_index.out.reference_fai)
        .combine(samtools_index.out.reference_gzi)
        .map { region_id, pops, vcfs, idxs, reference_fasta, reference_fai, reference_gzi ->
            tuple(region_id, pops, vcfs, idxs, reference_fasta, reference_fai, reference_gzi)
        }

    multipop_vcf = bcftools_merge(multipop_merge_in, 'raw')

    singlepop_vcf = pop_vcfs_ch.combine(
        npops.single_pop
        )
        .map { region_id, _pop, vcf, idx, _npops ->
        tuple(region_id, vcf, idx)
        }

    raw_vcfs = multipop_vcf
        .mix(singlepop_vcf)
    
    // extract some summary stats from raw vcf files
    raw_stats = raw_vcf_stats(raw_vcfs, 'raw_variants')
    // and combine these into one file
    combined_raw_stats = combine_raw_vcf_stats(
        raw_stats.ab_dp.collect(),
        raw_stats.qual_fmiss_maf.collect(),
        raw_stats.sample_stats.collect(),
        raw_stats.rec_counts.collect(),
        'raw_variants'
    )

    // If raw vcfs are to be stored, concatenate here
    if (params.store_raw_vcf == true) {
        sorted_raw_vcfs = raw_vcfs
            .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by region_id
            .map { sorted_list ->
                // Extract VCFs and indices in sorted order
                def vcfs = sorted_list.collect { it[1] }
                def idxs = sorted_list.collect { it[2] }
                tuple(vcfs, idxs)
            }
        // Concatenate raw VCF files
        bcftools_concat_raw(sorted_raw_vcfs, 'raw_variants')
    }

    // // ========================================================================= \\
    // // GATK HaplotypeCaller, per sample only 
    // // ======================================================================== \\

    // if (params.variant_caller == 'gatk_haplotypecaller') {
    //     // First we need to make this sequence dictionary that gatk requires
    //     seqdict = create_sequence_dict(ch_reference)
    //     // and add this one to the variant_ch
    //     variant_ch = variant_ch
    //         .combine(seqdict.dict)
    //         .map { sample_id, cram, crai, reference, index_files, chunk_id, interval, ref_dict ->
    //             return tuple(sample_id, cram, crai, reference, index_files + [ref_dict], chunk_id, interval)
    //         }
    //     haplotype_caller(variant_ch)
    //     sample_vcfs_ch = haplotype_caller.out.vcf
    // }
    
    // // ========================================================================= \\
    // // GATK HaplotypeCaller in ERC mode + CombineGVCFs + GenotypeGVCFs for joint genotyping
    // // ======================================================================== \\
    // if (params.variant_caller == 'gatk_joint') {
    //     // First we need to make this sequence dictionary that gatk requires
    //     seqdict = create_sequence_dict(ch_reference)
    //     // and add this one to the variant_ch
    //     variant_ch = variant_ch
    //         .combine(seqdict.dict)
    //         .map { sample_id, cram, crai, reference, index_files, chunk_id, interval, ref_dict ->
    //             return tuple(sample_id, cram, crai, reference, index_files + [ref_dict], chunk_id, interval)
    //         }
    //     haplotype_caller_erc(variant_ch)
    //     sample_vcfs_ch = haplotype_caller_erc.out.vcf

    //     // Extract a single GVCF and its index per sample.
    //     // Note: `haplotype_caller_erc` emits a file bundle via `*.vcf.gz*` (vcf + index files),
    //     // which would otherwise become nested lists after `groupTuple()`.
    //     per_sample_gvcf_ch = sample_vcfs_ch
    //         .map { row ->
    //             def region_id = row[0]
    //             def sample_id = row[1]
    //             def bundle    = row[2]
    //             def files     = (bundle instanceof List) ? bundle : [bundle]

    //             def gvcf = files.find { it.getName().endsWith('.vcf.gz') }
    //             def idx  = files.find { it.getName().endsWith('.tbi') } ?: files.find { it.getName().endsWith('.csi') }

    //             if( !gvcf ) throw new IllegalStateException("Missing .vcf.gz in bundle for ${sample_id} region ${region_id}: ${files}")
    //             if( !idx  ) throw new IllegalStateException("Missing index (.tbi/.csi) in bundle for ${sample_id} region ${region_id}: ${files}")

    //             tuple(region_id, sample_id, gvcf, idx)
    //         }

    //     // Barrier: emits one item per region only when all upstream tasks are done (unless `size:` is set).
    //     grouped_sample_vcfs_ch = per_sample_gvcf_ch
    //         .groupTuple(by: 0)
    //         .map { region_id, samples, gvcfs, idxs ->
    //             def zipped = [samples, gvcfs, idxs].transpose()
    //                 .sort { a, b -> a[0] <=> b[0] } // sort by sample_id
    //             def (samples_sorted, gvcfs_sorted, idxs_sorted) = zipped.transpose()
    //             tuple(region_id, samples_sorted, gvcfs_sorted, idxs_sorted)
    //         }

    //     // Attach the interval string for -L and then run CombineGVCFs.
    //     combine_in_ch = grouped_sample_vcfs_ch
    //         .join(refintervals_ch, by: 0)
    //         .combine(ch_reference)
    //         .combine(seqdict)
    //         .combine(samtools_index.out.reference_fai)
    //         .map { region_id, _samples, gvcfs, idxs, region, reference, seq_dict, ref_index ->
    //             tuple(region_id, region, reference, gvcfs, idxs, ref_index, seq_dict)
    //         }

    //     // run through combinegvcfs
    //     combined_ch = combine_gvcfs(combine_in_ch)
    //     // and directly into genotype
    //     combined_ch = combined_ch
    //         .combine(ch_reference)
    //         .combine(seqdict)
    //         .combine(samtools_index.out.reference_fai)
    //         .map { region_id, region, gvcf, gvcf_idx, reference, seq_dict, ref_index ->
    //             tuple(region_id, region, reference, gvcf, gvcf_idx, seq_dict, ref_index)
    //         }
    //     region_vcf_ch = genotype_gvcfs(combined_ch)

    // }

    // // ========================================================================= \\
    // // Freebayes
    // // ========================================================================= \\
    // if (params.variant_caller == 'freebayes') {
    //     freebayes(variant_ch)
    //     sample_vcfs_ch = freebayes.out.vcf
    // }

    // // Unless we're doing the gatk joint route, merge the per-sample vcf's into one per region here.
    // // to have a consistent output for the filtering steps regardless of caller choice.

    // if (params.variant_caller != 'gatk_joint') {
    //     per_sample_vcf_ch = sample_vcfs_ch
    //        .map { row ->
    //        def region_id = row[0]
    //        def sample_id = row[1]
    //        // handle both the case where the vcf and index are emitted as a bundle (list) or separately
    //        def rest  = (row.size() > 3) ? row[2..-1] : [ row[2] ]
    //        def files = rest.flatten().findAll { it != null }
    //        // now locate the vcf and index files within the emitted files
    //        def vcf = files.find { it.getName().endsWith('.vcf.gz') }
    //        def tbi = files.find { it.getName().endsWith('.tbi') } ?: null
    //        def csi = files.find { it.getName().endsWith('.csi') } ?: null
    //        // combine indices into a list
    //        def indices = [tbi, csi].findAll { it != null }
    //           tuple(region_id, sample_id, vcf, tbi, csi)
    //        }
        
    //     // now combine them per region
    //     // Now grouping produces 5-tuple: (region_id, [samples], [taxa], [vcfs], [idxs])
    //     bcftools_merge_in_ch = per_sample_vcf_ch
    //         .groupTuple(by: 0)
    //         .map { region_id, samples, vcfs, tbi, csi ->
    //             def zipped = [samples, vcfs, tbi, csi].transpose()
    //                 .sort { a, b -> a[0] <=> b[0] }  // sort by sample_id
    //             def (samples_sorted, vcfs_sorted, _tbi_sorted, csi_sorted) = zipped.transpose()
    //             tuple(region_id, samples_sorted, vcfs_sorted, csi_sorted, "raw")
    //         }

    //     bcftools_merge_in_ch
    //         .branch {
    //             region_id, samples, vcfs, csi, category ->
    //             multi:  samples.size() > 1
    //             single: true
    //         }
    //         .set { raw_merge_branch }

    //     single_raw_vcf_ch = raw_merge_branch.single
    //         .map { region_id, samples, vcfs, csi, category ->
    //             tuple(region_id, vcfs[0], csi[0])
    //         }

    //     region_vcf_ch = bcftools_merge(raw_merge_branch.multi).mix(single_raw_vcf_ch)
    // }

    // // run statistics on the raw vcfs for each region
    // raw_stats_ch = region_vcf_ch.
    //     map {
    //         region_id, vcf, idx ->
    //         tuple(region_id, vcf, idx, "raw")
    //     }
    // raw_stats = raw_vcf_stats(raw_stats_ch)
    
    // // Collect and combine all stat files for the combine_stats process
    // raw_stats.ab_dp.collect().toList()
    //     .combine(raw_stats.qual_fmiss_maf.collect().toList())
    //     .combine(raw_stats.rec_counts.collect().toList())
    //     .combine(raw_stats.sample_stats.collect().toList())
    //     .map { ab_dp_files, gt_files, rec_files, sample_files -> 
    //         return tuple(ab_dp_files, gt_files, sample_files, rec_files, "raw")
    //     }
    //     .set { combine_raw_stats_input }
    
    // raw_stats = combine_raw_vcf_stats(combine_raw_stats_input)
    
    // // if the raw vcf files are to be stored, concatenate now and output
    // if (params.store_raw_vcf) {
    //     // region_vcf_ch emits: (region_id, vcf, idx)
    //     all_region_vcfs_sorted_ch = region_vcf_ch
    //         .toSortedList { a, b -> a[0] <=> b[0] }   // sort by region_id
    //         .map { rows ->
    //             def vcfs = rows.collect { row -> row[1] }
    //             def idxs = rows.collect { row -> row[2] }
    //             tuple(vcfs, idxs, "raw")
    //         }

    //     // concatenate raw variants
    //     bcftools_concat(all_region_vcfs_sorted_ch)
    // }
    // // generate bedfiles with homozygous reference for each region and sample
    // invariants_in = sample_depth_cutoffs_autosomes
    //     .flatMap()
    //     .map { sample, mean_depth, min_depth, max_depth, category ->
    //         tuple(sample, min_depth, max_depth)
    //     }
    //     .combine(region_vcf_ch)
    //     .combine(sex_assignments, by: 0)
    //     .map() { sample, min_depth, max_depth, region_id, vcf, idx, x_z_depth, y_w_depth, depth_ratio, sex ->
    //             tuple(region_id, sample, min_depth, max_depth, vcf, idx, sex)
    //      }
    //      .combine(sex_linked_contigs.map { scaffold_list -> [scaffold_list]} )

    
    // homref = parse_invariants(invariants_in)
    // homref_combine = homref
    //     .groupTuple()
    //     .combine(reference_fai)

    // combine_homref(homref_combine)
    // // ========================================================================= \\
    // // ######################################################################### \\
    // //                              Filter genotypes
    // // ######################################################################### \\
    // // ========================================================================= \\


    // // we'll apply a bunch of separate filter steps depending on variant type and user parameters
    // // Finally, we'll output filtered SNPs, filtered indels and a bed file with monomorphic regions per sample
    
    // // first run vcf files through bcftools normalize
    raw_vcfs
        .combine(ch_reference)
        .combine(samtools_index.out.reference_fai)
        .combine(samtools_index.out.reference_gzi)
        .map { region_id, vcf, idx, reference, reference_fai, reference_gzi ->
            tuple(region_id, vcf, idx, reference, reference_fai, reference_gzi)
        }
        .set { bcftools_norm_in_ch }
    
    normalized_vcfs_ch = bcftools_norm(bcftools_norm_in_ch)

    // run through bcftools select to pull out SNPs, indels and invariant sites separately
    snp_ch = select_snps(normalized_vcfs_ch)
    indel_ch = select_indels(normalized_vcfs_ch)

    // check if custom filter critera are set, else go with the defaults depending on the variant caller
    def snp_filter_expr
    def indel_filter_expr
    if (params.snp_filter_expression) {
        snp_filter_expr = params.snp_filter_expression
        indel_filter_expr = params.indel_filter_expression
    }
    else {
        if (params.variant_caller == 'gatk_joint' || params.variant_caller == 'gatk_haplotypecaller') {
            snp_filter_expr = params.snp_filter_expression_gatk
            indel_filter_expr = params.indel_filter_expression_gatk
        }
        else if (params.variant_caller == 'freebayes') {
            snp_filter_expr = params.snp_filter_expression_freebayes
            indel_filter_expr = params.indel_filter_expression_freebayes
        }
        else if (params.variant_caller == 'bcftools') {
            snp_filter_expr = params.snp_filter_expression_bcftools
            indel_filter_expr = params.indel_filter_expression_bcftools
        }
    }

    snps_filtered_step1 = bcftools_filter_snps(snp_ch, snp_filter_expr, "snps")
    indels_filtered_step1 = bcftools_filter_indels(indel_ch, indel_filter_expr, "indels")

    // Prepare depth cutoff channel for filtering: reformat to put sample_id at the right position for joining
    depth_cutoffs_for_filter = depth_cutoffs
        .map {
            sample_id, min_dp, max_dp, sex_assignment, _depths ->
                tuple (min_dp, max_dp, sample_id, sex_assignment)
        }

    // Build SNP filtering input: combine intervals, callable regions, filtered SNPs, and depth cutoffs
    snps_filter_input = refintervals_ch
        .combine(callable_regions.out.callable)
        .combine(snps_filtered_step1, by: 0)
        .combine(depth_cutoffs_for_filter, by: 2)
        .map {
            sample_id, region_id, region, bedfile, vcf, idx, min_depth, max_depth, sex_assignment ->
                tuple(sample_id, region_id, region, bedfile, vcf, idx, min_depth, max_depth, sex_assignment)
        }

    ab_dp_filter_snps(snps_filter_input, sex_limited_contigs, 'snps')

    // Build indel filtering input: same structure as SNPs
    indels_filter_input = refintervals_ch
        .combine(callable_regions.out.callable)
        .combine(indels_filtered_step1, by: 0)
        .combine(depth_cutoffs_for_filter, by: 2)
        .map {
            sample_id, region_id, region, bedfile, vcf, idx, min_depth, max_depth, sex_assignment ->
                tuple(sample_id, region_id, region, bedfile, vcf, idx, min_depth, max_depth, sex_assignment)
        }

    ab_dp_filter_indels(indels_filter_input, sex_limited_contigs, 'indels')

    // merge samples together
    // Helper function to prepare variant channel for merging: group by region, sort samples, add reference files
    def prepare_merge_input(filtered_vcf_channel) {
        return filtered_vcf_channel
            .groupTuple(by: 0)
            .map {
                region_id, samples, vcfs, idxs ->
                    // Sort samples by ID to ensure consistent order for reproducible merging
                    def zipped = [samples, vcfs, idxs].transpose()
                        .sort { a, b -> a[0] <=> b[0] }
                    def (samples_sorted, vcfs_sorted, idxs_sorted) = zipped.transpose()
                    tuple(region_id, samples_sorted, vcfs_sorted, idxs_sorted)
            }
            .combine(samtools_index.out.reference_fasta)
            .combine(samtools_index.out.reference_fai)
            .combine(samtools_index.out.reference_gzi)
            .map { region_id, samples, vcfs, idxs, reference_fasta, reference_fai, reference_gzi ->
                tuple(region_id, samples, vcfs, idxs, reference_fasta, reference_fai, reference_gzi)
            }
    }

    // Prepare SNP and indel channels for merging with consistent format
    snps_to_merge = prepare_merge_input(ab_dp_filter_snps.out.ab_dp_filtered_vcf)
    indels_to_merge = prepare_merge_input(ab_dp_filter_indels.out.ab_dp_filtered_vcf)
    
    bcftools_merge_snps(snps_to_merge, 'snps')
    bcftools_merge_indels(indels_to_merge, 'indels')

    // remove sites where all samples are filtered out, or if a user specified fraction of missing genotypes is given
    bcftools_fmiss_maf_filtered_snps = bcftools_filter_fmiss_maf_snps(bcftools_merge_snps.out.vcf, 'snps')
    bcftools_fmiss_maf_filtered_indels = bcftools_filter_fmiss_maf_indels(bcftools_merge_indels.out.vcf, 'indels')
    

    // stats
    filtered_snp_stats_out = filtered_snp_stats(bcftools_fmiss_maf_filtered_snps, 'snps')
    filtered_indels_stats = filtered_indel_stats(bcftools_fmiss_maf_filtered_indels, 'indels')

    // combine stats
    combined_stats_snps = combine_filtered_snps_stats(
        filtered_snp_stats_out.ab_dp.collect(),
        filtered_snp_stats_out.qual_fmiss_maf.collect(),
        filtered_snp_stats_out.sample_stats.collect(),
        filtered_snp_stats_out.rec_counts.collect(),
        'filtered_snps'
    )

    combined_stats_indels = combine_filtered_indels_stats(
        filtered_indels_stats.ab_dp.collect(),
        filtered_indels_stats.qual_fmiss_maf.collect(),
        filtered_indels_stats.sample_stats.collect(),
        filtered_indels_stats.rec_counts.collect(),
        'filtered_indels'
    )

    // Helper function to sort VCF channel by region_id and extract sorted file lists
    def sort_and_extract_vcfs(vcf_channel) {
        return vcf_channel
            .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by region_id for deterministic output
            .map { sorted_list ->
                def vcfs = sorted_list.collect { it[1] }
                def idxs = sorted_list.collect { it[2] }
                tuple(vcfs, idxs)
            }
    }

    // Sort SNP and indel VCFs by region for concatenation in consistent order
    sorted_snps = sort_and_extract_vcfs(bcftools_merge_snps.out.vcf)
    sorted_indels = sort_and_extract_vcfs(bcftools_merge_indels.out.vcf)

    // Concatenate filtered SNPs and indels into genome-wide VCFs
    bcftools_concat_snps(sorted_snps, 'filtered_snps')
    bcftools_concat_indels(sorted_indels, 'filtered_indels')

    // finalize a bunch of different mask files, which may be useful downstream
    // Prepare input for mask finalization: combine regions, callable sites, and variant VCFs for ALL samples
    finalize_masks_input = refintervals_ch
        .combine(callable_regions.out.callable)
        .combine(bcftools_merge_snps.out.vcf, by: 0)
        .combine(bcftools_merge_indels.out.vcf, by: 0)
        .map {
            region_id, region, sample_id, callable_bed, snp_vcf, _snp_idx, indel_vcf, _indel_idx ->
                tuple(sample_id, region_id, region, callable_bed, snp_vcf, indel_vcf)
        }

    region_masks = finalize_masks(finalize_masks_input)
    
    // Extract and group by sample for final merging
    region_homrefs = region_masks.homref_invariants
        .groupTuple(by: 0)
    region_totalmask = region_masks.mappability_mask
        .groupTuple(by: 0)
    region_snpmask = region_masks.mappability_mask_snps
        .groupTuple(by: 0)

    // Concatenate regional masks into genome-wide files for output
    homrefs = combine_homref_invariants(region_homrefs, reference_fai, 'homref_invariants')
    total_mask = combine_mappability_masks(region_totalmask, reference_fai, 'total_mask')
    snp_mask = combine_mappability_masks_snps(region_snpmask, reference_fai, 'snp_mask')



    //======================================================================== \\
    // ######################################################################### \\
    //                             Generate sample output report                \\
    // ######################################################################### \\

    // sample_stats.view {
    //     row -> "sample stats: ${row}"
    // }
    // sample_depths.view {
    //     row -> "sample depths: ${row}"
    // }

    // If no downsampling was performed, just push NA values for these columns before sending it to the parse
     if (!params.downsample_bams) {
         sample_stats = sample_stats
            .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment ->
                tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, "NA", "NA", "NA", "NA", "NA")
            }
     }
    // remove the extra sex_assignment col from the downsampeld sample stats
    sample_stats = sample_stats
        .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp, _ds_ratio, _ds_sex_assignment ->
            tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp)
        }

    // run through a parser to fetch some summary statistics
    parse_summary_stats(pops
        .combine(
        homrefs.bedfile,
        by: 0
        )
        .combine(
            total_mask.bedfile,
            by: 0
        )
        .combine(
            snp_mask.bedfile,
            by: 0
        )
        .combine(
            sample_stats,
            by: 0
        )
        .combine(
            combined_stats_snps.combined_summary_statistics.toList(),
        )
        .combine(combined_stats_indels.combined_summary_statistics.toList()),
        sex_chrom_system
    )

    // and combine these
    combine_summary_tables(parse_summary_stats.out.summary_statistics.collect())
        

    

    // ========================================================================= \\
    // ######################################################################### \\
    //                             Publish outputs
    // ######################################################################### \\

    publish:
    // FastQC / MultiQC reports on raw and trimmed reads
    multiqc_rawreads_report   = multiqc_rawreads.out.report
    multiqc_cleanreads_report = multiqc_cleanreads.out.report

    // Per-sample CRAM QC and the aggregated MultiQC report
    // cramqc_reports      = cramqc_ch.qc_report
    // cram_multiqc_report = cram_multiqc.out.report
    qualimap_pre_dedup_reports = qualimap_pre_dedup.out.qualimap_report
    qualimap_reports = qualimap.out.qualimap_report
    qualimap_downsampled_reports = qualimap_downsampled.out.qualimap_report

    // Deduplicated CRAM files and duplicate-marking metrics
    cram_files   = bams_for_calling_ch
    cram_metrics = samtools_markdups.out.metrics

    // Reference genome
    reference_genome = bwa_index.out.reference

    // downsampled bam files
    downsampled_cram_files = bams_for_calling_ch
        .map { sample_id, cram, crai ->
            tuple(sample_id, cram, crai)
        }

    // // Variant statistics (raw calls)
    // raw_vcf_stats = combine_raw_vcf_stats.out.combined_summary_statistics

    // // Raw concatenated VCF (only written when params.store_raw_vcf is true)
    raw_vcf = bcftools_concat_raw.out.vcf
    filtered_snps = bcftools_concat_snps.out.vcf
    filtered_indels = bcftools_concat_indels.out.vcf
    // Raw coverage bed files (only for testing)
    //coverage_beds = parse_region_depths.out.sample_depth_beds

    // callable regions
    callable_regions = total_mask.bedfile
    snpable_regions = snp_mask.bedfile
    invariant_calls = homrefs.bedfile

    // filtered snp and indel stats
    filtered_snps_stats = combined_stats_snps.combined_summary_statistics
    filtered_indel_stats = combined_stats_indels.combined_summary_statistics

    raw_variant_stats = combined_raw_stats.combined_summary_statistics

    summary_statistics = combine_summary_tables.out.table

    damage_profiles = damage_profiler.out.damage_reports


}

output {

    multiqc_rawreads_report {
        path "00_reports/00_fastqc/00_raw_reads"
    }
    multiqc_cleanreads_report {
        path "00_reports/00_fastqc/01_clean_reads"
    }
    qualimap_pre_dedup_reports {
        path "00_reports/01_qualimap/00_qualimap_pre_markdups"
    }
    qualimap_reports {
        path "00_reports/01_qualimap/01_qualimap_post_markdups"
    }
    qualimap_downsampled_reports {
        path "00_reports/01_qualimap/02_qualimap_downsampled"
    }
    reference_genome {
        path "01_reference_genome"
    }
    cram_files {
        path "02_cramfiles"
    }
    downsampled_cram_files {
        path "02_cramfiles/downsampled"
    }
    cram_metrics {
        path "02_cramfiles"
    }
    // raw_vcf_stats {
    //     path "04_variant_stats"
    // }
    raw_vcf {
        path "03_genotypes/00_raw_variants"
        enabled params.store_raw_vcf
    }

    // coverage_beds {
    //     path "03_callable_regions/coverage_beds"
    // }
    filtered_snps_stats {
        path "00_reports/02_variantstats/01_filtered_snps"
    }
    filtered_indel_stats {
        path "00_reports/02_variantstats/02_filtered_indels"
    }
    raw_variant_stats {
        path "00_reports/02_variantstats/00_raw_variants"
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
        path "00_reports"
    }
    damage_profiles {
        path "00_reports/03_damage_profiles"
    }
}

/*
 * Completion handler
 */
workflow.onComplete {
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║                   Pipeline Complete                       ║
    ╠═══════════════════════════════════════════════════════════╣
    ║ Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}    ║
    ║ Duration  : ${workflow.duration}                          ║
    ║ Output    : ${params.outdir}                              ║
    ╚═══════════════════════════════════════════════════════════╝
    """.stripIndent()
}
