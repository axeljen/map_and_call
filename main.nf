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
include { clumpify_single } from './modules/bbmap/bbmap'
include { clumpify_paired } from './modules/bbmap/bbmap'
include { concat_reads } from './modules/concat_libs/concat_libs'
include { concat_collapsed } from './modules/concat_libs/concat_libs'

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
include { damage_profiler_rescale } from './modules/mapdamage/damageprofiler_rescale'
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
include { merge_historical_bams } from './modules/samtools/mergebams'
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
include { plot_variant_stats as plot_raw_stats } from './modules/variant_stats/plot_variant_stats'
include { plot_variant_stats as plot_snp_stats } from './modules/variant_stats/plot_variant_stats'
include { plot_variant_stats as plot_indel_stats } from './modules/variant_stats/plot_variant_stats'
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
 * Prepare variant channel for merging: group by region, sort samples, add reference files
 */
def prepare_merge_input(filtered_vcf_channel, reference_fasta, reference_fai, reference_gzi) {
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
        .combine(reference_fasta)
        .combine(reference_fai)
        .combine(reference_gzi)
        .map { rid, samples_s, vcfs_s, idxs_s, ref_fa, ref_fai, ref_gzi ->
            tuple(rid, samples_s, vcfs_s, idxs_s, ref_fa, ref_fai, ref_gzi)
        }
}

/*
 * Sort VCF channel by region_id and extract sorted file lists
 */
def sort_and_extract_vcfs(vcf_channel) {
    return vcf_channel
        .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by region_id for deterministic output
        .map { sorted_list ->
            def vcfs = sorted_list.collect { it[1] }
            def idxs = sorted_list.collect { it[2] }
            tuple(vcfs, idxs)
        }
}

/*
 * Main workflow
 */
workflow {
    main:
    
    // Setup sex chromosome system based on user parameters
    sex_config = setup_sex_chromosome_system()
    def sex_chrom_system = sex_config.sex_chrom_system
    def sex_linked_list = sex_config.sex_linked_list
    def sex_limited_list = sex_config.sex_limited_list
    def non_sex_limited_list = sex_config.non_sex_limited_list
    
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
    //                     1. PREPROCESSING QC: RAW READS
    // ═══════════════════════════════════════════════════════════════════════════════

    // FastQC on raw (untrimmed) reads
    ch_raw_reads_for_qc = ch_input.modern.map { sample_id, lane, _datatype, library, r1, r2 ->
            [sample_id, lane, library, [r1, r2], 'raw_reads']
        }
        // push historical reads through the same fastqc process
        .mix(ch_input.historical.map { sample_id, lane, _datatype, library, r1, r2 ->
            [sample_id, lane, library, [r1, r2], 'raw_reads']
        })
    fastqc_rawreads(ch_raw_reads_for_qc)
    // Collect all FastQC outputs and run MultiQC
    ch_multiqc_input = fastqc_rawreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'raw_reads') }

    multiqc_rawreads(ch_multiqc_input)

    // ═══════════════════════════════════════════════════════════════════════════════
    //                  2. READ PREPROCESSING & TRIMMING
    // ═══════════════════════════════════════════════════════════════════════════════

    // Trim modern reads with fastp (quality and adapter trimming)
    fastp(ch_input.modern)

    // Trim historical reads with AdapterRemoval (handles ancient DNA damage removal)
    adapterremoval(ch_input.historical)

    // now concatenate runs from the same library, prior to deduplication
    concat_pairs_in = fastp.out.reads
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, reads1, reads2 ->
            tuple(sample_id, library, datatypes[0], reads1, reads2)
        }
        // make sure that reads1 and reads2 are sorted in the same order
        .map { sample_id, library, datatype, reads1, reads2 ->
            def zipped = [reads1, reads2].transpose()
                .sort { a, b -> a[0].name <=> b[0].name }
            def (reads1_sorted, reads2_sorted) = zipped.transpose()
            tuple(sample_id, library, datatype, reads1_sorted, reads2_sorted)
        }
        .mix (adapterremoval.out.trimmed_pairs
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, reads1, reads2 ->
            tuple(sample_id, library, datatypes[0], reads1, reads2)
        }
        .map { sample_id, library, datatype, reads1, reads2 ->
            // sort reads in the same order as fastp output to ensure correct pairing during concatenation
            def zipped = [reads1, reads2].transpose()
                .sort { a, b -> a[0].name <=> b[0].name }
            def (reads1_sorted, reads2_sorted) = zipped.transpose()
            tuple(sample_id, library, datatype, reads1_sorted, reads2_sorted)
         }
        )
        // branch based on number of read pairs per library
        .branch { sample_id, library, datatype, reads1, reads2 ->
            single_lib: reads1.size() == 1
            multi_lib: reads1.size() > 1
        }

    concat_pairs = concat_reads(concat_pairs_in.multi_lib)
    
    // and then concatenate merged reads from adapterremoval
    concat_merged_in = adapterremoval.out.trimmed_collapsed
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, collapsed ->
            tuple(sample_id, library, datatypes[0], collapsed)
        }
        // sort collapsed reads to ensure consistent order for concatenation
        .map { sample_id, library, datatype, collapsed ->
            def sorted_collapsed = collapsed.sort { a, b -> a.name <=> b.name }
            tuple(sample_id, library, datatype, sorted_collapsed)
        }
        .branch { sample_id, library, datatype, collapsed ->
            single_lib: collapsed.size() == 1
            multi_lib: collapsed.size() > 1
        }
    concat_merged = concat_collapsed(concat_merged_in.multi_lib)

    if (params.premapping_dedup) {
        // if user wants to do pre-mapping deduplication, we can push the concatenated libraries through clumpify now
        // combine with samples that didn't need concatenation
        dedup_pairs_in = concat_pairs.reads_concat
            .mix(concat_pairs_in.single_lib)
        dedup_merged_in = concat_merged.collapsed_concat
            .mix(concat_merged_in.single_lib)

        // and push the concatenated libraries through clumpify
        clean_paired = clumpify_paired(dedup_pairs_in)
            
        clean_merged = clumpify_single(dedup_merged_in)
    } else {
        // if no pre-mapping deduplication, just use the concatenated libraries as the "clean" reads for mapping
        // single_lib items come from groupTuple and carry single-element lists — unwrap to plain paths
        clean_paired = concat_pairs.reads_concat
            .mix(concat_pairs_in.single_lib
                .map { sample_id, library, datatype, reads1, reads2 ->
                    tuple(sample_id, library, datatype, reads1[0], reads2[0])
                })
        clean_merged = concat_merged.collapsed_concat
            .mix(concat_merged_in.single_lib
                .map { sample_id, library, datatype, collapsed ->
                    tuple(sample_id, library, datatype, collapsed[0])
                })
    }


    // ═══════════════════════════════════════════════════════════════════════════════
    //                  3. PREPROCESSING QC: CLEAN READS
    // ═══════════════════════════════════════════════════════════════════════════════

    // FastQC on trimmed reads
    fastqc_cleanreads(clean_merged
        .map { sample_id, library, _datatype, collapsed -> tuple(sample_id, 'combined', library, [collapsed], 'clean')}
        .mix(clean_paired
        .map { sample_id, library, _datatype, reads1, reads2 -> tuple(sample_id, 'combined', library, [reads1, reads2], 'clean')}
        )
    )

    multiqc_cleanreads(fastqc_cleanreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'clean_reads') }
        )

    // ═══════════════════════════════════════════════════════════════════════════════
    //              4. REFERENCE GENOME PREPARATION & INDEXING
    // ═══════════════════════════════════════════════════════════════════════════════

    // Index reference genome and create intervals for parallel processing
    // check if the indices needed for bwa already exist, and if not, create them. 
    bwa_index_ch = bwa_index(ch_reference)
    faidx_and_chunks_ch = samtools_index(ch_reference)
    reference_fasta = faidx_and_chunks_ch.reference_fasta.first()
    reference_fai = faidx_and_chunks_ch.reference_fai.first()
    reference_gzi = faidx_and_chunks_ch.reference_gzi.first()
    
    // if user provides a list of scaffolds upon which to call variants, fetch this
    if (params.scaffold_list) {
        scaffolds_ch = channel.fromPath(params.scaffold_list, checkIfExists: true)
            .splitCsv(header: false)
            .map { row -> row[0] }
            .collect()
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

    // define a list of autosomes and potentially xz/yw chromosomes if given
    scaffolds = scaffolds_ch
            .flatten()
            .branch { scaffold ->
                autosomes: [params.y_scaffolds, params.w_scaffolds, params.x_scaffolds, params.z_scaffolds].flatten().contains(scaffold) == false
                sex_limited: [params.y_scaffolds, params.w_scaffolds].flatten().contains(scaffold)
                non_sex_limited: [params.x_scaffolds, params.z_scaffolds].flatten().contains(scaffold)
            }

    // now generate reference intervals - these will be used for parallelization in the variant calling 
    // and filtering step, so not used in the mapping/preprocessing
    reference_intervals = dochunks(samtools_index.out.reference_fai, params.chunk_size, scaffolds_ch.flatten())

    // Handle genomic intervals for variant calling
    refintervals_ch = reference_intervals
        .map { row -> row?.trim() }
        .flatMap { row -> row ? row.split(/\r?\n/) as List : [] }
        .filter { row -> row }
        .collect()
        .flatMap { intervals ->
            intervals.withIndex().collect { interval, idx -> tuple(idx + 1, interval) }
        }

    // ═══════════════════════════════════════════════════════════════════════════════
    //                    5. MAPPING: MODERN SAMPLES
    // ═══════════════════════════════════════════════════════════════════════════════

    // Map trimmed modern reads to reference
    modern_pairs_to_map = clean_paired
        .filter { _sample_id, _library, datatype, _reads1, _reads2 -> datatype == '1' }
        .combine(bwa_index_ch.reference)
    // Choose mapping approach based on user parameter
    if (params.mapper == 'bwa_mem') {
        rawbam_ch = bwa_mem(modern_pairs_to_map)
    }
    else {
        error "Unsupported mapper specified: ${params.mapper}. Currently only 'bwa_mem' is supported."
    }

    // ═══════════════════════════════════════════════════════════════════════════════
    //                    6. MAPPING: HISTORICAL SAMPLES
    // ═══════════════════════════════════════════════════════════════════════════════


    // Map paired-end reads from historical samples (optional)
    if (params.map_historical_pairs) {
        historical_pairs_to_map = clean_paired
            .filter { _sample_id, _library, datatype, _reads1, _reads2 -> datatype == '2' }
            .combine(bwa_index_ch.reference)
        historical_pairbams_ch = map_historical(historical_pairs_to_map)
    }
    else {
        // Create an empty channel if paired mapping is disabled
        historical_pairbams_ch = channel.empty()
    }
    
    historical_merged_to_map = clean_merged
        .filter { _sample_id, _library, datatype, _collapsed -> datatype == '2' }
        .combine(bwa_index_ch.reference)

    historical_bams_ch = map_merged(historical_merged_to_map)

    // if mapping historical paired end reads, merge within each library
    if (params.map_historical_pairs) {
        merge_historical = historical_pairbams_ch
            .mix(historical_bams_ch)
            .groupTuple(by: [0,1,2])
            .map { sample_id, library, datatype, bam_paths, _bam_indices ->
                return tuple(sample_id, library, datatype, bam_paths)
            }
        historical_bams_merged_ch = merge_historical_bams(merge_historical)
        historical_bams_ch = historical_bams_merged_ch
    }

    

    // ═══════════════════════════════════════════════════════════════════════════════
    //            7. BAM POST-PROCESSING: DEDUP (optional) and MERGING
    // ═══════════════════════════════════════════════════════════════════════════════

    // Select the correct historical bam channel
    //all_historical_bams = params.map_historical_pairs ? //historical_bams_merged_ch : historical_bams_ch

    // Mix all per-library bams [sample_id, library, datatype, bam, bai],
    // group by sample, merge if multiple libraries, then branch by datatype
    per_sample_bams = rawbam_ch
        .mix(historical_bams_ch)
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
    // post mapping dedup if needed
    if (params.postmapping_dedup) {
        dedup_bams = samtools_markdups(all_sample_bams.historical.mix(all_sample_bams.modern))
        all_sample_bams = dedup_bams.bam.branch { _sample_id, datatype, _bam, _bai ->
            modern:     datatype == '1'
            historical: datatype == '2'
        }
    }
        

    // ═══════════════════════════════════════════════════════════════════════════════
    //        8. HISTORICAL DNA: DAMAGE PROFILING & RESCALING
    // ═══════════════════════════════════════════════════════════════════════════════

    // run damageprofiler/damageprofiler_rescale on historical bams
    if (params.damageprofiler_rescale) {
        println "Running damage profiling and rescaling on historical samples..."
        rescaled_bams = damage_profiler_rescale(
            all_sample_bams.historical,
            ch_reference
        )
        final_bam_ch = all_sample_bams.modern
            .map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) }
            .mix(rescaled_bams.rescaled_bam)
    } else {
        println "Damage profiling and rescaling is disabled. Skipping this step and using original BAMs for downstream analyses."
        rescaled_bams = damage_profiler(
            all_sample_bams.historical.map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) },
        ch_reference
        )
        final_bam_ch = all_sample_bams.modern
            .map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) }
            .mix(all_sample_bams.historical.map { sample_id, _datatype, bam, bai -> tuple(sample_id, bam, bai) })
    }


    // ═══════════════════════════════════════════════════════════════════════════════
    //              9. MAPPING QC: BAM QUALITY ASSESSMENT
    // ═══════════════════════════════════════════════════════════════════════════════

    qualimap(final_bam_ch)

    // ═══════════════════════════════════════════════════════════════════════════════
    //       10. DEPTH PROFILING & SEX ASSIGNMENT FOR FILTERING
    // ═══════════════════════════════════════════════════════════════════════════════

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
        // add the sample bedfile to this one
        .combine(parse_region_depths.out.sample_depth_beds, by: 0)

    mappability_mask = callable_regions(depth_cutoffs, sex_limited_contigs, non_sex_limited_contigs, reference_fai)

    // Normalize sample_stats to consistent tuple size before variant calling
    // If no downsampling was performed, add NA values for downsampling columns
    if (!params.downsample_bams) {
        sample_stats = sample_stats
            .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment ->
                tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, "NA", "NA", "NA", "NA", "NA")
            }
    }
    // Remove the extra sex_assignment col from the downsampled sample stats
    sample_stats = sample_stats
        .map { sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp, _ds_ratio, _ds_sex_assignment ->
            tuple(sample_id, autosomal_dp, non_sex_limited_dp, sex_limited_dp, ratio, sex_assignment, ds_autosomal_dp, ds_non_sex_limited_dp, ds_sex_limited_dp)
        }

    // ═══════════════════════════════════════════════════════════════════════════════
    //           11-12. VARIANT CALLING: SAMPLE & POPULATION SETUP & VARIANT CALLING
    // ═══════════════════════════════════════════════════════════════════════════════

    if (!params.skip_variant_calling) {

    // Parse population assignments for joint genotyping (if provided)
    // if a popfile was provided, get this info, otherwise create a pops channel with each sample as its own population (i.e. independent calling)
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
    
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //             12. VARIANT CALLING
    // 
    // Supports four two different variant calling approaches:
    //   - bcftools mpileup: Fast, lightweight variant discovery
    //   - freebayes: Joint calling across all samples
    // Likely to be implemented soon:
    //   - gatk_haplotypecaller: Per-sample calling (commented out)
    //   - gatk_joint: Joint calling with GVCFs (commented out)
    // ═══════════════════════════════════════════════════════════════════════════════

    if ( !['freebayes','bcftools'].contains(params.variant_caller) ) {
        error "Unsupported variant caller specified: ${params.variant_caller}. Must be one of freebayes or bcftools."
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

    // Call variants using selected variant caller

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

    // Generate variant statistics PDF report for raw variants
    plot_raw_stats(combined_raw_stats.combined_summary_statistics.collect(), 'raw_variants')

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

    // Get filter expressions based on variant caller type or user custom expressions
    filter_expressions = get_filter_expressions(
        params.variant_caller,
        params.snp_filter_expression,      // custom SNP filter (or null)
        params.indel_filter_expression     // custom indel filter (or null)
    )

    snps_filtered_step1 = bcftools_filter_snps(snp_ch, filter_expressions.snp_filter_expr, "snps")
    indels_filtered_step1 = bcftools_filter_indels(indel_ch, filter_expressions.indel_filter_expr, "indels")

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
    // Prepare SNP and indel channels for merging with consistent format
    snps_to_merge = prepare_merge_input(ab_dp_filter_snps.out.ab_dp_filtered_vcf, reference_fasta, reference_fai, reference_gzi)
    indels_to_merge = prepare_merge_input(ab_dp_filter_indels.out.ab_dp_filtered_vcf, reference_fasta, reference_fai, reference_gzi)
    
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

    // Generate variant statistics PDF reports for filtered SNPs and indels
    plot_snp_stats(combined_stats_snps.combined_summary_statistics.collect(), 'filtered_snps')
    plot_indel_stats(combined_stats_indels.combined_summary_statistics.collect(), 'filtered_indels')

    // Helper function to sort VCF channel by region_id and extract sorted file lists
    // Sort SNP and indel VCFs by region for concatenation in consistent order
    sorted_snps = sort_and_extract_vcfs(bcftools_fmiss_maf_filtered_snps.vcf)
    sorted_indels = sort_and_extract_vcfs(bcftools_fmiss_maf_filtered_indels.vcf)

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
        

    } // End of variant calling section


    // ========================================================================= \\
    // ######################################################################### \\
    //                             Publish outputs
    // ######################################################################### \\

    publish:
    // FastQC / MultiQC reports on raw and trimmed reads
    multiqc_rawreads_report   = multiqc_rawreads.out.report
    multiqc_cleanreads_report = multiqc_cleanreads.out.report

    // clean reads
    clean_reads = clean_merged.mix(clean_paired)

    // Per-sample CRAM QC and the aggregated MultiQC report
    // cramqc_reports      = cramqc_ch.qc_report
    // cram_multiqc_report = cram_multiqc.out.report
    qualimap_reports = qualimap.out.qualimap_report
    qualimap_downsampled_reports = params.downsample_bams ? qualimap_downsampled.out.qualimap_report : channel.empty()

    // Deduplicated CRAM files and duplicate-marking metrics
    raw_crams = final_bam_ch
        
    cram_metrics = params.postmapping_dedup ? samtools_markdups.out.metrics : channel.empty()

    // Reference genome
    reference_genome = bwa_index.out.reference

    // downsampled bam files (or use non-downsampled if downsampling is disabled)
    downsampled_cram_files = params.downsample_bams ? samtools_downsample.out.downsampled_bam : bams_for_calling_ch

    // // Raw concatenated VCF (only written when params.store_raw_vcf is true)
    raw_vcf = bcftools_concat_raw.out.vcf
    filtered_snps = bcftools_concat_snps.out.vcf
    filtered_indels = bcftools_concat_indels.out.vcf

    // callable regions
    callable_regions = total_mask.bedfile
    snpable_regions = snp_mask.bedfile
    invariant_calls = homrefs.bedfile

    // raw variant stats
    raw_variant_stats = combined_raw_stats.combined_summary_statistics
    raw_vcf_stats_plot = plot_raw_stats.out.report

    // filtered snp and indel stats
    filtered_snps_stats = combined_stats_snps.combined_summary_statistics
    filtered_snps_stats_plot = plot_snp_stats.out.report
    filtered_indel_stats = combined_stats_indels.combined_summary_statistics
    filtered_indel_stats_plot = plot_indel_stats.out.report

    summary_statistics = combine_summary_tables.out.table

    damage_profiles = rescaled_bams.damage_reports
}

output {

    multiqc_rawreads_report {
        path "00_reports/00_fastqc/00_raw_reads"
    }
    multiqc_cleanreads_report {
        path "00_reports/00_fastqc/01_clean_reads"
    }
    qualimap_reports {
        path "00_reports/01_qualimap/01_qualimap_post_markdups"
    }
    qualimap_downsampled_reports {
        enabled params.downsample_bams
        path "00_reports/01_qualimap/02_qualimap_downsampled"
    }
    reference_genome {
        path "01_reference_genome"
    }
    raw_crams {
        path "02_cramfiles"
    }
    clean_reads {
        enabled params.store_cleanreads
        path "clean_reads"
    }
    downsampled_cram_files {
        enabled params.downsample_bams
        path "02_cramfiles/downsampled"
    }
    cram_metrics {
        enabled params.postmapping_dedup
        path "02_cramfiles"
    }
    // raw_vcf_stats {
    //     path "04_variant_stats"
    // }
    raw_vcf {
        path "03_genotypes/00_raw_variants"
        enabled params.store_raw_vcf && !params.skip_variant_calling
    }
    raw_variant_stats {
        path "00_reports/02_variantstats/00_raw_variants"
        enabled !params.skip_variant_calling
    }
    raw_vcf_stats_plot {
        path "00_reports/02_variantstats/00_raw_variants"
        enabled !params.skip_variant_calling
    }
    filtered_snps_stats {
        path "00_reports/02_variantstats/01_filtered_snps"
        enabled !params.skip_variant_calling
    }
    filtered_snps_stats_plot {
       path "00_reports/02_variantstats/01_filtered_snps"
        enabled !params.skip_variant_calling
    }
    filtered_indel_stats {
        path "00_reports/02_variantstats/02_filtered_indels"
        enabled !params.skip_variant_calling
    }
    filtered_indel_stats_plot {
        path "00_reports/02_variantstats/02_filtered_indels"
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
        path "00_reports"
        enabled !params.skip_variant_calling
    }
    damage_profiles {
        path "00_reports/03_damage_profiles"
    }
}

