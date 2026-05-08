#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: VARIANT_FILTERS
//
// Purpose: Normalize, filter, and merge variant calls; generate mask files
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { bcftools_norm } from '../modules/bcftools/bcftools_norm'
include { select_snps } from '../modules/bcftools/select_snps'
include { select_indels } from '../modules/bcftools/select_indels'
include { bcftools_filter as bcftools_filter_snps } from '../modules/bcftools/bcftools_filter'
include { bcftools_filter as bcftools_filter_indels } from '../modules/bcftools/bcftools_filter'
include { callability_filter as callability_filter_snps } from '../modules/bedtools/callability_filter'
include { callability_filter as callability_filter_indels } from '../modules/bedtools/callability_filter'
include { ab_filter as ab_filter_snps } from '../modules/custom_variant_filters/ab_filter'
include { ab_filter as ab_filter_indels } from '../modules/custom_variant_filters/ab_filter'
include { ab_dp_filter as ab_dp_filter_snps } from '../modules/custom_variant_filters/ab_dp_filter'
include { ab_dp_filter as ab_dp_filter_indels } from '../modules/custom_variant_filters/ab_dp_filter'
include { bcftools_merge as bcftools_merge_snps } from '../modules/bcftools/bcftools_merge'
include { bcftools_merge as bcftools_merge_indels } from '../modules/bcftools/bcftools_merge'
include { bcftools_filter_fmiss_maf as bcftools_filter_fmiss_maf_snps } from '../modules/bcftools/bcftools_filter_fmiss_maf'
include { bcftools_filter_fmiss_maf as bcftools_filter_fmiss_maf_indels } from '../modules/bcftools/bcftools_filter_fmiss_maf'
include { bcftools_concat as bcftools_concat_snps } from '../modules/bcftools/bcftools_concat'
include { bcftools_concat as bcftools_concat_indels } from '../modules/bcftools/bcftools_concat'
include { vcf_stats as filtered_snp_stats } from '../modules/variant_stats/vcf_stats_process'
include { vcf_stats as filtered_indel_stats } from '../modules/variant_stats/vcf_stats_process'
include { combine_stats as combine_filtered_snps_stats } from '../modules/variant_stats/combine_stats'
include { combine_stats as combine_filtered_indels_stats } from '../modules/variant_stats/combine_stats'
include { plot_variant_stats as plot_snp_stats } from '../modules/variant_stats/plot_variant_stats'
include { plot_variant_stats as plot_indel_stats } from '../modules/variant_stats/plot_variant_stats'
include { finalize_masks } from '../modules/bedtools/finalize_masks'
include { combine_bedfiles as combine_homref_invariants } from '../modules/bedtools/combine_bedfiles'
include { combine_bedfiles as combine_mappability_masks } from '../modules/bedtools/combine_bedfiles'
include { combine_bedfiles as combine_mappability_masks_snps } from '../modules/bedtools/combine_bedfiles'

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────
def get_filter_expressions(caller_type, custom_snp=null, custom_indel=null, 
                           snp_filter_expression_bcftools, snp_filter_expression_freebayes, snp_filter_expression_gatk,
                           indel_filter_expression_bcftools, indel_filter_expression_freebayes, indel_filter_expression_gatk) {
    // If custom expressions are provided, use those
    if (custom_snp && custom_indel) {
        return [snp_filter_expr: custom_snp, indel_filter_expr: custom_indel]
    }
    
    // Otherwise, return defaults based on variant caller
    if (caller_type in ['gatk_joint', 'gatk_haplotypecaller']) {
        return [snp_filter_expr: snp_filter_expression_gatk, indel_filter_expr: indel_filter_expression_gatk]
    }
    else if (caller_type == 'freebayes') {
        return [snp_filter_expr: snp_filter_expression_freebayes, indel_filter_expr: indel_filter_expression_freebayes]
    }
    else if (caller_type == 'bcftools') {
        return [snp_filter_expr: snp_filter_expression_bcftools, indel_filter_expr: indel_filter_expression_bcftools]
    }
    else {
        error "Unknown variant caller for filter expression: ${caller_type}"
    }
}

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

def sort_and_extract_vcfs(vcf_channel) {
    return vcf_channel
        .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by region_id for deterministic output
        .map { sorted_list ->
            def vcfs = sorted_list.collect { item -> item[1] }
            def idxs = sorted_list.collect { item -> item[2] }
            tuple(vcfs, idxs)
        }
}

workflow VARIANT_FILTERS {
    take:
    raw_vcfs               // tuple: [region_id, vcf, idx]
    ch_reference           // path: reference FASTA
    reference_fai          // path: reference .fai index
    reference_gzi          // path: reference .gzi index
    refintervals_ch        // tuple: [region_id, [intervals]]
    callable_regions       // tuple: [sample_id, bedfile]
    depth_cutoffs          // tuple: [sample_id, min_dp, max_dp, sex_assignment, depths]
    sex_limited_contigs    // channel value: sex-limited chromosomes
    variant_caller         // val: variant caller name
    snp_filter_expression  // val: custom SNP filter expression (or null)
    indel_filter_expression // val: custom indel filter expression (or null)
    snp_filter_expression_bcftools   // val: default bcftools SNP filter
    snp_filter_expression_freebayes  // val: default freebayes SNP filter
    snp_filter_expression_gatk       // val: default GATK SNP filter
    indel_filter_expression_bcftools // val: default bcftools indel filter
    indel_filter_expression_freebayes // val: default freebayes indel filter
    indel_filter_expression_gatk     // val: default GATK indel filter
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Normalize VCF files
    // ─────────────────────────────────────────────────────────────────────────────
    raw_vcfs
        .combine(ch_reference)
        .combine(reference_fai)
        .combine(reference_gzi)
        .map { region_id, vcf, idx, reference, ref_fai, ref_gzi ->
            tuple(region_id, vcf, idx, reference, ref_fai, ref_gzi)
        }
        .set { bcftools_norm_in_ch }
    
    bcftools_norm(bcftools_norm_in_ch)

    // ─────────────────────────────────────────────────────────────────────────────
    // Select SNPs and indels separately
    // ─────────────────────────────────────────────────────────────────────────────
    select_snps(bcftools_norm.out.vcf)
    select_indels(bcftools_norm.out.vcf)

    // ─────────────────────────────────────────────────────────────────────────────
    // Get filter expressions based on variant caller type
    // ─────────────────────────────────────────────────────────────────────────────
    filter_expressions = get_filter_expressions(
        variant_caller,
        snp_filter_expression,
        indel_filter_expression,
        snp_filter_expression_bcftools,
        snp_filter_expression_freebayes,
        snp_filter_expression_gatk,
        indel_filter_expression_bcftools,
        indel_filter_expression_freebayes,
        indel_filter_expression_gatk
    )

    // ─────────────────────────────────────────────────────────────────────────────
    // Apply initial filters
    // ─────────────────────────────────────────────────────────────────────────────
    bcftools_filter_snps(select_snps.out.vcf, filter_expressions.snp_filter_expr, "snps")
    bcftools_filter_indels(select_indels.out.vcf, filter_expressions.indel_filter_expr, "indels")

    refintervals_ch
        .combine(bcftools_filter_snps.out.vcf, by: 0)
        .combine(callable_regions)
    // run through callability filter
    callability_filter_snps(
        refintervals_ch
        .combine(bcftools_filter_snps.out.vcf, by: 0)
        .combine(callable_regions)
        )
    callability_filter_indels(
        refintervals_ch
        .combine(bcftools_filter_indels.out.vcf, by: 0)
        .combine(callable_regions)
        )
    
    // push through allele balance filter
    ab_filter_snps(
        callability_filter_snps.out.vcf
            .map { region_id, sample, vcf, csi -> tuple(region_id, sample, vcf, csi, 'snps') }
        )
    ab_filter_indels(
        callability_filter_indels.out.vcf
            .map { region_id, sample, vcf, csi -> tuple(region_id, sample, vcf, csi, 'indels') }
            )

    // // ─────────────────────────────────────────────────────────────────────────────
    // // Merge samples together per region
    // // ─────────────────────────────────────────────────────────────────────────────

    

    bcftools_merge_snps(ab_filter_snps.out.vcf
        .groupTuple(by: 0)
        .map {region_id, sample_ids, sample_vcfs, idxs -> tuple(region_id, sample_ids, sample_vcfs, idxs, 'snps_abfiltered')}
        // ensure that sample vcf files are always sorted in the same order so that we can concatenate later on
        .map { region_id, sample_ids, sample_vcfs, idxs, category ->
            def zipped = [sample_ids, sample_vcfs, idxs].transpose()
                .sort { a, b -> a[0] <=> b[0] } // sort by sample ID
            def (sample_ids_sorted, sample_vcfs_sorted, idxs_sorted) = zipped.transpose()
            tuple(region_id, sample_ids_sorted, sample_vcfs_sorted, idxs_sorted, category)
        }
    )
    bcftools_merge_indels(ab_filter_indels.out.vcf
        .groupTuple(by: 0)
        .map {region_id, sample_ids, sample_vcfs, idxs -> tuple(region_id, sample_ids, sample_vcfs, idxs, 'indels_abfiltered')}
        // ensure that sample vcf files are always sorted in the same order so that we can concatenate later on
        .map { region_id, sample_ids, sample_vcfs, idxs, category ->
            def zipped = [sample_ids, sample_vcfs, idxs].transpose()
                .sort { a, b -> a[0] <=> b[0] } // sort by sample ID
            def (sample_ids_sorted, sample_vcfs_sorted, idxs_sorted) = zipped.transpose()
            tuple(region_id, sample_ids_sorted, sample_vcfs_sorted, idxs_sorted, category)
        }
    )
    

    // // ─────────────────────────────────────────────────────────────────────────────
    // // Last set of popgen filters based on missingness and allele frequency
    // // ─────────────────────────────────────────────────────────────────────────────
    bcftools_fmiss_maf_filtered_snps = bcftools_filter_fmiss_maf_snps(bcftools_merge_snps.out.vcf, 'snps')
    bcftools_fmiss_maf_filtered_indels = bcftools_filter_fmiss_maf_indels(bcftools_merge_indels.out.vcf, 'indels')

    // // ─────────────────────────────────────────────────────────────────────────────
    // // Generate statistics for filtered variants
    // // ─────────────────────────────────────────────────────────────────────────────
    filtered_snp_stats_out = filtered_snp_stats(bcftools_fmiss_maf_filtered_snps, 'snps')
    filtered_indels_stats_out = filtered_indel_stats(bcftools_fmiss_maf_filtered_indels, 'indels')

    // // Combine statistics
    combined_stats_snps = combine_filtered_snps_stats(
        filtered_snp_stats_out.ab_dp.collect(),
        filtered_snp_stats_out.qual_fmiss_maf.collect(),
        filtered_snp_stats_out.sample_stats.collect(),
        filtered_snp_stats_out.rec_counts.collect(),
        'filtered_snps'
    )

    combined_stats_indels = combine_filtered_indels_stats(
        filtered_indels_stats_out.ab_dp.collect(),
        filtered_indels_stats_out.qual_fmiss_maf.collect(),
        filtered_indels_stats_out.sample_stats.collect(),
        filtered_indels_stats_out.rec_counts.collect(),
        'filtered_indels'
    )

    // // Generate variant statistics PDF reports
    plot_snp_stats(combined_stats_snps.combined_summary_statistics.collect(), 'filtered_snps')
    plot_indel_stats(combined_stats_indels.combined_summary_statistics.collect(), 'filtered_indels')

    // // ─────────────────────────────────────────────────────────────────────────────
    // // Sort and concatenate filtered VCFs
    // // ─────────────────────────────────────────────────────────────────────────────
    sorted_snps = sort_and_extract_vcfs(bcftools_fmiss_maf_filtered_snps.vcf)
    sorted_indels = sort_and_extract_vcfs(bcftools_fmiss_maf_filtered_indels.vcf)

    bcftools_concat_snps(sorted_snps, 'filtered_snps')
    bcftools_concat_indels(sorted_indels, 'filtered_indels')

    // // ─────────────────────────────────────────────────────────────────────────────
    // // Finalize mask files
    // // ─────────────────────────────────────────────────────────────────────────────
    finalize_masks_input = refintervals_ch
        .combine(callable_regions)
        .combine(bcftools_merge_snps.out.vcf, by: 0)
        .combine(bcftools_merge_indels.out.vcf, by: 0)
        .map {
            region_id, region, sample_id, callable_bed, snp_vcf, _snp_idx, indel_vcf, _indel_idx ->
                tuple(sample_id, region_id, region, callable_bed, snp_vcf, indel_vcf)
        }

    region_masks = finalize_masks(finalize_masks_input)
    
    // // Extract and group by sample for final merging
    region_homrefs = region_masks.homref_invariants
        .groupTuple(by: 0)
    region_totalmask = region_masks.mappability_mask
        .groupTuple(by: 0)
    region_snpmask = region_masks.mappability_mask_snps
        .groupTuple(by: 0)

    // // Concatenate regional masks into genome-wide files for output
    homrefs = combine_homref_invariants(region_homrefs, reference_fai, 'homref_invariants')
    total_mask = combine_mappability_masks(region_totalmask, reference_fai, 'mappability_mask')
    snp_mask = combine_mappability_masks_snps(region_snpmask, reference_fai, 'snp_mask')

    emit:
    filtered_snps = bcftools_concat_snps.out.vcf
    filtered_indels = bcftools_concat_indels.out.vcf
    callable_regions_bed = total_mask.bedfile
    snpable_regions_bed = snp_mask.bedfile
    invariant_calls_bed = homrefs.bedfile
    filtered_snps_stats = combined_stats_snps.combined_summary_statistics
    filtered_snps_stats_plot = plot_snp_stats.out.report
    filtered_indel_stats = combined_stats_indels.combined_summary_statistics
    filtered_indel_stats_plot = plot_indel_stats.out.report
}
