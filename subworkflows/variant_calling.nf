#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: VARIANT_CALLING
//
// Purpose: Variant calling using bcftools mpileup or freebayes
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { mpileup } from '../modules/bcftools/bcftools_mpileup'
include { freebayes } from '../modules/freebayes/freebayes'
include { bcftools_merge } from '../modules/bcftools/bcftools_merge'
include { bcftools_concat as bcftools_concat_raw } from '../modules/bcftools/bcftools_concat'
include { vcf_stats as raw_vcf_stats } from '../modules/variant_stats/vcf_stats_process'
include { combine_stats as combine_raw_vcf_stats } from '../modules/variant_stats/combine_stats'
include { plot_variant_stats as plot_raw_stats } from '../modules/variant_stats/plot_variant_stats'

workflow VARIANT_CALLING {
    take:
    bams_for_calling       // tuple: [sample_id, cram, crai]
    ch_reference           // path: reference FASTA
    reference_fai          // path: reference .fai index
    reference_gzi          // path: reference .gzi index
    bwa_index              // tuple: [reference, index_files]
    refintervals_ch        // tuple: [region_id, interval]
    variant_caller         // val: variant caller name ('bcftools' or 'freebayes')
    popfile                // path: population assignment file (optional, or null)
    sample_stats           // tuple: [sample_id, autosomal_dp, ...]
    store_raw_vcf          // val: boolean flag to store raw VCF
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Validate variant caller parameter
    // ─────────────────────────────────────────────────────────────────────────────
    if (!['freebayes','bcftools'].contains(variant_caller)) {
        error "Unsupported variant caller specified: ${variant_caller}. Must be one of freebayes or bcftools."
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Parse population assignments for joint genotyping
    // ─────────────────────────────────────────────────────────────────────────────
    if (popfile) {
        pops = channel.fromPath(popfile, checkIfExists: true)
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

    // ─────────────────────────────────────────────────────────────────────────────
    // Prepare BAM channels based on variant caller
    // ─────────────────────────────────────────────────────────────────────────────


    // ─────────────────────────────────────────────────────────────────────────────
    // Prepare reference bundle with all index files
    // ─────────────────────────────────────────────────────────────────────────────
    ref_index_ch = bwa_index
        .map { _reference, index_files -> index_files }
    ref_bundle_ch = ch_reference
        .combine(reference_fai)
        .combine(reference_gzi)
        .combine(ref_index_ch)
        .map { row ->
            def reference = row[0]
            def fai = row[1]
            def gzi = row[2]
            def bwa_indices = row[3..-1]
            tuple(reference, [fai, gzi] + bwa_indices)
        }

    varcall_ch = bams_for_calling
        .map { _sample_id, cram, crai -> tuple(cram, crai) }
        .combine(refintervals_ch)
        .groupTuple(by: [2,3])
         // Broadcast reference bundle to all variant calling tasks
        .combine(ref_bundle_ch)
        .map { bams, bais, region_id, regions, ref_genome, ref_indices ->
            tuple(bams, bais, ref_genome, ref_indices, region_id, regions)
        }
        .view { it -> println "Variant calling input before combining reference bundle: ${it}" }


    // ─────────────────────────────────────────────────────────────────────────────
    // Call variants using selected variant caller
    // ─────────────────────────────────────────────────────────────────────────────

    // Collect population assignments for freebayes
    popstring = pops
        .map { sample, pop -> "${sample}=${pop}" }
        .collect().toList()
    varcall_ch = varcall_ch
        .combine(popstring)
        //.view { it -> println "Variant calling input: ${it}" }

    if (variant_caller == 'bcftools') {
        mpileup(varcall_ch)
        raw_vcfs = mpileup.out.vcf
    }
    else if (variant_caller == 'freebayes') {
        freebayes(varcall_ch)
        raw_vcfs = freebayes.out.vcf
    }

//     // ─────────────────────────────────────────────────────────────────────────────
//     // Merge VCFs across populations if multiple populations exist
//     // ─────────────────────────────────────────────────────────────────────────────
//     npops = bams_for_calling_grouped.groupTuple(by: 0)
//         .map { population, _crams, _crais -> population }
//         .collect()
//         .map { populations -> populations.size() }
//         .branch { n ->
//              multiple_pops: n > 1
//              single_pop: n <= 1
//         }

//     multipop_merge_in = npops.multiple_pops.combine(pop_vcfs_ch.groupTuple(by: 0)
//             // Make sure that populations are always sorted in a consistent order
//             .map { region_id, populations, vcfs, idxs ->
//                 def zipped = [populations, vcfs, idxs].transpose()
//                     .sort { a, b -> a[0] <=> b[0] } // sort by population name
//                 def (populations_sorted, vcfs_sorted, idxs_sorted) = zipped.transpose()
//                 tuple(region_id, populations_sorted, vcfs_sorted, idxs_sorted)
//             }
//             .map { region_id, pop_list, vcf_list, idx_list ->
//                 tuple(region_id, pop_list, vcf_list, idx_list)
//             }
//         )
//         .map { _npops, region_id, pop_list, vcf_list, idx_list ->
//             tuple(region_id, pop_list, vcf_list, idx_list)
//         }
//         // Broadcast singleton reference files across all regions
//         .combine(reference_fai)
//         .combine(reference_gzi)
//         .combine(ch_reference)
//         .map { region_id, pop_list, vcf_list, idx_list, ref_fai, ref_gzi, ref_fa ->
//             tuple(region_id, pop_list, vcf_list, idx_list, ref_fa, ref_fai, ref_gzi)
//         }

// //    multipop_vcf = bcftools_merge(multipop_merge_in, 'raw')

//     singlepop_vcf = pop_vcfs_ch.combine(npops.single_pop)
//         .map { region_id, _pop, vcf, idx, _npops ->
//             tuple(region_id, vcf, idx)
//         }

//     raw_vcfs = singlepop_vcf

    // ─────────────────────────────────────────────────────────────────────────────
    // Extract summary statistics from raw VCF files
    // ─────────────────────────────────────────────────────────────────────────────
    raw_stats = raw_vcf_stats(raw_vcfs, 'raw_variants')
    
    // Combine statistics into one file
    combined_raw_stats = combine_raw_vcf_stats(
        raw_stats.ab_dp.collect(),
        raw_stats.qual_fmiss_maf.collect(),
        raw_stats.sample_stats.collect(),
        raw_stats.rec_counts.collect(),
        'raw_variants'
    )

    // Generate variant statistics PDF report for raw variants
    plot_raw_stats(combined_raw_stats.combined_summary_statistics.collect(), 'raw_variants')

    // ─────────────────────────────────────────────────────────────────────────────
    // Concatenate raw VCFs if requested
    // ─────────────────────────────────────────────────────────────────────────────
    if (store_raw_vcf == true) {
        sorted_raw_vcfs = raw_vcfs
            .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by region_id
            .map { sorted_list ->
                // Extract VCFs and indices in sorted order
                def vcfs = sorted_list.collect { item -> item[1] }
                def idxs = sorted_list.collect { item -> item[2] }
                tuple(vcfs, idxs)
            }
        // Concatenate raw VCF files
        bcftools_concat_raw(sorted_raw_vcfs, 'raw_variants')
        raw_vcf = bcftools_concat_raw.out.vcf
    } else {
        raw_vcf = channel.empty()
    }

    emit:
    raw_vcfs = raw_vcfs
    pops = pops
    raw_vcf = raw_vcf
    raw_variant_stats = combined_raw_stats.combined_summary_statistics
    raw_vcf_stats_plot = plot_raw_stats.out.report
}
