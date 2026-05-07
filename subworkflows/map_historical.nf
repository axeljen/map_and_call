#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: MAP_HISTORICAL
//
// Purpose: Map historical/ancient DNA sequencing reads to reference genome
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { bwa_mem as map_historical } from '../modules/bwa/bwa_mem'
include { bwa_mem_singlereads as map_merged } from '../modules/bwa/bwa_mem'
include { merge_historical_bams } from '../modules/samtools/mergebams'

workflow MAP_HISTORICAL {
    take:
    clean_paired           // tuple: [sample_id, library, data_type, reads1, reads2]
    clean_merged           // tuple: [sample_id, library, data_type, collapsed]
    bwa_index              // tuple: [reference, index_files]
    map_historical_pairs   // val: boolean flag to map historical paired-end reads
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Map paired-end reads from historical samples (optional)
    // ─────────────────────────────────────────────────────────────────────────────
    if (map_historical_pairs) {
        historical_pairs_to_map = clean_paired
            .filter { _sample_id, _library, datatype, _reads1, _reads2 -> datatype == '2' }
            .combine(bwa_index)
        historical_pairbams_ch = map_historical(historical_pairs_to_map)
    }
    else {
        // Create an empty channel if paired mapping is disabled
        historical_pairbams_ch = channel.empty()
    }
    
    // ─────────────────────────────────────────────────────────────────────────────
    // Map collapsed/merged reads from historical samples
    // ─────────────────────────────────────────────────────────────────────────────
    historical_merged_to_map = clean_merged
        .filter { _sample_id, _library, datatype, _collapsed -> datatype == '2' }
        .combine(bwa_index)

    historical_bams_ch = map_merged(historical_merged_to_map)

    // ─────────────────────────────────────────────────────────────────────────────
    // If mapping historical paired end reads, merge within each library
    // ─────────────────────────────────────────────────────────────────────────────
    if (map_historical_pairs) {
        merge_historical = historical_pairbams_ch
            .mix(historical_bams_ch)
            .groupTuple(by: [0,1,2])
            .map { sample_id, library, datatype, bam_paths, _bam_indices ->
                return tuple(sample_id, library, datatype, bam_paths)
            }
        historical_bams_merged_ch = merge_historical_bams(merge_historical)
        historical_bams_ch = historical_bams_merged_ch
    }

    emit:
    bam = historical_bams_ch
}
