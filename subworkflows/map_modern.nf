#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: MAP_MODERN
//
// Purpose: Map modern DNA sequencing reads to reference genome
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { bwa_mem } from '../modules/bwa/bwa_mem'

workflow MAP_MODERN {
    take:
    clean_paired           // tuple: [sample_id, library, data_type, reads1, reads2]
    bwa_index              // tuple: [reference, index_files]
    mapper                 // val: mapper name (e.g., 'bwa_mem')
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Filter for modern reads only and combine with BWA index
    // ─────────────────────────────────────────────────────────────────────────────
    modern_pairs_to_map = clean_paired
        .filter { _sample_id, _library, datatype, _reads1, _reads2 -> datatype == '1' }
        .combine(bwa_index)
    
    // ─────────────────────────────────────────────────────────────────────────────
    // Map reads using selected mapper
    // ─────────────────────────────────────────────────────────────────────────────
    if (mapper == 'bwa_mem') {
        rawbam_ch = bwa_mem(modern_pairs_to_map)
    }
    else {
        error "Unsupported mapper specified: ${mapper}. Currently only 'bwa_mem' is supported."
    }

    emit:
    bam = rawbam_ch
}
