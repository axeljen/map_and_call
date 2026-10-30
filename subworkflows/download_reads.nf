#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: DOWNLOAD_READS
//
// Purpose: Download DNA sequencing reads from SRA
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { sratools_fastq_dump } from '../modules/sratools/sratools_fastq_dump'

workflow DOWNLOAD_READS {
    take:
    sra_ids           // tuple: [sample_id, lane, data_type, library, sra_id]
    
    main:
    
    // download reads 
    reads = sratools_fastq_dump(sra_ids)


    emit:
    fastq_files = reads.reads
}
