#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: INDEX_REFERENCE
//
// Purpose: Index reference genome and generate genomic intervals for parallel processing
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { bwa_index } from '../modules/bwa/bwa_index'
include { samtools_index } from '../modules/samtools/index_reference'
include { dochunks } from '../modules/reference_intervals/dochunks'

workflow INDEX_REFERENCE {
    take:
    ch_reference           // path: reference FASTA file
    scaffold_list          // path: optional scaffold list file (or null)
    chunk_size             // val: chunk size for interval generation
    x_scaffolds            // list: X chromosome scaffold names
    y_scaffolds            // list: Y chromosome scaffold names
    z_scaffolds            // list: Z chromosome scaffold names
    w_scaffolds            // list: W chromosome scaffold names
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // Index reference genome with BWA and samtools
    // ─────────────────────────────────────────────────────────────────────────────
    bwa_index_ch = bwa_index(ch_reference)
    faidx_and_chunks_ch = samtools_index(ch_reference)
    
    reference_fasta = faidx_and_chunks_ch.reference_fasta.first()
    reference_fai = faidx_and_chunks_ch.reference_fai.first()
    reference_gzi = faidx_and_chunks_ch.reference_gzi.first()
    
    // ─────────────────────────────────────────────────────────────────────────────
    // Parse scaffold list or extract from reference index
    // ─────────────────────────────────────────────────────────────────────────────
    if (scaffold_list) {
        scaffolds_ch = channel.fromPath(scaffold_list, checkIfExists: true)
            .splitCsv(header: false)
            .map { row -> row[0] }
            .collect()
    }
    else {
        // Extract scaffold names from reference .fai index
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

    // ─────────────────────────────────────────────────────────────────────────────
    // Branch scaffolds into autosomes and sex chromosomes
    // ─────────────────────────────────────────────────────────────────────────────
    scaffolds = scaffolds_ch
        .flatten()
        .branch { scaffold ->
            autosomes: [y_scaffolds, w_scaffolds, x_scaffolds, z_scaffolds].flatten().contains(scaffold) == false
            sex_limited: [y_scaffolds, w_scaffolds].flatten().contains(scaffold)
            non_sex_limited: [x_scaffolds, z_scaffolds].flatten().contains(scaffold)
        }

    // ─────────────────────────────────────────────────────────────────────────────
    // Generate reference intervals for parallel processing
    // ─────────────────────────────────────────────────────────────────────────────
    reference_intervals = dochunks(
        samtools_index.out.reference_fai, 
        chunk_size, 
        scaffolds_ch.flatten()
    )

    // Format intervals into tuples with region IDs
    refintervals_ch = reference_intervals
        .map { row -> row?.trim() }
        .flatMap { row -> row ? row.split(/\r?\n/) as List : [] }
        .filter { row -> row }
        .collect()
        .flatMap { intervals ->
            intervals.withIndex().collect { interval, idx -> tuple(idx + 1, interval) }
        }

    emit:
    bwa_index          = bwa_index_ch.reference
    reference_fasta    = reference_fasta
    reference_fai      = reference_fai
    reference_gzi      = reference_gzi
    refintervals       = refintervals_ch
    scaffolds_autosomes = scaffolds.autosomes
    scaffolds_sex_limited = scaffolds.sex_limited
    scaffolds_non_sex_limited = scaffolds.non_sex_limited
}
