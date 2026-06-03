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

// ─────────────────────────────────────────────────────────────────────────────
// Helper functions
// ─────────────────────────────────────────────────────────────────────────────

/*
 * Check if BWA index files already exist for a reference genome
 */
def check_bwa_index_exists(reference_file) {
    def base = reference_file.toString()
    def required_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    def all_exist = required_extensions.every { ext -> 
        file("${base}${ext}").exists() 
    }
    if (all_exist) {
        println "INFO: BWA index files found for ${reference_file.name}, skipping indexing"
    }
    return all_exist
}

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
    // Check if BWA index files already exist and branch accordingly
    // ─────────────────────────────────────────────────────────────────────────────
    
    // Duplicate the reference channel using multiMap
    ch_ref_split = ch_reference
        .multiMap { ref ->
            for_bwa: ref
            for_samtools: ref
        }
    
    // Branch references based on whether BWA index exists
    ref_branches = ch_ref_split.for_bwa
        .map { ref ->
            def is_indexed = check_bwa_index_exists(ref)
            tuple(ref, is_indexed)
        }
        .branch {
            indexed: it[1] == true
            needs_indexing: it[1] == false
        }
    
    // ─────────────────────────────────────────────────────────────────────────────
    // Index reference genome with BWA (only if needed)
    // ─────────────────────────────────────────────────────────────────────────────
    // For already indexed references, create a channel matching bwa_index output format
    already_indexed = ref_branches.indexed
        .map { ref, is_indexed ->
            def index_files = [
                file("${ref}.amb"),
                file("${ref}.ann"),
                file("${ref}.bwt"),
                file("${ref}.pac"),
                file("${ref}.sa")
            ]
            tuple(ref, index_files)
        }
    
    // Run BWA index only on references that need it
    needs_indexing_refs = ref_branches.needs_indexing.map { ref, is_indexed -> ref }
    newly_indexed = bwa_index(needs_indexing_refs)
    
    // Combine both channels - one or the other will be empty
    bwa_index_ch = already_indexed.mix(newly_indexed.reference)
    
    // Index with samtools
    faidx_and_chunks_ch = samtools_index(ch_ref_split.for_samtools)
    //faidx_and_chunks_ch.reference_intervals.view()
    
    reference_fasta = faidx_and_chunks_ch.reference_fasta
    reference_fai = faidx_and_chunks_ch.reference_fai
    reference_gzi = faidx_and_chunks_ch.reference_gzi
    
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
        scaffolds_ch.flatten().toList()
    )

    // Format intervals into tuples with region IDs
    refintervals_ch = reference_intervals
        .flatMap { text -> text
                .trim()
                .split('\n')
                .toList()
                .withIndex()
                .collect { interval, i ->
                    tuple(i + 1, interval.split(',').toList())
                }
        }

    emit:
    bwa_index          = bwa_index_ch
    reference_fasta    = reference_fasta
    reference_fai      = reference_fai
    reference_gzi      = reference_gzi
    refintervals       = refintervals_ch
    scaffolds_autosomes = scaffolds.autosomes
    scaffolds_sex_limited = scaffolds.sex_limited
    scaffolds_non_sex_limited = scaffolds.non_sex_limited
}
