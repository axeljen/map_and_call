#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: PREPROCESS_MODERN
//
// Purpose: Quality control and preprocessing of modern DNA sequencing reads
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { fastqc as fastqc_rawreads } from '../modules/fastqc/fastqc_process'
include { fastqc as fastqc_cleanreads } from '../modules/fastqc/fastqc_process'
include { multiqc_fastqc as multiqc_rawreads } from '../modules/multiqc/multiqc_fastqc'
include { multiqc_fastqc as multiqc_cleanreads} from '../modules/multiqc/multiqc_fastqc'
include { fastp } from '../modules/fastp/trimming'
include { clumpify_paired } from '../modules/bbmap/bbmap'
include { concat_reads } from '../modules/concat_libs/concat_libs'

workflow PREPROCESS_MODERN {
    take:
    ch_input_modern        // tuple: [sample_id, lane, data_type, library, read_1, read_2]
    premapping_dedup       // val: boolean flag for pre-mapping deduplication
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // FastQC on raw (untrimmed) modern reads
    // ─────────────────────────────────────────────────────────────────────────────
    ch_raw_reads_for_qc = ch_input_modern
        .map { sample_id, lane, datatype, library, r1, r2 ->
            [sample_id, lane, library, [r1, r2], 'raw_reads']
        }
    
    fastqc_rawreads(ch_raw_reads_for_qc)
    
    // Collect all FastQC outputs and run MultiQC
    ch_multiqc_input = fastqc_rawreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'raw_reads_modern') }

    multiqc_rawreads(ch_multiqc_input)

    // ─────────────────────────────────────────────────────────────────────────────
    // Trim modern reads with fastp (quality and adapter trimming)
    // ─────────────────────────────────────────────────────────────────────────────
    fastp(ch_input_modern)

    // ─────────────────────────────────────────────────────────────────────────────
    // Concatenate reads from the same library
    // ─────────────────────────────────────────────────────────────────────────────
    concat_pairs_in = fastp.out.reads
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, reads1, reads2 ->
            tuple(sample_id, library, datatypes[0], reads1, reads2)
        }
        // Make sure that reads1 and reads2 are sorted in the same order
        .map { sample_id, library, datatype, reads1, reads2 ->
            def zipped = [reads1, reads2].transpose()
                .sort { a, b -> a[0].name <=> b[0].name }
            def (reads1_sorted, reads2_sorted) = zipped.transpose()
            tuple(sample_id, library, datatype, reads1_sorted, reads2_sorted)
        }
        .branch { sample_id, library, datatype, reads1, reads2 ->
            single_lib: reads1.size() == 1
            multi_lib: reads1.size() > 1
        }

    concat_pairs = concat_reads(concat_pairs_in.multi_lib)

    // ─────────────────────────────────────────────────────────────────────────────
    // Pre-mapping deduplication (optional)
    // ─────────────────────────────────────────────────────────────────────────────
    if (premapping_dedup) {
        // Combine with samples that didn't need concatenation
        dedup_pairs_in = concat_pairs.reads_concat
            .mix(concat_pairs_in.single_lib)

        // Push the concatenated libraries through clumpify
        clean_paired = clumpify_paired(dedup_pairs_in)
    } else {
        // If no pre-mapping deduplication, just use the concatenated libraries
        // single_lib items come from groupTuple and carry single-element lists — unwrap to plain paths
        clean_paired = concat_pairs.reads_concat
            .mix(concat_pairs_in.single_lib
                .map { sample_id, library, datatype, reads1, reads2 ->
                    tuple(sample_id, library, datatype, reads1[0], reads2[0])
                })
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // FastQC on clean (trimmed) reads
    // ─────────────────────────────────────────────────────────────────────────────
    fastqc_cleanreads(clean_paired
        .map { sample_id, library, _datatype, reads1, reads2 -> 
            tuple(sample_id, 'combined', library, [reads1, reads2], 'clean')
        }
    )

    multiqc_cleanreads(fastqc_cleanreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'clean_reads_modern') }
    )

    emit:
    clean_reads = clean_paired
    fastqc_raw = fastqc_rawreads.out.html
    fastqc_clean = fastqc_cleanreads.out.html
    multiqc_raw_report = multiqc_rawreads.out.report
    multiqc_clean_report = multiqc_cleanreads.out.report
}
