#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                   SUBWORKFLOW: PREPROCESS_HISTORICAL
//
// Purpose: Quality control and preprocessing of historical/ancient DNA sequencing reads
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { fastqc as fastqc_rawreads } from '../modules/fastqc/fastqc_process'
include { fastqc as fastqc_cleanreads } from '../modules/fastqc/fastqc_process'
include { multiqc_fastqc as multiqc_rawreads } from '../modules/multiqc/multiqc_fastqc'
include { multiqc_fastqc as multiqc_cleanreads} from '../modules/multiqc/multiqc_fastqc'
include { adapterremoval } from '../modules/adapterremoval/adapterremoval'
include { clumpify_single } from '../modules/bbmap/bbmap'
include { clumpify_paired } from '../modules/bbmap/bbmap'
include { concat_reads } from '../modules/concat_libs/concat_libs'
include { concat_collapsed } from '../modules/concat_libs/concat_libs'

workflow PREPROCESS_HISTORICAL {
    take:
    ch_input_historical    // tuple: [sample_id, lane, data_type, library, read_1, read_2]
    premapping_dedup       // val: boolean flag for pre-mapping deduplication
    
    main:
    // ─────────────────────────────────────────────────────────────────────────────
    // FastQC on raw (untrimmed) historical reads
    // ─────────────────────────────────────────────────────────────────────────────
    ch_raw_reads_for_qc = ch_input_historical
        .map { sample_id, lane, datatype, library, r1, r2 ->
            [sample_id, lane, library, [r1, r2], 'raw_reads']
        }
    
    fastqc_rawreads(ch_raw_reads_for_qc)
    
    // Collect all FastQC outputs and run MultiQC
    ch_multiqc_input = fastqc_rawreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'raw_reads_historical') }

    multiqc_rawreads(ch_multiqc_input)

    // ─────────────────────────────────────────────────────────────────────────────
    // Trim historical reads with AdapterRemoval (handles ancient DNA damage)
    // ─────────────────────────────────────────────────────────────────────────────
    adapterremoval(ch_input_historical)

    // ─────────────────────────────────────────────────────────────────────────────
    // Concatenate paired reads from the same library
    // ─────────────────────────────────────────────────────────────────────────────
    concat_pairs_in = adapterremoval.out.trimmed_pairs
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, reads1, reads2 ->
            tuple(sample_id, library, datatypes[0], reads1, reads2)
        }
        .map { sample_id, library, datatype, reads1, reads2 ->
            // Sort reads in the same order to ensure correct pairing during concatenation
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
    // Concatenate merged/collapsed reads from the same library
    // ─────────────────────────────────────────────────────────────────────────────
    concat_merged_in = adapterremoval.out.trimmed_collapsed
        .groupTuple(by: [0,3])
        .map {
            sample_id, _lanes, datatypes, library, collapsed ->
            tuple(sample_id, library, datatypes[0], collapsed)
        }
        // Sort collapsed reads to ensure consistent order for concatenation
        .map { sample_id, library, datatype, collapsed ->
            def sorted_collapsed = collapsed.sort { a, b -> a.name <=> b.name }
            tuple(sample_id, library, datatype, sorted_collapsed)
        }
        .branch { sample_id, library, datatype, collapsed ->
            single_lib: collapsed.size() == 1
            multi_lib: collapsed.size() > 1
        }
    
    concat_merged = concat_collapsed(concat_merged_in.multi_lib)

    // ─────────────────────────────────────────────────────────────────────────────
    // Pre-mapping deduplication (optional)
    // ─────────────────────────────────────────────────────────────────────────────
    if (premapping_dedup) {
        // Combine with samples that didn't need concatenation
        dedup_pairs_in = concat_pairs.reads_concat
            .mix(concat_pairs_in.single_lib)
        dedup_merged_in = concat_merged.collapsed_concat
            .mix(concat_merged_in.single_lib)

        // Push the concatenated libraries through clumpify
        clean_paired = clumpify_paired(dedup_pairs_in)
        clean_merged = clumpify_single(dedup_merged_in)
    } else {
        // If no pre-mapping deduplication, just use the concatenated libraries
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

    // ─────────────────────────────────────────────────────────────────────────────
    // FastQC on clean (trimmed) reads
    // ─────────────────────────────────────────────────────────────────────────────
    clean_merged
        .map { sample_id, library, _datatype, collapsed -> 
            tuple(sample_id, 'combined', library, [collapsed], 'clean')
        }
        .mix(clean_paired
            .map { sample_id, library, _datatype, reads1, reads2 -> 
                tuple(sample_id, 'combined', library, [reads1, reads2], 'clean')
            }
        )
    fastqc_cleanreads(clean_merged
        .map { sample_id, library, _datatype, collapsed -> 
            tuple(sample_id, 'combined', library, [collapsed], 'clean')
        }
        .mix(clean_paired
            .map { sample_id, library, _datatype, reads1, reads2 -> 
                tuple(sample_id, 'combined', library, [reads1, reads2], 'clean')
            }
        )
    )

    multiqc_cleanreads(fastqc_cleanreads.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, 'clean_reads_historical') }
    )

    emit:
    clean_paired = clean_paired
    clean_merged = clean_merged
    fastqc_rawreads = fastqc_rawreads.out.html
    fastqc_cleanreads = fastqc_cleanreads.out.html
    multiqc_raw_report = multiqc_rawreads.out.report
    multiqc_clean_report = multiqc_cleanreads.out.report
}
