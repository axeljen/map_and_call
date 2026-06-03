#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//             PREPROCESS READS WORKFLOW - Preprocessing Only
//
// Purpose: Run preprocessing (adapter trimming, quality filtering) and store clean reads
// Use case: Generate clean reads for downstream analysis with external tools
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Subworkflow imports
// ─────────────────────────────────────────────────────────────────────────────
include { PREPROCESS_MODERN } from './subworkflows/preprocess_modern'
include { PREPROCESS_HISTORICAL } from './subworkflows/preprocess_historical'

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/*
 * Parse input metadata CSV
 * Expected format: sample_id;data_type;library;read_1;read_2
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

// ═══════════════════════════════════════════════════════════════════════════════
//                           MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════════

workflow {
    main:
    
    println "Running PREPROCESS_READS workflow: preprocessing only, outputting clean reads to ${params.output_dir}/00_input_data/01_clean_reads"

    // Input channel from metadata file
    ch_input = parse_input(params.input)
        .branch {
            sample_id, lane, data_type, library, r1, r2 ->
            modern: data_type == '1'
            historical: data_type == '2'
        }
    
    // Run preprocessing subworkflows
    PREPROCESS_MODERN(
        ch_input.modern,
        params.premapping_dedup
    )

    PREPROCESS_HISTORICAL(
        ch_input.historical,
        params.premapping_dedup
    )

    // Combine multiqc reports
    multiqc_rawreads_report = PREPROCESS_MODERN.out.multiqc_raw_report
        .mix(PREPROCESS_HISTORICAL.out.multiqc_raw_report)
        .collect()
    
    multiqc_cleanreads_report = PREPROCESS_MODERN.out.multiqc_clean_report
        .mix(PREPROCESS_HISTORICAL.out.multiqc_clean_report)
        .collect()

    clean_reads = PREPROCESS_MODERN.out.clean_reads
        .mix(PREPROCESS_HISTORICAL.out.clean_paired)
        .mix(PREPROCESS_HISTORICAL.out.clean_merged)
    
    publish:
    fastqc_raw = PREPROCESS_MODERN.out.fastqc_raw
        .mix(PREPROCESS_HISTORICAL.out.fastqc_rawreads)
    fastqc_clean = PREPROCESS_MODERN.out.fastqc_clean
        .mix(PREPROCESS_HISTORICAL.out.fastqc_cleanreads)
    multiqc_rawreads_report = multiqc_rawreads_report
    multiqc_cleanreads_report = multiqc_cleanreads_report
    clean_reads = clean_reads
}

// ═══════════════════════════════════════════════════════════════════════════════
//                             OUTPUT DEFINITIONS
// ═══════════════════════════════════════════════════════════════════════════════

output {
    clean_reads {
        path "00_input_data/01_clean_reads"
    }
    fastqc_raw {
        enabled params.store_sample_fastqc
        path "01_reports/00_fastqc/00_raw_reads/sample_fastqc_reports"
    }
    fastqc_clean {
        enabled params.store_sample_fastqc
        path "01_reports/00_fastqc/01_clean_reads/sample_fastqc_reports"
    }
    multiqc_rawreads_report {
        path "01_reports/00_fastqc/00_raw_reads"
    }
    multiqc_cleanreads_report {
        path "01_reports/00_fastqc/01_clean_reads"
    }
}
