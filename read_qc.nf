#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//                  READ QC WORKFLOW - FastQC Only
//
// Purpose: Run FastQC on raw reads and optionally on clean reads
// Use case: Quick quality assessment before running the full pipeline
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Module imports
// ─────────────────────────────────────────────────────────────────────────────
include { fastqc } from './modules/fastqc/fastqc_process'
include { multiqc_fastqc } from './modules/multiqc/multiqc_fastqc'

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
    
    // Input channel from metadata file if an input file is given:
    if (params.input) {
    ch_input = parse_input(params.input)
        .map { sample_id, lane, _data_type, library, r1, r2 ->
            tuple(sample_id, lane, library, r1, r2, 'raw')
        }
    } else if (!params.reads_dir){
        println "No input file or reads directory provided. Please provide either --input <metadata.csv> or --reads_dir <path_to_reads>."
        System.exit(1)
    } else {
        println "No input file provided. Scanning reads directory for FASTQ files..."
        ch_input = channel
            .fromPath("${params.reads_dir}/*.{fq,fastq}{,.gz}", checkIfExists: true)
            .map {filename ->
                // just dummmy sample, lib, datatype, lane names using the basename of the file without extension for all fields. This won't matter downstream,
                def sample = filename.baseName.replaceFirst(/\.f(ast)?q$/, '')
                def lane = 'single_lane'
                def library = sample
                def datatype = 'unknown'
                tuple(sample, lane, library, filename, null, 'raw')
            }
    } 
    
    // Run FastQC on raw reads
    fastqc(ch_input)
    
    // Collect all FastQC outputs
    all_fastqc_outputs = fastqc.out.zip
        .map { _sample_id, _lane, _library, zips -> zips }
        .collect()
        .map { zips -> tuple(zips, params.name) }
    
    // Run MultiQC on all FastQC outputs
    multiqc_fastqc(all_fastqc_outputs)
    
    publish:
    fastqc_reports = fastqc.out.html
    multiqc_report = multiqc_fastqc.out.report
}

// ═══════════════════════════════════════════════════════════════════════════════
//                             OUTPUT DEFINITIONS
// ═══════════════════════════════════════════════════════════════════════════════

output {
    fastqc_reports {
        path "01_reports/00_fastqc"
    }
    multiqc_report {
        path "01_reports/00_fastqc"
    }
}
