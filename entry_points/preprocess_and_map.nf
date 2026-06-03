#!/usr/bin/env nextflow

// ═══════════════════════════════════════════════════════════════════════════════
//      PREPROCESS AND MAP WORKFLOW - Read preprocessing, mapping, and BAM processing
//
// Purpose: Run preprocessing, mapping, and BAM processing (skip variant calling)
// Use case: Generate BAM/CRAM files for downstream analysis, sex assignment, depth estimation
// ═══════════════════════════════════════════════════════════════════════════════

// ─────────────────────────────────────────────────────────────────────────────
// Subworkflow imports
// ─────────────────────────────────────────────────────────────────────────────
include { INDEX_REFERENCE } from '../subworkflows/index_reference'
include { PREPROCESS_MODERN } from '../subworkflows/preprocess_modern'
include { PREPROCESS_HISTORICAL } from '../subworkflows/preprocess_historical'
include { MAP_MODERN } from '../subworkflows/map_modern'
include { MAP_HISTORICAL } from '../subworkflows/map_historical'
include { PROCESS_BAMS } from '../subworkflows/process_bams'

// ─────────────────────────────────────────────────────────────────────────────
// Utility imports
// ─────────────────────────────────────────────────────────────────────────────
// import WorkflowUtils

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/*
 * Parse input metadata CSV
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
            
            def read_1_path = params.reads_dir ? "${params.reads_dir}/${row.read_1}" : row.read_1
            def read_2_path = params.reads_dir ? "${params.reads_dir}/${row.read_2}" : row.read_2
            
            def read_1 = file(read_1_path, checkIfExists: true)
            def read_2 = file(read_2_path, checkIfExists: true)
            return [sample_id, data_type, library, read_1, read_2]
        }
        .groupTuple(by: 0)
        .flatMap { sample_id, data_type, library, reads_1, reads_2 ->
            if (data_type.size() == 1) {
                return [[sample_id, 'single_lane', data_type[0], library[0], reads_1[0], reads_2[0]]]
            } else {
                return (0..<data_type.size()).collect { idx ->
                    [sample_id, "${idx + 1}", data_type[idx], library[idx], reads_1[idx], reads_2[idx]]
                }
            }
        }
}

/*
 * Setup sex chromosome system and identify sex-linked contigs
 * Returns map with: sex_chrom_system, sex_linked_list, sex_limited_list, non_sex_limited_list
 */
def setup_sex_chromosome_system() {
    def result = [
        sex_chrom_system: 'unknown',
        sex_linked_list: [],
        sex_limited_list: [],
        non_sex_limited_list: []
    ]
    
    if (params.x_scaffolds || params.y_scaffolds) {
        if (params.z_scaffolds || params.w_scaffolds) {
            error "Please only specify either X and/or Y, OR Z and/or W scaffolds, not both."
        }
        result.sex_chrom_system = 'XY'
        result.sex_linked_list = [params.x_scaffolds, params.y_scaffolds].flatten()
        result.sex_limited_list = [params.y_scaffolds].flatten()
        result.non_sex_limited_list = [params.x_scaffolds].flatten()
    }
    else if (params.z_scaffolds || params.w_scaffolds) {
        result.sex_chrom_system = 'ZW'
        result.sex_linked_list = [params.z_scaffolds, params.w_scaffolds].flatten()
        result.sex_limited_list = [params.w_scaffolds].flatten()
        result.non_sex_limited_list = [params.z_scaffolds].flatten()
    }
    
    return result
}

// ═══════════════════════════════════════════════════════════════════════════════
//                           MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════════

workflow {
    main:
    
    println "Running PREPROCESS_AND_MAP workflow: preprocessing, mapping, and BAM processing (no variant calling)"
    
    // Setup sex chromosome system
    sex_config = setup_sex_chromosome_system()
    def sex_chrom_system = sex_config.sex_chrom_system
    def sex_linked_list = sex_config.sex_linked_list
    def sex_limited_list = sex_config.sex_limited_list
    def non_sex_limited_list = sex_config.non_sex_limited_list
    
    println "Sex-limited contigs: ${sex_limited_list}"
    println "Non-sex-limited contigs: ${non_sex_limited_list}"
    println "Inferred sex chromosome system: ${sex_chrom_system}"
    
    // Create channels for sex-linked contigs
    sex_linked_contigs = channel.value(sex_linked_list)
    sex_limited_contigs = channel.value(sex_limited_list)
    non_sex_limited_contigs = channel.value(non_sex_limited_list)
    
    // Input channel from metadata file
    ch_input = parse_input(params.input)
        .branch {
            sample_id, lane, data_type, library, r1, r2 ->
            modern: data_type == '1'
            historical: data_type == '2'
        }
    
    // Reference genome (value channel - can be consumed multiple times)
    println "Reference genome: ${params.reference}"
    ch_reference = channel.fromPath(params.reference, checkIfExists: true)
    ch_reference = ch_reference.view { "Reference genome (channel): $it" }
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     INDEX REFERENCE
    // ═══════════════════════════════════════════════════════════════════════════════
    
    INDEX_REFERENCE(
        ch_reference,
        params.scaffold_list,
        params.chunk_size,
        params.x_scaffolds ?: [],
        params.y_scaffolds ?: [],
        params.z_scaffolds ?: [],
        params.w_scaffolds ?: []
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PREPROCESS READS
    // ═══════════════════════════════════════════════════════════════════════════════
    
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
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     MAP READS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    MAP_MODERN(
        PREPROCESS_MODERN.out.clean_reads,
        INDEX_REFERENCE.out.bwa_index,
        params.mapper
    )

    MAP_HISTORICAL(
        PREPROCESS_HISTORICAL.out.clean_paired,
        PREPROCESS_HISTORICAL.out.clean_merged,
        INDEX_REFERENCE.out.bwa_index,
        params.map_historical_pairs,
        params.historical_mapper
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PROCESS BAMS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    PROCESS_BAMS(
        MAP_MODERN.out.bam,
        MAP_HISTORICAL.out.bam,
        ch_reference,
        INDEX_REFERENCE.out.reference_fai,
        INDEX_REFERENCE.out.reference_gzi,
        INDEX_REFERENCE.out.refintervals,
        sex_limited_list,
        non_sex_limited_list,
        sex_limited_contigs,
        non_sex_limited_contigs,
        params.postmapping_dedup,
        params.damageprofiler_rescale,
        params.downsample_bams,
        params.downsample_bams_coverage,
        params.min_depth,
        params.max_depth,
        params.sex_assignment_lower_threshold,
        params.sex_assignment_upper_threshold
    )
    
    // ═══════════════════════════════════════════════════════════════════════════════
    //                     PUBLISH OUTPUTS
    // ═══════════════════════════════════════════════════════════════════════════════
    
    publish:
    fastqc_raw = PREPROCESS_MODERN.out.fastqc_raw
        .mix(PREPROCESS_HISTORICAL.out.fastqc_rawreads)
    fastqc_clean = PREPROCESS_MODERN.out.fastqc_clean
        .mix(PREPROCESS_HISTORICAL.out.fastqc_cleanreads)
    multiqc_rawreads_report = multiqc_rawreads_report
    multiqc_cleanreads_report = multiqc_cleanreads_report
    clean_reads = clean_reads
    reference_genome = INDEX_REFERENCE.out.bwa_index
    qualimap_reports = PROCESS_BAMS.out.qualimap_reports
    qualimap_downsampled_reports = params.downsample_bams ? PROCESS_BAMS.out.qualimap_downsampled_reports : channel.empty()
    bamfiles = PROCESS_BAMS.out.raw_bams
    bam_metrics = PROCESS_BAMS.out.cram_metrics
    cramfiles = params.store_crams ? PROCESS_BAMS.out.raw_crams : channel.empty()
    mapping_depths = PROCESS_BAMS.out.mapping_depths
    downsampled_cramfiles = (params.downsample_bams && params.store_crams) ? PROCESS_BAMS.out.downsampled_crams : channel.empty()
    downsampled_bamfiles = (params.downsample_bams && !params.store_crams) ? PROCESS_BAMS.out.final_bam : channel.empty()
    damage_profiles = PROCESS_BAMS.out.damage_reports
}

// ═══════════════════════════════════════════════════════════════════════════════
//                             OUTPUT DEFINITIONS
// ═══════════════════════════════════════════════════════════════════════════════

output {
    reference_genome {
        path "00_input_data/00_reference_genome"
    }
    clean_reads {
        enabled params.store_cleanreads
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
    qualimap_reports {
        path "01_reports/01_qualimap/01_qualimap_post_markdups"
    }
    qualimap_downsampled_reports {
        enabled params.downsample_bams
        path "01_reports/01_qualimap/02_qualimap_downsampled"
    }
    bamfiles {
        enabled !params.store_crams
        path "02_bamfiles"
    }
    mapping_depths {
        path "02_bamfiles/mapping_depths"
    }
    downsampled_bamfiles {
        enabled params.downsample_bams && !params.store_crams
        path "02_bamfiles/downsampled"
    }
    bam_metrics {
        enabled params.postmapping_dedup
        path "02_bamfiles/dedup_metrics"
    }
    cramfiles {
        enabled params.store_crams
        path "02_bamfiles"
    }
    downsampled_cramfiles {
        enabled params.downsample_bams && params.store_crams
        path "02_bamfiles/downsampled"
    }
    damage_profiles {
        path "01_reports/03_damage_profiles"
    }
}
