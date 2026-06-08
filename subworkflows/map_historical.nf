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
include { bwa_mem_singlereads as map_long_historical } from '../modules/bwa/bwa_mem'
include { bwa_mem_singlereads as map_merged } from '../modules/bwa/bwa_mem'
include { merge_historical_bams } from '../modules/samtools/mergebams'
include { merge_aln_mem as merge_short_and_long_read_bams } from '../modules/samtools/mergebams'
include { split_fq_by_length } from '../modules/seqkit/seqkit'
include { bwa_aln } from '../modules/bwa/bwa_aln'

workflow MAP_HISTORICAL {
    take:
    clean_paired           // tuple: [sample_id, library, data_type, reads1, reads2]
    clean_merged           // tuple: [sample_id, library, data_type, collapsed]
    bwa_index              // tuple: [reference, index_files]
    map_historical_pairs   // val: boolean flag to map historical paired-end reads
    historical_mapper     // val: either "split", "bwa_mem" or "bwa_aln"; if "split", historical reads will be split into those with lengths < 70 bp and >= 70 bp and mapped with bwa aln and bwa mem, respectively
    
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
    // if historical_mapper is set to "split", split reads into those shorter and longer than the short_reads_threshold and map with bwa aln and bwa mem, respectively; otherwise, map all reads with the specified mapper
    if (historical_mapper == "split") {
        historical_splitreads_ch = clean_merged
            .filter { _sample_id, _library, datatype, _collapsed -> datatype == '2' }
            .map { sample_id, library, datatype, collapsed ->
                return tuple(sample_id, library, datatype, collapsed)
            }
        split_fq_by_length(historical_splitreads_ch)
        // map shortreads with bwa aln
        bwa_aln(split_fq_by_length.out.shortreads
            .combine(bwa_index)
            )
        // longreads with mem
        map_long_historical(split_fq_by_length.out.longreads
            .combine(bwa_index)
            )
        // merge the resulting bams
        historical_bams_ch = merge_short_and_long_read_bams(
            bwa_aln.out.bam
                .mix(map_long_historical.out.bam)
                .groupTuple(by: [0,1,2])
                .map { sample_id, library, datatype, bam_paths, _bam_indices ->
                    return tuple(sample_id, library, datatype, bam_paths)
                }
        )
        historical_bams_ch.view {
            row -> "historical bams after merging short and long read mappings: ${row}"
        } 
    }
    else {
        // otherwise just prep the merged reads for mapping with the specified mapper
         historical_merged_to_map = clean_merged
            .filter { _sample_id, _library, datatype, _collapsed -> datatype == '2' }
            .combine(bwa_index)
        if (historical_mapper == "bwa_aln") {
            historical_bams_ch = bwa_aln(historical_merged_to_map)
        }
        else if (historical_mapper == "bwa_mem") {
            historical_bams_ch = map_merged(historical_merged_to_map)
        }
        else {
            error "Invalid option for historical_mapper: ${historical_mapper}. Must be either 'split', 'bwa_aln' or 'bwa_mem'."
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // If mapping historical paired end reads too, merge within each library
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
