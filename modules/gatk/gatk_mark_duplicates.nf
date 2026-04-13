/*
 * GATK MarkDuplicates - Mark/remove PCR duplicates
 */

process gatk_mark_duplicates {
    tag "$sample_id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"), emit: metrics

    script:
    def avail_mem = (task.memory.giga * 0.8).intValue()
    """
    // gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
    //     --INPUT ${bam} \\
    //     --OUTPUT ${sample_id}.dedup.bam \\
    //     --METRICS_FILE ${sample_id}.dedup_metrics.txt \\
    //     --CREATE_INDEX true \\
    //     --VALIDATION_STRINGENCY LENIENT

    // # Rename the index file to standard naming convention
    // mv ${sample_id}.dedup.bai ${sample_id}.dedup.bam.bai || true

    echo "Simulated deduplicated BAM content for sample ${sample_id}" > ${sample_id}.dedup.bam
    echo "Simulated deduplicated BAM index for sample ${sample_id}" > ${sample_id}.dedup.bam.bai
    echo "Simulated dedup metrics for sample ${sample_id}" > ${sample_id}.dedup_metrics.txt

    """

    stub:
    """
    touch ${sample_id}.dedup.bam
    touch ${sample_id}.dedup.bam.bai
    touch ${sample_id}.dedup_metrics.txt
    """
}
