process samtools_markdups {
    tag "$sample_id"
    label 'process_wide'
    conda "${moduleDir}/environment.yml"

    // publishDir "${params.outdir}/02_cramfiles", mode: 'copy'
    input:
    tuple val(sample_id), val(datatype), path(bam), path(bam_idx)

    output:
    tuple val(sample_id), val(datatype), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"), emit: metrics

    script:
    """
    samtools collate -@ ${params.samtools_md_threads} -O -u ${bam} | \
    samtools fixmate -@ ${params.samtools_md_threads} -m -u - - | \
    samtools sort -@ ${params.samtools_md_threads} -u - | \
    samtools markdup -@ ${params.samtools_md_threads} -r -f ${sample_id}.dedup_metrics.txt -O BAM - ${sample_id}.dedup.bam
    #--use-read-groups 
    samtools index ${sample_id}.dedup.bam

    """

    stub:
    """
    touch ${sample_id}.dedup.bam
    touch ${sample_id}.dedup.bam.bai
    touch ${sample_id}.dedup_metrics.txt
    """
}
